#----------------------------------------------------------------------------------
#telescope_output_differential_expression_analysis.R
#----------------------------------------------------------------------------------

#Karin Isaev
#Date: January 16th, 2020
#This script takes in individual files obtained by telescope
#for each individual BAM file from RNA-Seq FL TGL13 samples
#and conducts EdgeR differential expression analysis of these trancripts

#----------------------------------------------------------------------------------
#PACKAGES
#----------------------------------------------------------------------------------

date = Sys.Date()

library(data.table)
library(dplyr)
library(plyr)
library(reshape2)
library(edgeR)
library(tidyverse)
library(readxl)

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/concatenated_results") #or where ever the 136 tsv files are stored

#----------------------------------------------------------------------------------
#DATA
#----------------------------------------------------------------------------------

#The file below was created by script 003_A_get_ERV_counts_into_matrix.R
#It read in all 136 tsv files and put it together into one dataframe
#which we can use further
all_telescope = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/_2020-07-08_TELESCOPE_OUTPUT_WITH_SAMPLE_ANNOTATION.csv")

#information regarding each sample and which stage of disease and cluster they are part of
sample_info = as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.xlsx"))

#subset telescope results based on QC tiers
tiers=c("tier3", "tier2", "tier1")
#apply downstream code to each tier and then compare results
get_tier_summary = function(tier){
  if(tier == "tier3"){
    dat=sample_info
    print(length(unique(dat$SAMPLE_ID)))
    T3_all_telescope=as.data.table(filter(all_telescope, sample %in% dat$rna_seq_file_sample_ID))
    write.csv(T3_all_telescope, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/tiers", date,"T3_all_telescope.csv", sep="_"), quote=F, row.names=F)}
  if(tier == "tier2"){
    dat=sample_info
    z1 = which(dat$rRNAcontam > 40)
    z2 = which(dat$Uniquely.mapped < 10000000)
    z3 = which(dat$'Uniquely.mapped.reads..' < 50)
    z4 = which(dat$'X..of.reads.mapped.to.multiple.loci' > 20)
    total_lost = unique(c(z1, z2, z3, z4))
    dat=sample_info[-total_lost]
    print(length(unique(dat$SAMPLE_ID)))
    T2_all_telescope=as.data.table(filter(all_telescope, sample %in% dat$rna_seq_file_sample_ID))
    write.csv(T2_all_telescope, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/tiers", date,"T2_all_telescope.csv", sep="_"), quote=F, row.names=F)}
  if(tier == "tier1"){
    dat=sample_info
    z1 = which(dat$rRNAcontam > 40)
    z2 = which(dat$Uniquely.mapped < 10000000)
    z3 = which(dat$'Uniquely.mapped.reads..' < 70)
    z4 = which(dat$'X..of.reads.mapped.to.multiple.loci' > 20)
    total_lost = unique(c(z1, z2, z3, z4))
    dat=sample_info[-total_lost]
    print(length(unique(dat$SAMPLE_ID)))
    T1_all_telescope=as.data.table(filter(all_telescope, sample %in% dat$rna_seq_file_sample_ID))
    write.csv(T1_all_telescope, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/tiers", date,"T1_all_telescope.csv", sep="_"), quote=F, row.names=F)}
    return(dat)
    print("done tier analysis")
}
all_tiers = as.data.table(ldply(llply(tiers, get_tier_summary)))

T3_telescope = fread("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/tiers/T3_all_telescope.csv")
T2_telescope = fread("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/tiers/T2_all_telescope.csv")
T1_telescope = fread("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/tiers/T1_all_telescope.csv")

#re-arrange so that samples are in COLUMNS and transcript names are in one column
T3_telescope_vert = as.data.frame(dcast(T3_telescope, transcript ~ sample, value.var = "final_count"))
T2_telescope_vert = as.data.frame(dcast(T2_telescope, transcript ~ sample, value.var = "final_count"))
T1_telescope_vert = as.data.frame(dcast(T1_telescope, transcript ~ sample, value.var = "final_count"))


all_telescope = as.data.table(ldply(llply(results, get_res)))
#all_telescope = join(all_telescope, sample_info)
rownames(all_telescope_vertical) = all_telescope_vertical$transcript
all_telescope_vertical$transcript = NULL

#----------------------------------------------------------------------------------
#Differential expression analysis using EdgeR adapted from RNA-Seq code
#----------------------------------------------------------------------------------

x <- all_telescope_vertical
#for now replace all NAs with 0s so have full picture later reconsider if should actually remove those genes
#with NA median values
#meds = apply(x, 1, median)
#z = which(is.na(meds)) #13636 transcripts had at least one NA in one sample so these were removed
#x = x[-z,]
x[is.na(x)] <- 0

#get groups
sample_info = as.data.table(filter(sample_info, rna_seq_file_sample_ID %in% colnames(x), STAGE %in% c("ADVANCED", "LIMITED")))
z = which(colnames(x) %in% sample_info$rna_seq_file_sample_ID)
x = x[,z]

#reorder so same order of patients as in counts matrix
sample_info = sample_info[order(match(rna_seq_file_sample_ID, colnames(x)))]
group= sample_info$STAGE

#confirm everything is in right order
sample_info$rna_seq_file_sample_ID == colnames(x)

#create DGEList object
df <- DGEList(counts=x,group=group)
df_samples = as.data.frame(df$samples)
summary(df_samples$lib.size)

# Normalize.
# Normalization may not actually be required:
# Normalization issues arise only to the extent that technical factors
# have sample-specific effects.
df <- calcNormFactors(df)

# Filter out lowly expressed genes - feel free to change this
#this isn't necessarily best way to do it

sums = apply(df$counts, 1, sum)
z1 = which(sums > 1000)
z2 = which(sums < 25000)
keep=unique(names(sums)[c(z1,z2)])
df <- df[keep, keep.lib.sizes = FALSE]

# Write out expression matrix, annotated with gene symbols rather than ENSG identifiers
exprs.df <- data.frame(cpm(df, normalized.lib.sizes = FALSE, log = FALSE)) # Computes counts per million (CPM) values
exprs.df$ERV_ids <- row.names(exprs.df)

#---
# Set up comparisons----
# Only pairwise comparisons at this point, within same cell line
# ---

# Construct design matrix
design <- model.matrix(~0 + group, data = df$samples)
colnames(design) <- levels(df$samples$group)
design

# Estimate dispersion.
# Allowing gene-specific dispersion is necessary in order that differential
# expression is not driven by outliers. Therefore the tagwise dispersions are
# strongly recommended in model fitting and testing for differential expression.
# This method estimates common dispersion, trended dispersions and tagwise dispersions
# in one run and is recommended.
df <- estimateDisp(df, design)

# The dispersion estimates can be viewed in a BCV plot
#plotBCV(df)

# Fit generalized linear model
# Such a model is an extension of classical linear models to non-normally distributed response data.
fit <- glmQLFit(df, design)

# Make contrasts
my.contrasts <- makeContrasts(ADVANCED_LIMITED = ADVANCED-LIMITED,
                              levels = design)

# Set cut-offs for logFC and FDR
x <- 1 # logFCx
y <- 1 #keep all p-values for now can filter significant hits later

#get differential expression results summary from all contrasts
all_contrasts = colnames(my.contrasts)

get_res = function(contrast){
	print(contrast)
	contrast_dat <- glmQLFTest(fit, contrast = my.contrasts[,contrast])

	contrast_dat <- data.frame(topTags(contrast_dat, n = Inf)) %>%
  		mutate(ensembl_gene_id = row.names(.)) %>%
  			filter(abs(logFC) > x & FDR < y)
  	if(!(dim(contrast_dat)[1] == 0)){
	contrast_dat$contrast = contrast
	print("done")
	return(contrast_dat)
}
}

all_de_ervs = as.data.table(ldply(llply(all_contrasts, get_res)))
colnames(all_de_ervs)[6] = "transcript"
all_de_ervs = join(all_de_ervs,telescope_annotations)

#save results and upload to OneDrive
