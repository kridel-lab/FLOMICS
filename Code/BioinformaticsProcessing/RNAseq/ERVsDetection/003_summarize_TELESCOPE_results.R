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

setwd("/cluster/projects/kridelgroup/FLOMICS/TELESCOPE_ANALYSIS/concatenated_results") #or where ever the 147 tsv files are stored

#----------------------------------------------------------------------------------
#DATA
#----------------------------------------------------------------------------------

results = list.files(pattern=".tsv")

sample_info = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/sample_annotations_rcd6Nov2019.txt")
telescope_annotations = fread("/cluster/home/kisaev/data/telescope_hg19.gtf.bed")

#clean up 10th colunmn with ERV info
telescope_annotations = telescope_annotations %>% separate("V10", into=c("x", "xb", "gene_id", "xx", "transcript_id",
  "xxx", "locus", "xxxx", "repName", "xxxxx", "geneRegion", sep=" "))
telescope_annotations = telescope_annotations[,c(1:4, 12:13)]
colnames(telescope_annotations) = c("Chr", "Start", "End", "transcript", "ERV_ID", "locus")

#read-in all files and assemble into one data-table
get_res = function(file){
	f=fread(file)
	f$sample = unlist(strsplit(file, "Aligned"))[1]
	return(f)
}

all_telescope = as.data.table(ldply(llply(results, get_res)))
all_telescope = as.data.table(filter(all_telescope, !(transcript== "__no_feature")))
colnames(all_telescope)[12] = "SAMPLE_ID"

#summary of how many ERVs were evaluated for each sample
summ = as.data.table(table(all_telescope$transcript))
summ = summ[order(-N)]
#summ = as.data.table(filter(summ, N > 75)) #only keep those ERVs that have been detected in at least 55 samples
all_telescope = as.data.table(filter(all_telescope, transcript %in% summ$V1))

#re-arrange so that samples are in COLUMNS and transcript names are in one column
all_telescope_vertical = as.data.frame(dcast(all_telescope, transcript ~ SAMPLE_ID, value.var = "final_count"))
all_telescope = join(all_telescope, sample_info)

rownames(all_telescope_vertical) = all_telescope_vertical$transcript
all_telescope_vertical$transcript = NULL
write.csv(all_telescope, paste(date, "TELESCOPE_OUTPUT_WITH_SAMPLE_ANNOTATION.csv", sep="_"), quote=F, row.names=F)

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
sample_info = as.data.table(filter(sample_info, SAMPLE_ID %in% colnames(x), STAGE %in% c("ADVANCED", "LIMITED")))
z = which(colnames(x) %in% sample_info$SAMPLE_ID)
x = x[,z]

#reorder so same order of patients as in counts matrix
sample_info = sample_info[order(match(SAMPLE_ID, colnames(x)))]
group= sample_info$STAGE

#confirm everything is in right order
sample_info$SAMPLE_ID == colnames(x)

#create DGEList object
df <- DGEList(counts=x,group=group)
df_samples = as.data.frame(df$samples)
summary(df_samples$lib.size)

# Normalize.
# Normalization may not actually be required:
# Normalization issues arise only to the extent that technical factors
# have sample-specific effects.
df <- calcNormFactors(df)

# Examine the samples for outliers:
# Plot in which distances between samples correspond to leading biological
# coefficient of variation (BCV) between those samples.
pdf("MDS.pdf")
plotMDS(df)
dev.off()

# Filter out lowly expressed genes
sums = apply(df$counts, 1, sum)
z1 = which(sums > 1000)
z2 = which(sums < 25000)
keep=unique(names(sums)[c(z1,z2)])
df <- df[keep, keep.lib.sizes = FALSE]

# Write out expression matrix, annotated with gene symbols rather than ENSG identifiers
exprs.df <- data.frame(cpm(df, normalized.lib.sizes = FALSE, log = FALSE)) # Computes counts per million (CPM) values
exprs.df$ERV_ids <- row.names(exprs.df)
write.table(exprs.df, "CPM_exprs_matrix_norm_filt_ERVs_FL_109_cases.txt", sep = "\t", row.names = F, quote=F)
write.csv(exprs.df, "CPM_exprs_matrix_norm_filt_ERVs_FL_109_cases.csv", row.names = F, quote=F)

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
y <- 0.1 #FDR value

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

write.csv(all_de_ervs, file=paste(date, "Telescope_ERVs_differentially_ADVANCED_vs_LIMITED.csv", sep="_"), quote=F, row.names=F) #requierd columns 7, 1, 5
