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
#library(tidyverse)
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
tiers=c("tier1", "tier2", "tier3")
#apply downstream code to each tier and then compare results
get_telescope_tier_summary = function(tier){
  if(tier == "tier3"){
    dat=sample_info
    print(length(unique(dat$SAMPLE_ID)))}
  if(tier == "tier2"){
    dat=sample_info
    z1 = which(dat$rRNAcontam > 40)
    z2 = which(dat$Uniquely.mapped < 10000000)
    z3 = which(dat$'Uniquely.mapped.reads..' < 50)
    z4 = which(dat$'X..of.reads.mapped.to.multiple.loci' > 20)
    total_lost = unique(c(z1, z2, z3, z4))
    dat=sample_info[-total_lost]
    print(length(unique(dat$SAMPLE_ID)))}
  if(tier == "tier1"){
    dat=sample_info
    z1 = which(dat$rRNAcontam > 40)
    z2 = which(dat$Uniquely.mapped < 10000000)
    z3 = which(dat$'Uniquely.mapped.reads..' < 70)
    z4 = which(dat$'X..of.reads.mapped.to.multiple.loci' > 20)
    total_lost = unique(c(z1, z2, z3, z4))
    dat=sample_info[-total_lost]
    print(length(unique(dat$SAMPLE_ID)))}
    tier_telescope_filter=as.data.table(filter(all_telescope, sample %in% dat$rna_seq_file_sample_ID))
    print("done tier analysis")
    return(tier_telescope_filter)
}
alltiers_telescope = llply(tiers, get_telescope_tier_summary)

#length(unique(alltiers_telescope[[1]]$sample))=81, ncol=13352988
  #length(unique(alltiers_telescope[[1]]$transcript))=519257
#length(unique(alltiers_telescope[[2]]$sample))=104, nol=14967961
  #length(unique(alltiers_telescope[[1]]$transcript))=545956
#length(unique(alltiers_telescope[[3]]$sample))=132, nocl=18300436
  ##length(unique(alltiers_telescope[[1]]$transcript))=579881

vert_format=function(x){
  telescope_vert = as.data.frame(dcast(x, transcript ~ sample, value.var = "final_count"))
  rownames(telescope_vert) = telescope_vert$transcript
  telescope_vert$transcript = NULL
  return(telescope_vert)
}
vert_alltiers_telescope = llply(alltiers_telescope, vert_format)
#ncol(vert_alltiers_telescope[[1]])=81, nrow(vert_alltiers_telescope[[1]])=519257
#ncol(vert_alltiers_telescope[[2]])=104, nrow(vert_alltiers_telescope[[2]])=545956
#ncol(vert_alltiers_telescope[[3]])=132, nrow(vert_alltiers_telescope[[3]])=579881
#----------------------------------------------------------------------------------
#Differential expression analysis using EdgeR adapted from RNA-Seq code
#----------------------------------------------------------------------------------

#x <- all_telescope_vertical
#for now replace all NAs with 0s so have full picture later reconsider if should actually remove those genes
#with NA median values
#meds = apply(x, 1, median)
#z = which(is.na(meds)) #13636 transcripts had at least one NA in one sample so these were removed
#x = x[-z,]
#x[is.na(x)] <- 0

#get groups
#sample_info = as.data.table(filter(sample_info, rna_seq_file_sample_ID %in% colnames(x), STAGE %in% c("ADVANCED", "LIMITED")))
#z = which(colnames(x) %in% sample_info$rna_seq_file_sample_ID)
#x = x[,z]

#reorder so same order of patients as in counts matrix
#sample_info = sample_info[order(match(rna_seq_file_sample_ID, colnames(x)))]
#group= sample_info$STAGE

#confirm everything is in right order
#sample_info$rna_seq_file_sample_ID == colnames(x)

#format for EdgeR
filter_removeNA = function(file){
  x <- file
  x[is.na(x)] <- 0

  #don't use the same name as the actual dataframe here  - sample info could run
  #into trouble

  #sample_info = as.data.table(filter(sample_info, rna_seq_file_sample_ID %in% colnames(x), STAGE %in% c("ADVANCED", "LIMITED")))
  tier_sample_info = as.data.table(filter(sample_info, rna_seq_file_sample_ID %in% colnames(x), STAGE %in% c("ADVANCED", "LIMITED")))
  z = which(colnames(x) %in% tier_sample_info$rna_seq_file_sample_ID)
  x = x[,z]

  tier_sample_info = tier_sample_info[order(match(rna_seq_file_sample_ID, colnames(x)))]

  ##
  group=tier_sample_info$STAGE #group not really needed here because doesn't get used here
  ##

  print(tier_sample_info$rna_seq_file_sample_ID == colnames(x))
  return(x)
}
xtiers=llply(vert_alltiers_telescope,filter_removeNA)
#ncol(xtiers[[1]]) = 71, nrow(xtiers[[1]])=519257
#ncol(xtiers[[2]])= 94, nrow(xtiers[[2]])=545956
#ncol(xtiers[[3]]) = 121, nrow(xtiers[3])=579881


groups_on_tiers = function(file){
  x <- file
  x[is.na(x)] <- 0

  #don't use the same name as the actual dataframe here  - sample info could run
  #into trouble

  #sample_info = as.data.table(filter(sample_info, rna_seq_file_sample_ID %in% colnames(x), STAGE %in% c("ADVANCED", "LIMITED")))
  tier_sample_info = as.data.table(filter(sample_info, rna_seq_file_sample_ID %in% colnames(x), STAGE %in% c("ADVANCED", "LIMITED")))

  z = which(colnames(x) %in% tier_sample_info$rna_seq_file_sample_ID)
  x = x[,z]
  tier_sample_info = tier_sample_info[order(match(rna_seq_file_sample_ID, colnames(x)))]

  group=tier_sample_info$Cluster
  return(as.character(group))
}
group_tiers=llply(vert_alltiers_telescope,groups_on_tiers)


#length(group_tiers[[1]])=71
#length(group_tiers[[2]])=94
#length(group_tiers[[3]])=121
##################################################
#create DGEList object
#df <- DGEList(counts=x,group=group)
#df_samples = as.data.frame(df$samples)
#summary(df_samples$lib.size)

# Normalize.
# Normalization may not actually be required:
# Normalization issues arise only to the extent that technical factors
# have sample-specific effects.
#df <- calcNormFactors(df)

# Filter out lowly expressed genes - feel free to change this
#this isn't necessarily best way to do it

#sums = apply(df$counts, 1, sum)
#z1 = which(sums > 1000)
#z2 = which(sums < 25000)
#keep=unique(names(sums)[c(z1,z2)])
#df <- df[keep, keep.lib.sizes = FALSE]

#create DGEList object, normalize
tier_numbers=c(1,2,3)
make_DGEs = function(tier){
  df = DGEList(counts=xtiers[[tier]],
    group=group_tiers[[tier]])
  df_samples=as.data.frame(df$samples)
  print(summary(df_samples$lib.size))

 #normalize, filter lowly expressed genes
  df <- calcNormFactors(df)
  sums = apply(df$counts, 1, sum)
  z1 = which(sums > 500) #this can be changed`
  z2 = which(sums < 10000000) #this can be changed
  keep=unique(names(sums)[c(z1,z2)])
  df <- df[keep, keep.lib.sizes = FALSE]
  return(df)
  }
df_tiers=llply(tier_numbers,make_DGEs)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 815622 3718738 4160641 4031534 4421368 5490691
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 305584 2833351 3942428 3487440 4356880 7414654
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  67010 2170631 3767050 3233262 4246633 7414654

#nrow(df_tiers[[1]]$samples)=71,nrow(df_tiers[[3]]$counts)=519257
#nrow(df_tiers[[2]]$samples)=94,nrow(df_tiers[[3]]$counts)=545956
#nrow(df_tiers[[3]]$samples)=121, nrow(df_tiers[[3]]$counts)=579881

#Write out expression matrix, annotated with gene symbols rather than ENSG identifiers
  #we don't need this since no ENSG identifiers - do we still want a write out of CPM conversion?
#exprs.df <- data.frame(cpm(df, normalized.lib.sizes = FALSE, log = FALSE)) # Computes counts per million (CPM) values
#exprs.df$ERV_ids <- row.names(exprs.df) # don't need this now

convert_CPM=function(dftiers){
  exprs.df <- data.frame(cpm(dftiers, normalized.lib.sizes = FALSE, log = FALSE))
  return(exprs.df)
  ###should I write this out and to what directory?
}
exprs.df_tiers=llply(df_tiers,convert_CPM)

#---
# Set up comparisons----
# Only pairwise comparisons at this point, within same cell line
# ---

# Construct design matrix
#design <- model.matrix(~0 + group, data = df$samples)
#colnames(design) <- levels(df$samples$group)
#design

#construct design matrix
#tier_numbers=c(1,2,3)
design_matrix=function(tiernumbers){
    design <- model.matrix(~0 + group_tiers[[tiernumbers]],
      data=df_tiers[[tiernumbers]]$samples)
    colnames(design)=c("Cluster1", "Cluster2")
    return(design)
}
design_tiers=llply(tier_numbers,design_matrix)

# Estimate dispersion.
# Allowing gene-specific dispersion is necessary in order that differential
# expression is not driven by outliers. Therefore the tagwise dispersions are
# strongly recommended in model fitting and testing for differential expression.
# This method estimates common dispersion, trended dispersions and tagwise dispersions
# in one run and is recommended.
#df <- estimateDisp(df, design) #should we be using estimateTagwiseDisp(x)?

# The dispersion estimates can be viewed in a BCV plot
#plotBCV(df)

# Fit generalized linear model
# Such a model is an extension of classical linear models to non-normally distributed response data.
#fit <- glmQLFit(df, design)

#estimate dispersion
#tier_numbers=c(1,2,3)
get_disp_fit=function(tiernumbers){
  df_tiers[[tiernumbers]] <- estimateDisp(df_tiers[[tiernumbers]],
    design_tiers[[tiernumbers]])
  #plotBCV(df)
  fit <- glmQLFit(df_tiers[[tiernumbers]], design_tiers[[tiernumbers]])
  return(fit)}
fit_tiers=llply(tier_numbers,get_disp_fit)

# Make contrasts
#my.contrasts <- makeContrasts(ADVANCED_LIMITED = ADVANCED-LIMITED,
#                              levels = design)

#make contrasts
get_myconstrasts=function(tier_numbers){ #dont make input into function same
  #as the list you are applying function to!
  my.contrasts <- makeContrasts(Cluster1_Cluster1 = Cluster1-Cluster2,
                              levels = design_tiers[[tier_numbers]])
  return(my.contrasts)
  }
my.contrasts_tier=llply(tier_numbers,get_myconstrasts)

#get differential expression results summary from all contrasts
#all_contrasts = colnames(my.contrasts)

#get differential expression results summary from all contrasts
get_allcontrasts=function(tier_numbers){ #same here 
  all_contrasts = colnames(my.contrasts_tier[[tier_numbers]])
  return(all_contrasts)
}
all_contrasts_tier = llply(tier_numbers,get_allcontrasts)

# Set cut-offs for logFC and FDR
x <- 1 # logFCx
y <- 1 #keep all p-values for now can filter significant hits later

get_res = function(contrast){
	print(contrast)
	contrast_dat <- glmQLFTest(fit_tiers[[]], contrast = my.contrasts[[]][,1])

	contrast_dat <- data.frame(topTags(contrast_dat, n = Inf)) %>%
  		mutate(ensembl_gene_id = row.names(.)) %>%
  			filter(abs(logFC) > x & FDR < y)
  	if(!(dim(contrast_dat)[1] == 0)){
	contrast_dat$contrast = contrast
	print("done")
	return(contrast_dat)
}
}

all_de_ervs = as.data.table(llply(tier_numbers, get_res)))
colnames(all_de_ervs)[6] = "transcript"

#could just replace "mutate(ensembl_gene_id" with "mutate(transcript"
 #and exclude the above line?

#all_de_ervs = join(all_de_ervs,telescope_annotations)
#not sure what is meant by "telescope_annotations" above <-- this can be done
#to add back all the ERV coordinates and transcript info to main result file
#but not necessary!

#save results and upload to OneDrive
write.csv(all_de_ervs_cluster[[1]], paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope", date, "tier1_de_ervs_cluster.csv", sep="_"), quote=F, row.names=F)
write.csv(all_de_ervs_cluster[[2]], paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope", date, "tier2_de_ervs_cluster.csv", sep="_"), quote=F, row.names=F)
write.csv(all_de_ervs_cluster[[3]], paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope", date, "tier3_de_ervs_cluster.csv", sep="_"), quote=F, row.names=F)
