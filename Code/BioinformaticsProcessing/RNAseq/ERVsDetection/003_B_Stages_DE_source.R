#----------------------------------------------------------------------------------
#telescope_DE_analysis_stages_source.R
#----------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------
#FUNCTIONS
#----------------------------------------------------------------------------------
#make transcripts in rows, samlpes in columns

vert_format=function(x){
  telescope_vert = as.data.frame(dcast(x, transcript ~ sample, value.var = "final_count"))
  rownames(telescope_vert) = telescope_vert$transcript
  telescope_vert$transcript = NULL
  return(telescope_vert)
}
#vert_alltiers_telescope = llply(alltiers_telescope, vert_format)


#format for EdgeR,remove NAs
filter_removeNA = function(file){
  x <- file
  x[is.na(x)] <- 0
  tier_sample_info = as.data.table(filter(sample_info, rna_seq_file_sample_ID %in% colnames(x), STAGE %in% c("ADVANCED", "LIMITED")))
  z = which(colnames(x) %in% tier_sample_info$rna_seq_file_sample_ID)
  x = x[,z]
  tier_sample_info = tier_sample_info[order(match(rna_seq_file_sample_ID, colnames(x)))]
  print(tier_sample_info$rna_seq_file_sample_ID == colnames(x))
  return(x)
}
#xtiers=llply(vert_alltiers_telescope,filter_removeNA)

groups_on_tiers = function(file){
  x <- file
  x[is.na(x)] <- 0
  tier_sample_info = as.data.table(filter(sample_info, rna_seq_file_sample_ID %in% colnames(x), STAGE %in% c("ADVANCED", "LIMITED")))
  z = which(colnames(x) %in% tier_sample_info$rna_seq_file_sample_ID)
  x = x[,z]
  tier_sample_info = tier_sample_info[order(match(rna_seq_file_sample_ID, colnames(x)))]
  group=tier_sample_info$STAGE
  return(group)
}
#group_tiers=llply(vert_alltiers_telescope,groups_on_tiers)

#create DGEList object, normalize
#tier_numbers=c(1,2,3)
make_DGEs = function(whichtier){
  df = DGEList(counts=xtiers[[whichtier]], group=group_tiers[[whichtier]])
  df_samples=df$samples
  print(summary(df_samples$lib.size))

 #normalize, filter lowly expressed genes
  df <- calcNormFactors(df)
  sums = apply(df$counts, 1, sum)
  z1 = which(sums > 500)
  z2 = which(sums < 10000000)
  keep=unique(names(sums)[c(z1,z2)])
  df <- df[keep, keep.lib.sizes = FALSE]
  return(df)
  }
#df_tiers=llply(tier_numbers,make_DGEs)

#convert_CPM=function(dftiers){
#  exprs.df <- data.frame(cpm(dftiers, normalized.lib.sizes = FALSE, log = FALSE))
#  return(exprs.df)
#  ###should I write this out and to what directory?}
#exprs.df_tiers=llply(df_tiers,convert_CPM)

# Set up pairwise comparisons
#construct design matrix
#tier_numbers=c(1,2,3)
design_matrix=function(whichtier){
    design <- model.matrix(~0 + group_tiers[[whichtier]], data=df_tiers[[whichtier]]$samples)
    colnames(design)=levels(df_tiers[[whichtier]]$samples$group)
    return(design)
}
#design_tiers=llply(tier_numbers,design_matrix)

#estimate dispersion
#tier_numbers=c(1,2,3)
get_disp_fit=function(whichtier){
  df_disp_tiers <- estimateDisp(df_tiers[[whichtier]], design_tiers[[whichtier]])
  #plotBCV(df)
  fit <- glmQLFit(df_disp_tiers, design_tiers[[whichtier]])
  return(fit)
  }
#fit_tiers=llply(tier_numbers,get_disp_fit)

#make contrasts
get_myconstrasts=function(whichtier){
  my.contrasts <- makeContrasts(ADVANCED_LIMITED = ADVANCED-LIMITED, levels = design_tiers[[whichtier]])
  print(my.contrasts)
  return(my.contrasts)
  }
#my.contrasts_tier=llply(tier_numbers,get_myconstrasts)

#get differential expression results summary from all contrasts
get_allcontrasts=function(whichtier){
  all_contrasts = colnames(my.contrasts_tier[[whichtier]])
  print(all_contrasts)
  return(all_contrasts)
}
#all_contrasts_tier = llply(tier_numbers,get_allcontrasts)

get_res = function(whichtier){
	contrast_dat <- glmQLFTest(fit_tiers[[whichtier]], contrast = my.contrasts_tier[[1]][,"ADVANCED_LIMITED"])
	contrast_dat <- data.frame(topTags(contrast_dat, n = Inf)) %>%
  		mutate(transcript = row.names(.)) %>%
  			filter(abs(logFC) > x & FDR < y)
  	if(!(dim(contrast_dat)[1] == 0)){
	contrast_dat$contrast = "ADVANCED_LIMITED"
	print("done")
	return(contrast_dat)
}
}
#all_de_ervs_cluster = llply(tier_numbers, get_res)
#colnames(all_de_ervs)[6] = "transcript"

#all_de_ervs = join(all_de_ervs,telescope_annotations)
#not sure what is meant by "telescope_annotations" above <-- this can be done
#to add back all the ERV coordinates and transcript info to main result file
#but not necessary!

#save results and upload to OneDrive
