#----------------------------------------------------------------------------------
#002_mixcr_concatenate_all_clonotypes.R
#----------------------------------------------------------------------------------

#Sarah Russell
#Date: July 29th, 2020
#This script takes in individual "ALL" clonotype files obtained by mixcr
#for each paired reads from RNA-Seq FL TGL13 samples and concatenates into
#one file with sample annotations for stage, type and cluster.

#----------------------------------------------------------------------------------
#PACKAGES
#----------------------------------------------------------------------------------
date = Sys.Date()

library(data.table)
library(dplyr)
library(plyr)
library(tidyverse)
library(readxl)

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS")

#----------------------------------------------------------------------------------
#DATA & ANALYSIS
#----------------------------------------------------------------------------------

results = list.files(pattern="ALL.txt")
###read-in all files and assemble into one data-table
get_res = function(file){
	f=fread(file)
	f$sample = unlist(strsplit(file, ".clonotypes"))[1]
  print("done")
	return(f)
}
all_clones = as.data.table(ldply(llply(results, get_res)))

#> dim(all_clones)
#[1] 49246    36
#> length(unique(all_clones$sample))
#[1] 134

###information regarding each sample and which stage of disease and cluster they are part of
sample_info = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv")

###only keep ERVs from Telescope Annotation file in the main Telescope results
all_clones_filtered = as.data.table(filter(all_clones, sample %in% sample_info$rna_seq_file_sample_ID))

#> dim(all_clones_filtered)
#[1] 48292    36
#> length(unique(all_clones_filtered$sample))
#[1] 130
#132 samples total in sample_info file, but two samples failed to produce clonotypes = 130 samples total

###add stage, type, cluster annotations
for(i in 1:nrow(all_clones_filtered)){
	STAGE[i]=sample_info$STAGE[sample_info$rna_seq_file_sample_ID==all_clones_filtered$sample[i]]
	TYPE[i]=sample_info$TYPE[sample_info$rna_seq_file_sample_ID==all_clones_filtered$sample[i]]
	CLUSTER[i]=sample_info$Cluster[sample_info$rna_seq_file_sample_ID==all_clones_filtered$sample[i]]
	clones_info=cbind(all_clones_filtered,STAGE,TYPE,CLUSTER)
}
#could probably write above in function
dim(clones_info)
#[1] 48292    39
length(unique(clones_info$sample))
#[1] 130

###save final matrix for all MiXCR clonotypes on 132 samples
###seperate by tab since some commas exist within table elements
write.table(clones_info, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS/", date, "MIXCR_OUTPUT_WITH_SAMPLE_ANNOTATION.csv", sep="_"), sep="\t", quote=F, row.names=F)
