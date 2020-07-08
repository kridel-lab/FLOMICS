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

results = list.files(pattern=".tsv")
#read-in all files and assemble into one data-table
get_res = function(file){
	f=fread(file)
	f$sample = unlist(strsplit(file, "Aligned"))[1]
  print("done")
	return(f)
}

all_telescope = as.data.table(ldply(llply(results, get_res)))

#information regarding each sample and which stage of disease and cluster they are part of
sample_info = as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.xlsx"))

#telescope annotation file cleaned up - ERVs family only
telescope_annotations = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/hg19_clean_repeatmasker_annotations.bed")

#only keep ERVs from Telescope Annotation file in the main Telescope results
all_telescope = as.data.table(filter(all_telescope, transcript %in% telescope_annotations$V10))

#save final pre-differential expression analysis ERV count matrix for all 136 samples
write.csv(all_telescope, paste("/cluster/projects/kridelgroup/FLOMICS/DATA/", date, "TELESCOPE_OUTPUT_WITH_SAMPLE_ANNOTATION.csv", sep="_"), quote=F, row.names=F)
