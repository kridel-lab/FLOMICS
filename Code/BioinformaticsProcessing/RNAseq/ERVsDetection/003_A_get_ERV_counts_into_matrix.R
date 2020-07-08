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
saveRDS(all_telescope, file=paste("/cluster/projects/kridelgroup/FLOMICS/DATA/", date, "all_telescope_results_matrix.rds", sep="_"))
