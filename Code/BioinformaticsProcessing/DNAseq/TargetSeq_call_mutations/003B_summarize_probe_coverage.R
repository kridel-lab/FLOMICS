#----------------------------------------------------------------------
#Summarize probe coverage from picard - part two
#----------------------------------------------------------------------

#Karin Isaev
#Started August 5th 2020
#tested on R version 3.5.0

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#this script takes output from picard collect targeted pcr metrics
#merges it with sample information
#make summary csv file with coverage per gene per patient
#simple barplot summaries of coverage per gene/probe/region
#calculate differences in coverage between coding and non-coding regions

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)
date=Sys.Date()

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
"plyr",
"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)


#----------------------------------------------------------------------
#Load most up to date data matrix from previous script "all_res"
#----------------------------------------------------------------------

all_res = fread(list.files(pattern="picard_tools_coverage")[length(list.files(pattern="picard_tools_coverage"))])
f = fread(list.files(pattern="picard_tools_coverage")[length(list.files(pattern="picard_tools_coverage"))])

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#calculate mean coverage between coding and non-coding
all_res %>% group_by(type) %>% dplyr::summarize(mean_cov = mean(mean_coverage))

#plot barplot summary per gene coverage
