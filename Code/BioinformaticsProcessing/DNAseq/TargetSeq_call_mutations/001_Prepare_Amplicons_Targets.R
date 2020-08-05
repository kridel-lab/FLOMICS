#----------------------------------------------------------------------
#prepare intervals for picard collect metrics tool
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#Interval files provided by IDT and processed by Robert
#----------------------------------------------------------------------

#coding
c_amp = fread("coding_genes_probe_coordinates_n_1917.txt")
c_target = fread("coding_genes_target_coordinates_n_838.txt")

#noncoding
