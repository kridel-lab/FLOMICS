#----------------------------------------------------------------------
#prepare intervals for picard collect metrics tool
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
"plyr",
"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#Interval files provided by IDT and processed by Robert
#----------------------------------------------------------------------

#coding
c_amp = fread("coding_genes_probe_coordinates_n_1917.txt")
c_amp$V1 = sapply(c_amp$V1, function(x){unlist(strsplit(x, "chr"))[2]})
c_amp = c_amp[,c(1:4)]
c_amp$type = "coding"

c_target = fread("coding_genes_target_coordinates_n_838.txt")
c_target$V1 = sapply(c_target$V1, function(x){unlist(strsplit(x, "chr"))[2]})
c_target = c_target[,c(1:4)]
c_target$type = "coding"

#noncoding
nc_amp = fread("non_coding_probe_coordinates_n_1481.txt")
nc_amp = nc_amp[,c(3:5,1)]
nc_amp$type = "noncoding"
colnames(nc_amp) = colnames(c_amp)
#nc_amp$V1 = paste("chr", nc_amp$V1, sep="")
nc_target = fread("non_coding_target_coordinates_n_35.txt")
nc_target$V1 = sapply(nc_target$V1, function(x){unlist(strsplit(x, "chr"))[2]})
nc_target = nc_target[,c(1:4)]
nc_target$type = "noncoding"

#collect all amplicon intervals
all_amps = rbind(c_amp, nc_amp)
colnames(all_amps) = c("chr", "start", "stop", "region", "type")

#collect all target intervals
all_targets = rbind(c_target, nc_target)
colnames(all_targets) = c("chr", "start", "stop")

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------


#output files here:
"/cluster/projects/kridelgroup/FLOMICS/DATA/BC_TargetSeq_Calls"
"_target_coverage.txt"
