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
c_amp = c_amp[,c(1:3)]
c_amp$V1 = sapply(c_amp$V1, function(x){unlist(strsplit(x, "chr"))[2]})
c_target = fread("coding_genes_target_coordinates_n_838.txt")
c_target = c_target[,c(1:3)]
c_target$V1 = sapply(c_target$V1, function(x){unlist(strsplit(x, "chr"))[2]})

#noncoding
nc_amp = fread("non_coding_probe_coordinates_n_1481.txt")
nc_amp = nc_amp[,c(3:5)]
colnames(nc_amp) = colnames(c_amp)
#nc_amp$V1 = paste("chr", nc_amp$V1, sep="")
nc_target = fread("non_coding_target_coordinates_n_35.txt")
nc_target = nc_target[,c(1:3)]
nc_target$V1 = sapply(nc_target$V1, function(x){unlist(strsplit(x, "chr"))[2]})

#collect all amplicon intervals
all_amps = rbind(c_amp, nc_amp)
colnames(all_amps) = c("chr", "start", "stop")
write.table(all_amps, file="picard_tools_amps_input.bed", quote=F,
row.names=F, col.names=F, sep="\t")

#collect all target intervals
all_targets = rbind(c_target, nc_target)
colnames(all_targets) = c("chr", "start", "stop")
write.table(all_targets, file="picard_tools_targets_input.bed", quote=F,
row.names=F, col.names=F, sep="\t")
