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
probe_space = fread("capseq-space.probe_PLOSMED.bed", skip=1)
probe_space = probe_space[,c(1:3)]
colnames(probe_space) = c("chr", "start", "stop")
write.table(probe_space, file="PLOS_MED_probes_input_picard.bed", quote=F,
row.names=F, col.names=F, sep="\t")
