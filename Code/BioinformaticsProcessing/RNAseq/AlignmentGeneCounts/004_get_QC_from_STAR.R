#----------------------------------------------------------------------------------
#get QC from STAR
#----------------------------------------------------------------------------------

#Karin Isaev

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
library(hciR)
library(stringr)

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TGL_BAM_RNASEQ_sorted_FASTQ") #or where ever the 136 files are stored

#----------------------------------------------------------------------------------
#DATA
#----------------------------------------------------------------------------------

results = list.files(pattern="Log.final.out")
sample_info = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/sample_annotations_rcd6Nov2019.txt")

star_summary = as.data.table(read_STAR(path = ".", reshape = FALSE))
star_summary$sample = sapply(star_summary$sample, function(x){ unlist(strsplit(x, "Log"))[1]})
star_summary$SAMPLE_ID = star_summary$sample
colnames(star_summary)[1] = "rna_seq_file_sample_ID"

#dcast
star_summary_dcast <- as.data.table(dcast(star_summary, SAMPLE_ID + rna_seq_file_sample_ID ~ stat))

z=which(!(star_summary_dcast$SAMPLE_ID %in% sample_info$SAMPLE_ID))
star_summary_dcast$SAMPLE_ID[z] = paste(star_summary_dcast$SAMPLE_ID[z], "_T1",  sep="")

star_summary_dcast = merge(star_summary_dcast, sample_info, by="SAMPLE_ID")

write.table(star_summary_dcast, file=paste("/cluster/home/kisaev/data/FL_TGL_STAR_logQC", date, "summary_KI.txt", sep="_"), quote=F, row.names=F, sep=";")
