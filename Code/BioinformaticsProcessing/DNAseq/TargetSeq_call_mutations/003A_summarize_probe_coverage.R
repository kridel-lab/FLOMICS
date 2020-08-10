#----------------------------------------------------------------------
#Summarize probe coverage from picard
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
library(readxl)

#----------------------------------------------------------------------
#Load sample info
#----------------------------------------------------------------------

#Sample info (originally provided by BC team in first data upload)
samp_info = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/library_mapping_BC.csv")

#detailed sample info (provided by Anjali and Robert)
more_samp_info = as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/sample_annotations_rcd6Nov2019.xlsx"))

#----------------------------------------------------------------------
#Interval files provided by IDT and prepared by Robert
#----------------------------------------------------------------------

#coding
c_amp = fread("coding_genes_probe_coordinates_n_1917.txt")
c_amp$V1 = sapply(c_amp$V1, function(x){unlist(strsplit(x, "chr"))[2]})
c_amp = c_amp[,c(1:4)]
c_amp$type = "coding"

c_target = fread("coding_genes_target_coordinates_n_838.txt")
c_target$V1 = sapply(c_target$V1, function(x){unlist(strsplit(x, "chr"))[2]})
c_target$gene = sapply(c_target$V4, function(x){unlist(strsplit(x, "\\("  ))[2]})
c_target$gene = sapply(c_target$gene, function(x){unlist(strsplit(x, "\\)"  ))[1]})
c_target = c_target[,c(1:4,7)]
c_target$type = "coding"

#noncoding
nc_amp = fread("non_coding_probe_coordinates_n_1481.txt")
nc_amp = nc_amp[,c(3:5,1)]
nc_amp$type = "noncoding"
colnames(nc_amp) = colnames(c_amp)
#nc_amp$V1 = paste("chr", nc_amp$V1, sep="")
nc_target = fread("non_coding_target_coordinates_n_35.txt")
nc_target$V1 = sapply(nc_target$V1, function(x){unlist(strsplit(x, "chr"))[2]})
nc_target$gene = sapply(nc_target$V4, function(x){unlist(strsplit(x, "_"))[3]})
nc_target = nc_target[,c(1:4,7)]
nc_target$type = "noncoding"

#collect all amplicon intervals
all_amps = rbind(c_amp, nc_amp)
colnames(all_amps) = c("chr", "start", "stop", "region", "type")

#collect all target intervals
all_targets = rbind(c_target, nc_target)
colnames(all_targets) = c("chr", "start", "stop", "region", "gene","type")

#----------------------------------------------------------------------
#load in output files from picard tools
#----------------------------------------------------------------------

#all files with coverage per target
all_res = list.files("/cluster/projects/kridelgroup/FLOMICS/DATA/BC_TargetSeq_Calls", pattern="coverage")

#read in all files and append sample names
read_file = function(file_name){

  f = paste("/cluster/projects/kridelgroup/FLOMICS/DATA/BC_TargetSeq_Calls/", file_name, sep="")
  f = fread(f)
  sample = unlist(strsplit(file_name, "_"))[1]
  sample_check = unlist(strsplit(file_name, "-"))[1]

    if(!(sample_check == file_name)){
      sample = sample_check
    }

  f$Library = sample
  return(f)
}

all_res = as.data.table(ldply(llply(all_res, read_file)))
all_res$id = paste(all_res$chrom, all_res$end, sep="_")
all_targets$id = paste(all_targets$chr, all_targets$stop, sep="_")

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#merge information about target region to picard output file
all_res = merge(all_res, all_targets, by = "id")
colnames(all_res)[3] = "start_picard"
colnames(all_res)[18] = "start_target"

#merge with sample name
all_res = merge(all_res, samp_info, by="Library")

#merge with detaild sample information
colnames(more_samp_info)[1] = "External_identifier"
all_res = merge(all_res, more_samp_info, by = "External_identifier")
all_res$RNAseq_COMMENT=NULL

#save all_res
write.csv(all_res, file=paste(date,
  "picard_tools_coverage_summary_targets_DNA_sequencing.csv", sep="_"), quote=F, row.names=F)

#calculate mean coverage between coding and non-coding
all_res %>% group_by(type) %>% dplyr::summarize(mean_cov = mean(mean_coverage))
