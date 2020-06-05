#----------------------------------------------------------------------
#variants_003_read_in_VCFs_into_matrix.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
  "ggrepel", "stringr", "maftools", "VariantAnnotation")
lapply(packages, require, character.only = TRUE)
library(GenomicRanges)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
library(biomaRt)

date = Sys.Date()

print(date)
#args = commandArgs(trailingOnly = TRUE) #patient ID 
#index = args[1] 
#print(index) 

setwd("/cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder/MANTA_SVs")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?

#----------------------------------------------------------------------
#data 
#----------------------------------------------------------------------

all_SVs = list.files(pattern = "_2019-12-03")

read_f = function(filee){
  f = readRDS(filee)
  return(f)
}

all_SVs = as.data.table(ldply(llply(all_SVs, read_f)))
write.table(all_SVs, paste(date, "all_SVs_samples.txt", sep="_"), quote=F, row.names=F, sep="\t")