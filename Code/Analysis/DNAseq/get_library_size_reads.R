#----------------------------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/PROCESSING/library_sizes")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
"plyr",
"ggrepel", "stringr")
lapply(packages, require, character.only = TRUE)

#get data
all_lib_sizes = list.files(pattern="num_reads_in_bam.txt")

#read through them all
get_lib_size = function(file_lib){

  print(file_lib)

  f = fread(file_lib)
  num = f$V1[1]
  sample = unlist(strsplit(file_lib, "_"))[1]

  if(length(unlist(strsplit(sample, "-"))) > 1){
    sample=unlist(strsplit(sample, "-"))[1]
  }
  info = c(sample, num)
  return(info)
}

all_samples = as.data.table(ldply(llply(all_lib_sizes, get_lib_size)))
colnames(all_samples) = c("Library", "Library_Reads_Total")
write.table(all_samples, file="FL_Aug2020_131_samples_library_sizes.txt", quote=F, row.names=F, sep="\t")
