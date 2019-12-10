setwd("/Users/kisaev/FLOMICS_Anjali/FLOMICS_Anajli/Data/qc")

library(data.table)
library(plyr)
library(dplyr)

files=list.files(pattern=".captureKit.unique.tsv")
get_cov = function(file){
  print(file)
  f = fread(file)
  f=melt(f)
  f$pats = ""
  f$pats = sapply(as.character(f$variable), function(x){unlist(strsplit(x, "_"))[3]})
  f$type_measure = ""
  f$type_measure = sapply(as.character(f$variable), function(x){paste(unlist(strsplit(x, "_"))[1:2], collapse="_")})
  f$variable = NULL
  return(f)
}

all_covs = as.data.table(ldply(llply(files, get_cov)))
