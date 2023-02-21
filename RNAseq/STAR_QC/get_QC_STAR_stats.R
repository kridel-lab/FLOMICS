#----------------------------------------------------------------------------------
#get QC from STAR
#----------------------------------------------------------------------------------

#Karin Isaev created 
#Ting Liu updated

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
source("read_STAR.R")
source("read_sample_files.R")
source("extract_samples.R")
#library(hciR)
library(stringr)

setwd("/pwd/STAR_QC/log_final_out")
#----------------------------------------------------------------------------------
#DATA
#----------------------------------------------------------------------------------

results = list.files(pattern="_Log.final.out$")

star_summary = as.data.table(read_STAR(path = ".", pattern="Log.final.out", reshape = FALSE))
star_summary$sample = sapply(star_summary$sample, function(x){ unlist(strsplit(x, "Log"))[1]})

#dcast
star_summary_dcast <- as.data.table(dcast(star_summary, sample ~ stat))

write.table(star_summary_dcast, file=paste(date, "STAR_QC_summary_TL.txt", sep="_"), quote=F, row.names=F, sep="	")

