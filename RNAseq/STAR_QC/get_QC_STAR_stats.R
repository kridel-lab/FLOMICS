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
source("read_STAR.R")
source("read_sample_files.R")
source("extract_samples.R")
#library(hciR)
library(stringr)

setwd("/cluster/projects/kridelgroup/FLOMICS/2022_combined_RNA-Seq_analysis/STAR_QC/TGL_OICR_log_final_out")
#setwd("/cluster/projects/kridelgroup/FLOMICS/2022_combined_RNA-Seq_analysis/STAR_QC/E4402_log_final_out")
#setwd("/cluster/projects/kridelgroup/RNAseq_cell_lines/220131_A00827_0499_AHTY3MDRXY_Kridel_Michael/STAR_out/log_final_out")
#setwd("/cluster/projects/kridelgroup/RNAseq_cell_lines/211123_A00827_0459_BHNWW2DRXY_Kridel_Robert/STAR_log_final")
#----------------------------------------------------------------------------------
#DATA
#----------------------------------------------------------------------------------

results = list.files(pattern="_T1Log.final.out$")

star_summary = as.data.table(read_STAR(path = ".", pattern="Log.final.out", reshape = FALSE))
star_summary$sample = sapply(star_summary$sample, function(x){ unlist(strsplit(x, "Log"))[1]})

#dcast
star_summary_dcast <- as.data.table(dcast(star_summary, sample ~ stat))

write.table(star_summary_dcast, file=paste(date, "STAR_QC_summary_TL.txt", sep="_"), quote=F, row.names=F, sep="	")

