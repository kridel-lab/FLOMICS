#----------------------------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors = F)
options(scipen=999)#avoid scientific notation

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "tidyr",
"data.table", "plyr", "gridExtra", "limSolve", "xCell",
"tibble", "immunedeconv",
"ggrepel", "stringr", "maftools", "ggpubr", "readxl", "skimr",
 "edgeR", "annotables", "EnvStats")

lapply(packages, require, character.only = TRUE)

#source reset scripts
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/methNorSel.R")
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/eventScore.R")
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/methStatus.R")
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/FDRcal.R")
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/reset.R")
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/shared_functions.R")

#date
date=Sys.Date()

#getwd() --> FLOMICS teams folder
#cd /Users/kisaev/UHN/kridel-lab - Documents/FLOMICS

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#reset results
reset_res_sil = readRDS("Analysis-Files/RESET/silening_events_results.rds")
reset_res_enh = readRDS("Analysis-Files/RESET/enhancing_events_results.rds")

#get scores
sil_scores = as.data.table(reset_res_sil[[6]])
enh_scores = as.data.table(reset_res_enh[[6]])

sil_scores_fdr = as.data.table(reset_res_sil[[11]])
colnames(sil_scores_fdr)[1] = "Score"
sil_scores = merge(sil_scores, sil_scores_fdr, by="Score")

enh_scores_fdr = as.data.table(reset_res_enh[[11]])
colnames(enh_scores_fdr)[1] = "Score"
enh_scores = merge(enh_scores, enh_scores_fdr, by="Score")

sil_scores = sil_scores[order(-Score)]
enh_scores = enh_scores[order(-Score)]
