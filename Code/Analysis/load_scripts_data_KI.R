#----------------------------------------------------------------------
#Karin Isaev
#load libraries required for analysis
#load data stored locally for visualizations and some analysis
#that can be done locally without memory or time issues
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors = F)
#avoid scientific notation
options(scipen=999)
#load libraries
packages <- c("dplyr", "readr", "ggplot2", "tidyr",
"data.table", "plyr", "gridExtra",
"ggrepel", "stringr", "maftools", "ggpubr", "readxl", "skimr",
 "edgeR", "annotables", "EnvStats", "limSolve",
"xCell", "tibble", "immunedeconv" , "BisqueRNA", "Biobase", "Seurat")

lapply(packages, require, character.only = TRUE)

#load cibersort scripts
set_cibersort_binary("Analysis-Files/Immune-Deconvolution/CIBERSORT.R")
set_cibersort_mat("Analysis-Files/Immune-Deconvolution/LM22.txt")

#load RESET scripts
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/methNorSel.R")
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/eventScore.R")
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/methStatus.R")
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/FDRcal.R")
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/reset.R")
source("/Users/kisaev/github/FLOMICS/Code/Tools/RESET/R/shared_functions.R")

#date
date=Sys.Date()

#if do getwd() --> FLOMICS teams folder (starting point for all my local scripts)
#cd /Users/kisaev/UHN/kridel-lab - Documents/FLOMICS

#----------------------------------------------------------------------
#data filtering
#----------------------------------------------------------------------

#1. gene annotations
#UCSC gene classes - only protein coding genes
genes_class = as.data.table(grch37)
#genes_class = as.data.table(filter(genes_class, biotype == "protein_coding"))
#genes_class = as.data.table(filter(genes_class, !(is.na(entrez))))
genes_class = unique(genes_class[,c("ensgene", "symbol", "biotype")])
#keep only one ens id per gene name
z = which(duplicated(genes_class$symbol))
genes_class = genes_class[-z,]

#2. all FLOMICS samples included - load sample information++++++++++++++++++++++
all.samples.DNAseq.FLOMICS <- fread("metadata/sample_annotations_rcd6Nov2019.csv")

#3. sample info with rna-seq qc+++++++++++++++++++++++++++++++++++++++++++++++++
rnaseq_qc = fread("metadata/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")

#4. load cluster labels from SNF
snf = fread("/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/Cluster Labels/InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv")

#5. load in results from STAR containing raw counts+++++++++++++++++++++++++++++
exp <- fread("RNAseq/counts/2020-06-18STAR_quantmode_counts_matrix_FL_136_patients.txt")

print("ready for analysis!")
