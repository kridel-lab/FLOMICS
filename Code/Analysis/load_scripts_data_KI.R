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

#4. load in results from STAR containing raw counts+++++++++++++++++++++++++++++
exp <- fread("RNAseq/counts/2020-06-18STAR_quantmode_counts_matrix_FL_136_patients.txt")

#5. load in results from kallisto+++++++++++++++++++++++++++++++++++++++++++++++
tpm = fread("RNAseq/counts/2020-09-01_kallisto_gene_based_counts.txt", data.table=F)
colnames(tpm)[1] = "ensgene"
tpm = merge(tpm, genes_class, by = "ensgene")
tpm = as.data.frame(tpm)
rownames(tpm) = tpm$symbol
tpm$symbol = NULL
tpm$ensgene = NULL
#remove T2 samples from expression matrix
z = which(str_detect(colnames(tpm), "T2"))
tpm = tpm[,-z]
#add _T1 to sample names that don't have it because that's how the samples
#are listed in the methylation file
z = which(!(str_detect(colnames(tpm), "T1")) & !(str_detect(colnames(tpm), "DLC")) & !(str_detect(colnames(tpm), "RLN")))
colnames(tpm)[z] = paste(colnames(tpm)[z], "_T1", sep="")
#only include FL tumours in TPM matrix for analysis
z = which(str_detect(colnames(tpm), "LY_FL"))
tpm = tpm[,z] #122 FL samples
#exclude one patient that Anjali has filtered out
z = which(colnames(tpm) %in% rnaseq_qc$SAMPLE_ID)
tpm = tpm[,z] #121 FL samples
dim(tpm)
#[1] 17378   121

print("ready for analysis!")
