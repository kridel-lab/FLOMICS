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
#data filtering
#----------------------------------------------------------------------

#gene annotations
#UCSC gene classes - only protein coding genes
genes_class = as.data.table(grch37)
#genes_class = as.data.table(filter(genes_class, biotype == "protein_coding"))
#genes_class = as.data.table(filter(genes_class, !(is.na(entrez))))
genes_class = unique(genes_class[,c("ensgene", "symbol")])
#keep only one ens id per gene name
z = which(duplicated(genes_class$symbol))
genes_class = genes_class[-z,]

#methylation probes
meth = readRDS("methylation/2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.rds")

#load in results from kallisto
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
tpm = tpm[,z]

# All FLOMICS samples included - load sample information
all.samples.DNAseq.FLOMICS <- fread("metadata/sample_annotations_rcd6Nov2019.csv")

#sample info with rna-seq qc
rnaseq_qc = fread("metadata/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")

#load the probe list
load("Analysis-Files/RESET/promoter-probes-list.rdata")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

probes=promoter.probes.list
normal_matrix = meth[,c(which(colnames(meth) %in% c("LY_RLN_001","LY_RLN_002","LY_RLN_003",
"LY_RLN_004","LY_RLN_005")))]
tum_matrix = meth[,c(which(!(colnames(meth) %in% c("LY_RLN_001","LY_RLN_002","LY_RLN_003",
"LY_RLN_004","LY_RLN_005"))))]

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
