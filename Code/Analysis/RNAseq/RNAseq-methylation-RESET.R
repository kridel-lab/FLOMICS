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
genes_class = as.data.table(filter(genes_class, biotype == "protein_coding"))
genes_class = as.data.table(filter(genes_class, !(is.na(entrez))))
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

#prepare input for resest

#1. first function is to prepare TSS probes that are consistently silenced or enhanced
#across normal samples

probes=promoter.probes.list
normal_matrix = meth[,c(which(colnames(meth) %in% c("LY_RLN_001","LY_RLN_002","LY_RLN_003",
"LY_RLN_004","LY_RLN_005")))]

methNorSel_res = methNorSel(normal_matrix, probes)#output is a list of two

#2. define methylation matrix for tumour samples
tum_matrix = meth[,c(which(!(colnames(meth) %in% c("LY_RLN_001","LY_RLN_002","LY_RLN_003",
"LY_RLN_004","LY_RLN_005"))))]

#3. run main analysis using the reset function to get silencing/enhancing scores

#do everything using just the reset function
# INPUTS:
  # 1. selected normal-sample datasets (the probe should have either sil or enh probes condition as specified in 'methNorSel' function)
    ## rows: probe index
    ## columns: normal samples
  # 2. tumor methylome datasets
    ## columns: tumor samples
    ## probeIDs
  # 3. Transcriptome dataset of the tumor samples
    ## columns: tumor samples
    ## Gene name (as specified in the probe index)
  # 4. no.permutation for FDR calculation
  # 5. Methylation event type: either "sil" or "enh"

#reset=function(normal.db, meth.tumor, transcriptome,
#    methylation.event=c('sil','enh'),
#    FDR.permutation.no=100, seed=100)

reset_res_sil = reset(methNorSel_res[[1]], tum_matrix, tpm, FDR.permutation.no=100, "sil")
reset_res_enh = reset(methNorSel_res[[2]], tum_matrix, tpm, FDR.permutation.no=100, "enh")

saveRDS(reset_res_sil, "Analysis-Files/RESET/silening_events_results.rds")
saveRDS(reset_res_enh, "Analysis-Files/RESET/enhancing_events_results.rds")
