#-------------------------------------------------------------------------------
#RNA-seq-methylation-RESET.R
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors = F)
options(scipen=999) #avoid scientific notation

library(BisqueRNA)
library(Biobase)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "tidyr",
"data.table", "plyr", "gridExtra", "limSolve", "xCell",
"tibble", "immunedeconv",
"ggrepel", "stringr", "maftools", "ggpubr", "readxl", "skimr",
 "edgeR", "annotables", "EnvStats")

lapply(packages, require, character.only = TRUE)

#date
date=Sys.Date()

#getwd() --> FLOMICS teams folder
#cd /Users/kisaev/UHN/kridel-lab - Documents/FLOMICS

#-------------------------------------------------------------------------------
#data filtering
#-------------------------------------------------------------------------------

#gene annotations+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#UCSC gene classes - keep only protein coding genes
genes_class = as.data.table(grch37)
genes_class = as.data.table(filter(genes_class, biotype == "protein_coding"))
genes_class = as.data.table(filter(genes_class, !(is.na(entrez))))
genes_class = unique(genes_class[,c("ensgene", "symbol")])
#keep only one ens id per gene name
z = which(duplicated(genes_class$symbol))
genes_class = genes_class[-z,]

#sample info with rna-seq qc++++++++++++++++++++++++++++++++++++++++++++++++++++
rnaseq_qc = fread("metadata/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")

#methylation probes+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
meth = readRDS("methylation/2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.rds")
#only keep RLN and FL, remove dlbcl
z = which(str_detect(colnames(meth), "DLC"))
meth = meth[,-z] #155 FL + 5RN

#load in results from kallisto (gene expression)++++++++++++++++++++++++++++++++
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

# All FLOMICS samples included - load sample information+++++++++++++++++++++++++
all.samples.DNAseq.FLOMICS <- fread("metadata/sample_annotations_rcd6Nov2019.csv")


#-------------------------------------------------------------------------------
#Bisque analysis
#-------------------------------------------------------------------------------

#1. load bulk seq data into object
#Bulk RNA-seq data can be converted from a matrix
#(columns are samples, rows are genes) to an ExpressionSet as follows:
tpm_mat = as.matrix(tpm)
bulk.eset <- Biobase::ExpressionSet(assayData = tpm_mat)

#2. Single-cell data requires additional information in the ExpressionSet,
#specificially cell type labels and individual labels. Individual labels
#indicate which individual each cell originated from.

#To add this information, Biobase requires it to be stored in a data frame format.
#Assuming we have character vectors of cell type labels (cell.type.labels) and individual labels (individual.labels), we can convert
#scRNA-seq data (with counts also in matrix format) as follows:

#can automatically convert from seurat object
seurat.obj = readRDS("Analysis-Files/Seurat/combined_processed_snRNAseq_FL_seurat_object.rds")
sc.eset <- BisqueRNA::SeuratToExpressionSet(seurat.obj, delimiter="-", position=2, version="v3")
