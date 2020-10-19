#----------------------------------------------------------------------
#Karin Isaev
#Use RNA-seq count matrix (first normalize to TPM)
#Run through immune deconvolution tools to estimate fraction of immune
#cells present in each sample
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors = F)
#avoid scientific notation
options(scipen=999)
#load libraries
packages <- c("dplyr", "readr", "ggplot2", "tidyr",
"data.table", "plyr",
"ggrepel", "stringr", "maftools", "ggpubr", "readxl", "skimr",
 "edgeR", "annotables", "EnvStats")
library(gridExtra)
lapply(packages, require, character.only = TRUE)
library(limSolve)
library(xCell)
library(tibble)
library(immunedeconv) #<- main package with tools for immune deconvolution
set_cibersort_binary("Analysis-Files/Immune-Deconvolution/CIBERSORT.R")
set_cibersort_mat("Analysis-Files/Immune-Deconvolution/LM22.txt")

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

#load in results from kallisto
tpm = fread("RNAseq/counts/2020-09-01_kallisto_gene_based_counts.txt", data.table=F)
colnames(tpm)[1] = "ensgene"
tpm = merge(tpm, genes_class, by = "ensgene")
tpm = as.data.frame(tpm)
rownames(tpm) = tpm$symbol
tpm$symbol = NULL
tpm$ensgene = NULL

# All FLOMICS samples included - load sample information
all.samples.DNAseq.FLOMICS <- fread("metadata/sample_annotations_rcd6Nov2019.csv")

#sample info with rna-seq qc
rnaseq_qc = fread("metadata/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")

#load in results from STAR containing raw counts++++++++++++++++++++++++++++++++
exp <- fread("RNAseq/counts/2020-06-18STAR_quantmode_counts_matrix_FL_136_patients.txt")

# All FLOMICS samples included - load sample information++++++++++++++++++++++++
all.samples.DNAseq.FLOMICS <- fread("metadata/sample_annotations_rcd6Nov2019.csv")

bisque = readRDS("Analysis-Files/Seurat/bisque_decomposed_samples.rds")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#visualize results from bisque
immune_cells = as.data.frame(bisque)
immune_cells$cell_type = rownames(immune_cells)

#add tag based on whether sample is limited advanced and FL vs DLBCL
immune_cells = as.data.table(immune_cells)
immune_cells = melt((immune_cells))
colnames(immune_cells)[2] = "rna_seq_file_sample_ID"
immune_cells = merge(immune_cells, rnaseq_qc, by = "rna_seq_file_sample_ID")
immune_cells$Cluster = factor(immune_cells$Cluster)

head(immune_cells)

#plot distribution of each cell type frequency across disease and stages
#set up file for plotting
pdf(paste("Analysis-Files/Immune-Deconvolution/", date, "BISQUE_results.pdf", sep=""), width=15)
g1 = ggboxplot(filter(immune_cells, !(is.na(STAGE))), x="cell_type", y="value", fill="STAGE") +
theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
stat_compare_means(aes(group = STAGE), label = "p.signif")

g2 = ggboxplot(immune_cells, x="cell_type", y="value", fill="TYPE") +
theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
stat_compare_means(aes(group = TYPE), label = "p.signif")

g3 = ggboxplot(immune_cells, x="cell_type", y="value", fill="Cluster",palette = "jco") +
theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
stat_compare_means(aes(group = Cluster), label = "p.signif")

print(g1)
print(g2)
print(g3)

dev.off()
