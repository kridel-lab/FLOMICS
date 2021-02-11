#-------------------------------------------------------------------------------
#001_Seurat_tutorial.R
#karin isaev (using base code from RK)
#September 30th 2020
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors=F)

#make sure loading R/4.0.0 before running this script

library(dittoSeq)
library(scater)
library(loomR)
library(Seurat)
library(scRNAseq)
library(SingleCellExperiment)

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq")

#load functions to analyze seurat
source("/cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/doSeuratProc.R")
date=Sys.Date()

#output directory
output="/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/"

#load libraries
library(ellipsis)
library(tidyverse)
#library(splitstackshape)
packages <- c("readr", "data.table", "plyr",
	"stringr",
  "Seurat",
  "cowplot",
	"patchwork", "Biobase")

lapply(packages, require, character.only = TRUE)

#r=readRDS("combined_processed_seurat_object_rmFL277dim20.rds")
r=readRDS("/cluster/projects/kridelgroup/FLOMICS/DATA/2021-02-05_combined_processed_snRNAseq_FL_seurat_object.rds")

#-------------------------------------------------------------------------------
#purpose
#-------------------------------------------------------------------------------

#the purpose is obtain clusters from snRNAseq that are representative
#of cell types which we can then use to convolute what cell types are present
#in the bulk RNA-seq samples

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

combined = r

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------
genesall = c("CD79A", "BCL2", "CD19", "PTPRG", "BANK1",
"EBF1", "BACH2", "LMO2", "FCER2", "CD27", "ACSM3",
"BACH2", "KLHL6", "EBF1", "FCER2", "LMO2",
"CD200","MS4A1", "ARHGAP24", "IGHM",
"MKI67", "TOP2A", "BCL6", "TBC1D9",
"FCER2", "LMO2", "IL21R")

genesall = c("KIF18B", "MIR155HG", "TOP2A", "MKI67", "ARHGAP24",
"TBC1D9", "HMMR", "CD83", "TPX2", "MT-CO3", "UBE2C",
"MT-ND3", "IGSF3", "VAV3-AS1", "EPHB1", "SAMSN1",
"DEPDC1","PTPRG","SGOL2","TCF4","ASPM","POLQ","CD80","MT-ND4",
"KIF14","FANCA","MT-ATP6","BANK1","SMC4","CENPE","MT-ND1",
"RP11-789C2.1","ANKRD33B","NDC80","BMP7",
"NUSAP1","JADE3","PPM1E","MT-ND2","MT-CYB","PARVB",
"CENPF","AC023590.1","ECT2","SSTR2","CSGALNACT1",
"FAP","NR4A3","FCER2","MT-CO2","CD58","CADPS","EPDR1",
"DOCK10","KIAA0125","SLAMF1","BUB1B","CCL17","SGPP2","ARHGAP11B","IGLC2",
"GPM6A", "MIR4432HG", "ADAMTS6", "ARHGAP24", "SLC7A7", "KHDRBS2","APOLD1",
"MGLL","MIR181A2HG","EZH2","RFC3","CDCA2","RP11-669M16.1","C12orf42",
"HELLS","PECAM1","CTD-2034I4.1","ST6GALNAC3","ZNF804A","VAV3","BATF",
"GALNT14","CCL17","WDFY4","MEF2C","CD80","BCL2","KIAA2022","ARHGAP31","KIF11",
"LYPD6B","BUB1","LINC01572","NR6A1","SLCO5A1","IL21R","BACH2","C1orf186","FCRL5",
"HIP1","RNLS","TP63","KIF13A","NUF2","KLHL6","C2orf48","KYNU","RP11-624C23.1",
"LMO2","LYPD6B","UBE2E2","ADRBK2","ABTB2","RP11-475O6.1","CLSPN","DGKG","MARCKS",
"CCL22")
genesall = unique(genesall)

cells_t=c("B cells_0", "B cells_1", "B cells_2", "naive B or malignant B_9",
"proliferating B cell_11", "memory B cell_12", "Cluster 13")

DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)
combined <- ScaleData(combined, verbose = TRUE)

combined_t <- subset(combined, idents = cells_t)
subset.matrix <- combined_t[genesall, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest

# Annotating and ordering cells by some meaningful feature(s):
pdf(paste(date, "_" , "rmFL277dim20_Bcells_only_markers_dittoheatmap_AUC0.75.pdf", sep=""), height=20)
dittoHeatmap(subset.matrix, annot.by = c("seurat_clusters"),scaled.to.max = TRUE)
dev.off()

h = DoHeatmap(combined_t, features = genesall, assay="RNA", size=2)
pdf(paste(date, "_" , "rmFL277dim20_Bcells_only_markers_AUC0.75.pdf", sep=""), height=20)
print(h)
dev.off()
