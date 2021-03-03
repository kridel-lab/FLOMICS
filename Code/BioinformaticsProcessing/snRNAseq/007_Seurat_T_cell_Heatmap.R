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
output="/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/Feb2020/"

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
#r=readRDS("/cluster/projects/kridelgroup/FLOMICS/DATA/2021-02-05_combined_processed_snRNAseq_FL_seurat_object.rds")
r=readRDS(paste(output, "pc_genes_only_no_seurat_integrated_dim_20_2000_2021-02-23_samples_clusters.rds", sep=""))

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
genesall = c("CD3D", "CD3G", "CD4", "CD8A", "TCF7", "PTPRC", "TIGIT", "PDCD1", "TOX",
	"TOX2", "TNFSF8", "PTPN13", "ILR3", "BTLA", "CD200", "ICOS", "IL21", "CCL5", "GZMK", "GZMA",
	"PRDM1", "KLRG1", "TIGIT", "HAVCR2", "EOMES", "CTLA4", "TOX2", "FOXP3", "IL6R", "ICOS",
	"MKI67", "TCF7", "TOP2A")

cells_t=c(4, 5, 6, 8, 10, 16)

###
# T cells
###
#*P4 - CTLA4 - Treg
#P5
#P6
#P8
#P10
#*P16 - proliferating T cells


DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)
combined <- ScaleData(combined, verbose = TRUE)

combined_t <- subset(combined, idents = cells_t)
subset.matrix <- combined_t[genesall, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest

# Annotating and ordering cells by some meaningful feature(s):
pdf(paste(output, date, "_" , "T_cell_genes_dittoheatmap.pdf", sep=""), height=10, width=10)
dittoHeatmap(subset.matrix, annot.by = c("seurat_clusters"),scaled.to.max = TRUE)
dev.off()

h1 = DoHeatmap(combined_t, features = genesall, assay="RNA", size=2)
h2 = DoHeatmap(combined_t, features = genesall, assay="integrated", size=2)

pdf(paste(output, date, "_" , "T_cell_genes_heatmap_seurat.pdf", sep=""), height=10, width=10)
print(h1)
print(h2)
dev.off()
