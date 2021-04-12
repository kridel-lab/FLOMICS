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
#library(loomR)
library(Seurat)
library(scRNAseq)
library(SingleCellExperiment)
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(ReactomeGSA)

setwd("~/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April82021")

date=Sys.Date()

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

r=readRDS("pc_genes_only_no_seurat_integrated_dim_20_2000_2021-04-08_samples_clusters.rds")

mainDir="/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April82021"
subDir=paste("Figures_for_manuscript", date, sep="_")

dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir, subDir))

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
genesall = c("CD3D", "CD3G", "IL7R", "CD4", "CD8A", "TCF7", "PTPRC", "TIGIT", "PDCD1", "TOX",
             "TOX2", "TNFSF8", "PTPN13", "ILR3", "BTLA", "CD200", "ICOS", "IL21", "CCL5", "GZMK", "GZMA",
             "PRDM1", "KLRG1", "TIGIT", "HAVCR2", "EOMES", "CTLA4", "TOX2", "FOXP3", "IL6R", 
             "MKI67", "TCF7", "TOP2A", "IL2RA", "PLAC8", "KLF2", "CCL5", "GZMA", "NKG7", "CCL4", "CXCR5")

#FeaturePlot(subset.matrix, features = c("CD8A"), order=TRUE, label=TRUE)

cells_t=c(5, 6, 7, 8, 9, 15)

###
# T cells
###

DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)

#combined <- ScaleData(combined, verbose = TRUE)

combined_t <- subset(combined, idents = cells_t)
subset.matrix <- combined_t[genesall, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest

pdf("Figure6D.pdf", width=7, height=4)
dittoHeatmap(subset.matrix, annot.by = c("seurat_clusters"), scaled.to.max = TRUE)
dev.off()

pdf("T-cell_clusters_features.pdf", width=12, height=8)

dittoHeatmap(subset.matrix, annot.by = c("seurat_clusters"), scaled.to.max = TRUE)

FeaturePlot(subset.matrix, features = c("GZMK", "CCL5"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1)
FeaturePlot(subset.matrix, features = c("IL7R", "PLAC8"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1)
FeaturePlot(subset.matrix, features = c("PDCD1", "TOX"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1)
FeaturePlot(subset.matrix, features = c("IL2RA", "FOXP3"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1)
FeaturePlot(subset.matrix, features = c("TOP2A", "MKI67"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1)

dev.off()




