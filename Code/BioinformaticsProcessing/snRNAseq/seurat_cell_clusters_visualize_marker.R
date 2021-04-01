#-------------------------------------------------------------------------------
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

setwd("~/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April2021")

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

r=readRDS("pc_genes_only_no_seurat_integrated_dim_20_2000_2021-04-01_samples_clusters.rds")

#-------------------------------------------------------------------------------
#purpose
#-------------------------------------------------------------------------------

#visualize differences between B cell clusters

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

combined = r

#confirm clusters
DimPlot(combined, reduction = "umap", label = TRUE)

DefaultAssay(combined) <- "RNA"

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

FeaturePlot(combined, features = c("CCSER1", "KHDRBS2"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1)

FeaturePlot(combined, features = c("BCL2", "MS4A1"), order=TRUE, label=TRUE, blend.threshold=0.1, cols=c("grey", "thistle1", "steelblue", "red"))

# For performing differential expression after integration, we switch back to the original data
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
head(nk.markers)