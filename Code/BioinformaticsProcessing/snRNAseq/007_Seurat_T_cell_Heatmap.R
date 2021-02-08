#-------------------------------------------------------------------------------
#001_Seurat_tutorial.R
#karin isaev (using base code from RK)
#September 30th 2020
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors=F)

#make sure loading R/3.6.1 before running this script

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq")

#load functions to analyze seurat
source("/cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/doSeuratProc.R")

#output directory
output="/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/"

#load libraries
library(ellipsis)
library(tidyverse)
library(splitstackshape)
packages <- c("readr", "data.table", "plyr",
	"stringr",
  "Seurat",
  "cowplot",
	"patchwork", "Biobase")

lapply(packages, require, character.only = TRUE)

r=readRDS("combined_processed_seurat_object_rmFL277dim20.rds")

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

cells_t=c("Tfh cells_3", "CD8 T cells_4", "CD4 Treg cells_5", "naive T cells_7",
"memory T cells_8", "proliferating T cell_17")

DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)
combined <- ScaleData(combined, verbose = TRUE)

combined_t <- subset(combined, idents = cells_t)

dittoHeatmap(combined, genes,
    annot.by = c("ident", "sample"))

h = DoHeatmap(combined_t, features = genesall,assay="RNA", size=2)
pdf("rmFL277dim20_Tcells_only_markers.pdf")
print(h)
dev.off()
