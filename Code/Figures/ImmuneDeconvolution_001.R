#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors=F)

#make sure loading R/4.0.0 before running this script

#library(dittoSeq)
#library(scater)
#library(loomR)
library(Seurat)
#library(scRNAseq)
#library(SingleCellExperiment)
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(ReactomeGSA)
library(enrichR)
#load libraries
library(ellipsis)
library(tidyverse)
library(dittoSeq)

#library(splitstackshape)
packages <- c("readr", "data.table", "plyr",
	"stringr",
  "Seurat",
  "cowplot",
	"patchwork", "Biobase")

lapply(packages, require, character.only = TRUE)

setwd("~/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April82021")

date=Sys.Date()

r=readRDS("pc_genes_only_no_seurat_integrated_dim_20_2000_2021-04-08_samples_clusters.rds")

mainDir="/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April82021"
subDir=paste("Figures_for_manuscript", date, sep="_")

dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir, subDir))

cols = c(dittoColors(), dittoColors(1)[seq_len(7)])

#-------------------------------------------------------------------------------
#purpose
#-------------------------------------------------------------------------------

#make summary figures for UMAP

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

combined = r
cells_b = c(0, 1, 2, 4, 11, 12)

DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)

#only keep B cell clusters
combined_b <- subset(combined, idents = cells_b)

#plot
p1 <- DimPlot(combined, reduction = "umap", group.by = "sample")+
  theme(axis.line = element_line(colour = 'black', size = 1),
        text = element_text(size = 20), axis.text = element_text(size = 20))

p2 <- DimPlot(combined, reduction = "umap", label = TRUE, label.size=5, cols = cols)+
  theme(axis.line = element_line(colour = 'black', size = 1),
        text = element_text(size = 20), axis.text = element_text(size = 20))

#save 
pdf("Figure6A.pdf", width=5, height=5)
print(p1)
dev.off()

pdf("Figure6B.pdf", width=5, height=5)
print(p2)
dev.off()
