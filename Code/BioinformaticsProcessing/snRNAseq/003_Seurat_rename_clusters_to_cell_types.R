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

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA")

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

date=Sys.Date()

#-------------------------------------------------------------------------------
#purpose
#-------------------------------------------------------------------------------

#the purpose is obtain clusters from snRNAseq that are representative
#of cell types which we can then use to convolute what cell types are present
#in the bulk RNA-seq samples

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

combined = readRDS(paste(output, "seurat_integrated_dim_20_2000_2021-02-05_samples_clusters.rds", sep=""))

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

#I think using imsig_cells seems like it makes the most sense

#based on review of gene markers across the clusters
#label cell types

current.cluster.ids = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
new.cluster.ids = c("B cells_0", "B cells_1", "B cells_2", "Tfh cells_3", "CD8 T cells_4",
 					"CD4 Treg cells_5", "Cluster 6", "naive T cells_7", "memory T cells_8",
					"naive B or malignant B_9", "macrophage or monocyte_10", "proliferating B cell_11",
				"memory B cell_12", "Cluster 13", "stromal cells_14", "endothelial cells_15",
			"Cluster 16", "proliferating T cell_17", "macrophage or monocyte_18",
		"Cluster 19")

names(x = new.cluster.ids) <- levels(x = combined)
combined <- RenameIdents(object = combined, new.cluster.ids)

pdf(paste(output, date, "_", "seurat_integrated_samples_with_cell_types.pdf", sep=""), width=10, height=8)
DimPlot(object = combined, reduction = "umap", label = TRUE, pt.size = 0.5)+
theme(axis.line = element_line(colour = 'black', size = 1), text = element_text(size = 12), axis.text = element_text(size = 12))
dev.off()

#4. save seurat object for bisque analysis
saveRDS(combined, file=paste(date, "combined_processed_snRNAseq_FL_seurat_object.rds", sep="_"))
