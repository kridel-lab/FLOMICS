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

#-------------------------------------------------------------------------------
#purpose
#-------------------------------------------------------------------------------

#the purpose is obtain clusters from snRNAseq that are representative
#of cell types which we can then use to convolute what cell types are present
#in the bulk RNA-seq samples

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

combined = readRDS()

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

#I think using imsig_cells seems like it makes the most sense

#based on review of gene markers across the clusters
#label cell types

current.cluster.ids = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
new.cluster.ids = c("Malignant B cells", "Malignant B cells", "CD4 reg or Tfh cells_2",
					"CD8 cytotoxic T cells cells_3",
                    "Malignant B cells", "CD4 FTH cells_5",
										"stromal cells_6", "CD4 T cells_7", "CD4 reg or Tfh cells_2",
									"Macrophages_9", "Predicted Norm B cells_10", "potential_Monocytes_11", "Endothelial cells_12")
names(x = new.cluster.ids) <- levels(x = combined)
combined <- RenameIdents(object = combined, new.cluster.ids)

pdf(paste(output, "seurat_integrated_samples_cell_types_prelim.pdf", sep=""), width=18, height=12)
DimPlot(object = combined, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

#4. save seurat object for bisque analysis
saveRDS(combined, file="combined_processed_snRNAseq_FL_seurat_object.rds")
