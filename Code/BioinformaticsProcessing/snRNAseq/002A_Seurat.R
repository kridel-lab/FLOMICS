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

all_perms = list.files(output, pattern=".rds")
setwd(output)

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

get_marker_genes = function(dat){
	print(dat)
	combined = readRDS(dat)

	#1. Find cluster markers++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# find markers for every cluster compared to all remaining cells, report only the positive ones
	DefaultAssay(combined) <- "integrated"

	combined.markers <- FindAllMarkers(combined, only.pos = TRUE,
		 min.pct = 0.25, logfc.threshold = 0.25, test.use="roc")
	combined.markers = as.data.table(combined.markers)

	output_file = paste(unlist(strsplit(dat, ".rds")), "_markers.csv", sep="")

	write.csv(combined.markers, file=output_file,
	quote=F, row.names=F)

	#2. Visualize cluster markers+++++++++++++++++++++++++++++++++++++++++++++++++

	pdf(paste(unlist(strsplit(dat, ".rds")), "_some_markers.pdf", sep=""), width=18, height=12)

	genes=c("CD3D", "CD4", "CD8A", "NCAM1", "MS4A1", "PTPRC", "BCL6",
	"FOXP3", "PDCD1", "CD68", "CD163", "CXCR5", "CD69", "CD45RA", "CD45RO", "CCR7", "IL7R")

	v = VlnPlot(combined, features = genes)
	print(v)

	#overlay on UMAP clusters
	f = FeaturePlot(combined, features = genes,
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f)

	#v2 = VlnPlot(combined, features = c("CD79A", "CR2", "CD3D", "CCL5", "ICOS", "IGKC",
	#"IGLC2", "IGLC3", "CD19"))
	#print(v2)

	#overlay on UMAP clusters
	f2 = FeaturePlot(combined, features = c("CD79A", "CR2", "CD3D", "CCL5",
	"ICOS", "IGKC", "IGLC2", "IGLC3", "BCL2", "CD19"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f2)

	top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_diff)
	h = DoHeatmap(combined, features = top10$gene)
	print(h)

	#overlay on UMAP clusters
	f3= FeaturePlot(combined, features = c("STAT3", "CD27", "CCR5", "CCR4", "CXCR4", "LAG3",
	"CD19", "IGHM", "PTPRG", "TNFRSF18", "TNFRSF4", "FAS"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f3)

	f4 = FeaturePlot(combined, features = c("CD3G", "CD8A", "BANK1", "VIM", "LYZ",
	"ST8SIA1", "ICA1", "IL7R", "TRAC", "CTLA4", "IKZF2"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f4)

	dev.off()

	print("done")

}

get_marker_genes("seurat_integrated_dim_10_2000_samples_clusters.rds")
#llply(all_perms, get_marker_genes, .progress=".text")
