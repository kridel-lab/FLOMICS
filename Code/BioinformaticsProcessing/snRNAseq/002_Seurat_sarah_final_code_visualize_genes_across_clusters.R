#-------------------------------------------------------------------------------
#002B_Seurat_tutorial.R
#sarah russell (using base code from RK)
#September 30th 2020
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors=F)
date=Sys.Date()

#make sure loading R/4.0.0 before running this script

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA")

#load functions to analyze seurat
source("/cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/doSeuratProc.R")

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

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE) #patient ID
input = args[1]
print(input) #name of seurat object that should be evaluated
analysis_type=unlist(strsplit(unlist(strsplit(input, "/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/"))[2], ".rds"))

#testing
#input=paste(output, "pc_genes_only_yes_seurat_integrated_SCnorm_dim_20_2000_2021-02-17_samples_clusters.rds", sep="")

combined = readRDS(input)

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

get_more_markers = function(dat, analysis_type){

	combined = dat
	DefaultAssay(combined) <- "RNA"

	#1. Find cluster markers++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# find markers for every cluster compared to all remaining cells, report only the positive ones
	combined.markers <- FindAllMarkers(combined, only.pos = TRUE,
		 min.pct = 0.25, logfc.threshold = 0.25, test.use="roc")
	combined.markers = as.data.table(combined.markers)

	#2. Visualize cluster markers+++++++++++++++++++++++++++++++++++++++++++++++++

	combined <- NormalizeData(combined)

	pdf(paste(date, "_", analysis_type, ".pdf", sep=""),width=18, height=16)

	genes=c("CD3D","CD3G","CD4","CD8A","CCR7","SELL","TCF7","IL7R","TYMS","MKI67","PRF1","CCL5")
	genes2=c("GZMK","GZMB","GZMA","PRDM1","CD69","TNF","IFNG","PTPRC","CCL4","KLRG1","TIGIT","CTLA4")
	genes3=c("PDCD1","HAVCR2","LAG3","TOX","TOX2","TBX21","EOMES","CD274","FOXP3","ENTPD1","CXCL13","TNFSF8")
	genes4=c("IL21","PTPN13","IL6R","CXCR5","BTLA","CD200","BCL2","ICOS","TOP2A","CR2")

	genesall=c("CD3D","CD3G","CD4","CD8A","CCR7","SELL","TCF7","IL7R","TYMS","MKI67","PRF1","CCL5",
	"GZMK","GZMB","GZMA","PRDM1","CD69","TNF","IFNG","PTPRC","CCL4","KLRG1","TIGIT","CTLA4",
	"PDCD1","HAVCR2","LAG3","TOX","TOX2","TBX21","EOMES","CD274","FOXP3","ENTPD1","CXCL13","TNFSF8",
	"IL21","PTPN13","IL6R","CXCR5","BTLA","CD200","BCL2","ICOS","TOP2A","CR2")

	#overlay on UMAP clusters
	f = FeaturePlot(combined, features = c("NCAM1", "MS4A1", "PTPRC", "BCL6", "CD68", "CD163", "CXCR5",
					    "CD69", "CD45RA", "CD45RO", "CCR7", "ICOS"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f)

	#overlay on UMAP clusters
	f2 = FeaturePlot(combined, features = c("CD79A", "CR2", "CD3D", "CCL5",
	"ICOS", "IGKC", "IGLC2", "IGLC3", "BCL2", "CD19"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f2)

	#overlay on UMAP clusters
	f3= FeaturePlot(combined, features = c("STAT3", "CD27", "CCR5", "CCR4", "CXCR4", "LAG3",
	"CD19", "IGHM", "PTPRG", "TNFRSF18", "TNFRSF4", "FAS"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f3)

	#overlay on UMAP clusters
	f4 = FeaturePlot(combined, features = c("CD3G", "CD8A", "BANK1", "VIM", "LYZ",
	"ST8SIA1", "ICA1", "IL7R", "TRAC", "CTLA4", "IKZF2"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f4)

	#B cell genes
	f5 = FeaturePlot(combined, features = c("EBF1", "BACH2", "LMO2", "FCER2", "CD27", "ACSM3",
  "KLHL6", "BCL2", "CD200", "MS4A1", "ARHGAP24", "IGHM", "BANK1", "PTPRG",
	"MKI67", "TOP2A", "BCL6", "TBC1D9", "LMO2", "IL21R"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f5)

	#T cell genes
	f6 = FeaturePlot(combined, features = c("CD3D", "CD3G", "CD4", "CD8A", "TCF7", "PTPRC", "TIGIT", "PDCD1", "TOX",
		"TOX2", "TNFSF8", "PTPN13", "ILR3", "BTLA", "CD200", "ICOS", "IL21", "CCL5", "GZMK", "GZMA",
		"PRDM1", "KLRG1", "TIGIT", "HAVCR2", "EOMES", "CTLA4", "TOX2", "FOXP3", "IL6R", "ICOS",
		"MKI67", "TCF7", "TOP2A"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f6)

	#other cluster genes
	f7 = FeaturePlot(combined, features = c("CXCL12", "SDF1", "VCAM1", "CR2", "CD21", "CD36",
  "TFPI", "LDB2", "CD3G", "HAVCR2", "CD274", "ILR6", "CD68", "CSF1R", "IFNGR1",
	"SPARCL1", "PECAM1", "VCAM1", "CD36", "CUX2", "SLC7A11", "CD4"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f7)

  combined <- ScaleData(combined, verbose = FALSE)
	h = DoHeatmap(combined, features = genesall, assay="RNA")
	print(h)

	top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_diff)
	h3 = DoHeatmap(combined, features = top10$gene)
	print(h3)

	dev.off()
	print("done")

}

get_more_markers(combined, analysis_type)
