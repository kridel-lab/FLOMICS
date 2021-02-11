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

combined = readRDS(paste(output, "seurat_integrated_dim_20_2000_2021-02-05_samples_clusters.rds", sep=""))

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

get_more_markers = function(dat){
	combined = dat

	#1. Find cluster markers++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# find markers for every cluster compared to all remaining cells, report only the positive ones
	combined.markers <- FindAllMarkers(combined, only.pos = TRUE,
		 min.pct = 0.25, logfc.threshold = 0.25, test.use="roc")
	combined.markers = as.data.table(combined.markers)

	#output_file = paste(unlist(strsplit(dat, ".rds")), "_markers.csv", sep="")

	#write.csv(combined.markers, file=output_file,
	#quote=F, row.names=F)

	#2. Visualize cluster markers+++++++++++++++++++++++++++++++++++++++++++++++++

	DefaultAssay(combined) <- "RNA"
	combined <- NormalizeData(combined)

#	pdf(paste(date, unlist(strsplit(dat, ".rds")), "_gene_markers.pdf", sep="_"), width=18, height=12)

	pdf(paste(date, "dim20_no227_gene_markers.pdf", sep="_"),width=18, height=12)
	genes=c("CD3D","CD3G","CD4","CD8A","CCR7","SELL","TCF7","IL7R","TYMS","MKI67","PRF1","CCL5")
	#v = VlnPlot(combined, features = genes, pt.size = 0.0625)
	#print(v)

	genes2=c("GZMK","GZMB","GZMA","PRDM1","CD69","TNF","IFNG","PTPRC","CCL4","KLRG1","TIGIT","CTLA4")
	#v2 = VlnPlot(combined, features = genes2, pt.size = 0.0625)
	#print(v2)

	genes3=c("PDCD1","HAVCR2","LAG3","TOX","TOX2","TBX21","EOMES","CD274","FOXP3","ENTPD1","CXCL13","TNFSF8")
	#v3 = VlnPlot(combined, features = genes3, pt.size = 0.0625)
	#print(v3)

	genes4=c("IL21","PTPN13","IL6R","CXCR5","BTLA","CD200","BCL2","ICOS","TOP2A","CR2")
	#v4 = VlnPlot(combined, features = genes4, pt.size = 0.0625)
	#print(v4)

	genesall=c("CD3D","CD3G","CD4","CD8A","CCR7","SELL","TCF7","IL7R","TYMS","MKI67","PRF1","CCL5",
	"GZMK","GZMB","GZMA","PRDM1","CD69","TNF","IFNG","PTPRC","CCL4","KLRG1","TIGIT","CTLA4",
	"PDCD1","HAVCR2","LAG3","TOX","TOX2","TBX21","EOMES","CD274","FOXP3","ENTPD1","CXCL13","TNFSF8",
	"IL21","PTPN13","IL6R","CXCR5","BTLA","CD200","BCL2","ICOS","TOP2A","CR2")

	v5 = VlnPlot(combined, features = c("NCAM1", "MS4A1", "PTPRC", "BCL6", "CD68", "CD163", "CXCR5",
					    "CD69", "CD45RA", "CD45RO", "CCR7"))
	#print(v5)

	#overlay on UMAP clusters
	f = FeaturePlot(combined, features = c("NCAM1", "MS4A1", "PTPRC", "BCL6", "CD68", "CD163", "CXCR5",
					    "CD69", "CD45RA", "CD45RO", "CCR7", "ICOS"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f)

	v6 = VlnPlot(combined, features = c("CD79A", "CCL5", "IGKC", "IGLC2", "IGLC3", "CD19"))
	print(v6)

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

	f4 = FeaturePlot(combined, features = c("CD3G", "CD8A", "BANK1", "VIM", "LYZ",
	"ST8SIA1", "ICA1", "IL7R", "TRAC", "CTLA4", "IKZF2"),
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f4)

  combined <- ScaleData(combined, verbose = FALSE)
	h = DoHeatmap(combined, features = genesall, assay="RNA")
	print(h)

  #cell cycle markers
	#s.genes <- cc.genes$s.genes
	#g2m.genes <- cc.genes$g2m.genes

	#h2 = DoHeatmap(combined, features = s.genes,assay="RNA")
	#print(h2)

	#h3 = DoHeatmap(combined, features = g2m.genes,assay="RNA")
	#print(h3)

	top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_diff)
	h3 = DoHeatmap(combined, features = top10$gene)
	print(h3)

	dev.off()
	print("done")

}

get_more_markers(combined)
