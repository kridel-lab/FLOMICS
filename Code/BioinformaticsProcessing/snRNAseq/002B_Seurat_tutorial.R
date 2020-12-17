#-------------------------------------------------------------------------------
#002B_Seurat_tutorial.R
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
#data
#-------------------------------------------------------------------------------
all_perms = list.files(output, pattern=".rds")
setwd(output)


#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

get_more_markers = function(dat){
	print(dat)
	combined = readRDS(dat)

	#1. Find cluster markers++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# find markers for every cluster compared to all remaining cells, report only the positive ones
	DefaultAssay(combined) <- "RNA"
	combined <- NormalizeData(combined)

#	combined.markers <- FindAllMarkers(combined, only.pos = TRUE,
#		 min.pct = 0.25, logfc.threshold = 0.25, test.use="roc")
#	combined.markers = as.data.table(combined.markers)

#	output_file = paste(unlist(strsplit(dat, ".rds")), "_markers.csv", sep="")

#	write.csv(combined.markers, file=output_file,
#	quote=F, row.names=F)

	#2. Visualize cluster markers+++++++++++++++++++++++++++++++++++++++++++++++++

	pdf(paste(date, "dim30_samples_gene_markers.pdf", sep="_"), width=18, height=12)

	genes=c("CD3D","CD3G","CD4","CD8A","CCR7","SELL","TCF7","IL7R","TYMS","MKI67","PRF1","CCL5")
	v = VlnPlot(combined, features = genes, pt.size = 0.0625)
	print(v)

	genes2=c("GZMK","GZMB","GZMA","PRDM1","CD69","TNF","IFNG","PTPRC","CCL4","KLRG1","TIGIT","CTLA4")
	v2 = VlnPlot(combined, features = genes2, pt.size = 0.0625)
	print(v2)

	genes3=c("PDCD1","HAVCR2","LAG3","TOX","TOX2","TBX21","EOMES","CD274","FOXP3","ENTPD1","CXCL13","TNFSF8")
	v3 = VlnPlot(combined, features = genes3, pt.size = 0.0625)
	print(v3)

	genes4=c("IL21","PTPN13","IL6R","CXCR5","BTLA","CD200","BCL2","ICOS","TOP2A")
	v4 = VlnPlot(combined, features = genes4, pt.size = 0.0625)
	print(v4)

	genesall=c("CD3D","CD3G","CD4","CD8A","CCR7","SELL","TCF7","IL7R","TYMS","MKI67","PRF1","CCL5",
	"GZMK","GZMB","GZMA","PRDM1","CD69","TNF","IFNG","PTPRC","CCL4","KLRG1","TIGIT","CTLA4",
	"PDCD1","HAVCR2","LAG3","TOX","TOX2","TBX21","EOMES","CD274","FOXP3","ENTPD1","CXCL13","TNFSF8",
	"IL21","PTPN13","IL6R","CXCR5","BTLA","CD200","BCL2","ICOS","TOP2A")
	
  combined <- ScaleData(combined, verbose = FALSE)
	h = DoHeatmap(combined, features = genesall,assay="RNA")
	print(h)

  #cell cycle markers
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes

	h2 = DoHeatmap(combined, features = s.genes,assay="RNA")
	print(h2)

	h3 = DoHeatmap(combined, features = g2m.genes,assay="RNA")
	print(h3)
	

	dev.off()
	print("done")

}

get_more_markers("seurat_integrated_dim_30_2000_samples_clusters.rds")


#llply(all_perms, get_marker_genes, .progress=".text")
