#-------------------------------------------------------------------------------
#002B_Seurat_tutorial.R
#sarah russell
#December 2 2020
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors=F)

#R 4.0

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA")

#load functions to analyze seurat
source("/cluster/home/srussell/github/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/doSeuratProc.R")

#output directory
output="/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/sarah-seurat/final/"

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

#the purpose is to perform differential expression analysis between select
#seurat clusters (helping to differentiate cluster cell types)

#-------------------------------------------------------------------------------
#notes
#-------------------------------------------------------------------------------
#https://satijalab.org/seurat/v3.1/de_vignette.html
#bulk of Seurat’s differential expression features can be accessed through the FindMarkers function.
#As a default, Seurat performs differential expression based on the non-parameteric Wilcoxon rank sum test.
#To test for differential expression between two specific groups of cells, specify the ident.1 and ident.2 parameters.

##Prefilter features or cells to increase the speed of DE testing:
#Pre-filter features that are detected at <50% frequency ( min.pct = 0.5 )
#Pre-filter features that have less than a two-fold change ( logfc.threshold = log(2) )
#Pre-filter features whose detection percentages across the two groups are similar ( min.diff.pct = 0.25 )

# The ROC test returns the 'classification power' for any individual marker
#(ranging from 0 - random, to 1 - perfect). Though not a statistical test,
#it is often very useful for finding clean markers.

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

all_perms = list.files(output, pattern=".rds")
setwd(output)
clusters=readRDS("seurat_integrated_dim_10_2000_samples_clusters.rds")

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

get_degs = function(dat,clus1,clus2){
	combined = dat

	#1. Find cluster markers that are differentially expressed between clusters
	DefaultAssay(combined) <- "integrated"

	#find markers
	diffexp.markers <- FindMarkers(comb, ident.1 = clus1, ident.2 = clus2, min.pct = 0.25,
		 		logfc.threshold = log(2), test.use="roc")
	
	diffexp.markers$gene=rownames(diffexp.markers)
	diffexp.markers  = as.data.table(diffexp.markers)

	output_file = paste(clus1,"vs",clus2,"_diffexp_markers.csv", sep="")
	write.csv(diffexp.markers, file=output_file, quote=F, row.names=F)

	#2. Visualize cluster markers+++++++++++++++++++++++++++++++++++++++++++++++++
	pdf(paste("additional_markers.pdf", sep=""), width=18, height=12)

	genes=c("ICA1", "TOX2", "CTLA4", "ICOS", "TOX", "PTPRG", "BANK1",
	"IGHM", "CDK14", "BCL11A", "IKZF2", "MKI67", "TOP2A")

	v = VlnPlot(combined, features = genes)
	print(v)

	#overlay on UMAP clusters
	f = FeaturePlot(combined, features = genes,
	cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))
	print(f)

	dev.off()
	print("done")

}

get_degs(clusters,clus1="5",clus2="3")
get_degs(clusters,clus1="5",clus2="4")
get_degs(clusters,clus1="5",clus2="6")
get_degs(clusters,clus1="5",clus2="15")

get_degs(clusters,clus1="4",clus2="6")
get_degs(clusters,clus1="4",clus2="3")

get_degs(clusters,clus1="5",clus2="3")
get_degs(clusters,clus1="5",clus2="4")

###additional clusters
#get_degs(clusters,clus1="1",clus2="2")
		#no features pass logfc threshold
#get_degs(clusters,clus1="0",clus2="1")
		#no features pass logfc threshold
get_degs(clusters,clus1="0",clus2="2") #passed

get_degs(clusters,clus1="2",clus2="12")
get_degs(clusters,clus1="10",clus2="11")
get_degs(clusters,clus1="7",clus2="14")

