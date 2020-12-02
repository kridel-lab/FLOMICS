#-------------------------------------------------------------------------------
#002B_Seurat_tutorial.R
#sarah russell
#December 2 2020
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors=F)

#make sure loading R/3.6.1 before running this script

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

#the purpose is obtain differential expression analysis between a select
#group of seurat clusters to help differentiate clusters

#-------------------------------------------------------------------------------
#notes
#-------------------------------------------------------------------------------
#https://satijalab.org/seurat/v3.1/de_vignette.html
#bulk of Seurat’s differential expression features can be accessed through the FindMarkers function.
#As a default, Seurat performs differential expression based on the non-parameteric Wilcoxon rank sum test.
#To test for differential expression between two specific groups of cells, specify the ident.1 and ident.2 parameters.

## Find differentially expressed features between CD14+ and FCGR3A+ Monocytes
#monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono")

##Prefilter features or cells to increase the speed of DE testing
#features that are very infrequently detected in either group of cells, or
#features that are expressed at similar average levels, are unlikely to be differentially expressed

#Pre-filter features that are detected at <50% frequency
# min.pct = 0.5

#Pre-filter features that have less than a two-fold change
#logfc.threshold = log(2)

#Pre-filter features whose detection percentages across the two groups are similar
#min.diff.pct = 0.25

#there are alternative tests, test.use default = wilcox

# The ROC test returns the 'classification power' for any individual marker
#(ranging from 0 - random, to 1 - perfect). Though not a statistical test,
#it is often very useful for finding clean markers.
#ips.markers=find.markers(nbt,7,thresh.use = 2,test.use = "roc")

## Find markers that distinguish clusters 2 and 3 (markers with +avg_diff distinguish C2, markers with -avg_diff characterize C3)
# As noted in Pollen et al., these markers are consistent with different neuronal maturation states
#c2.markers=find.markers(nbt,2,3,thresh.use = 2,test.use = "roc")
#print(head(c2.markers,10))

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

all_perms = list.files(output, pattern=".rds")
setwd(output)

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

get_degs = function(dat,clus1,clus2){
	print(dat)
	combined = readRDS(dat)

	#1. Find cluster markers++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# find markers for every cluster compared to all remaining cells, report only the positive ones
	DefaultAssay(combined) <- "integrated"

	diffexp.markers <- FindMarkers(combined, ident.1 = clus1, ident.2 = clus2,
		 min.pct = 0.25, logfc.threshold = log(2), test.use="roc")
	diffexp.markers$gene=rownames(diffexp.markers)
	diffexp.markers  = as.data.table(diffexp.markers)
	#could also filter for adj p value

	output_file = paste(clus1,clus2,"diffexp_markers.csv", sep="_")
	write.csv(diffexp.markers, file=output_file, quote=F, row.names=F)

	#2. Visualize cluster markers+++++++++++++++++++++++++++++++++++++++++++++++++

	print("done")

}

get_degs("seurat_integrated_dim_10_2000_samples_clusters.rds",
				clus1="3",clus2="5")
