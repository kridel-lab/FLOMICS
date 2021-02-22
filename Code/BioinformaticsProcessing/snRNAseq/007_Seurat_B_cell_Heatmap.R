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

library(dittoSeq)
library(scater)
library(loomR)
library(Seurat)
library(scRNAseq)
library(SingleCellExperiment)

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq")

#load functions to analyze seurat
source("/cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/doSeuratProc.R")
date=Sys.Date()

#output directory
output="/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/Feb2020/"

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

#r=readRDS("combined_processed_seurat_object_rmFL277dim20.rds")
r=readRDS(paste(output, "pc_genes_only_yes_seurat_integrated_dim_20_2000_2021-02-19_samples_clusters.rds", sep=""))

#-------------------------------------------------------------------------------
#purpose
#-------------------------------------------------------------------------------

#visualize differences between B cell clusters

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

combined = r

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------
b_genes = c("CD79A", "BCL2", "CD19", "PTPRG", "BANK1",
"EBF1", "BACH2", "LMO2", "FCER2", "CD27", "ACSM3",
"BACH2", "KLHL6", "EBF1", "FCER2", "LMO2",
"CD200","MS4A1", "ARHGAP24", "IGHM",
"MKI67", "TOP2A", "BCL6", "TBC1D9",
"FCER2", "LMO2", "IL21R")

genesall = c("KIF18B", "MIR155HG", "TOP2A", "MKI67", "ARHGAP24",
"TBC1D9", "HMMR", "CD83", "TPX2", "MT-CO3", "UBE2C",
"MT-ND3", "IGSF3", "VAV3-AS1", "EPHB1", "SAMSN1",
"DEPDC1","PTPRG","SGOL2","TCF4","ASPM","POLQ","CD80","MT-ND4",
"KIF14","FANCA","MT-ATP6","BANK1","SMC4","CENPE","MT-ND1",
"RP11-789C2.1","ANKRD33B","NDC80","BMP7",
"NUSAP1","JADE3","PPM1E","MT-ND2","MT-CYB","PARVB",
"CENPF","AC023590.1","ECT2","SSTR2","CSGALNACT1",
"FAP","NR4A3","FCER2","MT-CO2","CD58","CADPS","EPDR1",
"DOCK10","KIAA0125","SLAMF1","BUB1B","CCL17","SGPP2","ARHGAP11B","IGLC2",
"GPM6A", "MIR4432HG", "ADAMTS6", "ARHGAP24", "SLC7A7", "KHDRBS2","APOLD1",
"MGLL","MIR181A2HG","EZH2","RFC3","CDCA2","RP11-669M16.1","C12orf42",
"HELLS","PECAM1","CTD-2034I4.1","ST6GALNAC3","ZNF804A","VAV3","BATF",
"GALNT14","CCL17","WDFY4","MEF2C","CD80","BCL2","KIAA2022","ARHGAP31","KIF11",
"LYPD6B","BUB1","LINC01572","NR6A1","SLCO5A1","IL21R","BACH2","C1orf186","FCRL5",
"HIP1","RNLS","TP63","KIF13A","NUF2","KLHL6","C2orf48","KYNU","RP11-624C23.1",
"LMO2","LYPD6B","UBE2E2","ADRBK2","ABTB2","RP11-475O6.1","CLSPN","DGKG","MARCKS",
"CCL22")
genesall = unique(genesall)

#get differentially expressed genes between different B cell clusters
get_degs = function(dat,clus1,clus2, clus3){
	seurat_obj = dat

	#1. Find cluster markers that are differentially expressed between clusters
	DefaultAssay(seurat_obj) <- "RNA"
	seurat_obj <- NormalizeData(seurat_obj)

	#find markers
	diffexp.markers <- FindMarkers(seurat_obj, ident.1 = clus1,
		ident.2 = c(clus2, clus3), only.pos=TRUE, min.pct = 0.3, logfc.threshold = 0.25)

	diffexp.markers$gene=rownames(diffexp.markers)
	diffexp.markers  = as.data.table(diffexp.markers)
	diffexp.markers$clust1 = clus1
	diffexp.markers$clust2 = paste(clus2, clus3, sep="_")
	print(head(diffexp.markers))
	return(diffexp.markers)
	print("done")

}

c0vs1 = get_degs(combined, 0, 1, 2)
c0vs4 = get_degs(combined, 0, 4, 10)
c0vs13 = get_degs(combined, 0, 12, 13)
c0vs19 = get_degs(combined, 0, 14, 19)
c1vs2 = get_degs(combined, 1, 0, 2)
c1vs4 = get_degs(combined, 1, 4, 10)
c2vs4 = get_degs(combined, 2, 0, 1)
c2vs10 = get_degs(combined, 2, 4, 10)

all_genes = rbind(c0vs1, c0vs4, c0vs13, c0vs19,
c1vs2, c1vs4, c2vs4, c2vs10)
integrated_genes = rownames(combined)
genes_b_plot_int = unique(filter(all_genes, pct.1 > 0.6, pct.2 < 0.4, avg_logFC > 0.5, gene %in% integrated_genes)$gene)
#genes_b_plot = unique(filter(all_genes, (clust1 == 0 & clust2==1) | (clust1 == 0 & clust2==2) | (clust1 == 1 & clust2==2), avg_logFC > 0.5)$gene)

#cells_t=c("B cells_0", "B cells_1", "B cells_2", "naive B or malignant B_9",
#"proliferating B cell_11", "memory B cell_12", "Cluster 13")

cells_t = c(0, 1, 2, 4, 10, 13, 12, 19, 14)

DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)

#only keep B cell clusters
combined_t <- subset(combined, idents = cells_t)
#only keep differentially expressed genes related to different B cells
subset.matrix <- combined_t[genes_b_plot_int, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
#scale the data
combined_scaled <- ScaleData(subset.matrix, verbose = TRUE)

#make some feature plots for the genes of interest in B cells
pdf(paste(output, date, "_" , "B_cell_genes_featureplot.pdf", sep=""), height=10, width=10)
FeaturePlot(combined, features = genes_b_plot_int, label=TRUE, cols=c("yellow", "purple"))
DotPlot(combined_t, features = genes_b_plot_int) + RotatedAxis()
DoHeatmap(subset(combined_t, downsample = 100), features = genes_b_plot_int, size = 3, assay="integrated")
dev.off()

# Annotating and ordering cells by some meaningful feature(s):
pdf(paste(output, date, "_" , "rmFL277dim20_Bcells_only_markers_dittoheatmap_AUC0.75.pdf", sep=""), height=12, width=15)
dittoHeatmap(combined_scaled, heatmap.colors.max.scaled=colorRamps::blue2yellow(25), scaled.to.max=TRUE, annot.by = c("seurat_clusters"), assay="integrated")
dev.off()

h = DoHeatmap(combined_t, features = genes_b_plot, assay="RNA", size=2)
pdf(paste(date, "_" , "rmFL277dim20_Bcells_only_markers_AUC0.75.pdf", sep=""), height=20)
print(h)
dev.off()
