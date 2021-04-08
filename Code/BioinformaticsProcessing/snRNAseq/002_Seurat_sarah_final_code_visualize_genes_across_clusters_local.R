#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors=F)

#make sure loading R/4.0.0 before running this script

library(dittoSeq)
library(scater)
#library(loomR)
library(Seurat)
library(scRNAseq)
library(SingleCellExperiment)
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(ReactomeGSA)

setwd("~/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April82021")

date=Sys.Date()

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

r=readRDS("pc_genes_only_no_seurat_integrated_dim_20_2000_2021-04-08_samples_clusters.rds")

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

combined = r
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)

#1. Find cluster markers++++++++++++++++++++++++++++++++++++++++++++++++++++++

# find markers for every cluster compared to all remaining cells, report only the positive ones
#using integrated seurat object

#RUN ONCE 
#combined.markers <- FindAllMarkers(combined, only.pos = TRUE,
#	 min.pct = 0.4, logfc.threshold = 0.3)
#combined.markers = as.data.table(combined.markers)

#save results 
#saveRDS(combined.markers, file="combined_markers_all_clusters.rds")

combined.markers = readRDS("combined_markers_all_clusters.rds")

top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
combined_markers <- combined[top10$gene, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest

pdf("marker_genes_top_10_heatmap.pdf", height=20, width=10)
dittoHeatmap(combined_markers, annot.by = c("seurat_clusters"), scaled.to.max = TRUE)
dev.off()

#save markers in file
markers = as.data.table(combined.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC))
file_name=paste(date, "_samples_clusters_markers_all_integrated_object.txt", sep="")
write.table(markers, file_name, row.names=F, sep=";", quote=F)

#2. Visualize cluster markers+++++++++++++++++++++++++++++++++++++++++++++++++

pdf(paste(date, "_", "seurat_cluster_markers", ".pdf", sep=""),width=18, height=20)

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
					    "CD69", "CD45RA", "CD45RO", "CCR7", "ICOS"), order=TRUE, label=TRUE,
cols=c("grey", "thistle1", "steelblue", "red"))
print(f)

igh = FeaturePlot(combined, features = c("IGHM", "IGHD", "IGKC", "IGLC2", "IGLC3", "BCL2", "CD19"),order=TRUE, label=TRUE,
                 cols=c("grey", "thistle1", "steelblue", "red"))
print(igh)

#overlay on UMAP clusters
f2 = FeaturePlot(combined, features = c("CD79A", "CR2", "CD3D", "CCL5",
"ICOS", "IGKC", "IGLC2", "IGLC3", "BCL2", "CD19"),order=TRUE, label=TRUE,
cols=c("grey", "thistle1", "steelblue", "red"))
print(f2)

#overlay on UMAP clusters
f3= FeaturePlot(combined, features = c("STAT3", "CD27", "CCR5", "CCR4", "CXCR4", "LAG3",
"CD19", "IGHM", "PTPRG", "TNFRSF18", "TNFRSF4", "FAS"),order=TRUE, label=TRUE,
cols=c("grey", "thistle1", "steelblue", "red"))
print(f3)

#overlay on UMAP clusters
f4 = FeaturePlot(combined, features = c("CD3G", "CD8A", "BANK1", "VIM", "LYZ",
"ST8SIA1", "ICA1", "IL7R", "TRAC", "CTLA4", "IKZF2"),order=TRUE, label=TRUE,
cols=c("grey", "thistle1", "steelblue", "red"))
print(f4)

#B cell genes
f5 = FeaturePlot(combined, features = c("EBF1", "BACH2", "LMO2", "FCER2", "CD27", "ACSM3",
"KLHL6", "BCL2", "CD200", "MS4A1", "ARHGAP24", "IGHM", "BANK1", "PTPRG",
	"MKI67", "TOP2A", "BCL6", "TBC1D9", "LMO2", "IL21R", "CD44", "PDL2", "CD80", "NT5E"),order=TRUE, label=TRUE,
cols=c("grey", "thistle1", "steelblue", "red"))
print(f5)

#more B cell genes from https://www.nature.com/articles/s41590-018-0181-4#Sec32
f5b = FeaturePlot(combined, features = c("NOTCH2","JAM3", "PRDM1", "SDC1", "BCL6", "AICDA", "IRF4",
"MS4A1", "PAX5", "HLA-DRA", "POLH", "P2RY8", "GNA13", "SDC1"),order=TRUE, label=TRUE,
cols=c("grey", "thistle1", "steelblue", "red"))
print(f5b)

#T cell genes
f6 = FeaturePlot(combined, features = c("CD3D", "CD3G", "CD4", "CD8A", "TCF7", "PTPRC", "TIGIT", "PDCD1", "TOX",
	"TOX2", "TNFSF8", "PTPN13", "ILR3", "BTLA", "CD200", "ICOS", "IL21", "CCL5", "GZMK", "GZMA",
	"PRDM1", "KLRG1", "TIGIT", "HAVCR2", "EOMES", "CTLA4", "TOX2", "FOXP3", "IL6R", "ICOS",
	"MKI67", "TCF7", "TOP2A"),order=TRUE, label=TRUE,
cols=c("grey", "thistle1", "steelblue", "red"))
print(f6)

#more T cell genes based on plot from https://www.nature.com/articles/s41556-020-0532-x#Fig2
f7 = FeaturePlot(combined, features = c("IL7R", "CD4", "PLAC8", "KLF2","GZMK", "CCL5",
"GZMA", "NKG7", "CCL4", "CD8A", "PDCD1", "TOX", "TOX2", "CD200", "CXCR5", "ICOS",
"IL2RA", "FOXP3"),order=TRUE, label=TRUE,cols=c("grey", "thistle1", "steelblue", "red"))
print(f7)

#other cluster genes
f8 = FeaturePlot(combined, features = c("CXCL12", "SDF1", "VCAM1", "CR2", "CD21", "CD36",
"TFPI", "LDB2", "CD3G", "HAVCR2", "CD274", "ILR6", "CD68", "CSF1R", "IFNGR1",
"SPARCL1", "PECAM1", "VCAM1", "CD36", "CUX2", "SLC7A11", "CD4"),order=TRUE, label=TRUE,
cols=c("grey", "thistle1", "steelblue", "red"))
print(f8)

dev.off()

pdf("IGHM_IGHD.pdf", width=15, height=6)
FeaturePlot(combined, features = c("IGHM", "IGHD"), order=TRUE, label=TRUE, 
            blend.threshold=0.6, blend=TRUE)
dev.off()

pdf("BCL2_BCL6.pdf", width=15, height=6)
FeaturePlot(combined, features = c("BCL2", "BCL6"), order=TRUE, label=TRUE, 
            blend.threshold=0.7, blend=TRUE)
dev.off()

