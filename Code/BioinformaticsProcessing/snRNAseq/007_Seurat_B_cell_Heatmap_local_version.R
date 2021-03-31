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
#library(loomR)
library(Seurat)
library(scRNAseq)
library(SingleCellExperiment)
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(ReactomeGSA)

setwd("~/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/Feb232021")

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

r=readRDS("pc_genes_only_no_seurat_integrated_dim_20_2000_2021-02-23_samples_clusters.rds")

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

genes_check = c("NOTCH2", "JAM3", "PRDM1", "IRF4", "MS4A1", "PAX5",
                "HLA-DRA", "BCL6", "POLH", "P2RY8", "GNA13", "AICDA",
                "LMO2")

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
	diffexp.markers = diffexp.markers[order(-avg_log2FC)]

	genes = filter(diffexp.markers, p_val_adj < 0.05)$gene
	gs = as.data.table(gprofiler2::gost(genes, exclude_iea=TRUE, ordered_query=TRUE))

	print(gs$result.term_name)

	return(diffexp.markers)
	print("done")

}

###
# B cells
###
#P0
#*P1
#*P2
#*P7
#*P9
#P12 - proliferating B cells
#P13

c0vs1 = get_degs(combined, 0, 1, 2)
c0vs9 = get_degs(combined, 0, 9, 7)
c0vs13 = get_degs(combined, 0, 9, 13)

c1vs2 = get_degs(combined, 1, 0, 2)
c1vs9 = get_degs(combined, 1, 9, 7)

c2vs1 = get_degs(combined, 2, 0, 1)
c2vs9 = get_degs(combined, 2, 9, 7)

c9vs1 = get_degs(combined, 9, 0, 1)
c9vs2 = get_degs(combined, 9, 1, 2)

c13c1 = get_degs(combined, 13, 0, 1)
c13c2 = get_degs(combined, 13, 1, 2)
c13c9 = get_degs(combined, 13, 1, 9)

c7c1 = get_degs(combined, 7, 0, 1)
c7c2 = get_degs(combined, 7, 1, 2)
c7c9 = get_degs(combined, 7, 1, 9)

all_genes = rbind(c0vs1, c0vs9, c0vs13, c1vs2, c1vs9, c2vs1, c2vs9, c9vs1, c9vs2,
                  c13c1, c13c2, c13c9, c7c1, c7c2,  c7c9)
integrated_genes = rownames(combined)
genes_b_plot_int = unique(filter(all_genes, pct.1 > 0.5, pct.2 < 0.5, avg_log2FC > 0.5, gene %in% integrated_genes)$gene)

cells_b = c(0, 1, 2, 7, 9, 12, 13)

DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)

#only keep B cell clusters
combined_b <- subset(combined, idents = cells_b)
#only keep differentially expressed genes related to different B cells
subset.matrix <- combined_b[genes_b_plot_int, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
#scale the data
combined_scaled <- ScaleData(subset.matrix, verbose = TRUE)

gsva_result <- analyse_sc_clusters(combined_b, verbose = TRUE, use_interactors=F)
pathway_expression <- pathways(gsva_result)
pathway_expression$X12.Seurat = NULL

# find the maximum differently expressed pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max_difference$diff <- max_difference$max - max_difference$min

# sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

head(max_difference)
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1])

#naive B cell markers
pdf("IGHD_IGHM_coexpression_Bcells.pdf", width=10, height=6)
cells_wexp = WhichCells(object = combined_b, expression = IGHD > 1 & IGHM >1)
FeaturePlot(combined_b, features = c("IGHD", "IGHM"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1)
FeaturePlot(combined_b, features = c("IGHD", "IGHM"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1, cells=cells_wexp)
dev.off()

pdf("naive_Bcell_marker_diff_exp_analysis.pdf", width=10, height=6)
cells_wexp = WhichCells(object = combined_b, expression = CCSER1 > 1 & KHDRBS2 >1)
FeaturePlot(combined_b, features = c("CCSER1", "KHDRBS2"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1)
FeaturePlot(combined_b, features = c("CCSER1", "KHDRBS2"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1, cells=cells_wexp)
dev.off()

pdf("markers_of_follicular_B_cells.pdf", width=10, height=6)
cells_wexp = WhichCells(object = combined_b, expression = PAX5 > 3 & 'HLA-DRA' >3)
FeaturePlot(combined_b, features = c("PAX5", "HLA-DRA"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.6)
FeaturePlot(combined_b, features = c("PAX5", "HLA-DRA"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.6, cells=cells_wexp)
dev.off()

pdf("HLA_genes_B_cells.pdf", width=10, height=6)
#cells_wexp = WhichCells(object = combined_b, expression = "HLA-A" > 1 & "HLA-B" >1)
FeaturePlot(combined_b, features = c("HLA-B", "HLA-A"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.6)
#FeaturePlot(combined_b, features = c("HLA-B", "HLA-A"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.6, cells=cells_wexp)
dev.off()

FeaturePlot(combined_b, features = c("CD80"), order=TRUE, label=TRUE, blend.threshold=0.1, cols=c("grey", "thistle1", "steelblue", "red"))

FeaturePlot(combined_b, features = c2vs1$gene[1:12], cols=c("grey", "thistle1", "steelblue", "red"),
            order=TRUE, min.cutoff='q15', label=TRUE, ncol=3)
ggsave("Cluster2_top12_enriched_Bcell_genes_featureplots.pdf", width=12, height=15)

FeaturePlot(combined_b, features = c9vs2$gene[1:12], cols=c("grey", "thistle1", "steelblue", "red"),
            order=TRUE, min.cutoff='q15', label=TRUE, ncol=3)
ggsave("Cluster9_top12_enriched_Bcell_genes_featureplots.pdf", width=12, height=15)

FeaturePlot(combined_b, features = genes_check, cols=c("grey", "thistle1", "steelblue", "red"),
            order=TRUE, min.cutoff='q15', label=TRUE, ncol=3)
ggsave("Bcell_genes_featureplots.pdf", width=12, height=15)

DotPlot(combined_b, features = genes_check, idents = cells_b) + RotatedAxis()
ggsave("Bcell_genes_dotplot.pdf")

DZ_genes=c("AICDA", "CXCR4", "GCSAM", "CD27", "SEMA4B", "MKI67", "EZH2", "CCNB1", "BACH2", "AURKA", "RAD51", "POLH")
DotPlot(combined_b, features = DZ_genes, idents = cells_b) + RotatedAxis()
ggsave("Bcell_DZ_dotplot.pdf")

FeaturePlot(combined_b, features = DZ_genes, cols=c("grey", "thistle1", "steelblue", "red"),
            order=TRUE, min.cutoff='q15', label=TRUE, ncol=3)
ggsave("Bcell_DZ_genes_featureplots.pdf", width=12, height=15)

LZ_genes=c("SEMA7A", "BCL2A1", "CD83", "PRDM1", "SLAMF1", "IRF4", "LY75", "CD40", "PTPN6")

FeaturePlot(combined_b, features = LZ_genes, cols=c("grey", "thistle1", "steelblue", "red"),
            order=TRUE, min.cutoff='q15', label=TRUE, ncol=3)
ggsave("Bcell_LZ_genes_featureplots.pdf", width=12, height=15)

DotPlot(combined_b, features = LZ_genes, idents = cells_b) + RotatedAxis()
ggsave("Bcell_LZ_dotplot.pdf")

FeaturePlot(combined_b, features = "LRMP", cols=c("grey", "thistle1", "steelblue", "red"),
            order=TRUE, min.cutoff='q15', label=TRUE, ncol=3)

#make some feature plots for the genes of interest in B cells
pdf(paste(date, "_" , "B_cell_genes_featureplot.pdf", sep=""), height=10, width=10)
DotPlot(combined_b, features = genes_b_plot_int) + RotatedAxis()
dittoHeatmap(combined_b, genes = genes_b_plot_int, scaled.to.max=TRUE, assay="integrated")
DoHeatmap(subset(combined_b, downsample = 100), features = genes_b_plot_int, size = 3, assay="integrated")
dev.off()

# Annotating and ordering cells by some meaningful feature(s):
pdf(paste(date, "_" , "B_cell_genes_ditto_heatmap_scaled.pdf", sep=""), height=20, width=20)
FeaturePlot(combined, features = genes_b_plot_int, cols=c("yellow", "purple"))
dittoHeatmap(combined_scaled, scaled.to.max=TRUE, assay="integrated")
dittoHeatmap(combined_scaled, scaled.to.max=TRUE, assay="RNA")
dev.off()

h = DoHeatmap(combined_b, features = genes_b_plot_int, assay="integrated", size=2)
pdf(paste(date, "_" , "B_cell_genes_ditto_heatmap_seurat.pdf", sep=""), height=6)
print(h)
dev.off()
