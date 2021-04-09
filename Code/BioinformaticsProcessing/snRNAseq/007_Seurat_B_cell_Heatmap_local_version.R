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
library("xlsx")
library(tidyr)

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
#purpose
#-------------------------------------------------------------------------------

#visualize differences between B cell clusters

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

combined = r
DimPlot(combined)
dev.off()

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

normalized_combined = combined
DefaultAssay(normalized_combined) <- "RNA"
normalized_combined <- NormalizeData(normalized_combined)


#get differentially expressed genes between different B cell clusters
get_degs = function(clus1,clus2, clus3){
	
  print(paste(clus1, clus2, clus3))
  seurat_obj = normalized_combined

	#1. Find cluster markers that are differentially expressed between clusters

	#find markers
	diffexp.markers <- FindMarkers(seurat_obj, ident.1 = clus1,
		ident.2 = c(clus2, clus3), only.pos=TRUE, min.pct = 0.3, logfc.threshold = 0.25)

	diffexp.markers$gene=rownames(diffexp.markers)
	diffexp.markers  = as.data.table(diffexp.markers)
	diffexp.markers$clust1 = clus1
	diffexp.markers$clust2 = paste(clus2, clus3, sep="_")
#	diffexp.markers = filter(diffexp.markers, p_val_adj < 0.05)
	diffexp.markers = diffexp.markers[order(-avg_log2FC)]
	
	genes = filter(diffexp.markers, p_val_adj < 0.05)$gene
	gs = as.data.table(gprofiler2::gost(genes, exclude_iea=TRUE, ordered_query=TRUE, 
	                                 user_threshold = 0.01, sources=c("REAC")))

	print(head(gs$result.term_name))
  print(head(diffexp.markers))
	return(as.data.table(diffexp.markers))
	print("done")

}

###
# B cells
###

#P0 - 
#P1 - 
#P2 - 
#P4 - 
#P11 -  
#P12 - proliferating

#####################################
#RUN once
#####################################


cells_b = c(0, 1, 2, 4, 11, 12)

#all possible comparisons 
comps = as.data.table(crossing(var1 = cells_b, var2 = cells_b, var3 = cells_b)) %>% 
    filter(!(var1 == var2), !(var1==var3), !(var2==var3), !(var1 == 12))
comps

#all_genes = mapply(get_degs, comps$var1, comps$var2, comps$var3, SIMPLIFY=FALSE)
#all_genes_combined <- do.call("rbind", all_genes)
#tail(all_genes_combined)
#all_genes = as.data.table(ldply(all_genes))

#saveRDS(all_genes, file="B_cell_clusters_diff_exp_genes.rds")

all_genes = readRDS("B_cell_clusters_diff_exp_genes.rds")

all_genes = filter(all_genes, p_val_adj < 0.05)
write.xlsx(all_genes, file, file = "Bcell_marker_genes_across_clusters.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

integrated_genes = rownames(combined)
genes_b_plot_int = unique(filter(all_genes, pct.1 > 0.5, pct.2 < 0.5, avg_log2FC > 0.5, gene %in% integrated_genes)$gene)

DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)

#only keep B cell clusters
combined_b <- subset(combined, idents = cells_b)
#only keep differentially expressed genes related to different B cells
subset.matrix <- combined_b[genes_b_plot_int, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
#scale the data
combined_scaled <- ScaleData(subset.matrix, verbose = TRUE)

#naive B cell markers
pdf("IGHD_IGHM_coexpression_Bcells.pdf", width=10, height=6)
cells_wexp = WhichCells(object = combined_b, expression = IGHD > 1 & IGHM >1)
FeaturePlot(combined_b, features = c("IGHD", "IGHM"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1)
FeaturePlot(combined_b, features = c("IGHD", "IGHM"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.1, cells=cells_wexp)
dev.off()

pdf("BCL2_BCL6.pdf", width=8, height=6)
FeaturePlot(combined_b, features = c("BCL2", "BCL6"), order=TRUE, label=TRUE, 
            blend.threshold=0.7, cols=c("grey", "thistle1", "steelblue", "red"))
dev.off()

FeaturePlot(combined_b, features = c("CD83", "CXCR4"), blend = TRUE, order=TRUE, label=TRUE, blend.threshold=0.6)
ggsave("LZ_DZ_CD83_CXCR4.pdf")

FeaturePlot(combined_b, features = filter(all_genes, clust1==1)[1:12]$gene, cols=c("grey", "thistle1", "steelblue", "red"),
            order=TRUE, min.cutoff='q15', label=TRUE, ncol=3)
ggsave("Cluster1_top12_enriched_Bcell_genes_featureplots.pdf", width=12, height=15)

FeaturePlot(combined_b, features = filter(all_genes, clust1==2)[1:12]$gene, cols=c("grey", "thistle1", "steelblue", "red"),
            order=TRUE, min.cutoff='q15', label=TRUE, ncol=3)
ggsave("Cluster2_top12_enriched_Bcell_genes_featureplots.pdf", width=12, height=15)

DZ_genes=c("AICDA", "CXCR4", "GCSAM", "CD27", "SEMA4B", "MKI67", "EZH2", "CCNB1", "BACH2", "AURKA", "RAD51", "POLH")
DotPlot(combined_b, features = DZ_genes, idents = cells_b, cols=c("green", "red"), cluster.idents=TRUE) + RotatedAxis()
ggsave("Bcell_DZ_dotplot.pdf", width=7, height=5)

FeaturePlot(combined_b, features = DZ_genes, cols=c("grey", "thistle1", "steelblue", "red"),
            order=TRUE, min.cutoff='q15', label=TRUE, ncol=3)
ggsave("Bcell_DZ_genes_featureplots.pdf", width=12, height=15)

LZ_genes=c("SEMA7A", "BCL2A1", "CD83", "PRDM1", "SLAMF1", "IRF4", "LY75", "CD40", "PTPN6")

FeaturePlot(combined_b, features = LZ_genes, cols=c("grey", "thistle1", "steelblue", "red"),
            order=TRUE, min.cutoff='q15', label=TRUE, ncol=3)
ggsave("Bcell_LZ_genes_featureplots.pdf", width=12, height=15)

DotPlot(combined_b, features = LZ_genes, idents = cells_b, cols=c("green", "red"), cluster.idents=TRUE) + RotatedAxis()
ggsave("Bcell_LZ_dotplot.pdf", width=7, height=5)

#make some feature plots for the genes of interest in B cells
pdf(paste(date, "_" , "B_cell_genes_dotplot.pdf", sep=""), height=7, width=15)
DotPlot(combined_b, features = genes_b_plot_int, cols=c("green", "red"), cluster.idents=TRUE) + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf(paste(date, "_" , "B_cell_genes_heatmap.pdf", sep=""), height=8, width=8)
DoHeatmap(combined_b, features = genes_b_plot_int, size = 3, assay="integrated")
dev.off()



