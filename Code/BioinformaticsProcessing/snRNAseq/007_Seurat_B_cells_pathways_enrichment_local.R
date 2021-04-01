#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors=F)

#make sure loading R/4.0.0 before running this script

#library(dittoSeq)
#library(scater)
#library(loomR)
library(Seurat)
#library(scRNAseq)
#library(SingleCellExperiment)
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(ReactomeGSA)
library(enrichR)
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

setwd("~/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/Feb232021")

date=Sys.Date()

r=readRDS("pc_genes_only_no_seurat_integrated_dim_20_2000_2021-02-23_samples_clusters.rds")

mainDir="/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat"
subDir=paste("B_cells_pathways", date, sep="_")

dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir, subDir))

#-------------------------------------------------------------------------------
#purpose
#-------------------------------------------------------------------------------

#visualize differences between B cell clusters

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

combined = r
cells_b = c(0, 1, 2, 7, 9, 12, 13)

DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)

#only keep B cell clusters
combined_b <- subset(combined, idents = cells_b)
gsva_result <- analyse_sc_clusters(combined_b, verbose = TRUE, use_interactors=F)

#-------------------------------------------------------------------------------
#functions
#-------------------------------------------------------------------------------

#1. GSVA

gsva_analysis = function(dat, clust1, clust2){

	pathway_expression <- pathways(dat)
	colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))

	# find the maximum differently expressed pathway
	get_vals = function(x){
		values <- as.numeric(x[2:length(x)])
		clust1_val = x[which(names(x) == clust1)]
		clust2_val = x[which(names(x) == clust2)]
    return(data.frame(name = x[1], min = min(values),
		max = max(values),
		mean = mean(values),
		clust1=as.numeric(clust1_val),
		clust2=as.numeric(clust2_val)))
	}

	max_difference <- do.call(rbind, apply(pathway_expression, 1, get_vals))
	max_difference = subset(max_difference, clust1 >0, clust2 >0)
	max_difference$diff <- max_difference$clust1 - max_difference$clust2

	# sort based on the difference
	max_difference <- max_difference[order(abs(max_difference$diff), decreasing = T), ]

	my_palette <- colorRampPalette(c("purple", "grey", "orange"))(n = 299)

	head(max_difference)
	p = plot_gsva_heatmap(dat,
		scale = "row",
		#margins = c(6,30),
		col=my_palette,
		truncate_names=TRUE,
		#dendrogram = "col",
		pathway_ids = rownames(max_difference)[1:20],
		#lwid=c(0.1,4),
		key = FALSE, main = paste(clust1, "vs", clust2, "top 20 pathways"))

	print(p)
}

#2. DEenrichRPlot (Seurat + EnrichR)

enrichr_analysis = function(dat, clust1, clust2, db){

	paths = DEenrichRPlot(dat, ident.1 = clust1,
		ident.2=clust2, balanced=TRUE,
		enrich.database = db, max.genes = 75,
		logfc.threshold = 0.5, num.pathway=20)
	print(paths)

}

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

#1. get GSVA pathways comparing different clusters

pdf("GSVA_analysis_cluster1_vs_cluster2.pdf", width=10)
gsva_analysis(gsva_result, "X1", "X2") #X1 = seurat cluster 1
dev.off()

pdf("GSVA_analysis_cluster1_vs_cluster9.pdf", width=10)
gsva_analysis(gsva_result, "X1", "X9")
dev.off()

pdf("GSVA_analysis_cluster2_vs_cluster9.pdf", width=10)
gsva_analysis(gsva_result, "X2", "X9")
dev.off()

pdf("GSVA_analysis_cluster7_vs_cluster13.pdf", width=10)
gsva_analysis(gsva_result, "X7", "X13")
dev.off()

#2. Analyze pathways using EnrichR comparing different clusters

paths_db = c("MSigDB_Hallmark_2020",
"WikiPathways_2019_Human",
"KEGG_2019_Human",
"MGI_Mammalian_Phenotype_Level_4_2019")

get_db_enrichr_plots = function(db_type, clust1, clust2){
	print(db_type)
	print(clust1)
	print(clust2)
	enrichr_analysis(combined_b, clust1, clust2, db_type)
}

pdf("EnrichR_analysis_cluster1_vs_cluster2.pdf", width=12, height=8)
llply(paths_db, get_db_enrichr_plots, 1, 2)
dev.off()
