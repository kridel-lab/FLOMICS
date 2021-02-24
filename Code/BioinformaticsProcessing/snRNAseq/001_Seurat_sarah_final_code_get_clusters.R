#-------------------------------------------------------------------------------
#001_Seurat.R
#karin isaev (using base code from RK)
#September 30th 2020
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors=F)
library(tidyverse)

#make sure loading R/4.0.0 before running this script

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA")

#load functions to analyze seurat
source("/cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/doSeuratProc.R")

#output directory
output="/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/Feb2020/"

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "tidyr", "data.table", "plyr",
	"stringr",
  "Seurat",
  "cowplot",
	"patchwork")

lapply(packages, require, character.only = TRUE)

library(annotables)

date=Sys.Date()

args = commandArgs(trailingOnly = TRUE) #patient ID
input = args[1]
print(input) #whether only protein coding genes should be included or not
norm_type = args[2]
print(norm_type)

#genes = as.data.table(grch38)
genes = as.data.table(grch37)
pc_genes = unique(filter(genes, biotype == "protein_coding")$symbol)

#-------------------------------------------------------------------------------
#purpose
#-------------------------------------------------------------------------------

#the purpose is obtain clusters from snRNAseq that are representative
#of cell types which we can then use to convolute what cell types are present
#in the bulk RNA-seq samples

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

#1. read in data++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Following vignette in https://urldefense.com/v3/__https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html__;!!CjcC7IQ!ZRlLMV8nCj5QEeXhD5snLFSleb1Xd8O3zRcRHiR-f4VmextcM9CV5IZg961roEutzHzYWQ$ [satijalab[.]org] latest updated version
# and this one for integration of multiple seurat objects:
# https://urldefense.com/v3/__https://satijalab.org/seurat/v3.0/immune_alignment.html__;!!CjcC7IQ!ZRlLMV8nCj5QEeXhD5snLFSleb1Xd8O3zRcRHiR-f4VmextcM9CV5IZg961roEsQLA-62A$ [satijalab[.]org]

#most recent folder locations
#LY_FL_062_T1: snRNAseq/191218_A00827_0099_AHMW73DMXX_Robert_Kridel/Robert_Kridel__LY_FL_062_T1
#LY_FL_064_T1: snRNAseq/200420_A00827_0152_AHT2YJDMXX_Kridel_Robert/Kridel_Robert__LY_FL_064_T1
#LY_FL_076_T1: snRNAseq/200420_A00827_0152_AHT2YJDMXX_Kridel_Robert/Kridel_Robert__LY_FL_076_T1
#LY_FL_227_T1.: snRNAseq/200317_A00827_0147_BH75CGDSXY_Kridel_Robert/Kridel_Robert__LY_FL_227_T1

data_dir_FL062 <- "snRNAseq/191218_A00827_0099_AHMW73DMXX_Robert_Kridel/Robert_Kridel__LY_FL_062_T1/filtered_feature_bc_matrix"
data_dir_FL064 <- "snRNAseq/200420_A00827_0152_AHT2YJDMXX_Kridel_Robert/Kridel_Robert__LY_FL_064_T1/filtered_feature_bc_matrix"
data_dir_FL076 <- "snRNAseq/200420_A00827_0152_AHT2YJDMXX_Kridel_Robert/Kridel_Robert__LY_FL_076_T1/filtered_feature_bc_matrix"
#data_dir_FL227 <- "snRNAseq/200317_A00827_0147_BH75CGDSXY_Kridel_Robert/Kridel_Robert__LY_FL_227_T1/filtered_feature_bc_matrix"

list.files(data_dir_FL062)
list.files(data_dir_FL064)
list.files(data_dir_FL076)

#list.files(data_dir_FL227) #"barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"

expression_matrix_FL062 <- Read10X(data.dir = data_dir_FL062)
expression_matrix_FL064 <- Read10X(data.dir = data_dir_FL064)
expression_matrix_FL076 <- Read10X(data.dir = data_dir_FL076)
#expression_matrix_FL227 <- Read10X(data.dir = data_dir_FL227)

#2. Setup the Seurat objects and normalize++++++++++++++++++++++++++++++++++++++

exp_matrices = list(expression_matrix_FL062, expression_matrix_FL064, expression_matrix_FL076) #, expression_matrix_FL227)
names(exp_matrices) = c("FL062", "FL064", "FL076") #, "FL227")
samps = c("FL062", "FL064", "FL076") #, "FL227")

pdf(paste(output, "pc_genes_only_", input, "_", norm_type, "_", "seurat_objects_qc_vln_plots.pdf", sep=""), width=16, height=8)
all_objects = mapply(doSeuratProc, exp_matrices, samps, mito_rm="yes", nc_rm=input, norm_type=norm_type)
dev.off()

# plot variable features with and without labels (before sample integration)
pdf(paste(output, "pc_genes_only_", input, "_", norm_type, "_", "seurat_top10_genes_per_sample.pdf", sep=""), width=16, height=8)
for(i in 1:3){
  print(i)
  # Identify the 10 most highly variable genes in the first sample of the list
  top10 <- head(VariableFeatures(all_objects[[i]]), 20)
  plot1 <- VariableFeaturePlot(all_objects[[i]]) + ggtitle(all_objects[[i]]@meta.data$sample[1])
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  p = plot1 + plot2
  print(p)
}
dev.off()

#3. Set up anchors++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#dat = all_objects
#dim = 10
#anch_features = 3000

get_integrated_obj = function(dat, dim, anch_features, norm_method_used){

	set.seed(100)
	print(dim)
	print(anch_features)

	if(norm_method_used == "SC"){

		# Select the most variable features to use for integration
		integ_features <- SelectIntegrationFeatures(object.list = dat,
                                            nfeatures = 2000)

		# Prepare the SCT list object for integration
		dat <- PrepSCTIntegration(object.list = dat,
			 anchor.features = integ_features)

		# Find best buddies - can take a while to run
		integ_anchors <- FindIntegrationAnchors(object.list = dat,
			                                         normalization.method = "SCT",
			                                         anchor.features = integ_features)

		# Integrate across conditions
		seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

		# Run PCA
		seurat_integrated <- RunPCA(object = seurat_integrated)

		# Run UMAP
		seurat_integrated <- RunUMAP(seurat_integrated,
		                             dims = 1:dim, reduction = "pca")

		# Explore heatmap of PCs
		pcs_plot = DimHeatmap(seurat_integrated,
			           dims = 1:dim,
			           cells = 500,
			           balanced = TRUE)

		# Determine the K-nearest neighbor graph
		seurat_integrated <- FindNeighbors(object = seurat_integrated,
				                                dims = 1:dim)

		# Determine the clusters for various resolutions
		seurat_integrated <- FindClusters(object = seurat_integrated, resolution = c(0.4, 0.5, 0.6, 0.8, 1.0, 1.4))

		# Assign identity of clusters
		Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
		# Plot the UMAP
		umap_integrated_res = DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

		p3 = DimPlot(seurat_integrated,
							label = TRUE,
							split.by = "sample")  + NoLegend()

		saveRDS(seurat_integrated, file = paste(output, "pc_genes_only_", input, "_", "seurat_integrated_SCnorm_dim_", dim , "_", anch_features, "_", date, "_samples_clusters.rds", sep=""))
	}

	if(!(norm_method_used == "SC")){

	anchors <- FindIntegrationAnchors(object.list = dat, dims = 1:dim,
		anchor.features = anch_features)

	combined <- IntegrateData(anchorset = anchors, dims = 1:dim)

	# specify that we will perform downstream analysis on the corrected data note that the original
	# unmodified data still resides in the 'RNA' assay
	DefaultAssay(combined) <- "integrated"

	#4. Run clustering++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# Run the standard workflow for visualization and clustering
	combined <- ScaleData(combined, verbose = FALSE)
	combined <- RunPCA(combined, verbose = FALSE)

	# Examine and visualize PCA results a few different ways
	print(combined[["pca"]], dims = 1:20, nfeatures = 5)
	pdf(paste(output, "pc_genes_only_", input, "_", "seurat_integrated_dim_dimensionality_exploration", dim , "_", anch_features, "_", date, "_samples_clusters.pdf", sep=""))
	v = VizDimLoadings(combined, dims = 1:2, reduction = "pca")
	d = DimHeatmap(combined, dims = 1:dim, cells = 500, balanced = TRUE)
	el = ElbowPlot(combined)
	print(v)
	print(d)
	print(el)
	dev.off()

	# t-SNE and Clustering
	combined <- FindNeighbors(combined, reduction = "pca", dims = 1:dim)
	combined <- FindClusters(combined, resolution = 0.3)
	combined <- RunUMAP(combined, reduction = "pca", dims = 1:dim)

	pdf(paste(output, "pc_genes_only_", input, "_", "seurat_integrated_dim_", dim , "_", anch_features, "_", date, "_samples_clusters.pdf", sep=""), width=13, height=6)
	p1 <- DimPlot(combined, reduction = "umap", group.by = "sample")+
	theme(axis.line = element_line(colour = 'black', size = 1), text = element_text(size = 20), axis.text = element_text(size = 20))
	p2 <- DimPlot(combined, reduction = "umap", label = TRUE, label.size=6)+
	theme(axis.line = element_line(colour = 'black', size = 1), text = element_text(size = 20), axis.text = element_text(size = 20))
	p3 = DimPlot(combined,
        label = TRUE,
        split.by = "sample")  + NoLegend()
	print(p1 + p2)
	print(p3)
	dev.off()

	saveRDS(combined, file = paste(output, "pc_genes_only_", input, "_", "seurat_integrated_dim_", dim , "_", anch_features, "_", date, "_samples_clusters.rds", sep=""))
	print("finished this analysis")

}
} #end function

get_integrated_obj(all_objects, dim=20, 2000, norm_type)

sessionInfo()
