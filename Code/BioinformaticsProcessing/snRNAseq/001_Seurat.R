#-------------------------------------------------------------------------------
#001_Seurat.R
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
packages <- c("dplyr", "readr", "ggplot2", "tidyr", "data.table", "plyr",
	"stringr",
  "Seurat",
  "cowplot",
	"patchwork")

lapply(packages, require, character.only = TRUE)

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

# Following vignette in https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html latest updated version
# and this one for integration of multiple seurat objects:
# https://satijalab.org/seurat/v3.0/immune_alignment.html

#most recent folder locations
#LY_FL_062_T1: snRNAseq/191218_A00827_0099_AHMW73DMXX_Robert_Kridel/Robert_Kridel__LY_FL_062_T1
#LY_FL_064_T1: snRNAseq/200420_A00827_0152_AHT2YJDMXX_Kridel_Robert/Kridel_Robert__LY_FL_064_T1
#LY_FL_076_T1: snRNAseq/200420_A00827_0152_AHT2YJDMXX_Kridel_Robert/Kridel_Robert__LY_FL_076_T1
#LY_FL_227_T1.: snRNAseq/200317_A00827_0147_BH75CGDSXY_Kridel_Robert/Kridel_Robert__LY_FL_227_T1

data_dir_FL062 <- "snRNAseq/191218_A00827_0099_AHMW73DMXX_Robert_Kridel/Robert_Kridel__LY_FL_062_T1/filtered_feature_bc_matrix"
data_dir_FL064 <- "snRNAseq/200420_A00827_0152_AHT2YJDMXX_Kridel_Robert/Kridel_Robert__LY_FL_064_T1/filtered_feature_bc_matrix"
data_dir_FL076 <- "snRNAseq/200420_A00827_0152_AHT2YJDMXX_Kridel_Robert/Kridel_Robert__LY_FL_076_T1/filtered_feature_bc_matrix"
data_dir_FL227 <- "snRNAseq/200317_A00827_0147_BH75CGDSXY_Kridel_Robert/Kridel_Robert__LY_FL_227_T1/filtered_feature_bc_matrix"

list.files(data_dir_FL062)
list.files(data_dir_FL064)
list.files(data_dir_FL076)
list.files(data_dir_FL227) #"barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"

expression_matrix_FL062 <- Read10X(data.dir = data_dir_FL062)
expression_matrix_FL064 <- Read10X(data.dir = data_dir_FL064)
expression_matrix_FL076 <- Read10X(data.dir = data_dir_FL076)
expression_matrix_FL227 <- Read10X(data.dir = data_dir_FL227)

#2. Setup the Seurat objects and normalize++++++++++++++++++++++++++++++++++++++

exp_matrices = list(expression_matrix_FL062, expression_matrix_FL064, expression_matrix_FL076, expression_matrix_FL227)
names(exp_matrices) = c("FL062", "FL064", "FL076", "FL227")
samps = c("FL062", "FL064", "FL076", "FL227")

pdf(paste(output, "seurat_objects_qc_vln_plots.pdf", sep=""), width=16, height=8)
all_objects = mapply(doSeuratProc, exp_matrices, samps)
dev.off()

# plot variable features with and without labels (before sample integration)
pdf(paste(output, "seurat_top10_genes_per_sample.pdf", sep=""), width=16, height=8)
for(i in 1:4){
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

get_integrated_obj = function(dat, dim, anch_features){

	anchors <- FindIntegrationAnchors(object.list = dat, dims = 1:dim, anchor.features = anch_features)

	#We then pass these anchors to the IntegrateData function, which returns a Seurat object.
	#The returned object will contain a new Assay, which holds an integrated
	#(or 'batch-corrected') expression matrix for all cells, enabling them to be jointly analyzed.

	combined <- IntegrateData(anchorset = anchors, dims = 1:dim)

	#After running IntegrateData, the Seurat object will contain a new Assay with
	#the integrated expression matrix. Note that the original (uncorrected values)
	#are still stored in the object in the other assay, so you can switch back and forth.

	#We can then use this new integrated matrix for downstream analysis and
	#visualization. Here we scale the integrated data, run PCA, and
	#visualize the results with UMAP. The integrated datasets cluster by cell type, instead of by technology.

	# switch to integrated assay. The variable features of this assay are automatically
	# set during IntegrateData

	DefaultAssay(combined) <- "integrated"

	#4. Run clustering++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# Run the standard workflow for visualization and clustering
	combined <- ScaleData(combined, verbose = FALSE)
	combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
	# t-SNE and Clustering
	combined <- RunUMAP(combined, reduction = "pca", dims = 1:dim)
	combined <- FindNeighbors(combined, reduction = "pca", dims = 1:dim)
	combined <- FindClusters(combined, resolution = 0.5)

	head(Idents(combined), 5)

	pdf(paste(output, "seurat_integrated_dim_", dim , anch_features, "_samples_clusters.pdf", sep=""), width=16, height=8)
	p1 <- DimPlot(combined, reduction = "umap", group.by = "sample")
	p2 <- DimPlot(combined, reduction = "umap", label = TRUE)
	p1 + p2
	dev.off()

	saveRDS(combined, file = paste(output, "seurat_integrated_dim_", dim , anch_features,  "_samples_clusters.rds", sep=""))

}

get_integrated_obj(all_objects, 10, 3000)
get_integrated_obj(all_objects, 20, 3000)
get_integrated_obj(all_objects, 30, 3000)

get_integrated_obj(all_objects, 10, 2000)
get_integrated_obj(all_objects, 20, 2000)
get_integrated_obj(all_objects, 30, 2000)
