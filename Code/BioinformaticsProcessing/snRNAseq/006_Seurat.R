#-------------------------------------------------------------------------------
#006_Seurat.R
#sarah russell
#Dec 9, 2020
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors=F)
date=Sys.Date()

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
	packages <- c("dplyr", "readr", "ggplot2", "tidyr", "data.table", "plyr",
		"stringr",
	  "Seurat",
	  "cowplot",
		"patchwork","Biobase")

lapply(packages, require, character.only = TRUE)

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

setwd(output)
#2. Setup the Seurat objects and normalize++++++++++++++++++++++++++++++++++++++
exp_matrices = list(expression_matrix_FL062, expression_matrix_FL064, expression_matrix_FL076, expression_matrix_FL227)
names(exp_matrices) = c("FL062", "FL064", "FL076", "FL227")
samps = c("FL062", "FL064", "FL076", "FL227")

pdf(paste(output, "seurat_objects_qc_vln_plots.pdf", sep=""), width=16, height=8)
all_objects = mapply(doSeuratProc, exp_matrices, samps)
dev.off()

## A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#3. Set up anchors++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#dat = all_objects
#dim = 10
#anch_features = 3000

get_integrated_obj = function(dat, dim, anch_features){

	set.seed(100)
	print(dim)
	print(anch_features)

	anchors <- FindIntegrationAnchors(object.list = all_objects, dims = 1:dim, anchor.features = anch_features)

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

	# Assign Cell-Cycle Scores based on its expression of G2/M and S phase markers
	combined <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

	# Visualize the distribution of cell cycle markers across
	pdf(paste(date,"cell_cycle_markers.pdf", sep="_"), width=16, height=8)
	r=RidgePlot(combined, features = c( "TOP2A", "MKI67"), ncol = 1)
	print(r)

	#these two not found in integrated assay, will take from "RNA" assay
	r2=RidgePlot(combined, features = c( "PCNA", "MCM6"), ncol = 1)
	print(r2)

	#Running a PCA on cell cycle genes
	combined <- RunPCA(combined, features = c(s.genes, g2m.genes))
	d=DimPlot(combined)
	print(d)

	# Regress out cell cycle scores during data scaling
	combined <- ScaleData(combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(combined))

	# Now, a PCA on the variable genes no longer returns components associated with cell cycle
	combined <- RunPCA(combined, features = c(s.genes, g2m.genes))
	d2=DimPlot(combined)
	print(d2)

	#assuming I return back to these steps?
	combined <- RunUMAP(combined, reduction = "pca", dims = 1:dim)
	combined <- FindNeighbors(combined, reduction = "pca", dims = 1:dim)
	combined <- FindClusters(combined, resolution = 0.5)


	p1 <- DimPlot(combined, reduction = "umap", group.by = "sample")
	p2 <- DimPlot(combined, reduction = "umap", label = TRUE)
	print(p1 + p2)
	dev.off()

	saveRDS(combined, file = "seurat_integrated_cycle_regression_samples_clusters.rds")
	print("finished this analysis")


}

#get_integrated_obj(all_objects, 10, 3000)
#get_integrated_obj(all_objects, 20, 3000)
#get_integrated_obj(all_objects, 30, 3000)

get_integrated_obj(all_objects, 10, 2000)
