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

all_objects = mapply(doSeuratProc, exp_matrices, samps)

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

anchors <- FindIntegrationAnchors(object.list = all_objects, dims = 1:30)

#We then pass these anchors to the IntegrateData function, which returns a Seurat object.
#The returned object will contain a new Assay, which holds an integrated
#(or 'batch-corrected') expression matrix for all cells, enabling them to be jointly analyzed.

combined <- IntegrateData(anchorset = anchors, dims = 1:30)

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

combined <- FindNeighbors(combined, dims = 1:10)
combined <- FindClusters(combined, resolution = 0.5)

head(Idents(combined), 5)

combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:10)
combined <- RunUMAP(combined, dims = 1:10)

pdf(paste(output, "seurat_integrated_samples_clusters.pdf", sep=""), width=16, height=8)
p1 <- DimPlot(combined, reduction = "umap", group.by = "sample")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE)
p1 + p2
dev.off()

saveRDS(combined, file = "combined_processed_snRNAseq_FL.rds")




















# t-SNE and Clustering
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

DefaultAssay(combined) <- "RNA"
markers <- FindConservedMarkers(combined, ident.1 = 1, grouping.var = "sample", verbose = FALSE)
head(markers)

pdf(paste(output, "seurat_integrated_sample_feature_plot.pdf", sep=""), width=16, height=12)
FeaturePlot(combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "IL4R",
                                          "CCL2", "AC023590.1"), min.cutoff = "q9")
dev.off()

object <- merge(x = object_FL062, y = c(object_FL064, object_FL076, object_FL227),
                add.cell.ids = c("FL062", "FL064", "FL076", "FL227"), project = "FL")

# Standard pre-processing workflow
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
head(object@meta.data, 5)
#VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))
object <- subset(object,
                 subset = nFeature_RNA > 200 &
                          nFeature_RNA < 2500 &
                          percent.mt < 5)
# Normalizing the data
object <- NormalizeData(object,
                        normalization.method = "LogNormalize",
                        scale.factor = 10000)

# Identification of highly variable features (feature selection)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)

# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(object)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))

# Scaling the data
all.genes <- rownames(object)
object <- ScaleData(object, features = all.genes)

# Use SCTransform instead of NormalizeData to remove technical variation while retaining biological heterogeneity
# object <- SCTransform(object, vars.to.regress = "percent.mt", verbose = FALSE)
# Perform linear dimensional reduction

object <- RunPCA(object, features = VariableFeatures(object = object))
print(object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object, dims = 1:2, reduction = "pca")
DimPlot(object, reduction = "pca")
DimHeatmap(object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the 'dimensionality' of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time

object <- JackStraw(object, num.replicate = 100)
object <- ScoreJackStraw(object, dims = 1:20)
JackStrawPlot(object, dims = 1:15)
ElbowPlot(object)
# Cluster the cells
object <- FindNeighbors(object, dims = 1:10)
object <- FindClusters(object, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(object), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
object <- RunUMAP(object, dims = 1:10)
object <- RunTSNE(object, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(object, reduction = "umap")
DimPlot(object, reduction = "tsne")
# Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 1
cluster1.markers <- FindMarkers(object, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(object, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# find all markers of cluster 3
cluster3.markers <- FindMarkers(object, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
object.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(object, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(object, features = c("BCL2", "PAX5", "MS4A1"))
# you can plot raw counts as well
VlnPlot(object, features = c("BCL2", "PAX5", "MS4A1"), slot = "counts", log = TRUE)
FeaturePlot(object, features = c("MS4A1", "BCL2", "CD3E", "CD4",
                                 "CD8", "CXCR5", "BTLA",
                                 "B2M", "TMSB4X", "CD74",
                                 "CD163", "CD68", "CR2",
                                 "IKZF2", "MIR4435-2HG", "IL7R"))
# cluster5.markers <- FindMarkers(object, ident.1 = 5, ident.2 = c(1, 3), min.pct = 0.25)
top10 <- object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object, features = top10$gene) + NoLegend()
# Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(object)
object <- RenameIdents(object, new.cluster.ids)
DimPlot(object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(object, file = "../output/pbmc3k_final.rds")
