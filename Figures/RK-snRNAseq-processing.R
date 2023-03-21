#-------------------------------------------------------------------------------
# RK_seurat_analysis_script.R
# Robert Kridel
# December 2022
#-------------------------------------------------------------------------------

setwd("~/github/FLOMICS/")
# R 4.1.0
# Seurat 4.2.0

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors = F)
packages <- c("dplyr", "readr", "ggplot2", "tidyr", "data.table", "plyr",
	            "stringr", "Seurat", "tidyverse", "annotables", "clustree", 
              "cowplot", "patchwork", "DoubletFinder", "ggpubr")
lapply(packages, require, character.only = TRUE)

#-------------------------------------------------------------------------------
#set up environment
#-------------------------------------------------------------------------------

# output directory
output = "/Users/robert/github/FLOMICS/Analysis-Files/Seurat/December2022/"

#snRNA data
input_data = "~/snRNAseq/"

#date
date = Sys.Date()

#-------------------------------------------------------------------------------
# purpose
#-------------------------------------------------------------------------------

# the purpose is obtain clusters from snRNAseq that are representative
# of cell types which we can then use to convolute what cell types are present
# in the bulk RNA-seq samples

#-------------------------------------------------------------------------------
# analysis
#-------------------------------------------------------------------------------

# 1. read in data++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# The input files were aligned to hg19

# Locations
data_dir_FL062 <- paste(input_data, "/", "LY_FL_062_T1/data", sep="")
data_dir_FL064 <- paste(input_data, "/", "LY_FL_064_T1/data", sep="")
data_dir_FL076 <- paste(input_data, "/", "LY_FL_076_T1/data", sep="")

# List files in data directories
list.files(data_dir_FL062)
list.files(data_dir_FL064)
list.files(data_dir_FL076)

# Read in expression matrices
expression_matrix_FL062 <- Read10X(data.dir = data_dir_FL062)
expression_matrix_FL064 <- Read10X(data.dir = data_dir_FL064)
expression_matrix_FL076 <- Read10X(data.dir = data_dir_FL076)

# Create Seurat objects
FL062 <- CreateSeuratObject(counts = expression_matrix_FL062, project = "FL062", min.cells = 10, min.features = 400)
print(FL062)

FL064 <- CreateSeuratObject(counts = expression_matrix_FL064, project = "FL064", min.cells = 10, min.features = 400)
print(FL064)

FL076 <- CreateSeuratObject(counts = expression_matrix_FL076, project = "FL076", min.cells = 10, min.features = 400)
print(FL076)

# Remove expression matrices and free up memory
rm(expression_matrix_FL062, expression_matrix_FL064, expression_matrix_FL076, expression_matrix_FL227)
gc()

# Display count matrix and the metatada
as.data.frame(FL062@assays$RNA@counts[1:10, 1:2])
head(FL062@meta.data, 10)

as.data.frame(FL064@assays$RNA@counts[1:10, 1:2])
head(FL064@meta.data, 10)

as.data.frame(FL076@assays$RNA@counts[1:10, 1:2])
head(FL076@meta.data, 10)

# 2. QC++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# calculate mitochondrial percentage
FL062[["percent.mt"]] <- PercentageFeatureSet(FL062, pattern = "^MT-")
print(summary(FL062[["percent.mt"]]))
FL064[["percent.mt"]] <- PercentageFeatureSet(FL064, pattern = "^MT-")
print(summary(FL064[["percent.mt"]]))
FL076[["percent.mt"]] <- PercentageFeatureSet(FL076, pattern = "^MT-")
print(summary(FL076[["percent.mt"]]))

# calculate ribosomal percentage
FL062[["percent.ribo"]] <- PercentageFeatureSet(FL062, pattern = "^RP[SL]")
print(summary(FL062[["percent.ribo"]]))
FL064[["percent.ribo"]] <- PercentageFeatureSet(FL064, pattern = "^RP[SL]")
print(summary(FL064[["percent.ribo"]]))
FL076[["percent.ribo"]] <- PercentageFeatureSet(FL076, pattern = "^RP[SL]")
print(summary(FL076[["percent.ribo"]]))

p1 <- VlnPlot(FL062, features = c("nFeature_RNA", "percent.mt", "percent.ribo"), ncol = 3)
p2 <- VlnPlot(FL064, features = c("nFeature_RNA", "percent.mt", "percent.ribo"), ncol = 3)
p3 <- VlnPlot(FL076, features = c("nFeature_RNA", "percent.mt", "percent.ribo"), ncol = 3)
patchwork::wrap_plots(p1, p2, p3, ncol = 1)

p1 <- FeatureScatter(FL062, "nCount_RNA", "nFeature_RNA")
p2 <- FeatureScatter(FL064, "nCount_RNA", "nFeature_RNA")
p3 <- FeatureScatter(FL076, "nCount_RNA", "nFeature_RNA")
patchwork::wrap_plots(p1, p2, p3, ncol = 4)

# remove cells with high mitochondrial percentage
FL062 <- subset(FL062, subset = percent.mt < 5)
FL064 <- subset(FL064, subset = percent.mt < 5)
FL076 <- subset(FL076, subset = percent.mt < 5)

# remove cells with high ribo percentage
FL062 <- subset(FL062, subset = percent.ribo < 5)
FL064 <- subset(FL064, subset = percent.ribo < 5)
FL076 <- subset(FL076, subset = percent.ribo < 5)

# Identify genes that are highly expressed
par(mar = c(4, 8, 2, 1))
C <- FL062@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
boxplot(t(as.matrix(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

par(mar = c(4, 8, 2, 1))
C <- FL064@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
boxplot(t(as.matrix(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

par(mar = c(4, 8, 2, 1))
C <- FL076@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
boxplot(t(as.matrix(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

# MALAT1 gene outlier, represents ~3-5% of all reads -> remove
FL062 <- FL062[!grepl("MALAT1", rownames(FL062)), ]
FL064 <- FL064[!grepl("MALAT1", rownames(FL064)), ]
FL076 <- FL076[!grepl("MALAT1", rownames(FL076)), ]

# Sex: all 3 samples from F patients
# Remove X and Y chromosome genes
VlnPlot(FL062, features = c("XIST"))
VlnPlot(FL064, features = c("XIST"))
VlnPlot(FL076, features = c("XIST"))

genes <- as.data.table(grch37)
chrX.gene <- genes %>% filter(chr == "X") %>% .$symbol %>% unique()
chrY.gene <- genes %>% filter(chr == "Y") %>% .$symbol %>% unique()

counts <- GetAssayData(FL062, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(chrX.gene, chrY.gene))),]
FL062 <- subset(FL062, features = rownames(counts))

counts <- GetAssayData(FL064, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(chrX.gene, chrY.gene))),]
FL064 <- subset(FL064, features = rownames(counts))

counts <- GetAssayData(FL076, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(chrX.gene, chrY.gene))),]
FL076 <- subset(FL076, features = rownames(counts))

# Identify doublets
FL062_d <- NormalizeData(FL062)
FL062_d = FindVariableFeatures(FL062_d, verbose = TRUE)
FL062_d = ScaleData(FL062_d, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = TRUE)
FL062_d = RunPCA(FL062_d, verbose = TRUE, npcs = 20)
FL062_d = RunUMAP(FL062_d, dims = 1:10, verbose = TRUE)
nExp <- round(ncol(FL062_d) * 0.04)  # expect 4% doublets
FL062_d <- doubletFinder_v3(FL062_d, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
DF.name = colnames(FL062_d@meta.data)[grepl("DF.classification", colnames(FL062_d@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(FL062_d, group.by = "orig.ident") + NoAxes(),
                   DimPlot(FL062_d, group.by = DF.name) + NoAxes())
VlnPlot(FL062_d, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
FL062_singlets <- FL062_d@meta.data %>% filter(.[[7]]  == "Singlet") %>% row.names()
FL062_doublets <- FL062_d@meta.data %>% filter(.[[7]] == "Doublet") %>% row.names()
FL062 <- FL062[, FL062_singlets]

FL064_d <- NormalizeData(FL064)
FL064_d = FindVariableFeatures(FL064_d, verbose = TRUE)
FL064_d = ScaleData(FL064_d, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = TRUE)
FL064_d = RunPCA(FL064_d, verbose = TRUE, npcs = 20)
FL064_d = RunUMAP(FL064_d, dims = 1:10, verbose = TRUE)
nExp <- round(ncol(FL064_d) * 0.04)  # expect 4% doublets
FL064_d <- doubletFinder_v3(FL064_d, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
DF.name = colnames(FL064_d@meta.data)[grepl("DF.classification", colnames(FL064_d@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(FL064_d, group.by = "orig.ident") + NoAxes(),
                   DimPlot(FL064_d, group.by = DF.name) + NoAxes())
VlnPlot(FL064_d, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
FL064_singlets <- FL064_d@meta.data %>% filter(.[[7]] == "Singlet") %>% row.names()
FL064_doublets <- FL064_d@meta.data %>% filter(.[[7]] == "Doublet") %>% row.names()
FL064 <- FL064[, FL064_singlets]

FL076_d <- NormalizeData(FL076)
FL076_d = FindVariableFeatures(FL076_d, verbose = TRUE)
FL076_d = ScaleData(FL076_d, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = TRUE)
FL076_d = RunPCA(FL076_d, verbose = TRUE, npcs = 20)
FL076_d = RunUMAP(FL076_d, dims = 1:10, verbose = TRUE)
nExp <- round(ncol(FL076_d) * 0.04)  # expect 4% doublets
FL076_d <- doubletFinder_v3(FL076_d, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
DF.name = colnames(FL076_d@meta.data)[grepl("DF.classification", colnames(FL076_d@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(FL076_d, group.by = "orig.ident") + NoAxes(),
                   DimPlot(FL076_d, group.by = DF.name) + NoAxes())
VlnPlot(FL076_d, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
FL076_singlets <- FL076_d@meta.data %>% filter(.[[7]] == "Singlet") %>% row.names()
FL076_doublets <- FL076_d@meta.data %>% filter(.[[7]] == "Doublet") %>% row.names()
FL076 <- FL076[, FL076_singlets]

saveRDS(FL062, paste0(output, "FL062_qc.rds"))
saveRDS(FL064, paste0(output, "FL064_qc.rds"))
saveRDS(FL076, paste0(output, "FL076_qc.rds"))

# 3. Integration ++++++++++++++++++++++++++++++++++++++

# FL062 <- readRDS("Analysis-Files/Seurat/December2022/FL062_qc.rds")
# FL064 <- readRDS("Analysis-Files/Seurat/December2022/FL064_qc.rds")
# FL076 <- readRDS("Analysis-Files/Seurat/December2022/FL076_qc.rds")

# SCTransform: Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().

exp = list(FL062, FL064, FL076) 
names(exp) = c("FL062", "FL064", "FL076")

exp <- lapply(X = exp, FUN = function(x) {
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo"), verbose = TRUE)
})

features <- SelectIntegrationFeatures(object.list = exp, nfeatures = 2000)
exp <- PrepSCTIntegration(object.list = exp, anchor.features = features)
#takes 12min to run:
exp.anchors <- FindIntegrationAnchors(object.list = exp, normalization.method = "SCT", anchor.features = features)
exp.combined.sct <- IntegrateData(anchorset = exp.anchors, normalization.method = "SCT")
exp.combined.sct <- RunPCA(exp.combined.sct, verbose = TRUE)
exp.combined.sct <- RunUMAP(exp.combined.sct, reduction = "pca", dims = 1:20)

# Examine and visualize PCA results a few different ways
print(exp.combined.sct[["pca"]], dims = 1:20, nfeatures = 5)
VizDimLoadings(exp.combined.sct, dims = 1:2, reduction = "pca")
DimHeatmap(exp.combined.sct, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(exp.combined.sct)

pdf(paste0(output, date, "_seurat_umap_3_samples_in_one.pdf"), width = 4, height = 5)
p <- DimPlot(exp.combined.sct, reduction = "umap", label = FALSE, repel = TRUE) + NoLegend()
p
dev.off()

pdf(paste0(output, date, "_seurat_umap_3_samples_separated.pdf"), width = 10, height = 5)
p <- DimPlot(exp.combined.sct, reduction = "umap", label = FALSE, repel = TRUE, split.by = "orig.ident") + NoLegend()
p
dev.off()

# Determine the K-nearest neighbor graph
exp.combined.sct <- FindNeighbors(object = exp.combined.sct, dims = 1:20)

# Determine the clusters for various resolutions
exp.combined.sct <- FindClusters(object = exp.combined.sct, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0))

pdf(paste0(output, date, "_clustree_3_samples.pdf"), width = 10, height = 10)
clustree(exp.combined.sct)
dev.off()

# Assign identity of clusters
Idents(object = exp.combined.sct) <- "integrated_snn_res.0.3"

# Plot the UMAP
pdf(paste0(output, date, "_seurat_umap_clusters_3_samples_0.3_snn_res.pdf"), width = 10, height = 3.5)
p1 <- DimPlot(exp.combined.sct, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
p2 <- DimPlot(exp.combined.sct, reduction = "umap", label = FALSE, repel = TRUE, split.by = "orig.ident")
p1 + p2 + plot_layout(widths = c(1, 3))
dev.off()

saveRDS(exp.combined.sct, paste0(output, "FL_qc_integrated_clusters_3_samples_0.3_snn_res.rds"))
