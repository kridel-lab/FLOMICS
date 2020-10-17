#-------------------------------------------------------------------------------
#001_Seurat_tutorial.R
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
library(ellipsis)
library(tidyverse)
library(splitstackshape)
packages <- c("readr", "data.table", "plyr",
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
#data
#-------------------------------------------------------------------------------

#processed data obtained for 4 samples
combined = readRDS("combined_processed_snRNAseq_FL.rds")

#gene marker versus cell type data
cells = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/PanglaoDB_markers_27_Mar_2020.tsv")
cells = as.data.table(filter(cells, species %in% c("Mm Hs", "Hs"), organ=="Immune system"))
colnames(cells)[c(2,3,5, 6, 9)] = c("gene", "cell", "ubiquitousness_index", "product", "germlayer")
cells = as.data.table(cells %>% select(gene, cell, ubiquitousness_index, product, sensitivity_human, specificity_human) %>%
	filter(sensitivity_human > 0.1, specificity_human > 0.1))

imsig_cells = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/imsig_gene_cells.csv")

cell_markers = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Human_cell_markers_HRBMUS_EDU.txt")
cell_markers = as.data.table(filter(cell_markers, tissueType=="Blood"))
cell_markers = cSplit(cell_markers, "cellMarker", ",", direction = "long")
cell_markers = cell_markers %>% select(tissueType, cellType, cellName, cellMarker)
colnames(cell_markers)[4] = "gene"

cell_markers_single_cell = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Single_cell_markers_HRBMU_EDU.txt")
cell_markers_single_cell = as.data.table(filter(cell_markers_single_cell, tissueType=="Blood"))
cell_markers_single_cell = cSplit(cell_markers_single_cell, "proteinName", ",", direction = "long")
cell_markers_single_cell = cell_markers_single_cell %>% select(tissueType, cellType, cellName, proteinName)
colnames(cell_markers_single_cell)[4] = "gene"

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

#1. Find cluster markers++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(combined) <- "integrated"

combined.markers <- FindAllMarkers(combined, only.pos = TRUE,
	 min.pct = 0.25, logfc.threshold = 0.25, test.use="roc")

#visualize markers across clusters

pdf(paste(output, "seurat_integrated_samples_clusters_wMarkers.pdf", sep=""), width=18, height=12)

VlnPlot(combined, features = c("CD79A", "CD3D", "CCL5", "ICOS", "IGKC", "IGLC2", "IGLC3"))

#overlay on UMAP clusters
FeaturePlot(combined, features = c("CD79A", "CD3D", "CCL5", "ICOS", "IGKC", "IGLC2", "IGLC3", "BCL2", "BCL6"),
cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))

top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_diff)
DoHeatmap(combined, features = top10$gene)
dev.off()

#2. Link markers to cell types++++++++++++++++++++++++++++++++++++++++++++++++++

clusters_genes = as.data.table(table(combined.markers$gene, combined.markers$cluster)) #569 unique genes left
colnames(clusters_genes)=c("gene", "cluster", "N")
clusters_genes = as.data.table(filter(clusters_genes, N >0))
clusters_cells_1 = merge(clusters_genes, cell_markers, by="gene")
clusters_cells_2 = merge(clusters_genes, cells, by="gene")
clusters_cells_3 = merge(clusters_genes, imsig_cells, by="gene")
clusters_cells_4 = merge(clusters_genes, cell_markers_single_cell, by="gene")

#I think using imsig_cells seems like it makes the most sense

#3. Label clusters based on marker genes++++++++++++++++++++++++++++++++++++++++

for(i in 0:12){
	print(filter(clusters_cells_3, cluster==i))
}

current.cluster.ids = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
new.cluster.ids = c("B cells_0", "cluster1", "T cells_2", "Tcells_3",
                    "B cells_4", "T cells_5", "Macrophages_6", "T cells_7", "cluster_8",
									"Macrophages_9", "B cells_10", "Monocytes_11", "Monocytes_12")
names(x = new.cluster.ids) <- levels(x = combined)
combined <- RenameIdents(object = combined, new.cluster.ids)

pdf(paste(output, "seurat_integrated_samples_cell_types_prelim.pdf", sep=""), width=18, height=12)
DimPlot(object = combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

#4. save seurat object for bisque analysis
saveRDS(combined, file="combined_processed_snRNAseq_FL_seurat_object.rds")
