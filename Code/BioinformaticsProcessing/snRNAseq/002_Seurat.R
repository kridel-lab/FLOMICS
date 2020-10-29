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
	"patchwork", "Biobase")

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

#data from nature methods paper for CellAssign
fl_combined_wov = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/41592_2019_529_MOESM4_ESM_HGSC_FL_combined_edited.csv")
fl_cells = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/41592_2019_529_MOESM4_ESM_FL_light_chain.csv")
fl_lightchain = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/41592_2019_529_MOESM4_ESM_FL_celltype.csv")

fl_combined_wov = melt(fl_combined_wov)
fl_combined_wov$data = "fl_combined_wov"

fl_cells = melt(fl_cells)
fl_cells$data = "fl_cells"

fl_lightchain = melt(fl_lightchain)
fl_lightchain$data = "fl_lightchain"

all_cellassign = rbind(fl_lightchain, fl_combined_wov, fl_cells)
colnames(all_cellassign)[1:2] = c("gene", "cell")
all_cellassign = as.data.table(filter(all_cellassign, value==1))

#data from christian steidl paper (manually put into spreadsheet)
steidl = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/FL_cells_types_aoki_cancerdiscovery.csv")

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

#1. Find cluster markers++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(combined) <- "integrated"

combined.markers <- FindAllMarkers(combined, only.pos = TRUE,
	 min.pct = 0.25, logfc.threshold = 0.25, test.use="roc")
combined.markers = as.data.table(combined.markers)
write.csv(combined.markers, file="seurat_roc_cluster_marker_genes_minPCT_25_logFCthresh_25.csv",
quote=F, row.names=F)

#additional markers

pdf(paste(output, "seurat_integrated_samples_clusters_wAdditional_Markers.pdf", sep=""), width=18, height=12)

genes=c("CD3D", "CD4", "CD8A", "NCAM1", "MS4A1", "PTPRC", "BCL6",
"FOXP3", "PDCD1", "CD68", "CD163", "CXCR5", "CD69", "CD45RA", "CD45RO", "CCR7", "IL7R")

VlnPlot(combined, features = genes)

#overlay on UMAP clusters
FeaturePlot(combined, features = genes,
cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))

dev.off()


#visualize markers across clusters

pdf(paste(output, "seurat_integrated_samples_clusters_wMarkers.pdf", sep=""), width=18, height=12)

VlnPlot(combined, features = c("CD79A", "CR2", "CD3D", "CCL5", "ICOS", "IGKC", "IGLC2", "IGLC3", "CD19"))

#overlay on UMAP clusters
FeaturePlot(combined, features = c("CD79A", "CR2", "CD3D", "CCL5", "ICOS", "IGKC", "IGLC2", "IGLC3", "BCL2", "CD19"),
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
clusters_cells_5 = merge(clusters_genes, all_cellassign, by="gene")
clusters_cells_6 = merge(clusters_genes, steidl, by="gene")

#based on some results above plot expression of those genes across clusters
pdf(paste(output, "seurat_integrated_samples_clusters_wMarkers_B_cells.pdf", sep=""), width=18, height=12)

VlnPlot(combined, features = c("BCL2"))
#overlay on UMAP clusters
FeaturePlot(combined, features = c("IGKC", "IGLC2", "IGLC3", "BCL2", "STAT3", "CD27", "CD19", "IGHM"),
cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))

dev.off()

#diff B cell genes
pdf(paste(output, "seurat_integrated_samples_clusters_wMarkers_cellAssign.pdf", sep=""), width=18, height=12)
#overlay on UMAP clusters
FeaturePlot(combined, features = c("CD3G", "CD8A", "BANK1", "VIM", "LYZ", "ST8SIA1", "ICA1", "IL7R", "TRAC", "CTLA4", "IKZF2"),
cols=c("antiquewhite", "cadetblue3", "chartreuse3", "red"))

dev.off()

#I think using imsig_cells seems like it makes the most sense

#3. Label clusters based on marker genes++++++++++++++++++++++++++++++++++++++++

for(i in 0:12){
	print(filter(clusters_cells_5, cluster==i))
}

#combine cluster 2 and 8

current.cluster.ids = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
new.cluster.ids = c("Malignant B cells", "Malignant B cells", "CD4 reg or Tfh cells_2",
					"CD8 cytotoxic T cells cells_3",
                    "Malignant B cells", "CD4 FTH cells_5",
										"stromal cells_6", "CD4 T cells_7", "CD4 reg or Tfh cells_2",
									"Macrophages_9", "Predicted Norm B cells_10", "potential_Monocytes_11", "Endothelial cells_12")
names(x = new.cluster.ids) <- levels(x = combined)
combined <- RenameIdents(object = combined, new.cluster.ids)

pdf(paste(output, "seurat_integrated_samples_cell_types_prelim.pdf", sep=""), width=18, height=12)
DimPlot(object = combined, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

#4. save seurat object for bisque analysis
saveRDS(combined, file="combined_processed_snRNAseq_FL_seurat_object.rds")
