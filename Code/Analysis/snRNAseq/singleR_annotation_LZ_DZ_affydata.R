### SingleR using normal hematopoietic references for annotating blood cancer cell clusters
### Andy Zeng
### Jan 13, 2021

#RUN in R/4.0.0

library(data.table)
library(tidyverse)
library(Seurat)
library(SingleR)
library(ggpubr)

date=Sys.Date()

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/JD_AZ_stem_cells")
mypal = c("#E5DFD9","#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

# Load the expression matrix from normal hematopoiesis
exprmat_CB <- fread("Expression_matrix_LZ_DZ.txt")
exprmat_CB$V1=NULL
exprmat_CB$ids=exprmat_CB$hgnc_symbol
exprmat_CB$hgnc_symbol = NULL
colnames(exprmat_CB)[-1] = sapply(colnames(exprmat_CB)[-1], function(x){unlist(strsplit(x, ".CEL"))})

anno <- tibble('Sample' = colnames(exprmat_CB)[-1]) %>%
  separate(Sample, c("Patient", "CellType", "Sample"), sep='_', remove=FALSE, extra='merge')

# Filter expression matrix to keep the celltypes of interest
genes <- exprmat_CB$ids
exprmat_CB <- exprmat_CB %>%
  select(colnames(exprmat_CB)[-1]) %>%
  data.matrix()
rownames(exprmat_CB) <- genes

#get Seurat object and normalize it
output="/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/April2021/"
SeuratObject=readRDS(paste(output, "pc_genes_only_no_seurat_integrated_dim_20_2000_2021-04-08_samples_clusters.rds", sep=""))

DefaultAssay(SeuratObject) <- "RNA"
SeuratObject <- NormalizeData(SeuratObject)

cells_b = c(0, 1, 2, 4, 11, 12)
SeuratObject <- subset(SeuratObject, idents = cells_b)

# Run scores at single cell level. Using the raw count matrix
# I use both raw counts and normalized counts and see which results make sense,
#sometimes one signature will dominate depending on normalization
pred_CBheme_singleR <- SingleR(test = GetAssayData(SeuratObject, assay="RNA"),
    ref = exprmat_CB, labels = anno$CellType, method='single')
pred_CBheme_singleR

# Add SingleR Labels to Seurat object
SeuratObject$SingleR.label <- pred_CBheme_singleR$labels
#pdf(paste(date, "_", "SingleR_normalized_seurat_clusters_labels_from_LZ_DZ_genes.pdf", sep=""))

pdf("Figure6C_umap.pdf", width=5, height=5)
DimPlot(SeuratObject, group.by = 'SingleR.label', label = FALSE, cols=c("darkgoldenrod1", "cyan4")) + #+ NoLegend()
theme(axis.line = element_line(colour = 'black', size = 1),
text = element_text(size = 20), axis.text = element_text(size = 20))
dev.off()

# Add individual SingleR scores to Seurat metadata
  # This way you can inspect your own clusters for enrichment of each cell type
  # This is probably more important than the labels
rownames(pred_CBheme_singleR$scores) <- rownames(pred_CBheme_singleR)
colnames(pred_CBheme_singleR$scores) <- paste0(colnames(pred_CBheme_singleR$scores), '_score_SingleR')
SeuratObject <- AddMetaData(SeuratObject, as.data.frame(pred_CBheme_singleR$scores))

#get summary of singleR clusters versus seurat clusters
cell_types = as.data.table((table(SeuratObject@meta.data$SingleR.label, SeuratObject@meta.data$seurat_clusters)))
colnames(cell_types)[1:3]=c("singleR_label", "seurat_cluster", "singleR_cells_in_seurat_cluster")
cell_types = as.data.table(filter(cell_types, singleR_cells_in_seurat_cluster >0))
cell_types = cell_types[order(-singleR_cells_in_seurat_cluster)]
seurat=as.data.table(table(SeuratObject@meta.data$seurat_clusters))
colnames(seurat)=c("seurat_cluster", "tota_seurat_cluster_cells")
cell_types=merge(seurat, cell_types, by="seurat_cluster")
cell_types
cell_types$seurat_explained_by_singleR = cell_types$singleR_cells_in_seurat_cluster/cell_types$tota_seurat_cluster_cells
cell_types
write.csv(cell_types, file="SingleR_normalized_seurat_clusters_labels_from_LZ_DZ_genes_summary_data.csv", quote=F, row.names=F)

#barplot summary of cell_types
cell_types$seurat_cluster = as.numeric(cell_types$seurat_cluster)
cell_types = cell_types[order(seurat_cluster)]
cell_types$seurat_cluster = factor(cell_types$seurat_cluster, levels=unique(cell_types$seurat_cluster))

#pdf(paste(date, "_", "SingleR_normalized_seurat_clusters_labels_from_LZ_DZ_genes_barplot.pdf", sep=""))
pdf("Figure6C_barplot.pdf", width=5, height=5)
ggbarplot(cell_types, x="seurat_cluster", y="seurat_explained_by_singleR",
fill="singleR_label", palette=c("darkgoldenrod1", "cyan4"))+
theme(axis.line = element_line(colour = 'black', size = 1),
       text = element_text(size = 20), axis.text = element_text(size = 20))
dev.off()
