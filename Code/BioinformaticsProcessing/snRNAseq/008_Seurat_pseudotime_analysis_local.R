library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(data.table)
library(plyr)
library(cowplot)
library(dplyr)
library(dittoSeq)

set.seed(1234)

setwd("~/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April82021")

date = Sys.Date()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#DATA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#seurat object from our single cell data
r=readRDS("pc_genes_only_no_seurat_integrated_dim_20_2000_2021-04-08_samples_clusters.rds")

#define B cells
cells_b = c(0, 1, 2, 4, 11)
DefaultAssay(r) <- "RNA"
r <- NormalizeData(r)

mainDir="/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April82021"
subDir=paste("Pseudotime_analysis", date, sep="_")

dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir, subDir))

cols = c(dittoColors(), dittoColors(1)[seq_len(7)])

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, cluster_name){

  cell_ids <- which(colData(cds)[, "seurat_clusters"] == cluster_name)

  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}

get_pseudotime_estimates = function(root_cells){

    cell_start = root_cells
    r.cds <- order_cells(r.cds, root_pr_nodes=get_earliest_principal_node(r.cds, cell_start))

   # plot trajectories colored by pseudotime
   p = plot_cells(cds = r.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE, label_branch_points=FALSE,  graph_label_size=3) + ggtitle(paste(cell_start, "as root"))

  #add final pseudotime to main seurat object
  r <- AddMetaData(
    object = r,
    metadata = r.cds@principal_graph_aux@listData$UMAP$pseudotime,
    col.name = "Pseudotime"
  )

  f = FeaturePlot(r, "Pseudotime", pt.size = 0.1, label = TRUE) & scale_color_viridis_c(option="plasma")
  p = plot_grid(p, f, labels = c('A', 'B'), label_size = 12)
  print(p)

  #get median pseudotime for each cluster when this one is set as root
  med_pseuds = as.data.table(r[[c("seurat_clusters", "Pseudotime")]]) %>%
    filter(seurat_clusters %in% cells_b) %>% group_by(seurat_clusters) %>%
    dplyr::summarize(median =median(Pseudotime))
  med_pseuds = as.data.table(med_pseuds)
  med_pseuds$root = root_cells
  return(med_pseuds)
}

get_pseudotime_umap_figure = function(root_cells){
  
  cell_start = root_cells
  r.cds <- order_cells(r.cds, root_pr_nodes=get_earliest_principal_node(r.cds, cell_start))
  
  #add final pseudotime to main seurat object
  r <- AddMetaData(
    object = r,
    metadata = r.cds@principal_graph_aux@listData$UMAP$pseudotime,
    col.name = "Pseudotime"
  )
  
  f = FeaturePlot(r, "Pseudotime", pt.size = 0.1, label = TRUE) + scale_color_viridis_c(option="plasma")+
    theme(axis.line = element_line(colour = 'black', size = 1),
          text = element_text(size = 20), axis.text = element_text(size = 20))
  
  print(f)
  
}

get_pseudotime_vs_cluster = function(root_cells){
  
  cell_start = root_cells
  r.cds <- order_cells(r.cds, root_pr_nodes=get_earliest_principal_node(r.cds, cell_start))
  
  #add final pseudotime to main seurat object
  r <- AddMetaData(
    object = r,
    metadata = r.cds@principal_graph_aux@listData$UMAP$pseudotime,
    col.name = "Pseudotime"
  )
  
  #get pseudotime for each cluster when this one is set as root and gene expression for gene of interest
  med_pseuds = as.data.table(r[[c("seurat_clusters", "Pseudotime")]])
  med_pseuds$seurat_clusters = as.character(med_pseuds$seurat_clusters)
  med_pseuds = filter(med_pseuds, seurat_clusters %in% cells_b, !(Pseudotime == "Inf")) 
  medians = as.data.table(med_pseuds %>% dplyr::group_by(seurat_clusters) %>% dplyr::summarize(median = median(Pseudotime)))
  medians = medians[order(median)]
  
  f = ggboxplot(med_pseuds, x="seurat_clusters", y="Pseudotime", order=medians$seurat_clusters, color="seurat_clusters",
                palette = cols[as.numeric(medians$seurat_clusters)+1], legend="none")+
    theme(axis.line = element_line(colour = 'black', size = 1),
          text = element_text(size = 20), axis.text = element_text(size = 20)) + xlab("Cell Population")
  
  print(f)
  
}

get_gene_vs_pseudotime = function(gene, root_cells){
  
  #test
  #gene="BCL2"
  
  cell_start = root_cells
  r.cds <- order_cells(r.cds, root_pr_nodes=get_earliest_principal_node(r.cds, cell_start))
  
  # plot trajectories colored by pseudotime
  p = plot_cells(cds = r.cds,
                 color_cells_by = "pseudotime",
                 show_trajectory_graph = TRUE, label_branch_points=FALSE,  graph_label_size=3) + ggtitle(paste(cell_start, "as root"))
  
  #add final pseudotime to main seurat object
  r <- AddMetaData(
    object = r,
    metadata = r.cds@principal_graph_aux@listData$UMAP$pseudotime,
    col.name = "Pseudotime"
  )
  
  #gene expression 
  rna_exp = as.numeric(r$RNA[which(rownames(r$RNA) == gene)])
  
  #get pseudotime for each cluster when this one is set as root and gene expression for gene of interest
  med_pseuds = as.data.table(r[[c("seurat_clusters", "Pseudotime")]])
  med_pseuds$gene_exp = rna_exp
  med_pseuds$seurat_clusters = as.character(med_pseuds$seurat_clusters)
  med_pseuds = filter(med_pseuds, seurat_clusters %in% cells_b, !(Pseudotime == "Inf")) 
  plot = ggscatter(med_pseuds, x="Pseudotime", y="gene_exp", color="seurat_clusters")+geom_smooth(span = 0.3)+
    ggtitle(paste(gene, "vs Pseudotime"))+ylab("Gene Expression")
  print(plot)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ANALYSIS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#convert to CellDataSet object
r.cds <- as.cell_data_set(r)
#create partitions for pseudotime analysis
r.cds <- cluster_cells(cds = r.cds, reduction_method = "UMAP")
#learn graph
r.cds <- learn_graph(r.cds, use_partition = TRUE)
#check partitions created
#dev.off()

pdf("Estimated_partitions_by_Monocle_for_pseudotime.pdf")
plot_cells(r.cds, color_cells_by = "partition")
dev.off()

pdf("Bcells_pseudotime_estimation_all_possible_cells_as_roots.pdf", width=10, height=8)
all_res = as.data.table(ldply(llply(cells_b, get_pseudotime_estimates)))
dev.off()

#summary
ggline(all_res, "root", "median", color="seurat_clusters") + ylab("Median Pseudotime") +xlab("Cluster set as root") + theme_bw()
ggsave("Median_pseudotimes_vs_cell_type_diff_roots.pdf")

#get pseudotime vs gene expression values
pdf("Bcells_pseudotime_root_is_11_vs_gene_expression.pdf", width=8, height=6)
get_gene_vs_pseudotime("BACH2", 11)
get_gene_vs_pseudotime("CXCR4", 11)
get_gene_vs_pseudotime("CD83", 11)
get_gene_vs_pseudotime("LMO2", 11)
get_gene_vs_pseudotime("PAX5", 11)
get_gene_vs_pseudotime("HLA-DRA", 11)
get_gene_vs_pseudotime("BCL6", 11)
get_gene_vs_pseudotime("BIRC3", 11)
get_gene_vs_pseudotime("GS1-410F4.2", 11)
get_gene_vs_pseudotime("VAV3", 11)
get_gene_vs_pseudotime("ADRBK2", 11)
dev.off()

pdf("Bcells_pseudotime_root_is_4_vs_gene_expression.pdf", width=8, height=6)
get_gene_vs_pseudotime("BACH2", 4)
get_gene_vs_pseudotime("CXCR4", 4)
get_gene_vs_pseudotime("CD83", 4)
get_gene_vs_pseudotime("LMO2", 4)
get_gene_vs_pseudotime("PAX5", 4)
get_gene_vs_pseudotime("HLA-DRA", 4)
get_gene_vs_pseudotime("BCL6", 4)
get_gene_vs_pseudotime("BIRC3", 4)
get_gene_vs_pseudotime("GS1-410F4.2", 4)
get_gene_vs_pseudotime("VAV3", 4)
get_gene_vs_pseudotime("ADRBK2", 4)
dev.off()

#save plots for manuscript figure 6
mainDir="/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April82021"
subDir=paste("Figures_for_manuscript", date, sep="_")
dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir, subDir))

#UMAP with pseudotime on it
pdf("Figure6_pseudotime_umap.pdf", width=5, height=5)
get_pseudotime_umap_figure(11) #11 set as root
dev.off()

#Clusters versus pseudotime
pdf("Figure6_pseudotime_boxplots_vs_clusters.pdf", width=5, height=5)
get_pseudotime_vs_cluster(11)
dev.off()
