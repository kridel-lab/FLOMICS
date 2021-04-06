library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(ggpubr)

set.seed(1234)

setwd("~/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April2021")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#DATA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#seurat object from our single cell data
r=readRDS("pc_genes_only_no_seurat_integrated_dim_20_2000_2021-04-01_samples_clusters.rds")

#confirm UMAP plot
DimPlot(r, label = TRUE)
dev.off()

#define B cells
cells_b = c(12, 2, 0, 13, 6, 16, 1, 10)
DefaultAssay(r) <- "RNA"
r <- NormalizeData(r)

mainDir="/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April2021"
subDir=paste("Pseudotime_analysis", date, sep="_")

dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir, subDir))

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

get_gene_vs_pseudotime = function(gene){

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

pdf("Estimated_partitions_by_Monocle_for_pseudotime.pdf")
plot_cells(r.cds, color_cells_by = "partition")
dev.off()

pdf("Bcells_pseudotime_estimation_all_possible_cells_as_roots.pdf", width=10, height=8)
all_res = as.data.table(ldply(llply(cells_b, get_pseudotime_estimates)))
dev.off()

#summary
ggline(all_res, "root", "median", color="seurat_clusters") + ylab("Median Pseudotime") +xlab("Cluster set as root") + theme_bw()
ggsave("Median_pseudotimes_vs_cell_type_diff_roots.pdf")
