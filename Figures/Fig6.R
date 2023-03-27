#-------------------------------------------------------------------------------
# This script reads in the processed snRNAseq data and produces Fig 6 panels
#-------------------------------------------------------------------------------

setwd("~/your working directory/FLOMICS/")
# R 4.1.0
# Seurat 4.2.0

#--
# Load packages and set up environment
#--

options(stringsAsFactors = FALSE)
packages <- c("dplyr", "ggplot2", "Seurat", "patchwork", "annotables", "IOBR",
 "BisqueRNA")
lapply(packages, require, character.only = TRUE)

#snRNA data
input_data <- "~/snRNAseq/"

#date
date <- Sys.Date()

#--
# Read Seurat objet in
#--
#exp.combined.sct
exp_combined_sct <- readRDS(
  "FL_qc_integrated_clusters_3_samples_0.3_snn_res.rds")

#--
# Generate Fig 6A
#--

pdf(paste0("img/",  date, " _seurat_umap_clusters_3_samples_0.3_snn_res.pdf"),
 width = 10, height = 3.5)
p1 <- DimPlot(exp.combined.sct,
 reduction = "umap", label = TRUE, repel = TRUE) +
 NoLegend()
p2 <- DimPlot(exp.combined.sct,
 reduction = "umap", label = FALSE, repel = TRUE, split.by = "orig.ident")
p1 + p2 + plot_layout(widths = c(1, 3))
dev.off()

#--
# Generate Fig 6B
#--

DefaultAssay(exp_combined_sct) <- "RNA"
exp_combined_sct <- NormalizeData(exp_combined_sct)
b_cell_exp_combined_sct <- subset(exp_combined_sct, idents = c(0, 1, 8, 11))

# Read in Attaf GC and memory genes
attaf_gc_up <- read.csv("Attaf_GC_up.csv") %>% filter(avg_logFC > 1) %>% .$Gene
attaf_mem_up <- read.csv("Attaf_Mem_up.csv") %>%
 filter(avg_logFC < -1) %>%
  .$Gene

b_cell_exp_combined_sct <- AddModuleScore(b_cell_exp_combined_sct,
 features = list(attaf_gc_up), name = "Attaf_GC_up")
b_cell_exp_combined_sct <- AddModuleScore(b_cell_exp_combined_sct,
 features = list(attaf_mem_up), name = "Attaf_Mem_up")

# Plot scores
pdf(paste0("img/",  date, " _Attaf scores.pdf"), width = 3, height = 4)
p1 <- FeaturePlot(b_cell_exp_combined_sct, features = "Attaf_GC_up1",
 label = TRUE, repel = TRUE,
   cols = c("grey", "thistle1", "steelblue", "red")) +
xlim(-6, 5) +
ylim(-10, 3) +
  ggtitle("Expression of GC-like genes") +
  theme(plot.title = element_text(size = 7),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.text = element_text(size = 7))
p2 <- FeaturePlot(B.cell.exp.combined.sct, features = "Attaf_Mem_up1",
 label = TRUE, repel = TRUE,
  cols = c("grey", "thistle1", "steelblue", "red")) +
  xlim(-6, 5) +
  ylim(-10, 3) +
  ggtitle("Expression of Mem-like genes") +
  theme(plot.title = element_text(size = 7),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.text = element_text(size = 7))
p1 + p2 + plot_layout(ncol = 1)
dev.off()

#--
# Generate Fig 6C
#--

# Genes defining given populations
proliferation <- c("CENPE", "CENPF", "EZH2", "POLQ", "TOP2A", "MKI67")
pop1_genes <- c("BCL6", "KLHL6",  "RGS13", "PIP4K2A")
pop0_genes <- c("BANK1", "CXCR4", "FCMR", "FOXP1", "CCND3", "TXNIP")
tfh_educated <- c("CD83", "IL4I1", "IL21R", "MIR155HG", "TNFAIP3", "TRAF1")

marker_genes <- c(proliferation, pop1_genes, pop0_genes, tfh_educated)

levels(b_cell_exp_combined_sct) <- c("8", "1", "0", "11")

pdf(paste0("img/",  date, " dotplot_3_samples_0.3_Attaf_dotplot.pdf"),
 width = 3.5, height = 3.75)
DotPlot(b_cell_exp_combined_sct, features = rev(marker_genes),
        cols = c("#4575b4", "#d73027"), cluster.idents = FALSE, dot.scale = 3) +
  RotatedAxis() +
  coord_flip() +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, face = "italic"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.margin = margin(-10, -10, -10, -10)) +
  ylab("Population")
dev.off()

#--
# Prepare eset for Bisque
#--
seurat_object <- exp_combined_sct
delimiter <- "_"
position <- 2

get_cell_names <- function(obj) base::colnames(obj)
get_ident <- function(obj) Seurat::Idents(object = obj)
get_raw_data <- function(obj) Seurat::GetAssayData(object = obj[["RNA"]], slot = "counts")

individual_ids <- base::sapply(base::strsplit(get.cell.names(seurat_object),
 delimiter), `[[`, position)

base::names(individual_ids) <- get.cell.names(seurat_object)
individual_ids <- base::factor(individual_ids)
n_individuals <- base::length(base::levels(individual_ids))
base::message(base::sprintf("Split sample names by \"%s\"", delimiter),
              base::sprintf(" and checked position %i.", position),
              base::sprintf(" Found %i individuals.", n.individuals))
base::message(base::sprintf("Example: \"%s\" corresponds to individual \"%s\".",
                            get.cell.names(seurat_object)[1],
                             individual_ids[1]))
sample_ids <- base::names(get.ident(seurat_object))
sc_pheno <- base::data.frame(check.names = FALSE, check.rows = FALSE,
                             stringsAsFactors = FALSE,
                             row.names = sample_ids,
                             SubjectName = individual_ids,
                             cellType = get.ident(seurat_object))
sc_meta <- base::data.frame(labelDescription = base::c("SubjectName",
 "cellType"), row.names = base::c("SubjectName", "cellType"))
sc_pdata <- methods::new("AnnotatedDataFrame", data = sc_pheno,
 varMetadata = sc_meta)

sc_data <- base::as.matrix(get.raw.data(seurat_object)[, sample.ids,
 drop = FALSE])
sc_eset <- Biobase::ExpressionSet(assayData = sc.data, phenoData = sc.pdata)

#--
# Run Bisque
#--

# Read in genetic clusters
clusters <- read.csv("GMM_Cluster_Labels_flexmix_clusters.csv") %>%
  select(SAMPLE_ID, ClusterAIC, cohort)

# Read in RNAseq QC data
passing_rnaseq_qc <- read.table("metadata_passed_290.txt", sep = "\t",
 header = TRUE)

# Prepare expression matrix
# bulk gene expression matrix is not normalized
## as the analysis is based on proportion of given populations/sample
bulk_matrix <- read.table("2022-10-25_df_counts_adj.txt", sep = " ",
 header = TRUE)
bulk_matrix <- bulk_matrix[, passing_rnaseq_qc$id]

bulk_matrix$ensgene <- row.names(bulk_matrix)
bulk_matrix <- bulk_matrix %>%
 dplyr::inner_join(grch37[, c("ensgene", "symbol", "chr")],
  by = c("ensgene")) %>%
  filter(!chr %in% c("X", "Y"))
bulk_matrix <- remove_duplicate_genes(eset = bulk_matrix,
 column_of_symbol = "symbol", method = "mean")
summary(duplicated(rownames(bulk_matrix)))
bulk_matrix <- data.matrix(bulk_matrix)
bulk_matrix <- bulk_matrix[, !(colnames(bulk_matrix) %in% c("ensgene", "chr"))]
z <- which(colnames(bulk_matrix) %in%
 c("LY_FL_062_T1", "LY_FL_064_T1", "LY_FL_076_T1"))
colnames(bulk_matrix)[z] <- c(1, 2, 3)

bulk_eset <- Biobase::ExpressionSet(assayData = bulk_matrix)

# Run Bisque
res <- ReferenceBasedDecomposition(bulk_eset, sc_eset, markers = NULL,
 use.overlap = TRUE)
ref_based_estimates <- res$bulk_props

my_comparisons_p0_1 <- list(
  c("C5", "C1"), c("C5", "C2"), c("C5", "C3"), c("C5", "C4"))

p0 <- ref_based_estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(0)) %>%
  mutate(Var1 = paste0("Population ", Var1)) %>%
  select(Population = Var1, SAMPLE_ID = Var2, value) %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  mutate(Cluster = factor(
    ClusterAIC, levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  ggboxplot("Cluster", "value", color = "Cluster", ncol = 5, lwd = 0.3,
   add = "jitter", add.params = list(size = 0.1, jitter = 0.2)) +
  theme_bw() +
  stat_compare_means(size = 2, label.x.npc = 0.15, label.y = 0.97) +
  stat_compare_means(comparisons = my_comparisons_p0_1, method = "wilcox.test",
   size = 2, label.y = c(0.6, 0.68, 0.76, 0.84)) +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7)) +
  ylab("Estimatd proportion") +
  ylim(-0.001, 1.0)

p1 <- ref.based.estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(1)) %>%
  mutate(Var1 = paste0("Population ", Var1)) %>%
  select(Population = Var1, SAMPLE_ID = Var2, value) %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  mutate(Cluster = factor(
    ClusterAIC, levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  ggboxplot("Cluster", "value", color = "Cluster", ncol = 5, lwd = 0.3,
   add = "jitter", add.params = list(size = 0.1, jitter = 0.2)) +
  theme_bw() +
  stat_compare_means(size = 2, label.x.npc = 0.15, label.y = 0.58) +
  stat_compare_means(comparisons = my_comparisons_p0_1, method = "wilcox.test",
   size = 2, label.y = c(0.32, 0.37, 0.42, 0.47)) +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7)) +
  ylab("Estimatd proportion") +
  ylim(-0.001, 0.6)

my_comparisons_p8 <- list(
  c("C3", "C1"), c("C3", "C2"), c("C3", "C4"), c("C3", "C5"))

p8 <- ref.based.estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(8)) %>%
  mutate(Var1 = paste0("Population ", Var1)) %>%
  select(Population = Var1, SAMPLE_ID = Var2, value) %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  mutate(Cluster = factor(
    ClusterAIC, levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  ggboxplot("Cluster", "value", color = "Cluster", ncol = 5, lwd = 0.3,
   add = "jitter", add.params = list(size = 0.1, jitter = 0.2)) +
  theme_bw() +
  stat_compare_means(size = 2, label.x.npc = 0.15, label.y = 0.58) +
  stat_compare_means(comparisons = my_comparisons_p8, method = "wilcox.test",
   size = 2, label.y = c(0.32, 0.37, 0.42, 0.47)) +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7)) +
  ylab("Estimatd proportion") +
  ylim(-0.001, 0.6)

p11 <- ref.based.estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(11)) %>%
  mutate(Var1 = paste0("Population ", Var1)) %>%
  select(Population = Var1, SAMPLE_ID = Var2, value) %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  mutate(Cluster = factor(
    ClusterAIC, levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  ggboxplot("Cluster", "value", color = "Cluster", ncol = 5, lwd = 0.3,
   add = "jitter", add.params = list(size = 0.1, jitter = 0.2)) +
  theme_bw() +
  stat_compare_means(size = 2, label.x.npc = 0.15, label.y = 0.58) +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7)) +
  ylab("Estimatd proportion") +
  ylim(-0.001, 0.6)

g <- gridExtra::arrangeGrob(p0, p1, p8, p11, nrow = 1)
ggsave(file = paste0("img/",
  date, " _bisque_deconvolution_not_norm_full_matrix_adj_B.pdf"),
   g, width = 16, height = 5, units = "cm")
