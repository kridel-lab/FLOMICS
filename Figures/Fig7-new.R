#-------------------------------------------------------------------------------
# This script reads in the processed snRNAseq data and produces Fig 7 panels
# Author: Robert Kridel
#-------------------------------------------------------------------------------

setwd("~/your working directory/FLOMICS/")
# R 4.1.0
# Seurat 4.2.0

#--
# Load packages and set up environment
#--

options(stringsAsFactors = FALSE)
packages <- c("dplyr", "ggplot2", "Seurat", "patchwork", "annotables", "IOBR",
 "BisqueRNA", "tidyr", "data.table")
lapply(packages, require, character.only = TRUE)

#date
date = Sys.Date()

#--
# Read genetic clusters in
#--

clusters <- read.csv("GMM_Cluster_Labels_flexmix_clusters.csv") %>%
  select(SAMPLE_ID, ClusterAIC) %>%
  mutate(ClusterAIC = factor(ClusterAIC,
   levels = c("CS", "TT", "GM", "Q", "AR")))

#--
# Read Ecotyper results in
#--

ecotyper.abundance <- read.table(file = "Ecotype_Abundance.txt",
 sep = "\t", header = TRUE)
b.cell.states.abundance <- read.table(file = "B.cells_Cell_State_Abundance.txt",
 sep = "\t", header = TRUE)

#--
# Read Seurat object in
#--

exp.combined.sct <- readRDS("FL_qc_integrated_clusters_3_samples_0.3_snn_res.rds")

#--
# Read RNAseq in
#--

passing.RNAseq.QC <- read.table("metadata_all_22_23.txt",sep = "\t", header = TRUE) %>%
  filter(!sample_id %in% c(paste("LY_DLC_00", seq(1:9), sep = ""), "LY_DLC_010")) %>%
  filter(qc_tier2 == "Y")

bulk.matrix <- read.table("RNA_seq_analysis_combined_2022_2023_2023-05-01_df_counts_adj.txt",
 sep = " ", header = TRUE)
bulk.matrix <- bulk.matrix[, passing.RNAseq.QC$sample_id]

#--
# Generate Figures 7A-C
#--

b.cell.states.abundance.stacked.barplot <- b.cell.states.abundance %>%
  pivot_longer(!ID, names_to = "B_cell_state", values_to = "abundance") %>%
  inner_join(clusters, by = c("ID" = "SAMPLE_ID")) %>%
  ggplot(aes(fill = B_cell_state, y = abundance, x = ClusterAIC)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("#EB7D5B", "#FED23F",
   "#B5D33D", "#6CA2EA", "#442288")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.margin = margin(-10, -10, -10, -10)) +
  labs(fill = "B cell states") +
  ylab("Relative mean abundance") +
  xlab("Subtype")

ecotyper.abundance.stacked.barplot <- ecotyper.abundance %>%
  pivot_longer(!ID, names_to = "Ecotype", values_to = "abundance") %>%
  inner_join(clusters, by = c("ID" = "SAMPLE_ID")) %>%
  ggplot(aes(fill = Ecotype, y = abundance, x = ClusterAIC)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("#3E2883", "#76A2E4",
   "#8E549F", "#458833", "#BBD058", "#FFFD61",
    "#F8D25D", "#F3A83B", "#EC5428")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.margin = margin(-10, -10, -10, -10)) +
  labs(fill = "Ecotypes") +
  ylab("Relative mean abundance") +
  xlab("Subtype")

ecotyper.abundance.boxplots <- ecotyper.abundance %>%
  pivot_longer(!ID, names_to = "Ecotype", values_to = "abundance") %>%
  inner_join(clusters, by = c("ID" = "SAMPLE_ID")) %>%
  filter(Ecotype %in% c("LE1", "LE8")) %>%
  ggboxplot("ClusterAIC", "abundance", color = "ClusterAIC",
   facet.by = "Ecotype", ncol = 1,
            add = "jitter", add.params = list(size = 0.05, jitter = 0.2)) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        strip.text.x = element_text(size = 7, hjust = 0, vjust = -1.5)) +
  stat_compare_means(method = "kruskal.test", size = 2,
   label.y.npc = 0.9, label.x.npc = 0.1) +
  ylab("Abundance") + xlab("Subtype")

g <- gridExtra::arrangeGrob(b.cell.states.abundance.stacked.barplot,
                            ecotyper.abundance.stacked.barplot,
                            ecotyper.abundance.boxplots,
                            ncol = 3)
ggsave(file = paste0("img/",  date, " ecotyper.pdf"), g,
 width = 16, height = 7, units = "cm")
    
#--
# Generate Figure 7D
#--

pdf(paste0("img/",  date, " _seurat_umap_clusters_3_samples_0.3_snn_res.pdf"),
 width = 10, height = 3.5)
p1 <- DimPlot(exp.combined.sct, reduction = "umap",
 label = TRUE, repel = TRUE) +
 NoLegend()
p2 <- DimPlot(exp.combined.sct, reduction = "umap",
 label = FALSE, repel = TRUE, split.by = "orig.ident")
p1 + p2 + plot_layout(widths = c(1, 3))
dev.off()

#--
# Run Bisque and generate Figure 7E
#--

# Prepare single nuclei dataset for Bisque

seurat.object <- exp.combined.sct
delimiter <- "_"
position <- 2

get.cell.names <- function(obj) base::colnames(obj)
get.ident <- function(obj) Seurat::Idents(object = obj)
get.raw.data <- function(obj) Seurat::GetAssayData(object = obj[["RNA"]], slot = "counts")

individual.ids <- base::sapply(base::strsplit(get.cell.names(seurat.object), delimiter), `[[`, position)

base::names(individual.ids) <- get.cell.names(seurat.object)
individual.ids <- base::factor(individual.ids)
n.individuals <- base::length(base::levels(individual.ids))
base::message(base::sprintf("Split sample names by \"%s\"", delimiter),
              base::sprintf(" and checked position %i.", position),
              base::sprintf(" Found %i individuals.", n.individuals))
base::message(base::sprintf("Example: \"%s\" corresponds to individual \"%s\".",
                            get.cell.names(seurat.object)[1], individual.ids[1]))
sample.ids <- base::names(get.ident(seurat.object))
sc.pheno <- base::data.frame(check.names = FALSE, check.rows = FALSE,
                             stringsAsFactors = FALSE,
                             row.names = sample.ids,
                             SubjectName = individual.ids,
                             cellType = get.ident(seurat.object))
sc.meta <- base::data.frame(labelDescription = base::c("SubjectName", "cellType"),
                            row.names = base::c("SubjectName", "cellType"))
sc.pdata <- methods::new("AnnotatedDataFrame",
 data = sc.pheno, varMetadata = sc.meta)

sc.data <- base::as.matrix(get.raw.data(seurat.object)[, sample.ids, drop = FALSE])
sc.eset <- Biobase::ExpressionSet(assayData = sc.data, phenoData = sc.pdata)

# Prepare bulk RNAseq for Bisque

bulk.matrix$ensgene <- row.names(bulk.matrix)
bulk.matrix <- bulk.matrix %>%
 dplyr::inner_join(grch37[, c("ensgene", "symbol", "chr")],
  by = c("ensgene")) %>%
 filter(!chr %in% c("X", "Y"))
bulk.matrix <- remove_duplicate_genes(eset = bulk.matrix,
 column_of_symbol = "symbol", method = "mean")
summary(duplicated(rownames(bulk.matrix)))
bulk.matrix <- data.matrix(bulk.matrix)
bulk.matrix <- bulk.matrix[, !(colnames(bulk.matrix) %in% c("ensgene", "chr"))]
z = which(colnames(bulk.matrix) %in% c("LY_FL_062_T1", "LY_FL_064_T1", "LY_FL_076_T1"))
colnames(bulk.matrix)[z] <- c(1, 2, 3)

bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix)

# Run Bisque

res <- ReferenceBasedDecomposition(bulk.eset, sc.eset, markers = NULL,
 use.overlap = TRUE)
ref.based.estimates <- res$bulk.props
ref.based.estimates <- data.frame(ref.based.estimates)

# Identify cased with zero P0 and P1

df <- t(ref.based.estimates)
colnames(df) <- paste0("P", colnames(df))
df <- data.frame(df)
cases_noP0_noP1 <- df %>%
  filter(P0 == 0 & P1 == 0) %>%
  rownames(.)

# Remove these cases from estimated proportions

ref.based.estimates$Population <- paste0("Population ",
 row.names(ref.based.estimates))
ref.based.estimates <- ref.based.estimates[ , !names(ref.based.estimates) %in% cases_noP0_noP1] 

# Generate plots

my_comparisons_p0_1 <- list(c("AR", "CS"), c("AR", "TT"),
 c("AR", "GM"), c("AR", "Q"))
my_comparisons_p8 <- list(c("GM", "CS"), c("GM", "TT"),
 c("GM", "Q"), c("GM", "AR"))

p0 <- ref.based.estimates %>%
  pivot_longer(!Population, names_to = "SAMPLE_ID",
   values_to = "estimated_proportion") %>%
  filter(Population == "Population 0") %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  ggboxplot("ClusterAIC", "estimated_proportion",
   color = "ClusterAIC", ncol = 5, lwd = 0.3, add = "jitter",
    add.params = list(size = 0.1, jitter = 0.2)) +
  theme_bw() +
  stat_compare_means(size = 2, label.x.npc = 0.15, label.y = 0.97) +
  stat_compare_means(comparisons = my_comparisons_p0_1,
   method = "wilcox.test", size = 2, label.y = c(0.6, 0.68, 0.76, 0.84)) +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
         axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
         axis.title.x = element_text(size = 7),
        plot.title = element_text(size = 7, vjust = -1.5)) +
  ggtitle("Population 0") +
  xlab("Subtype") +
  ylab("Estimated proportion") +
  ylim(-0.001, 1.0)

p1 <- ref.based.estimates %>%
  pivot_longer(!Population, names_to = "SAMPLE_ID",
   values_to = "estimated_proportion") %>%
  filter(Population == "Population 1") %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  ggboxplot("ClusterAIC", "estimated_proportion",
   color = "ClusterAIC", ncol = 5, lwd = 0.3, add = "jitter",
    add.params = list(size = 0.1, jitter = 0.2)) +
  theme_bw() +
  stat_compare_means(size = 2, label.x.npc = 0.15, label.y = 0.58) +
  stat_compare_means(comparisons = my_comparisons_p0_1,
   method = "wilcox.test", size = 2, label.y = c(0.32, 0.37, 0.42, 0.47)) +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
         axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
         axis.title.x = element_text(size = 7),
        plot.title = element_text(size = 7, vjust = -1.5)) +
  ggtitle("Population 1") +
  xlab("Subtype") +
  ylab("Estimated proportion") +
  ylim(-0.001, 0.6)

p8 <- ref.based.estimates %>% 
  pivot_longer(!Population, names_to = "SAMPLE_ID",
   values_to = "estimated_proportion") %>%
  filter(Population == "Population 8") %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  ggboxplot("ClusterAIC", "estimated_proportion",
   color = "ClusterAIC", ncol = 5, lwd = 0.3, add = "jitter",
    add.params = list(size = 0.1, jitter = 0.2)) +
  theme_bw() +
  stat_compare_means(size = 2, label.x.npc = 0.15, label.y = 0.58) +
  stat_compare_means(comparisons = my_comparisons_p8,
   method = "wilcox.test", size = 2, label.y = c(0.32, 0.37, 0.42, 0.47)) +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
         axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
         axis.title.x = element_text(size = 7),
        plot.title = element_text(size = 7, vjust = -1.5)) +
  ggtitle("Population 8") +
  xlab("Subtype") +
  ylab("Estimated proportion") +
  ylim(-0.001, 0.6)

p11 <- ref.based.estimates %>%
  pivot_longer(!Population, names_to = "SAMPLE_ID",
   values_to = "estimated_proportion") %>%
  filter(Population == "Population 11") %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  ggboxplot("ClusterAIC", "estimated_proportion",
   color = "ClusterAIC", ncol = 5, lwd = 0.3, add = "jitter",
    add.params = list(size = 0.1, jitter = 0.2)) +
  theme_bw() +
  stat_compare_means(size = 2, label.x.npc = 0.15, label.y = 0.58) +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
         axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
         axis.title.x = element_text(size = 7),
        plot.title = element_text(size = 7, vjust = -1.5)) +
  ggtitle("Population 11") +
  xlab("Subtype") +
  ylab("Estimated proportion") +
  ylim(-0.001, 0.6)

g <- gridExtra::arrangeGrob(p0, p1, p8, p11, nrow = 1)
ggsave(file = paste0("img/",  date, " _bisque_deconvolution_not_norm_full_matrix_adj_B.pdf"), g,
 width = 16, height = 5.2, units = "cm")

#--
# Generate Fig 7F
#--

DefaultAssay(exp.combined.sct) <- "RNA"
exp.combined.sct <- NormalizeData(exp.combined.sct)
B.cell.exp.combined.sct <- subset(exp.combined.sct, idents = c(0, 1, 8, 11))

diffexp.markers <- FindAllMarkers(B.cell.exp.combined.sct, only.pos = TRUE,
 min.pct = 0.1,logfc.threshold = 0.25, return.thresh = 0.01)
diffexp.markers <- as.data.table(diffexp.markers)
diffexp.markers <- diffexp.markers %>%
  filter(!grepl("RP", gene)) %>%
  filter(!grepl("MT-", gene))

top20 <- diffexp.markers %>%
 group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
   .$gene %>%
    unique()
pdf(paste0("img/", date, " dotplot_3_samples_0.3_all.marker_genes.pdf"),
 width = 17/2.54, height = 5/2.54)
DotPlot(B.cell.exp.combined.sct, features = top20,
        cols = c("#91bfdb", "#fc8d59"),
        cluster.idents = FALSE, dot.scale = 3) + RotatedAxis() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 5, angle = 90, face = "italic",
        vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 5),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.margin = margin(-10, -10, -10, -10)) +
  xlab("Gene") + ylab("Population")
dev.off()
