#-------------------------------------------------------------------------------
# This script reads in the processed snRNAseq data and
# produces Fig 7 panels
#-------------------------------------------------------------------------------

setwd("~/your working directory/FLOMICS/")
# R 4.1.0
# Seurat 4.2.0

#--
# Load packages and set up environment
#--

options(stringsAsFactors = FALSE)
packages <- c("dplyr", "ggplot2", "ggpubr", "Seurat",
 "patchwork", "annotables", "IOBR", "BisqueRNA")
lapply(packages, require, character.only = TRUE)

#snRNA data
input_data <- "~/snRNAseq/"

#date
date <- Sys.Date()

#--
# Panels A - C
#--

clusters <- read.csv("GMM_Cluster_Labels_flexmix_clusters.csv") %>%
  select(SAMPLE_ID, ClusterAIC, cohort)

ecotyper.abundance <- read.table(file = "Ecotype_Abundance.txt",
 sep = "\t", header = TRUE)

ecotyper.abundance.boxplots <- ecotyper.abundance %>%
 reshape2::melt() %>%
  select(SAMPLE_ID = ID, Ecotype = variable, value) %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  mutate(ClusterAIC = factor(ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
  ggboxplot("ClusterAIC", "value", color = "ClusterAIC", facet.by = "Ecotype", ncol = 9,
            add = "jitter", add.params = list(size = 0.05, jitter = 0.2)) +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7)) + ylab("Proportion")
ecotyper.abundance.boxplots + stat_compare_means(method = "kruskal.test")

ecotyper.abundance.boxplots <- ecotyper.abundance %>%
 reshape2::melt() %>%
  select(SAMPLE_ID = ID, Ecotype = variable, value) %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  filter(Ecotype %in% c("LE1", "LE8")) %>%
  mutate(ClusterAIC = factor(ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
  ggboxplot("ClusterAIC", "value", color = "ClusterAIC", facet.by = "Ecotype", ncol = 1,
            add = "jitter", add.params = list(size = 0.05, jitter = 0.2)) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        strip.text.x = element_text(size = 7, hjust = 0, vjust = -1.5)) +
  stat_compare_means(method = "kruskal.test", size = 2, label.y.npc = 0.9, label.x.npc = 0.1) +
  ylab("Abundance") + xlab("Subtype")
ecotyper.abundance.boxplots

ecotyper.abundance.stacked.barplot <- ecotyper.abundance %>%
 reshape2::melt() %>%
  group_by(ID, variable) %>%
  summarize(mean = mean(value)) %>%
  inner_join(clusters, by = c("ID" = "SAMPLE_ID")) %>%
  mutate(ClusterAIC = factor(ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
  ggplot(aes(fill = variable, y = mean, x = ClusterAIC)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("#3E2883", "#76A2E4", "#8E549F", "#458833", "#BBD058",
                               "#FFFD61", "#F8D25D", "#F3A83B", "#EC5428")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.margin = margin(-10,-10,-10,-10)) +
  labs(fill = "Ecotypes") + ylab("Relative mean abundance") + xlab("Subtype")
ecotyper.abundance.stacked.barplot

#--
# B cell states
#--

b.cell.states.assignment <- read.table(file = "B.cells_Cell_State_Assignment.txt", sep = "\t", header = TRUE)

b.cell.states.assignment %>%
  inner_join(clusters, by = c("ID" = "SAMPLE_ID")) %>%
  count(ClusterAIC, Cell.State) %>%
  group_by(ClusterAIC) %>% 
  mutate(prop = prop.table(n)) %>%
  print(n = 100)

b.cell.states.abundance <- read.table(file = "B.cells_Cell_State_Abundance.txt", sep = "\t", header = TRUE)

p <- b.cell.states.abundance %>% reshape2::melt() %>%
  select(SAMPLE_ID = ID, Ecotype = variable, value) %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  mutate(ClusterAIC = factor(ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
  ggboxplot("ClusterAIC", "value", color = "ClusterAIC", facet.by = "Ecotype", ncol = 5,
            add = "jitter", add.params = list(size = 0.05, jitter = 0.2)) +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7)) + ylab("Proportion")
p + stat_compare_means(method = "kruskal.test")

b.cell.states.abundance.stacked.barplot <- b.cell.states.abundance %>%
 reshape2::melt() %>%
  group_by(ID, variable) %>%
  summarize(mean = mean(value)) %>%
  inner_join(clusters, by = c("ID" = "SAMPLE_ID")) %>%
  mutate(ClusterAIC = factor(ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
  ggplot(aes(fill = variable, y = mean, x = ClusterAIC)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("#EB7D5B", "#FED23F", "#B5D33D", "#6CA2EA", "#442288")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.margin = margin(-10,-10,-10,-10)) +
  labs(fill = "B cell states") + ylab("Relative mean abundance") + xlab("Subtype")
b.cell.states.abundance.stacked.barplot

#--
# Save figure
#--

g <- gridExtra::arrangeGrob(b.cell.states.abundance.stacked.barplot,
                            ecotyper.abundance.stacked.barplot,
                            ecotyper.abundance.boxplots,
                            ncol = 3)
ggsave(file = paste0("img/",  date, " ecotyper.pdf"), g, width = 16, height = 7, units = "cm")


#--
# Read Seurat objet in
#--
#exp.combined.sct
exp_combined_sct <- readRDS("FL_qc_integrated_clusters_3_samples_0.3_snn_res.rds")

#--
# Generate Panel D
#--

pdf(paste0("img/",  date, " _seurat_umap_clusters_3_samples_0.3_snn_res.pdf"),
 width = 10, height = 3.5)
p1 <- DimPlot(exp_combined_sct,
 reduction = "umap", label = TRUE, repel = TRUE) +
 NoLegend()
p2 <- DimPlot(exp_combined_sct,
 reduction = "umap", label = FALSE, repel = TRUE, split.by = "orig.ident")
p1 + p2 + plot_layout(widths = c(1, 3))
dev.off()

#--
# populations
#--

DefaultAssay(exp_combined_sct) <- "RNA"
exp_combined_sct <- NormalizeData(exp_combined_sct)
b_cell_exp_combined_sct <- subset(exp_combined_sct, idents = c(0, 1, 8, 11))

# Read in Attaf GC and memory genes
attaf_gc_up <- read.csv("Attaf_GC_up.csv") %>%
 filter(avg_logFC > 1) %>% .$Gene
attaf_mem_up <- read.csv("Attaf_Mem_up.csv") %>%
 filter(avg_logFC < -1) %>%
  .$Gene

b_cell_exp_combined_sct <- AddModuleScore(b_cell_exp_combined_sct,
 features = list(attaf_gc_up), name = "Attaf_GC_up")
b_cell_exp_combined_sct <- AddModuleScore(b_cell_exp_combined_sct,
 features = list(attaf_mem_up), name = "Attaf_Mem_up")

# Genes defining given populations
proliferation <- c("CENPE", "CENPF", "EZH2", "POLQ", "TOP2A", "MKI67")
pop1_genes <- c("BCL6", "KLHL6",  "RGS13", "PIP4K2A")
pop0_genes <- c("BANK1", "CXCR4", "FCMR", "FOXP1", "CCND3", "TXNIP")
tfh_educated <- c("CD83", "IL4I1", "IL21R", "MIR155HG", "TNFAIP3", "TRAF1")

marker_genes <- c(proliferation, pop1_genes, pop0_genes, tfh_educated)

levels(b_cell_exp_combined_sct) <- c("8", "1", "0", "11")

#--
# Prepare eset for Bisque
#--
seurat_object <- exp_combined_sct
delimiter <- "_"
position <- 2

get_cell_names <- function(obj) base::colnames(obj)
get_ident <- function(obj) Seurat::Idents(object = obj)
get_raw_data <- function(obj) Seurat::GetAssayData(object = obj[["RNA"]], slot = "counts")

individual_ids <- base::sapply(base::strsplit(get_cell_names(seurat_object),
 delimiter), `[[`, position)

base::names(individual_ids) <- get_cell_names(seurat_object)
individual_ids <- base::factor(individual_ids)
n_individuals <- base::length(base::levels(individual_ids))
base::message(base::sprintf("Split sample names by \"%s\"", delimiter),
              base::sprintf(" and checked position %i.", position),
              base::sprintf(" Found %i individuals.", n_individuals ))
base::message(base::sprintf("Example: \"%s\" corresponds to individual \"%s\".",
                            get_cell_names(seurat_object)[1],
                             individual_ids[1]))
sample_ids <- base::names(get_ident(seurat_object))
sc_pheno <- base::data.frame(check.names = FALSE, check.rows = FALSE,
                             stringsAsFactors = FALSE,
                             row.names = sample_ids,
                             SubjectName = individual_ids,
                             cellType = get_ident(seurat_object))
sc_meta <- base::data.frame(labelDescription = base::c("SubjectName",
 "cellType"), row.names = base::c("SubjectName", "cellType"))
sc_pdata <- methods::new("AnnotatedDataFrame", data = sc_pheno,
 varMetadata = sc_meta)

sc_data <- base::as.matrix(get_raw_data(seurat_object)[, sample_ids,
 drop = FALSE])
sc_eset <- Biobase::ExpressionSet(assayData = sc_data, phenoData = sc_pdata)

#--
# Run Bisque
#--

# Read in genetic clusters
clusters <- read.csv("GMM_Subtype_Labels_flexmix_clusters.csv") %>%
  select(SAMPLE_ID, ClusterAIC, cohort)

# Read in RNAseq QC data
passing_rnaseq_qc <-
  read.table("metadata_all_22_23.txt", "\t", header = TRUE) %>%
  dplyr::filter(qc_tier2 == "Y")

# Prepare expression matrix
# bulk gene expression matrix is not normalized
## as the analysis is based on proportion of given populations/sample
bulk_matrix <- read.table("RNA_seq_analysis_combined_2022_2023_2023-05-01_df_counts_adj_filtered.txt", sep = " ",
 header = TRUE)
bulk_matrix <- bulk_matrix[, passing_rnaseq_qc$sample_id]

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
ref_based_estimates <- res$bulk.props

my_comparisons_p0_1 <- list(
  c("AR", "CS"), c("AR", "TT"), c("AR", "GM"), c("AR", "Q"))

#--
# Generate panel E
#--

#learning what is the maximum value, to inform the limits of the y-axis
p0_test <- ref_based_estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(0))
max(p0_test$value) # 0.5645319

p0 <- ref_based_estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(0)) %>%
  mutate(Var1 = paste0("Population ", Var1)) %>%
  select(Population = Var1, SAMPLE_ID = Var2, value) %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  mutate(Cluster = factor(
    ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
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
  ylim(-0.001, 1.0) +
  xlab("Subtype")
p0 <- p0 + labs(subtitle = "Population 0")

#learning what is the maximum value, to inform the limits of the y-axis
p1_test <- ref_based_estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(1))
max(p1_test$value) # 0.3209829

p1 <- ref_based_estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(1)) %>%
  mutate(Var1 = paste0("Population ", Var1)) %>%
  select(Population = Var1, SAMPLE_ID = Var2, value) %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  mutate(Cluster = factor(
    ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
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
  ylim(-0.001, 0.6) +
  xlab("Subtype")
p1 <- p1 + labs(subtitle = "Population 1") 

my_comparisons_p8 <- list(
  c("GM", "CS"), c("GM", "TT"), c("GM", "Q"), c("GM", "AR"))

#learning what is the maximum value, to inform the limits of the y-axis
p8_test <- ref_based_estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(8))
max(p8_test$value) #0.6720876

p8 <- ref_based_estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(8)) %>%
  mutate(Var1 = paste0("Population ", Var1)) %>%
  select(Population = Var1, SAMPLE_ID = Var2, value) %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  mutate(Cluster = factor(
    ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
  ggboxplot("Cluster", "value", color = "Cluster", ncol = 5, lwd = 0.3,
   add = "jitter", add.params = list(size = 0.1, jitter = 0.2)) +
  theme_bw() +
  stat_compare_means(size = 2, label.x.npc = 0.15, label.y = 0.97) +
  stat_compare_means(comparisons = my_comparisons_p8, method = "wilcox.test",
   size = 2, label.y = c(0.70, 0.77, 0.735, 0.84)) +
  theme(legend.position = "none", strip.background = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7)) +
  ylab("Estimatd proportion") +
  ylim(-0.001, 1) +
  xlab("Subtype")
p8 <- p8 + labs(subtitle = "Population 8")

#learning what is the maximum value, to inform the limits of the y-axis
p11_test <- ref_based_estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(11))
max(p11_test$value) # 0.1040954

p11 <- ref_based_estimates %>%
  reshape2::melt() %>%
  filter(Var1 %in% c(11)) %>%
  mutate(Var1 = paste0("Population ", Var1)) %>%
  select(Population = Var1, SAMPLE_ID = Var2, value) %>%
  inner_join(clusters, by = "SAMPLE_ID") %>%
  mutate(Cluster = factor(
    ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
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
  ylim(-0.001, 0.6) +
  xlab("Subtype")

p11 <- p11 + labs(subtitle = "Population 11")

g <- gridExtra::arrangeGrob(p0, p1, p8, p11, nrow = 1)
ggsave(file = paste0(
  date, " _bisque_deconvolution_not_norm_full_matrix_adj_B.pdf"),
   g, width = 16, height = 5, units = "cm")

#--
# generate panel F
#--

