#---
# This script plots Figures 2B-E
# Author: Robert Kridel
#---

packages <- c("dplyr", "ggplot2", "survival", "ggpubr", "stringr", "tidyr",
 "gridExtra", "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg19", "magrittr",
  "data.table")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

setwd("~/your working directory/FLOMICS/")

#--
# Read in data
#--

# Read in gene panels and define common genes
genes_plosmed <- read.csv(
  "DNAseq/Final_target_files/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
  filter(PLOS_MED_PANEL == "YES") %>%
   .$Gene.Name
genes_uhn <- read.csv(
  "DNAseq/Final_target_files/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
  filter(FLOMICS_PANEL == "YES") %>%
   .$Gene.Name
genes_common <- intersect(genes_PLOSMED, genes_UHN)

# Read in list of all mutations
capseq_t1 <- read.csv("filtered_MUTECT2_calls.csv") %>%
  mutate(TIMEPOINT = str_sub(SAMPLE_ID, -2, -1)) %>%
  filter(Hugo_Symbol %in% genes_common) %>%
  filter(SAMPLE_ID != "LY_FL_571_T1") %>%
  filter(!SAMPLE_ID %in% c("LY_FL_1135_T1", "LY_FL_1156_T1")) %>%
  filter(TIMEPOINT == "T1")

capseq_t1_maf <- capseq_t1 %>%
  dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position,
   Reference_Allele, Tumor_Seq_Allele2, Protein_Change = AAChange.ensGene,
   Variant_Classification, Variant_Type, Tumor_Sample_Barcode = SAMPLE_ID,
   VAF = Var_Freq)

# Read in mutation matrix
capseq_t1_mat <- read.table("gene_vs_sample_SNV_matrix.txt",
                            sep = ";", header = TRUE)
capseq_t1_mat <- capseq_t1_mat[, !(names(capseq_T1_mat) %in%
 c("LY_FL_1135_T1", "LY_FL_1156_T1"))]

# Read in cluster information
clusters <- read.csv("GMM_Cluster_Labels_flexmix_clusters.csv")

#--
# Identify mutations within SHM motifs
#--

hg19_genome <- BSgenome.Hsapiens.UCSC.hg19
seqnames(hg19_genome) <- gsub("chr", "", seqnames(hg19.genome))

rgyw_motif <- DNAString("RGYW")
wrcy_motif <- DNAString("WRCY")

snvs_rgyw_gr <- capseq_t1_maf %>%
  dplyr::select(chr = Chromosome, pos = Start_Position,
   ref = Reference_Allele) %>%
  filter(ref == "G") %>%
  distinct(chr, pos, ref) %$%
  GRanges(chr, IRanges(start = pos - 1, end = pos + 2),
   id = str_c(chr, ":", pos))

snvs_wrcy_gr <- capseq_t1_maf %>%
  dplyr::select(chr = Chromosome, pos = Start_Position,
   ref = Reference_Allele) %>%
  filter(ref == "C") %>%
  distinct(chr, pos, ref) %$%
  GRanges(chr, IRanges(start = pos - 2, end = pos + 1),
   id = str_c(chr, ":", pos))

snvs_rgyw_seq <- getSeq(hg19.genome, snvs_rgyw_gr)
names(snvs_rgyw_seq) <- mcols(snvs_rgyw_gr)[, "id"]

snvs_wrcy_seq <- getSeq(hg19.genome, snvs_wrcy_gr)
names(snvs_wrcy_seq) <- mcols(snvs_wrcy_gr)[, "id"]

snvs_rgyw_seq_df <- data.table(posID = names(snvs_rgyw_seq),
 RGYWMotif = as.character(snvs_rgyw_seq))
snvs_wrcy_seq_df <- data.table(posID = names(snvs_wrcy_seq),
 WRCYMotif = as.character(snvs_wrcy_seq))

# Check if any of the SNVs overlap with a SHM motifs
snvs_rgyw_vmatch <- vmatchPattern(RGYW.motif, snvs_rgyw_seq, fixed = FALSE)
snvs_wrcy_vmatch <- vmatchPattern(WRCY.motif, snvs_wrcy_seq, fixed = FALSE)

capseq_t1_maf_motif <- capseq_t1_maf %>%
  mutate(posID = str_c(Chromosome, ":", Start_Position)) %>%
  mutate(motif = ifelse(posID %in% names(unlist(snvs_rgyw_vmatch)), "RGYW",
                        ifelse(posID %in% names(unlist(snvs_wrcy_vmatch)),
                         "WRCY", "None"))) %>%
  left_join(snvs_rgyw_seq_df) %>%
  left_join(snvs_wrcy_seq_df) %>%
  mutate(motif = factor(motif, levels = c("None", "RGYW", "WRCY")))

capseq_t1_maf_motif <- capseq_t1_maf_motif %>%
 arrange(Tumor_Sample_Barcode, Chromosome, Start_Position)

# Percentage SHM per cluster
capseq_t1_maf_motif %>%
  dplyr::select(Hugo_Symbol, SAMPLE_ID = Tumor_Sample_Barcode, motif) %>%
  left_join(clusters) %>%
  group_by(ClusterAIC, motif) %>%
  summarize(count = n()) %>%
  pivot_wider(names_from = motif, values_from = count) %>%
  dplyr::mutate(SHM = RGYW + WRCY) %>%
  dplyr::select(None, SHM) %>%
  mutate(perc.SHM = SHM / (None + SHM))

#--
# Generate Fig 2 panels B-E
#--

gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors <- gg_color_hue(5)

my_comparisons <- list(c("C3", "C1"), c("C3", "C2"),
 c("C3", "C4"), c("C3", "C5"))

p1 <- capseq_t1_maf %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  left_join(clusters[, c("SAMPLE_ID", "ClusterAIC")],
   by = c("Tumor_Sample_Barcode" = "SAMPLE_ID")) %>%
  mutate(ClusterAIC = factor(ClusterAIC,
   levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  ggboxplot("ClusterAIC", "count", color = "ClusterAIC",
   ncol = 5, lwd = 0.3, add = "jitter",
    add.params = list(size = 0.001, jitter = 0.2)) +
  stat_compare_means(size = 2, label.x.npc = 0.15) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     size = 2, label.y = c(25, 29, 33, 37)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "none") +
  ylim(-4, 45) + ylab("Number of mutations/sample")

p2 <- capseq_t1_mat %>%
  tibble::rownames_to_column(var = "gene") %>%
  pivot_longer(-gene, names_to = "SAMPLE_ID", values_to = "mutation") %>%
  filter(mutation == 1) %>%
  group_by(SAMPLE_ID) %>%
  summarize(count = n()) %>%
  left_join(clusters[, c("SAMPLE_ID", "ClusterAIC")]) %>%
  mutate(ClusterAIC = factor(ClusterAIC,
   levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  ggboxplot("ClusterAIC", "count", color = "ClusterAIC",
   ncol = 5, lwd = 0.3, add = "jitter",
    add.params = list(size = 0.001, jitter = 0.2)) +
  stat_compare_means(size = 2, label.x.npc = 0.15) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     size = 2, label.y = c(15, 18, 21, 24)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "none") +
  ylim(-4, 28) + ylab("Number of mutated genes/sample")

p3 <- capseq_t1_maf_motif %>%
  filter(!motif == "None") %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  left_join(clusters[, c("SAMPLE_ID", "ClusterAIC")],
   by = c("Tumor_Sample_Barcode" = "SAMPLE_ID")) %>%
  dplyr::select(SAMPLE_ID = Tumor_Sample_Barcode, count, ClusterAIC) %>%
  right_join(clusters[, c("SAMPLE_ID", "ClusterAIC")]) %>%
  replace(is.na(.), 0) %>%
  mutate(ClusterAIC = factor(ClusterAIC,
   levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  ggboxplot("ClusterAIC", "count", color = "ClusterAIC",
   ncol = 5, lwd = 0.3, add = "jitter",
    add.params = list(size = 0.001, jitter = 0.2)) +
  stat_compare_means(size = 2, label.x.npc = 0.15) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     size = 2, label.y = c(7.5, 9, 10.5, 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "none") +
  ylim(-1.6, 14) + ylab("Number of SHM mutations/sample")

p4 <- data.frame(ClusterAIC = c("C1", "C2", "C3", "C4", "C5"),
                 Stability = c(0.619217706, 0.585277076, 0.865382864,
                  1.005451774, 1.062214523)) %>%
  ggplot(aes(x = ClusterAIC, y = Stability, fill = ClusterAIC)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "none")

g <- arrangeGrob(p1, p2, p3, p4, nrow = 1)
ggsave(file = paste0("img/",  date, " number mutations per cluster.pdf"),
 g, width = 16, height = 5, units = "cm")
