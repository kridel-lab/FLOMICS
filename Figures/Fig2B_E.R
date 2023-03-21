#---
# This script plots Figures 2B-E
# Author: Robert Kridel
#---

packages <- c("dplyr", "ggplot2", "survival", "ggpubr", "stringr", "tidyr", "gridExtra",
              "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg19", "magrittr", "data.table")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

setwd("~/github/FLOMICS/")

#--
# Read in data
#--

# Read in gene panels and define common genes
genes_PLOSMED <- read.csv("DNAseq/Final_target_files/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
  filter(PLOS_MED_PANEL == "YES") %>% .$Gene.Name # n=86
genes_UHN <- read.csv("DNAseq/Final_target_files/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
  filter(FLOMICS_PANEL == "YES") %>% .$Gene.Name # n=71
genes_common <- intersect(genes_PLOSMED, genes_UHN) # n=57

# Read in list of all mutations
capseq_T1 <- read.csv("DNAseq/CAPSEQ/Mutect2/all-cohorts/2022-05-24_AllCohorts_E2408-TO_filtered_MUTECT2_calls_BCL2_coding_and_noncoding.csv") %>%
  mutate(TIMEPOINT = str_sub(SAMPLE_ID, -2, -1)) %>%
  filter(Hugo_Symbol %in% genes_common) %>%
  filter(SAMPLE_ID != "LY_FL_571_T1") %>%
  filter(!SAMPLE_ID %in% c("LY_FL_1135_T1", "LY_FL_1156_T1")) %>%
  filter(TIMEPOINT == "T1")

capseq_T1_maf <- capseq_T1 %>%
  dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Protein_Change = AAChange.ensGene,
                Variant_Classification, Variant_Type, Tumor_Sample_Barcode = SAMPLE_ID, VAF = Var_Freq)

# Read in mutation matrix
capseq_T1_mat <- read.table("DNAseq/CAPSEQ/SNV_Cluster_Analysis/2022-05-15_UHN_E4402_PLOSMED-T-O_E2408-T-O_codingBCL2only/2022-05-15_gene_vs_sample_SNV_matrix.txt",
                            sep = ";", header = TRUE)
capseq_T1_mat <- capseq_T1_mat[ , !(names(capseq_T1_mat) %in% c("LY_FL_1135_T1", "LY_FL_1156_T1"))]

# Read in cluster information
clusters <- read.csv("DNAseq/CAPSEQ/SNV_Cluster_Analysis/2022-12-20_UHN_E4402_PLOSMED-T-O_E2408-T-O_codingBCL2only_713-samples/2023-01-19_most_confident_cluster_run_5_clusters/2023-01-19_14_GMM_Cluster_Labels_flexmix_clusters.csv")

#--
# Identify mutations within SHM motifs
#--

hg19.genome <- BSgenome.Hsapiens.UCSC.hg19
seqnames(hg19.genome) <- gsub("chr", "", seqnames(hg19.genome))

RGYW.motif <- DNAString("RGYW")
WRCY.motif <- DNAString("WRCY")

snvs.RGYW.gr <- capseq_T1_maf %>%
  dplyr::select(chr = Chromosome, pos = Start_Position, ref = Reference_Allele) %>%
  filter(ref == "G") %>%
  distinct(chr, pos, ref) %$%
  GRanges(chr, IRanges(start = pos - 1, end = pos + 2), id = str_c(chr, ":", pos))

snvs.WRCY.gr <- capseq_T1_maf %>%
  dplyr::select(chr = Chromosome, pos = Start_Position, ref = Reference_Allele) %>%
  filter(ref == "C") %>%
  distinct(chr, pos, ref) %$%
  GRanges(chr, IRanges(start = pos - 2, end = pos + 1), id = str_c(chr, ":", pos))

snvs.RGYW.seq <- getSeq(hg19.genome, snvs.RGYW.gr)
names(snvs.RGYW.seq) <- mcols(snvs.RGYW.gr)[, "id"]

snvs.WRCY.seq <- getSeq(hg19.genome, snvs.WRCY.gr)
names(snvs.WRCY.seq) <- mcols(snvs.WRCY.gr)[, "id"]

snvs.RGYW.seq.df <- data.table(posID = names(snvs.RGYW.seq), RGYWMotif = as.character(snvs.RGYW.seq))
snvs.WRCY.seq.df <- data.table(posID = names(snvs.WRCY.seq), WRCYMotif = as.character(snvs.WRCY.seq))

# Check if any of the SNVs overlap with a SHM motifs
snvs.RGYW.vmatch <- vmatchPattern(RGYW.motif, snvs.RGYW.seq, fixed = FALSE)
snvs.WRCY.vmatch <- vmatchPattern(WRCY.motif, snvs.WRCY.seq, fixed = FALSE)

capseq_T1_maf_motif <- capseq_T1_maf %>%
  mutate(posID = str_c(Chromosome, ":", Start_Position)) %>%
  mutate(motif = ifelse(posID %in% names(unlist(snvs.RGYW.vmatch)), "RGYW", 
                        ifelse(posID %in% names(unlist(snvs.WRCY.vmatch)), "WRCY", "None"))) %>%
  left_join(snvs.RGYW.seq.df) %>%
  left_join(snvs.WRCY.seq.df) %>%
  mutate(motif = factor(motif, levels = c("None", "RGYW", "WRCY")))

capseq_T1_maf_motif <- capseq_T1_maf_motif %>% arrange(Tumor_Sample_Barcode, Chromosome, Start_Position)

# Percentage SHM per cluster
capseq_T1_maf_motif %>%
  dplyr::select(Hugo_Symbol, SAMPLE_ID = Tumor_Sample_Barcode, motif) %>%
  left_join(clusters) %>%
  group_by(ClusterAIC, motif) %>%
  summarize(count = n()) %>%
  pivot_wider(names_from = motif, values_from = count) %>%
  dplyr::mutate(SHM = RGYW+WRCY) %>%
  dplyr::select(None, SHM) %>%
  mutate(perc.SHM = SHM/(None+SHM))

#--
# Generate Fig 2 panels B-E
#--

gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors <- gg_color_hue(5)

my_comparisons <- list(c("C3", "C1"), c("C3", "C2"), c("C3", "C4"), c("C3", "C5"))

p1 <- capseq_T1_maf %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  left_join(clusters[,c("SAMPLE_ID", "ClusterAIC")], by = c("Tumor_Sample_Barcode" = "SAMPLE_ID")) %>%
  mutate(ClusterAIC = factor(ClusterAIC, levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  # ggpubr::ggviolin(x = "ClusterAIC", y = "count", draw_quantiles = 0.5, #add = "jitter", add.params = list(size = 0.01, jitter = 0.2),
                   # color = "ClusterAIC",  palette = colors) +
  ggboxplot("ClusterAIC", "count", color = "ClusterAIC", ncol = 5, lwd = 0.3, add = "jitter", add.params = list(size = 0.001, jitter = 0.2)) +
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
  ylim(-4,45) + ylab("Number of mutations/sample")

p2 <- capseq_T1_mat %>%
  tibble::rownames_to_column(var = "gene") %>%
  pivot_longer(-gene, names_to = "SAMPLE_ID", values_to = "mutation") %>%
  filter(mutation == 1) %>%
  group_by(SAMPLE_ID) %>%
  summarize(count = n()) %>%
  left_join(clusters[,c("SAMPLE_ID", "ClusterAIC")]) %>%
  mutate(ClusterAIC = factor(ClusterAIC, levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  # ggviolin(x = "ClusterAIC", y = "count",  draw_quantiles = 0.5, #add = "jitter", add.params = list(size = 0.05, jitter = 0.2),
                   # color = "ClusterAIC",  palette = colors) +
  ggboxplot("ClusterAIC", "count", color = "ClusterAIC", ncol = 5, lwd = 0.3, add = "jitter", add.params = list(size = 0.001, jitter = 0.2)) +
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
  ylim(-4,28) + ylab("Number of mutated genes/sample")

p3 <- capseq_T1_maf_motif %>%
  filter(!motif == "None") %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  left_join(clusters[,c("SAMPLE_ID", "ClusterAIC")], by = c("Tumor_Sample_Barcode" = "SAMPLE_ID")) %>%
  dplyr::select(SAMPLE_ID = Tumor_Sample_Barcode, count, ClusterAIC) %>%
  right_join(clusters[,c("SAMPLE_ID", "ClusterAIC")]) %>%
  replace(is.na(.), 0) %>%
  mutate(ClusterAIC = factor(ClusterAIC, levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  # ggviolin(x = "ClusterAIC", y = "count",  draw_quantiles = 0.5, color = "ClusterAIC",  palette = colors) +
  ggboxplot("ClusterAIC", "count", color = "ClusterAIC", ncol = 5, lwd = 0.3, add = "jitter", add.params = list(size = 0.001, jitter = 0.2)) +
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
  ylim(-1.6,14) + ylab("Number of SHM mutations/sample")

p4 <- data.frame(ClusterAIC = c("C1", "C2", "C3", "C4", "C5"),
                 Stability = c(0.619217706, 0.585277076, 0.865382864, 1.005451774, 1.062214523)) %>%
  ggplot(aes(x = ClusterAIC, y = Stability, fill = ClusterAIC)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "none")

g <- arrangeGrob(p1, p2, p3, p4, nrow = 1)
ggsave(file = paste0("img/",  date, " number mutations per cluster.pdf"), g, width = 16, height = 5, units = "cm")

