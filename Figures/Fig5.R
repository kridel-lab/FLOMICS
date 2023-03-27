#---
# This script analyzes TRUST4 data and generated Figure 5
# Author: Robert Kridel
#---

packages <- c("dplyr", "readr", "ggpubr")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

setwd("~/your working directory/FLOMICS/")

# Read in RNAseq QC
rnaseq_samples_passing_qc <-
  read.table("metadata_passed_290.txt", "\t", header = TRUE) %>% .$id

# Read in clusters
clusters <- read.csv("GMM_Cluster_Labels_flexmix_clusters.csv") %>%
  select(SAMPLE_ID, ClusterAIC, cohort)

# Read in sample annotation
sample_annotation <- read.table("metadata/20221228_sample_annotations.txt",
 header = TRUE, sep = "\t")
fl_cases <- sample_annotation %>% filter(TYPE == "FL") %>% .$SAMPLE_ID

# Read in clinical data
clin <- read.csv("metadata/20220525_clinical_data.csv") %>%
  mutate(LY_FL_ID = paste0(LY_FL_ID, "_T1")) %>%
  select(LY_FL_ID, ANN_ARBOR_STAGE)

# Read in TRUST4 output
airr_tgl_ocir <- read_tsv("TRUST4/TGL_OICR/IMGT_airr_merged.tsv") %>%
  mutate(sample_id = ifelse(sample_id %in%
   c("LY_FL_234", "LY_FL_235",
   "LY_FL_236", "LY_FL_237", "LY_FL_238", "LY_FL_246", "LY_FL_247",
    "LY_FL_249", "LY_FL_252", "LY_FL_253", "LY_FL_254", "LY_FL_255",
     "LY_FL_256", "LY_FL_257", "LY_FL_258", "LY_FL_261"),
      paste0(sample_id, "_T1"), sample_id))

airr_e4402 <- read_tsv("TRUST4/E4402/IMGT_airr_merged.tsv")

# Combine and retain only those cases for which RNAseq data passed QC
airr <- rbind(airr_tgl_ocir, airr_e4402) %>%
 filter(sample_id %in% rnaseq_samples_passing_qc)

# Number of cases
length(unique(airr$sample_id))

#--
# Identify most common isotype vs clusters
#--

clonotypes_igh_fl <- airr %>%
  filter(sample_id %in% fl_cases) %>%
  filter(grepl("IGH", v_call))

clonotypes_igh_fl_tot <- clonotypes_igh_fl %>%
  group_by(sample_id) %>%
  summarize(count.tot = sum(consensus_count))

summary(clonotypes_igh_fl_tot$count.tot)

clonotypes_igh_fl_vdjgrouped <- clonotypes_igh_fl %>%
  mutate(v_call = sub("\\*.*", "", v_call)) %>%
  mutate(d_call = substr(d_call, 1, 5)) %>%
  mutate(j_call = substr(j_call, 1, 5)) %>%
  mutate(c_call = sub("\\*.*", "", c_call)) %>%
  group_by(sample_id, v_call, d_call, j_call, c_call) %>%
  summarize(count = sum(consensus_count))
dim(clonotypes_igh_fl_vdjgrouped)

clonotypes_igh_fl_vdjgrouped_dominant_all <- clonotypes_igh_fl_vdjgrouped %>%
  group_by(sample_id) %>%
  slice_max(n = 1, order_by = count)
summary(clonotypes_igh_fl_vdjgrouped_dominant_alll$count)

clonotypes_igh_fl_vdjgrouped_dominant <- clonotypes_igh_fl_vdjgrouped %>%
  group_by(sample_id) %>%
  slice_max(n = 1, order_by = count) %>%
  filter(!is.na(c_call)) %>%
  left_join(clonotypes.IGH.FL_tot) %>%
  mutate(rel.prop.dom = count / count.tot)
dim(clonotypes_igh_fl_vdjgrouped_dominant)

# Determine whether either IgG or IgM expression is enriched in genetic clusters
df <- clonotypes_igh_fl_vdjgrouped_dominant %>%
  mutate(isotype = c_call) %>%
  left_join(clusters[, c("SAMPLE_ID", "ClusterAIC")],
   by = c("sample_id" = "SAMPLE_ID")) %>%
  filter(!is.na(ClusterAIC)) %>%
  mutate(isotype = ifelse(isotype == "IGHM", "IGHM", "non-IGHM")) %>%
  group_by(ClusterAIC, isotype) %>%
  dplyr::summarize(n = n()) %>%
  tidyr::spread(ClusterAIC, n) %>%
  replace(is.na(.), 0)

b1 <- binom.test(as.numeric(df[2, "C1"]), colSums(df[, "C1"]), p = 0.5,
 alternative = c("two.sided"), conf.level = 0.95)$p.value
b1 <- round(b1, 3)
b2 <- binom.test(as.numeric(df[2, "C2"]), colSums(df[, "C2"]), p = 0.5,
 alternative = c("two.sided"), conf.level = 0.95)$p.value
b2 <- round(b2, 3)
b3 <- binom.test(as.numeric(df[2, "C3"]), colSums(df[, "C3"]), p = 0.5,
 alternative = c("two.sided"), conf.level = 0.95)$p.value
b3 <- round(b3, 3)
b4 <- binom.test(as.numeric(df[2, "C4"]), colSums(df[, "C4"]), p = 0.5,
 alternative = c("two.sided"), conf.level = 0.95)$p.value
b4 <- round(b4, 3)
b5 <- binom.test(as.numeric(df[2, "C5"]), colSums(df[, "C5"]), p = 0.5,
 alternative = c("two.sided"), conf.level = 0.95)$p.value
b5 <- formatC(b5, format = "e", digits = 2)

#--
# Plot the numbers of IgG or IgM expressing cases by genetic cluster
#--

# First, define annotation data frame:
ann_text <- data.frame(ClusterAIC = c("C1", "C2", "C3", "C4", "C5"),
                       isotype = rep(c("IGHG"), 5),
                       label = paste0("P = ", c(b1, b2, b3, b4, b5)),
                       n = rep(c(42), 5))

p1 <- clonotypes_igh_fl_vdjgrouped_dominant %>%
  mutate(isotype = c_call) %>%
  left_join(clusters[, c("SAMPLE_ID", "ClusterAIC")],
   by = c("sample_id" = "SAMPLE_ID")) %>%
  filter(!is.na(ClusterAIC)) %>%
  mutate(isotype = ifelse(isotype %in%
                     c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), "IGHG",
                   ifelse(isotype %in% c("IGHA1", "IGHA2"), "IGHA",
                   ifelse(isotype == "IGHM", "IGHM", isotype)))) %>%
  group_by(ClusterAIC, isotype) %>%
  dplyr::summarize(n = n()) %>%
  mutate(isotype = factor(isotype, levels = c("IGHM", "IGHG", "IGHA"))) %>%
  ggplot(aes(x = isotype, y = n, fill = isotype, label = n)) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  facet_wrap(~ClusterAIC, ncol = 5) +
  geom_text(data = ann_text, size = 2, aes(label = label, x = isotype, y = n)) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  ylab("Number of samples with given isotype") +
  theme(axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 5),
        axis.title.x = element_blank(), legend.position = "none",
        strip.background = element_blank(),
        plot.subtitle = element_blank()) +
  ggtitle(element_blank()) +
  scale_fill_manual(values = c("#a6dba0", "#c2a5cf", "#9970ab"))
p1

#--
# Plot the V gene identity in IgM and IgG-expressing cases
#--

df2 <- clonotypes_igh_fl %>%
  mutate(v_call = sub("\\*.*", "", v_call)) %>%
  mutate(d_call = substr(d_call, 1, 5)) %>%
  mutate(j_call = substr(j_call, 1, 5)) %>%
  mutate(c_call = sub("\\*.*", "", c_call)) %>%
  group_by(sample_id, v_call, d_call, j_call, c_call, v_identity) %>%
  summarize(count = sum(consensus_count)) %>%
  group_by(sample_id) %>%
  slice_max(n = 1, order_by = count) %>%
  filter(!is.na(c_call)) %>%
  left_join(clonotypes.IGH.FL_tot) %>%
  mutate(rel.prop.dom = count / count.tot) %>%
  left_join(clusters[, c("SAMPLE_ID", "ClusterAIC")],
   by = c("sample_id" = "SAMPLE_ID")) %>%
  filter(ClusterAIC == factor(ClusterAIC,
   levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  filter(!is.na(ClusterAIC)) %>%
  data.frame()

df2 %>% mutate(c_call = ifelse(c_call %in%
                         c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), "IGHG",
                        ifelse(c_call %in%
                         c("IGHA1", "IGHA2"), "IGHA", c_call))) %>%
  filter(c_call %in% c("IGHM", "IGHG")) %>%
  group_by(c_call) %>%
   summarize(mean(v_identity))

p2 <- df2 %>% mutate(c_call = ifelse(c_call %in%
                               c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), "IGHG",
                              ifelse(c_call %in%
                               c("IGHA1", "IGHA2"), "IGHA", c_call))) %>%
  filter(c_call %in% c("IGHM", "IGHG")) %>%
  ggboxplot("c_call", "v_identity", color = "c_call",
   palette = c("#a6dba0", "#c2a5cf"), ncol = 5, lwd = 0.5, add = "jitter",
    add.params = list(size = 0.1, jitter = 0.2)) +
  stat_compare_means(method = "wilcox.test", size = 2, label.x.npc = 0.25,
   label.y.npc = 0.95) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "none") +
  ylab("V gene identity (%)") + xlab("Isotype") + ylim(70, 110)

p2

#--
# Plot the V gene identity
#--

my_comparisons <- list(
  c("C3", "C1"), c("C3", "C2"), c("C3", "C4"), c("C3", "C5"))

p3 <- df2 %>% mutate(ClusterAIC = factor(ClusterAIC,
 levels = c("C1", "C2", "C3", "C4", "C5"))) %>%
  select(Cluster = ClusterAIC, v_identity) %>%
  ggboxplot("Cluster", "v_identity", color = "Cluster", ncol = 5, lwd = 0.3,
   add = "jitter", add.params = list(size = 0.1, jitter = 0.2)) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     size = 2, label.y = c(100, 103, 106, 109)) +
  stat_compare_means(method = "kruskal.test",
                     size = 2,
                     label.x.npc = 0.01, label.y.npc = 0.01) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7),
   axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
         axis.text.y = element_text(size = 7),
        legend.position = "none") +
  ylab("V gene identity (%)")

g <- gridExtra::grid.arrange(p1, gridExtra::arrangeGrob(p2, p3,
 ncol = 2, widths = c(0.8, 1)), ncol = 1)
ggsave(file = paste0("img/",  date, " Fig5.pdf"), g, width = 14,
 height = 12, units = "cm")
