#---
# This script analyzes TRUST4 data and generated Figure 6 - new
# Author: Robert Kridel & Victoria Shelton
#---

packages <- c("dplyr", "readr", "ggpubr", "ijtiff", "png")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

setwd("~/your working directory/TRUST4/")

# Read in RNAseq QC
rnaseq_samples_passing_qc <-
  read.table("metadata_all_22_23.txt", "\t", header = TRUE) %>%
  dplyr::filter(squence_core != "TGL_23") %>%
  dplyr::filter(qc_tier2 == "Y") %>% #290 samples
  .$sample_id #selecting only sample names

# Read in clusters
clusters <- read.csv("GMM_Subtype_Labels_flexmix_clusters.csv") %>%
  select(SAMPLE_ID, ClusterAIC, cohort)

# Read in sample annotation
sample_annotation <- read.table("sample_annotations.txt",
 header = TRUE, sep = "\t")
fl_cases <- sample_annotation %>% filter(TYPE == "FL") %>% .$SAMPLE_ID

# Read in clinical data
clin <- read.csv("20220525_clinical_data.csv") %>%
  mutate(LY_FL_ID = paste0(LY_FL_ID, "_T1")) %>%
  select(LY_FL_ID, ANN_ARBOR_STAGE)

# Read in TRUST4 output
airr_tgl_ocir <- read_tsv("TGL_OICR_IMGT_airr_merged.tsv") %>%
  mutate(sample_id = ifelse(sample_id %in%
   c("LY_FL_234", "LY_FL_235",
   "LY_FL_236", "LY_FL_237", "LY_FL_238", "LY_FL_246", "LY_FL_247",
    "LY_FL_249", "LY_FL_252", "LY_FL_253", "LY_FL_254", "LY_FL_255",
     "LY_FL_256", "LY_FL_257", "LY_FL_258", "LY_FL_261"),
      paste0(sample_id, "_T1"), sample_id))
length(unique(airr_tgl_ocir$sample_id)) #153 samples

airr_e4402 <- read_tsv("E4402_IMGT_airr_merged.tsv")
length(unique(airr_e4402$sample_id)) #209 samples

# Combine and retain only those cases for which RNAseq data passed QC
airr <- rbind(airr_tgl_ocir, airr_e4402) %>%
 filter(sample_id %in% rnaseq_samples_passing_qc)

# Number of cases
length(unique(airr$sample_id)) #290 samples in total

# Read in IHC information
IHC_results <- read_csv("Isotype_analysis_IHC.csv")

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
summary(clonotypes_igh_fl_vdjgrouped_dominant_all$count)

clonotypes_igh_fl_vdjgrouped_dominant <- clonotypes_igh_fl_vdjgrouped %>%
  group_by(sample_id) %>%
  slice_max(n = 1, order_by = count) %>%
  filter(!is.na(c_call)) %>%
  left_join(clonotypes_igh_fl_tot) %>%
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
paste0("CS samples: ", sum(df$CS), " ; ", "TT samples: ", sum(df$TT), " ; ", "GM samples: ", sum(df$GM), " ; ", "Q samples: ", sum(df$Q), " ; ","AR samples: ", sum(df$AR))
 #CS samples: 73 ; TT samples: 39 ; GM samples: 45 ; Q samples: 42 ; AR samples: 34"
paste0("total samples: ", sum(df$AR, df$CS, df$GM, df$Q, df$TT)) #233 samples


b1 <- binom.test(as.numeric(df[2, "CS"]), colSums(df[, "CS"]), p = 0.5,
 alternative = c("two.sided"), conf.level = 0.95)$p.value
b1 <- round(b1, 3)
b2 <- binom.test(as.numeric(df[2, "TT"]), colSums(df[, "TT"]), p = 0.5,
 alternative = c("two.sided"), conf.level = 0.95)$p.value
b2 <- round(b2, 3)
b3 <- binom.test(as.numeric(df[2, "GM"]), colSums(df[, "GM"]), p = 0.5,
 alternative = c("two.sided"), conf.level = 0.95)$p.value
b3 <- round(b3, 3)
b4 <- binom.test(as.numeric(df[2, "Q"]), colSums(df[, "Q"]), p = 0.5,
 alternative = c("two.sided"), conf.level = 0.95)$p.value
b4 <- round(b4, 3)
b5 <- binom.test(as.numeric(df[2, "AR"]), colSums(df[, "AR"]), p = 0.5,
 alternative = c("two.sided"), conf.level = 0.95)$p.value
b5 <- formatC(b5, format = "e", digits = 2)

#--
# Panel A - Plot the numbers of IgG or IgM expressing cases by genetic cluster
#--

# First, define annotation data frame:
ann_text <- data.frame(ClusterAIC = c("CS", "TT", "GM", "Q", "AR"),
                       isotype = rep(c("IGHG"), 5),
                       label = paste0("P = ", c(b1, b2, b3, b4, b5)),
                       n = rep(c(42), 5))
ann_text$ClusterAIC <- factor(ann_text$ClusterAIC,
 levels = c("CS", "TT", "GM", "Q", "AR"))
subtype_order <- c("CS", "TT", "GM", "Q", "AR")

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
  arrange(factor(ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR")),
    factor(isotype, levels = c("IGHM", "IGHG", "IGHA"))) #%>%
p1$ClusterAIC <- factor(p1$ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))

p1 <- p1 %>% ggplot(aes(x = isotype, y = n, fill = isotype, label = n)) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  facet_wrap(~ClusterAIC, ncol = 5) +
  geom_text(data = ann_text, size = 2,
   aes(label = label, x = isotype, y = n + 4)) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  ylab("Number of samples with given isotype") +
  theme(axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 5),
        axis.title.x = element_blank(), legend.position = "none",
        strip.background = element_blank(),
        plot.subtitle = element_blank(),
        strip.text.x = element_text(hjust = 0)) +
  ggtitle(element_blank()) +
  scale_fill_manual(values = c("#a6dba0", "#c2a5cf", "#9970ab"))

#--
# Panels B & C - reading in images
#--
img1 <- readPNG("IGM EAB2953 A6 100x.png")
img2 <- readPNG("IGM EAB2953 A6 200x.png")

im_A <- ggplot() +
  background_image(img1) +
  # This ensures that the image leaves some space at the edges
  theme(plot.margin = margin(t = 1, l = 0.6,
                             r = 0.2,
                             b = 1.5, unit = "cm"))

im_B <- ggplot() + background_image(img2) +
  theme(plot.margin = margin(t = 1,
                             l = 0.2,
                             r = 0.4, b = 1.5, unit = "cm"))

#--
# Panel D - Plot the IHC expression and non-expression
#--

IHC_results_E2408 <- IHC_results %>%
  mutate(OTHER_ID = as.character(OTHER_ID)) %>%
  select(-SAMPLE_ID) %>%
  filter(INSTITUTION == "E2408") %>%
  left_join(sample_annotation[, c("SAMPLE_ID", "OTHER_ID")],
   by = "OTHER_ID") %>%
  mutate(TIME_POINT = substr(SAMPLE_ID, 12, 13)) %>%
  filter(TIME_POINT == "T1") %>%
  select(SAMPLE_ID, OTHER_ID, INSTITUTION, IGM, IGG, COMMENT)

IHC_results <- IHC_results %>%
  filter(INSTITUTION == "UHN") %>%
  rbind(IHC_results_E2408) %>%
  mutate(IGM = ifelse(IGM == "-", "NEG",
                      ifelse(IGM == "+?", "UNDETERMINED",
                             ifelse(IGM == "+", "POS",
                                    ifelse(IGM == "++",
                                     "POS_STRONG", NA))))) %>%
  mutate(IGG = ifelse(IGG == "-", "NEG",
                      ifelse(IGG == "+?", "UNDETERMINED",
                             ifelse(IGG == "+", "POS", NA)))) %>%
  left_join(clonotypes_igh_fl_vdjgrouped_dominant,
   by = c("SAMPLE_ID" = "sample_id")) %>%
  left_join(clusters[, c("SAMPLE_ID", "ClusterAIC")], by = c("SAMPLE_ID"))

nrow(IHC_results)
IHC_results_undetermined <- IHC_results %>%
  dplyr::filter(IGM == "UNDETERMINED") %>%
  nrow()
#11 samples are undetermined

table(IHC_results$IGM, useNA = "ifany")

table(IHC_results$IGG, useNA = "ifany")

table(IHC_results$c_call, useNA = "ifany")

table(IHC_results$ClusterAIC, useNA = "ifany")

table("cluster" = IHC_results$ClusterAIC, "IGM" = IHC_results$IGM)

table("cluster" = IHC_results$ClusterAIC, "IGG" = IHC_results$IGG)

table("IHC" = IHC_results$c_call, "RNAseq" = IHC_results$IGM, useNA = "ifany")

table("IHC" = IHC_results$c_call, "RNAseq" = IHC_results$IGG, useNA = "ifany")

df_IHC <- IHC_results %>%
  filter(!is.na(ClusterAIC)) %>%
  filter(IGM != "UNDETERMINED") %>%
  mutate(IGM = ifelse(IGM == "POS_STRONG", "POS", IGM)) %>%
  group_by(ClusterAIC) %>%
  dplyr::summarize(n = n()) %>%
  tidyr::spread(ClusterAIC, n) %>%
  replace(is.na(.), 0)

ann_text_ihc <- data.frame(ClusterAIC = c("CS", "TT", "GM", "Q", "AR"),
                           total = paste0("N=",
                            c(df_IHC$CS, df_IHC$TT,
                             df_IHC$GM, df_IHC$Q, df_IHC$AR)),
                           IGHM = rep("NEG", 5),
                           n = rep(c(1.01), 5))
ann_text_ihc$ClusterAIC <- factor(ann_text_ihc$ClusterAIC,
 levels = c("CS", "TT", "GM", "Q", "AR"))
subtype_order <- c("CS", "TT", "GM", "Q", "AR")

#filled barplolt for pos vs neg status of ighm ihc
p2 <-
  IHC_results %>%
  filter(!is.na(ClusterAIC)) %>%
  filter(IGM != "UNDETERMINED") %>%
  mutate(IGM = ifelse(IGM == "POS_STRONG", "POS", IGM)) %>%
  dplyr::rename(IGHM = IGM) %>%
  group_by(ClusterAIC, IGHM) %>%
  dplyr::mutate(IGM = factor(IGHM, levels = c("NEG", "POS"))) %>%
  mutate(ClusterAIC = factor(ClusterAIC,
   levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
  dplyr::summarize(n = n()) %>%
  ggplot(aes(x = ClusterAIC, y = n, fill = IGHM, label = n)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(data = ann_text_ihc,
   aes(label = total, x = ClusterAIC, y = 1.0, label = label),
    vjust = -0.5, size = 1.4) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  ylab("Proportion") +
  xlab("Subtype") +
  theme(axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.margin = margin(-10, -10, -10, -10)) +
  ggtitle(element_blank()) +
  scale_fill_manual(values = c("#DCDCDC", "#a6dba0"))

p2 <- p2 + theme(plot.margin = margin(t = 0.6, l = 0.00,
                                      r = 0.2,
                                      b = 0.8, unit = "cm"))

#--
# Panel E - Plot the V gene identity in IgM and IgG-expressing cases
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
  left_join(clonotypes_igh_fl_tot) %>%
  mutate(rel.prop.dom = count / count.tot) %>%
  left_join(clusters[, c("SAMPLE_ID", "ClusterAIC")],
   by = c("sample_id" = "SAMPLE_ID")) %>%
  filter(ClusterAIC == factor(ClusterAIC,
   levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
  filter(!is.na(ClusterAIC)) %>%
  data.frame()

df2 %>% mutate(c_call = ifelse(c_call %in%
                         c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), "IGHG",
                        ifelse(c_call %in%
                         c("IGHA1", "IGHA2"), "IGHA", c_call))) %>%
  filter(c_call %in% c("IGHM", "IGHG")) %>%
  group_by(c_call) %>%
   summarize(mean(v_identity))

p3 <- df2 %>% mutate(c_call = ifelse(c_call %in%
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

#--
# Panel F - Plot the V gene identity
#--

my_comparisons <- list(
  c("GM", "CS"), c("GM", "TT"), c("GM", "Q"), c("GM", "AR"))

p4 <- df2 %>% mutate(ClusterAIC = factor(ClusterAIC,
 levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
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
  ylab("V gene identity (%)") +
  xlab("Subtype")

#----
# Save Figure 6 plots
#----

img_plot <- gridExtra::grid.arrange(im_A, im_B, p2,
 widths = c(1, 1, 0.95), ncol = 3, nrow = 1)

gg4 <- gridExtra::grid.arrange(p1,
                               gridExtra::arrangeGrob(img_plot, widths = 0.8),
                               gridExtra::arrangeGrob(p3, p4,
                                ncol = 2, heights = c(1), widths = c(0.8, 1)),
                               ncol = 1, heights = c(1, 1, 1))

dir <- "/figures/"

ggsave(file = paste0(dir, date, " Fig6.pdf"), gg4, width = 20,
       height = 20, units = "cm")
