#---
# FLOMICS manuscript Figure 7 
# Author: Robert Kridel & Victoria Shelton
#---

packages <- c("dplyr", "readr", "ggpubr", "ijtiff", "png", "forcats",
 "patchwork", "magick")
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
  dplyr::select(SAMPLE_ID, ClusterAIC, cohort)

# Read in sample annotation
sample_annotation <- read_excel("sample_annotations.xlsx")
fl_cases <- sample_annotation %>% 
    filter(TYPE == "FL") %>% 
    .$SAMPLE_ID

# Read in clinical data
clin <- read.csv("~/20231111_clinical_data.csv") %>%
  dplyr::mutate(LY_FL_ID = paste0(LY_FL_ID, "_T1")) %>%
  dplyr::select(LY_FL_ID, ANN_ARBOR_STAGE)

# Read in TRUST4 output
airr_tgl_ocir <- read_tsv("TGL_OICR_IMGT_airr_merged.tsv") %>%
  mutate(sample_id = ifelse(sample_id %in%
   c("LY_FL_234", "LY_FL_235",
   "LY_FL_236", "LY_FL_237", "LY_FL_238", "LY_FL_246", "LY_FL_247",
    "LY_FL_249", "LY_FL_252", "LY_FL_253", "LY_FL_254", "LY_FL_255",
     "LY_FL_256", "LY_FL_257", "LY_FL_258", "LY_FL_261"),
      paste0(sample_id, "_T1"), sample_id))

airr_e4402 <- read_tsv("E4402_IMGT_airr_merged.tsv")

# Combine and retain only those cases for which RNAseq data passed QC
airr <- rbind(airr_tgl_ocir, airr_e4402) %>%
 filter(sample_id %in% rnaseq_samples_passing_qc)

# Read in IHC information
ihc_results <- read_csv("~/Isotype_analysis_IHC.csv")

#--
# Identify most common isotype vs clusters
#--
clonotypes_igh_fl <- airr %>%
  filter(sample_id %in% fl_cases) %>%
  filter(grepl("IGH", v_call))

clonotypes_igh_fl_tot <- clonotypes_igh_fl %>%
  group_by(sample_id) %>%
  dplyr::summarize(count.tot = sum(consensus_count))

clonotypes_igh_fl_vdjgrouped <- clonotypes_igh_fl %>%
  mutate(v_call = sub("\\*.*", "", v_call)) %>%
  mutate(d_call = substr(d_call, 1, 5)) %>%
  mutate(j_call = substr(j_call, 1, 5)) %>%
  mutate(c_call = sub("\\*.*", "", c_call)) %>%
  group_by(sample_id, v_call, d_call, j_call, c_call) %>%
  dplyr::summarize(count = sum(consensus_count))

clonotypes_igh_fl_vdjgrouped_dominant_all <- clonotypes_igh_fl_vdjgrouped %>%
  group_by(sample_id) %>%
  slice_max(n = 1, order_by = count)

clonotypes_igh_fl_vdjgrouped_dominant <- clonotypes_igh_fl_vdjgrouped %>%
  group_by(sample_id) %>%
  slice_max(n = 1, order_by = count) %>%
  filter(!is.na(c_call)) %>%
  left_join(clonotypes_igh_fl_tot) %>%
  mutate(rel.prop.dom = count / count.tot)

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

#--
# Panel 7A 
#--

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

# annotation data frame:
ann_text <- data.frame(ClusterAIC = c("CS", "TT", "GM", "Q", "AR"),
                       isotype = rep(c("IGHG"), 5),
                       label = paste0("p = ", c(b1, b2, b3, b4, b5)),
                       n = rep(c(46), 5))
ann_text$ClusterAIC <- factor(ann_text$ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))
subtype_order <- c("CS", "TT", "GM", "Q", "AR")

plota <- clonotypes_igh_fl_vdjgrouped_dominant %>%
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
  arrange(factor(ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR")), factor(isotype, levels = c("IGHM", "IGHG", "IGHA"))) #%>%
plota$ClusterAIC <- factor(p1$ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))

plota <- plota %>% ggplot(aes(x = isotype, y = n, fill = isotype, label = n)) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  facet_wrap(~ClusterAIC, ncol = 5) +
  geom_text(data = ann_text, size = 2, aes(label = label, x = isotype, y = n-3)) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  ylab("# samples with given isotype") +
  theme(axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 5),
        axis.title.x = element_blank(), legend.position = "none",
        strip.background = element_blank(),
        plot.subtitle = element_blank(),
        strip.text.x = element_text(hjust = 0)) +
  ggtitle(element_blank()) +
  scale_fill_manual(values = c("#a6dba0", "#c2a5cf", "#9970ab"))

plota <- plota + theme(plot.margin = margin(t=0.0, l=0.2,
                                      r=0.2,
                                      b=0.8, unit = "cm"))


#--
# Panel 7B 
#--

ihc_results_e2408 <- ihc_results %>%
  mutate(OTHER_ID = as.character(OTHER_ID)) %>%
  dplyr::select(-SAMPLE_ID) %>%
  filter(INSTITUTION == "E2408") %>%
  left_join(sample_annotation[,c("SAMPLE_ID", "OTHER_ID")], by = "OTHER_ID") %>%
  mutate(TIME_POINT = substr(SAMPLE_ID, 12, 13)) %>%
  filter(TIME_POINT == "T1") %>%
  dplyr::select(SAMPLE_ID, OTHER_ID, INSTITUTION, IGM, IGG, COMMENT)

ihc_results <- ihc_results %>%
  filter(INSTITUTION == "UHN") %>%
  rbind(IHC_results_E2408) %>%
  mutate(IGM = ifelse(IGM == "-", "NEG",
                      ifelse(IGM == "+?", "UNDETERMINED",
                             ifelse(IGM == "+", "POS",
                                    ifelse(IGM == "++", "POS_STRONG", NA))))) %>%
  mutate(IGG = ifelse(IGG == "-", "NEG",
                      ifelse(IGG == "+?", "UNDETERMINED",
                             ifelse(IGG == "+", "POS", NA)))) %>%
  left_join(clonotypes_igh_fl_vdjgrouped_dominant, by = c("SAMPLE_ID" = "sample_id")) %>%
  left_join(clusters[,c("SAMPLE_ID", "ClusterAIC")], by = c("SAMPLE_ID"))

ihc_results_undetermined <- ihc_results %>%
  dplyr::filter(IGM == "UNDETERMINED") %>%
  nrow()

df_ihc <- ihc_results %>%
  filter(!is.na(ClusterAIC)) %>%
  #filter(IGM != "UNDETERMINED") %>% #going to include UNDETERMINED counts as well now
  mutate(IGM = ifelse(IGM == "POS_STRONG", "POS", IGM)) %>%
  group_by(ClusterAIC) %>%
  dplyr::summarize(n = n()) %>%
  tidyr::spread(ClusterAIC, n) %>%
  replace(is.na(.), 0)

plotb <- 
  ihc_results %>%
  filter(!is.na(ClusterAIC)) %>% #89 samples
  mutate(IGM = ifelse(IGM == "POS_STRONG", "POS", IGM)) %>%
  dplyr::rename(IGHM = IGM) %>%
  group_by(ClusterAIC, IGHM) %>%
  dplyr::mutate(IGM = factor(IGHM, levels = c("UNDETERMINED", "NEG", "POS"))) %>%
  mutate(ClusterAIC = factor(ClusterAIC, levels = c("CS", "TT", "GM", "Q", "AR"))) %>%
  dplyr::summarize(n = n())

plotb.other <- plotb %>%
  dplyr::filter(!ClusterAIC == 'AR') %>%
  dplyr::mutate(ClusterAIC = 'Other') %>%
  group_by(IGHM) %>% 
  dplyr::transmute(Total=sum(n)) %>%
  unique(.) %>%
  dplyr::mutate(ClusterAIC = 'Other') %>%
  dplyr::rename(n = Total)
plotb.ar <- plotb %>%
  dplyr::filter(ClusterAIC == 'AR')
plotb.bind <- rbind(plotb.ar, plotb.other) %>%
  dplyr::mutate(IGHM = factor(IGHM, levels = c("UNDETERMINED", "NEG", "POS")))
plotb.bind$IGHM <- relevel(plotb.bind$IGHM, "UNDETERMINED")

ann_text_ihc_bind <- data.frame(ClusterAIC = c("Other", "AR"),
                           total = paste0("N=", c(sum(df_ihc[,2:5]), sum(df_ihc$AR))),
                           IGHM = rep("NEG", 2),
                           n = rep(c(1.01), 2))
ann_text_ihc_bind$ClusterAIC <- factor(ann_text_ihc_bind$ClusterAIC, levels = c("Other", "AR"))
subtype_order <- c("Other", "AR")

plotb.1 <- plotb.bind %>% 
  ggplot(aes (x = factor(ClusterAIC, levels=c("Other", "AR")), y = n, fill = factor(IGHM, levels=c("UNDETERMINED", "NEG", "POS")), label = n))+
  geom_bar(position = "fill", stat = "identity") +
  geom_text(data = ann_text_ihc_bind, aes(label = total, x = ClusterAIC, y = 1.0), vjust = -0.5, size = 1.4) + #, label = label
  theme_bw() +
  theme(text = element_text(size = 10)) +
  ylab("Proportion") +
  xlab("Subtype") +
  theme(axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 5),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.margin = margin(-10, -10, -10, -10)) +
  ggtitle(element_blank()) +
  scale_fill_manual(breaks=c("UNDETERMINED", "NEG", "POS"), 
                    values=c("UNDETERMINED"="#f4f0ec", "NEG"="#DCDCDC", "POS"="#a6dba0"),
                    name = "IGHM")


#--
# Panel 7C (200x magnification)
#-- 

img <- image_read("IGM EAB2953 A6 200x_smaller.svg")

plotc <- ggplot() + background_image(img) + 
  theme(plot.margin = margin(t=0.5, 
                             l=0.4, 
                             r=0.7, b=0.9, unit = "cm"))


#--
# Panel 7D - Plot the V gene identity in IgM and IgG-expressing cases
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

plotd <- df2 %>% mutate(c_call = ifelse(c_call %in%
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
  ylab("V gene identity (%)") + xlab("Isotype") + ylim(70, 113)

plotd <- plotd + theme(plot.margin = margin(t=0.6, l=0.3,
                                          r=0.2,
                                          b=0.5, unit = "cm"))


#--
# Panel 7E - Plot the V gene identity
#--

my_comparisons <- list(
  c("GM", "CS"), c("GM", "TT"), c("GM", "Q"), c("GM", "AR"))

plote <- df2 %>% mutate(ClusterAIC = factor(ClusterAIC,
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
  xlab("Subtype") +
  ylim(70,113)

plote <- plote + theme(plot.margin = margin(t=0.6, l=0.3,
                                      r=0.2,
                                      b=0.5, unit = "cm"))
                                      

bar_plots <- gridExtra::grid.arrange(plota, plotb, ncol = 2, nrow = 1, widths = c(1.9, 1))
Vgene_plots <- gridExtra::arrangeGrob(plot2, plotd, plote, ncol = 3, heights = c(1), widths = c(1, 1, 1))
gg4 <- gridExtra::grid.arrange(bar_plots,
                               Vgene_plots,
                               ncol = 1, heights=c(1,1))

dir <- "~/figures/"

ggsave(file = paste0(dir, date, "_Fig7.pdf"), gg4, width = 18,
       height = 13, units = "cm")

