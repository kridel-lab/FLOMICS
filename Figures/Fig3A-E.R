#---
# This script generates plots for Figures 4, panels A to E
# Author: Robert Kridel
#---

packages <- c("dplyr", "ggplot2", "survival", "survminer")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

setwd("~/your working directory/FLOMICS/")

#--
# Read in data
#--

# Read in clinical data
clinical <- read.csv("metadata/20220525_clinical_data.csv") %>%
  mutate(GRADE_1.2_3a = ifelse(GRADE_1.2_3a == "", NA, GRADE_1.2_3a)) %>%
  mutate(TTP = as.numeric(TTP)) %>%
  mutate(survival_cohort = ifelse(LIMITED_STAGE_PROJECT_INCLUDE == "YES",
   "LIMITED",
                           ifelse(PRIM_TX_CAT %in% c("BR", "R-CHOP",
                            "R-CHOP + RAD", "R-CVP"), "ADVANCED_R_CHEMO", NA)))

# Read in clusters
clusters <- read.csv("GMM_Cluster_Labels_flexmix_clusters.csv")

# Combine clinical data with clusters
clusters_clinical <- clusters %>%
  mutate(PATIENT_ID = stringr::str_remove(SAMPLE_ID, "_T1")) %>%
  left_join(clinical, by = c("PATIENT_ID" = "LY_FL_ID"))

# Redefine this data frame, but just for R-chemo treated,
## advanced stage patients
clusters_clinical_r_chemo <- clusters_clinical %>%
  filter(survival_cohort == "ADVANCED_R_CHEMO" & ANN_ARBOR_STAGE %in% c(3, 4))

#--
# Generate Figure 3, panels A to C
#--

p_grade <- clusters_clinical %>%
  group_by(ClusterAIC, GRADE_1.2_3a) %>%
  filter(!is.na(GRADE_1.2_3a) & !GRADE_1.2_3a == "") %>%
  filter(!is.na(ClusterAIC)) %>%
  mutate(GRADE_1.2_3a = factor(GRADE_1.2_3a, levels = c("3A", "1 to 2"))) %>%
  summarize(count = n()) %>%
  mutate(percent = (count / sum(count)),
   label = ifelse(GRADE_1.2_3a == "1 to 2", paste0("N=", sum(count)), "")) %>%
  ggplot(aes(x = ClusterAIC, y = count, fill = GRADE_1.2_3a)) +
  geom_col(position = "fill") +
  geom_text(aes(y = 1.0, label = label), vjust = -0.5, size = 1) +
  scale_fill_manual(values = c("#f1a340", "#998ec3")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.margin = margin(-10, -10, -10, -10)) +
  guides(fill = guide_legend(title = "Grade")) + ylab("Proportion")

p_14_18 <- clusters_clinical %>%
  group_by(ClusterAIC, T_14_18) %>%
  filter(!is.na(T_14_18)) %>%
  filter(!is.na(ClusterAIC)) %>%
  mutate(T_14_18 = factor(T_14_18, levels = c("0", "1"))) %>%
  summarize(count = n()) %>%
  mutate(T_14_18 = ifelse(T_14_18 == "0", "no",
                          ifelse(T_14_18 == "1", "yes", T_14_18))) %>%
  mutate(T_14_18 = factor(T_14_18, levels = c("no", "yes"))) %>%
  mutate(percent = (count / sum(count)),
   label = ifelse(T_14_18 == "yes", paste0("N=", sum(count)), "")) %>%
  ggplot(aes(x = ClusterAIC, y = count, fill = T_14_18)) +
  geom_col(position = "fill") +
  geom_text(aes(y = 1.0, label = label), vjust = -0.5, size = 1) +
  scale_fill_manual(values = c("#ef8a62", "#999999")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.margin = margin(-10, -10, -10, -10)) +
  guides(fill = guide_legend(title = "t(14;18)")) + ylab("Proportion")
p_14_18

p_stage <- clusters_clinical %>%
  group_by(ClusterAIC, ANN_ARBOR_STAGE) %>%
  filter(!is.na(ANN_ARBOR_STAGE)) %>%
  filter(!is.na(ClusterAIC)) %>%
  mutate(ANN_ARBOR_STAGE = factor(ANN_ARBOR_STAGE,
   levels = c("4", "3", "2", "1"))) %>%
  summarize(count = n()) %>%
  mutate(percent = (count / sum(count)),
   label = ifelse(ANN_ARBOR_STAGE == "4", paste0("N=", sum(count)), "")) %>%
  ggplot(aes(x = ClusterAIC, y = count, fill = ANN_ARBOR_STAGE)) +
  geom_text(aes(y = 1.0, label = label), vjust = -0.5, size = 1) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("#d6604d", "#f4a582", "#92c5de", "#4393c3")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.margin = margin(-10, -10, -10, -10)) +
  guides(fill = guide_legend(title = "Ann Arbor stage")) + ylab("Proportion")
p_stage

g <- gridExtra::arrangeGrob(p_grade, p_14_18, p_stage, nrow = 1,
 widths = c(1, 1, 1.17))
ggsave(file = paste0("img/",  date, " Fig3A_C.pdf"), g, width = 16, height = 5,
 units = "cm")

#--
# Generate Figure 3, panels D-E
#--

surv_plots <- list()

fit <- survfit(Surv(TTP, CODE_TTP) ~ ClusterAIC, data = clusters_clinical)
names(fit$strata) <- gsub("ClusterAIC=", "", names(fit$strata))
surv_plots[[1]] <- ggsurvplot(fit, data = clusters_clinical,
                              xlim = c(0, 10), xlab = "Time (years)",
                               break.time.by = 2,
                              size = 0.5, censor.size = 2.5,
                              pval = TRUE, pval.size = 3,
                               pval.coord = c(0.5, 0.05),
                              risk.table = TRUE, fontsize = 2, legend = "none",
                              title = "TTP - all patients",
                               ggtheme = theme_classic2(base_size = 7))

fit <- survfit(Surv(TTP, CODE_TTP) ~ ClusterAIC,
 data = clusters_clinical_R_CHEMO)
names(fit$strata) <- gsub("ClusterAIC=", "", names(fit$strata))
surv_plots[[2]] <- ggsurvplot(fit, data = clusters_clinical_R_CHEMO,
                              xlim = c(0, 10), xlab = "Time (years)",
                               break.time.by = 2,
                              size = 0.5, censor.size = 2.5,
                              pval = TRUE, pval.size = 3,
                               pval.coord = c(0.5, 0.05),
                              risk.table = TRUE, fontsize = 2, legend = "none",
                              title = "TTP - immunochemotherapy",
                               ggtheme = theme_classic2(base_size = 7))

res <- arrange_ggsurvplots(surv_plots, print = TRUE, ncol = 2, nrow = 1,
 risk.table.height = 0.25)
dev.off()
ggsave(paste0("img/",  date, " Fig3D_E.pdf"), res, width = 15, height = 10,
 units = "cm")
