
#---
# This script reads in T2 variant (capseq) data and
## performs clustering prediction
# Authors: Robert Kridel + Victoria Shelton
# Last modified: March 21st, 2023 
#---

#--
# Load in packages
#--
packages <- c("dplyr", "ggplot2", "flexmix", "tidyverse", "gridExtra")
lapply(packages, require, character.only = TRUE)

packages2 <- c("magrittr", "cluster", "cluster.datasets", "cowplot", "NbClust",
 "clValid", "ggfortify", "clustree", "dendextend", "factoextra", "FactoMineR",
  "corrplot", "GGally", "ggiraphExtra", "knitr", "kableExtra", "taRifx",
   "survminer", "tibble", "tidyr", "ggalluvial", "survival", "stringr")
lapply(packages2, require, character.only = TRUE)

#--
# Setup working space
#--

date <- Sys.Date()
setwd("~/your working directory/T2_Predictions/")

#--
# Assign T2 samples to clusters
#--

#list those samples determined to not be T2 and will be removed from the dataset
t2_remove <- c("LY_FL_460_T2", "LY_FL_176_T2", "LY_FL_470_T2", "LY_FL_399_T2",
 "LY_FL_449_T2", "LY_FL_432_T2")

capseq_plosmed_t2 <- read.table("PLOSMED_exonic_filtered_MUTECT2_calls.txt",
 sep = ";", header = TRUE) %>%
  mutate(TIMEPOINT = substr(SAMPLE_ID, 11, 12)) %>%
  filter(TIMEPOINT == "T2") %>%
  filter(TYPE == "TFL") %>%
  filter(!SAMPLE_ID %in% T2_remove) %>%
  dplyr::distinct(SAMPLE_ID, Hugo_Symbol) %>%
  mutate(var = 1) %>%
  reshape2::dcast(Hugo_Symbol ~ SAMPLE_ID) %>%
  replace(is.na(.), 0)

plosmed_t2_genes <- as.data.frame(capseq_plosmed_t2$Hugo_Symbol) %>%
  dplyr::rename(Hugo_Symbol = `capseq_plosmed_t2$Hugo_Symbol`) #88 genes

genes_common_t2 <- as.data.frame(genes_common$Gene.Name) %>%
  dplyr::rename(Hugo_Symbol = `genes_common$Gene.Name`)

genes_to_add_plosmed_t2 <- dplyr::setdiff(genes_common_t2, plosmed_t2_genes)
#0 genes to add
dplyr::setdiff(genes_common$Gene.Name,
 as.data.frame(capseq_plosmed_t2$Hugo_Symbol))
# no genes within common genes targeted that do not have variants

capseq_plosmed_t2 <- capseq_plosmed_t2 %>%
  filter(Hugo_Symbol %in% genes_common$Gene.Name)

capseq_plosmed_t2 <- as.matrix(capseq_plosmed_t2)
dim(capseq_plosmed_t2) #146 samples, 57 genes (housed in the first column)

capseq_uhn_t2 <- read.table("UHN_exonic_filtered_MUTECT2_calls.txt", sep = ";",
 header = TRUE) %>%
  mutate(TIMEPOINT = substr(Library, 11, 12)) %>%
  filter(TIMEPOINT == "T2") %>%
  filter(!Library %in% T2_remove) %>%
  dplyr::distinct(Library, Hugo_Symbol) %>%
  mutate(var = 1) %>%
  reshape2::dcast(Hugo_Symbol ~ Library) %>%
  replace(is.na(.), 0)
uhn_t2_genes <- as.data.frame(capseq_uhn_t2$Hugo_Symbol) %>%
  dplyr::rename(Hugo_Symbol = `capseq_uhn_t2$Hugo_Symbol`)
genes_common_t2 <- as.data.frame(genes_common$Gene.Name) %>%
  dplyr::rename(Hugo_Symbol = `genes_common$Gene.Name`)
genes_to_add_uhn_t2 <- dplyr::setdiff(genes_common_t2, uhn_t2_genes)
#genes within common genes that do not have variants
genes_to_remove_uhn_t2 <- dplyr::setdiff(uhn_t2_genes, genes_common_t2)
#gene within UHN that is targeted but not common

capseq_uhn_t2_2 <- capseq_uhn_t2 %>%
  filter(!Hugo_Symbol %in% genes_to_remove_uhn_t2) #12 genes
capseq_uhn_t2_2 <- bind_rows(capseq_uhn_t2_2, genes_to_add_uhn_t2) #57 genes now
capseq_uhn_t2_2[is.na(capseq_uhn_t2_2)] <- 0
capseq_uhn_t2_2 <- as.matrix(capseq_uhn_t2_2)
dim(capseq_uhn_t2_2) #57 genes

df_t2 <- full_join(as.data.frame(capseq_plosmed_t2),
 as.data.frame(capseq_uhn_t2_2))
dim(df_t2) #148 samples, 57 genes
remove_dup <- c("LY_FL_1156_T2", "LY_FL_1135_T2")
remove_dup %in% colnames(df_t2)

row.names(df_t2) <- df_t2$Hugo_Symbol
df_t2$Hugo_Symbol <- NULL
df_t2_2 <- as.matrix(df_t2)
df_t2_3 <- t(df_t2_2) %>% data.frame() %>% mutate(SAMPLE_ID = row.names(.))
df_t2_4 <- cbind(SAMPLE_ID = df_t2_3$SAMPLE_ID,
 df_t2_3[, c(1:(ncol(df_t2_3) - 1))])
row.names(df_t2_4) <- NULL
df_t2_4 <- df_t2_4[, -1]
df_t2_4 <- sapply(df_t2_4, as.numeric)
dim(df_t2_4) #148 57
df_t2_archive <- df_t2_4

row.names(df_t2_3) <- NULL
df_t2_all <- df_t2_3
df_t2_df <- df_t2_all[, -58]
df_t2_df <- sapply(df_t2_df, as.numeric)
dim(df_t2_df) #148 57

#--
# Make cluster predictions
#--

predict_clusters <- function(model, newdata) {
  muts_df <<- newdata
  probs <- flexmix::clusters(model, newdata = data.frame(newdata))
  probs
}

predict_posteriors <- function(model, newdata) {
  muts_df_pos <<- newdata
  probs <- flexmix::posterior(model, newdata = data.frame(newdata))
  probs
}

#--
# Store cluster predictions
#--

#import aic model
aic <- readRDS("most_cofident_clustering_run.Rdata")

predictions <- predict_clusters(aic, df_t2_4)
pos_pred <- predict_posteriors(aic, df_t2_4)

df_t2_4 <- cbind(LY_FL_ID = row.names(df_t2_3), predictions)
df_t2_4 <- data.frame(df_t2_4)
df_t2_4$PATIENT_ID <- substr(df_t2_4$LY_FL_ID, 1, 9)
df_t2_4 <- df_t2_4 %>%
 select(PATIENT_ID, T2_prediction = predictions) %>%
 dplyr::rename(LY_FL_ID = PATIENT_ID)

df_t2_all$ClusterAIC <- as.factor(paste0("C", predictions))

genes <- colnames(df_t2_df)
criteria <- "AIC"
i <- "T2"
df_t2_aic <- heatmap_mutation_extended(df = df_t2_all,
 genes, group_var = "ClusterAIC", y_order = "group", idcol = "SAMPLE_ID")
grid.arrange(df_t2_aic)

ggsave(file = paste0(date, "_t2_heatmap_aic.png"),
 df_t2_aic, width = 25, height = 25, units = "cm")

aic_assignments <- read.csv("GMM_Cluster_Labels_flexmix_clusters.csv") %>%
  dplyr::select(SAMPLE_ID, ClusterAIC, cohort) %>%
  dplyr::rename(LY_FL_ID = SAMPLE_ID) %>%
  mutate(LY_FL_ID = trimws(LY_FL_ID, which = "right", whitespace = "_T1"))

muts_all_t1_t2_inner <- aic_assignments %>%
  inner_join(df_t2_4, by = c("LY_FL_ID")) %>%
  filter(!is.na(T2_prediction)) #116 samples in total
dim(muts_all_t1_t2_inner)
print("T2 prediction count table")
table(muts_all_t1_t2_inner$ClusterAIC, muts_all_t1_t2_inner$T2_prediction)

print("T2 prediction proportion table")
round(prop.table(table(muts_all_t1_t2_inner$ClusterAIC,
 muts_all_t1_t2_inner$T2_prediction), margin = 1), 2) * 100

print("T2 prediction Chi-squared")
chisq.test(muts_all_T1_T2_inner$ClusterAIC, muts_all_T1_T2_inner$cohort)

muts_all_t1_t2_inner <- muts_all_t1_t2_inner %>%
  dplyr::rename(T1.AIC = ClusterAIC, T2.AIC = T2_prediction) %>%
  select(1, 2, 4, 3) %>%
  mutate(T1.AIC = str_replace(T1.AIC, "C1", "1")) %>%
  mutate(T1.AIC = str_replace(T1.AIC, "C2", "2")) %>%
  mutate(T1.AIC = str_replace(T1.AIC, "C3", "3")) %>%
  mutate(T1.AIC = str_replace(T1.AIC, "C4", "4")) %>%
  mutate(T1.AIC = str_replace(T1.AIC, "C5", "5"))

muts_all_t1_t2 <- aic_assignments %>%
  full_join(df_t2_4, by = c("LY_FL_ID"))
dim(muts_all_t1_t2) #745 samples

muts_all_t1_t2 <- muts_all_t1_t2 %>%
  dplyr::rename(T1.AIC = ClusterAIC, T2.AIC = T2_prediction) %>%
  select(1, 2, 4, 3) %>%
  mutate(T1.AIC = str_replace(T1.AIC, "C1", "1")) %>%
  mutate(T1.AIC = str_replace(T1.AIC, "C2", "2")) %>%
  mutate(T1.AIC = str_replace(T1.AIC, "C3", "3")) %>%
  mutate(T1.AIC = str_replace(T1.AIC, "C4", "4")) %>%
  mutate(T1.AIC = str_replace(T1.AIC, "C5", "5"))

dir <- "/T2_Predictions/"
write.csv(df_t2_4,
 file = paste0(dir, date, "_T2_Cluster_Labels_flexmix_clusters.csv"),
  row.names = FALSE)

write.csv(muts_all_t1_t2,
 file = paste0(dir, date, "_T1vsT2_mat_fulljoin.csv"), row.names = FALSE)

write.csv(muts_all_t1_t2_inner,
 file = paste0(dir, date, "_T1vsT2_mat.csv"), row.names = FALSE)
#116 samples

# read in table again
muts_all_t1_t2_inner <- read.csv("/T2_Predictions/T1vsT2_mat.csv")

### creating the chart one more time
table(muts_all_t1_t2_inner$T1.AIC, muts_all_t1_t2_inner$T2.AIC)

#--
# Making alluvial plot, Fig 2F
#--

t1_t2_alluvial <- muts_all_t1_t2_inner %>%
  group_by(T1.AIC, T2.AIC) %>%
  dplyr::summarize(count = n()) %>%
  ggplot(aes(axis1 = T1.AIC, axis2 = T2.AIC, y = count)) +
  scale_x_discrete(limits = c("FL", "tFL"), expand = c(.2, .05)) +
  xlab("Demographic") +
  ylab("Number of Samples") +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 10)) +
  geom_alluvium(aes(fill = T1.AIC)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme(legend.position = "none", text = element_text(size = 24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(),
                axis.line = element_line(colour = "black"))

ggsave(file = paste0(date, "_T2_prediction_aic.png"),
 width = 25, height = 25, units = "cm")


#--
# Making survival plot, Fig 2G
#--

clinical_data <- read.csv("/sample annotation/20220525_clinical_data.csv") %>%
  dplyr::select(LY_FL_ID, CODE_TTP, CODE_PFS, CODE_TRANSF, CODE_OS, CODE_DSS,
   TTT)
t1_aic_assignments <- read.csv("GMM_Cluster_Labels_flexmix_clusters.csv") %>%
  dplyr::select(SAMPLE_ID, ClusterAIC, cohort) %>%
  dplyr::rename(LY_FL_ID = SAMPLE_ID, T1.AIC = ClusterAIC) %>%
  mutate(LY_FL_ID = trimws(LY_FL_ID, which = "right", whitespace = "_T1"))
t2_aic_assignments <- read.csv("T2_Cluster_Labels_flexmix_clusters.csv") %>%
  dplyr::rename(T2.AIC = T2_prediction) %>%
  mutate(T2.AIC = ifelse(T2.AIC == "1", "C1",
                         ifelse(T2.AIC == "2", "C2",
                                ifelse(T2.AIC == "3", "C3",
                                       ifelse(T2.AIC == "4", "C4",
                                              ifelse(T2.AIC == "5", "C5",
                                               T2.AIC))))))

t2_remove <- c("LY_FL_460", "LY_FL_176", "LY_FL_470", "LY_FL_399", "LY_FL_449",
 "LY_FL_432")

t1_t2_comp <- inner_join(t1_aic_assignments, t2_aic_assignments) %>%
  mutate(T1_equals_T2 = ifelse(T1.AIC == T2.AIC, 1,
                               ifelse(T1.AIC != T2.AIC, 0, NaN))) %>%
  relocate(cohort, .after = last_col()) %>%
  filter(!LY_FL_ID %in% t2_remove) # 116 samples
dim(t1_t2_comp)

write.csv(t1_t2_comp,
 file = paste0(dir, date, "_T1vsT2_mat.csv"), row.names = FALSE)

t1_t2_comp_full_join <- full_join(t1_aic_assignments, t2_aic_assignments) %>%
  mutate(T1_equals_T2 = ifelse(T1.AIC == T2.AIC, 1,
                               ifelse(T1.AIC != T2.AIC, 0, NaN))) %>%
  relocate(cohort, .after = last_col())
dim(t1_t2_comp_full_join)

setdiff(muts_all_t1_t2$LY_FL_ID, t1_t2_comp_full_join$LY_FL_ID)
#no difference between these two lists

write.csv(t1_t2_comp_full_join,
 file = paste0(dir, date, "_T1vsT2_mat_fulljoin.csv"), row.names = FALSE)

t1_t2_comp_ttt <- left_join(t1_t2_comp, clinical_data) %>%
  relocate(cohort, .after = last_col())
dim(t1_t2_comp_ttt)

t1_t2_comp_ttt_new <- t1_t2_comp_ttt
t1_t2_comp_ttt_new$T1_equals_T2 <- factor(t1_t2_comp_ttt_new$T1_equals_T2,
                         levels = c("1", "0"))

km_t1_t2_fit_transf <- survival::survfit(Surv(TTT, CODE_TRANSF) ~ T1_equals_T2,
 data = t1_t2_comp_ttt_new)

ttt_fit_transf <- ggsurvplot(km_t1_t2_fit_transf, data = t1_t2_comp_ttt,
                             fun = "event",
                             conf.int = FALSE,
                             pval = TRUE,
                             pval.coord = c(0, 0.94),
                             pval.method = TRUE,
                             pval.method.coord = c(0, 1),
                             xlim = c(0, 10),
                             break.x.by = 2,
                             xlab = "Time (years)",
                             legend.labs = c("Same Subtype",
                              "Switched Subtype"),
                             legend.title = "",
                             legend = c(0.7, 0.2),
                             pval.size = 8,
                             ggtheme = theme_classic2(base_size = 24),
                             font.title = c(19, "bold"),
                             title = "Time To Transformation"
                             )


pdf("ttt_fit_trasnf.pdf",
    width = 8, onefile = FALSE)
ttt_fit_trasnf
dev.off()
