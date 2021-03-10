
library(dplyr)
library(ggpubr)

TRA <- read.csv("Analysis-Files/Mixcr/Arash_script_output/divstats_TRA.csv") %>%
  mutate(clonotype = "TRA") %>%
  mutate(X = gsub("\\..*", "", X)) %>%
  rename(SAMPLE_ID = X) %>%
  mutate(tmp = nchar(SAMPLE_ID)) %>%
  mutate(SAMPLE_ID = ifelse(tmp == 9, paste0(SAMPLE_ID, "_T1"), SAMPLE_ID)) %>%
  mutate(tmp = nchar(SAMPLE_ID)) %>%
  mutate(TIME_POINT = ifelse(tmp == 12, substr(SAMPLE_ID, 11, 12), NA)) %>%
  select(-tmp) %>%
  filter(TIME_POINT %in% c(NA, "T1"))

TRB <- read.csv("Analysis-Files/Mixcr/Arash_script_output/divstats_TRB.csv") %>%
  mutate(clonotype = "TRB") %>%
  mutate(X = gsub("\\..*", "", X)) %>%
  rename(SAMPLE_ID = X) %>%
  mutate(tmp = nchar(SAMPLE_ID)) %>%
  mutate(SAMPLE_ID = ifelse(tmp == 9, paste0(SAMPLE_ID, "_T1"), SAMPLE_ID)) %>%
  mutate(tmp = nchar(SAMPLE_ID)) %>%
  mutate(TIME_POINT = ifelse(tmp == 12, substr(SAMPLE_ID, 11, 12), NA)) %>%
  select(-tmp) %>%
  filter(TIME_POINT %in% c(NA, "T1"))

clin <- read.csv("metadata/clinical_data_rcd11Aug2020.csv") %>%
  mutate(LY_FL_ID = paste0(LY_FL_ID, "_T1")) %>%
  mutate(POD24 = ifelse(CODE_TTP == 1 & TTP < 2, "YES", "NO")) %>%
  mutate(POD24_stage = paste0("POD24_", POD24, "_", TYPE))

SNF.clust <- read.csv("Cluster Labels/InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv") %>%
  select(ID, SNFClust = SNFClust10Feb2021)

sample.annotation <- read.csv("metadata/sample_annotations_rcd6Nov2019.csv")

TRA <- TRA %>%
  left_join(sample.annotation[, c("SAMPLE_ID", "TYPE", "STAGE")]) %>%
  left_join(SNF.clust, by = c("SAMPLE_ID" = "ID")) %>%
  filter(Sample_Coverage > 0.75)

TRB <- TRB %>%
  left_join(sample.annotation[, c("SAMPLE_ID", "TYPE", "STAGE")]) %>%
  left_join(SNF.clust, by = c("SAMPLE_ID" = "ID")) %>%
  filter(Sample_Coverage > 0.75)

##
# Plot Shannon diversity
##
p <- TRA %>% filter(TYPE %in% c("DLBCL", "FL")) %>%
  ggboxplot(x = "TYPE", y = "observed_Shannon", color = "TYPE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- TRB %>% filter(TYPE %in% c("DLBCL", "FL")) %>%
  ggboxplot(x = "TYPE", y = "observed_Shannon", color = "TYPE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- TRA %>% filter(STAGE %in% c("LIMITED", "ADVANCED")) %>%
  ggboxplot(x = "STAGE", y = "observed_Shannon", color = "STAGE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- TRB %>% filter(STAGE %in% c("LIMITED", "ADVANCED")) %>%
  ggboxplot(x = "STAGE", y = "observed_Shannon", color = "STAGE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- TRA %>% filter(SNFClust %in% c(1, 2)) %>%
  ggboxplot(x = "SNFClust", y = "observed_Shannon", color = "SNFClust", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- TRB %>% filter(SNFClust %in% c(1, 2)) %>%
  ggboxplot(x = "SNFClust", y = "observed_Shannon", color = "SNFClust", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

