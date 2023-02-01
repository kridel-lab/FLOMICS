
#---
# This script reads in mutation (capseq) data and performs clustering
# Authors: Robert Kridel + Victoria Shelton
# Last modified: January 17, 2022
# Tested on R/4.1.0
#---

#--
# Load in packages
#--
packages <- c("dplyr", "ggplot2", "flexmix", "tidyverse", "gridExtra")
lapply(packages, require, character.only = TRUE)

packages2 <- c("magrittr", "cluster", "cluster.datasets", "cowplot",
               "ggfortify", "dendextend", "factoextra", "FactoMineR",
               "GGally", "ggiraphExtra", "knitr", "kableExtra", "taRifx")
lapply(packages2, require, character.only = TRUE)

#--
# Setup working space
#--

date <- Sys.Date()
setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/")
dir <- "/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/"

#----------------------------------
# Read in data
#----------------------------------

### Read in gene panels
genes_PLOSMED <- read.csv("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/PROBES/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
  filter(PLOS_MED_PANEL == "YES") %>% dplyr::select(Gene.Name) # n=86
print("this is genes_PLOSMED")
print(genes_PLOSMED)

genes_UHN <- read.csv("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/PROBES/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
  filter(FLOMICS_PANEL == "YES") %>% dplyr::select(Gene.Name) # n=71
print("this is genes_UHN")
print(genes_UHN)

### Determine the common genes between the two utilized gene panels
genes_common <- intersect(genes_PLOSMED, genes_UHN) # n=57
print("these are the common genes between genes_PLOSMED and genes_UHN")
print(genes_common)

### Read in UHN filtered exonic SNV calls
capseq_UHN <- read.table("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/pre-processed-BC_CAPSEQ/2021-10-22_BC_CAPSEQ_exonic_filtered_MUTECT2_calls_with_splice.txt", sep = ";", header = T) %>%
  dplyr::rename(External_identifier = Library) %>%
  mutate(TIMEPOINT = substr(External_identifier, 11, 12)) %>%
  filter(TIMEPOINT != "T2") %>%
  filter(External_identifier != "LY_FL_571_T1") %>% # LY_FL_571_T1 and LY_FL_572_T1 from same pt
  dplyr::distinct(External_identifier, Hugo_Symbol) %>%
  mutate(var = 1) %>%
  reshape2::dcast(Hugo_Symbol ~ External_identifier) %>%
  dplyr::mutate(LY_FL_136_T1 = NA, LY_FL_179_T1 = NA, LY_FL_524_T1 = NA) # 3 samples did not return SNV calls from Mutect2 analysis
print("created capseq_UHN")
#print(capseq_UHN)
#print(capseq_UHN$Hugo_Symbol)
#print(genes_common$Gene.Name)
dplyr::setdiff(genes_common$Gene.Name, capseq_UHN$Hugo_Symbol) # CD83 within common genes targeted that are not mutated
capseq_UHN <- rbind(capseq_UHN, c("CD83", rep(NA, ncol(capseq_UHN)-1))) %>% replace(is.na(.), 0)
row.names(capseq_UHN) <- capseq_UHN$Hugo_Symbol
capseq_UHN$Hugo_Symbol <- NULL
capseq_UHN_mat <- as.matrix(capseq_UHN)

### Read in E4402 filtered exonic SNV calls
capseq_E4402 <- read.table("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/E4402/2021-10-17_E4402_exonic_filtered_MUTECT2_calls.txt", sep = ";", header = T) %>%
  mutate(PATIENT = substr(SAMPLE_ID, 1, 9)) %>%
  dplyr::distinct(SAMPLE_ID, Hugo_Symbol) %>%
  mutate(var = 1) %>%
  reshape2::dcast(Hugo_Symbol ~ SAMPLE_ID) %>%
  dplyr::mutate(LY_FL_1026_T1 = NA) %>% # 1 sample did not return SNV calls from Mutect2 analysis
  replace(is.na(.), 0)
print("created capseq_E4402")
dplyr::setdiff(genes_common$Gene.Name, capseq_E4402$Hugo_Symbol) # no genes within common genes are not mutated within the E4402 cohort
row.names(capseq_E4402) <- capseq_E4402$Hugo_Symbol
capseq_E4402$Hugo_Symbol <- NULL
capseq_E4402_mat <- as.matrix(capseq_E4402)

### Read in E2408 filtered exonic SNV calls
capseq_E2408 <- read.table("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/E2408-Tumor-Only/2021-11-29_E2408-Tumor-Only_exonic_filtered_MUTECT2_calls.txt", sep = ";", header = T) %>% # this is the path for the E2408-Tumor-Only analysis
#capseq_E2408 <- read.table("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/E2408-Tumor-Normal/MUTECT2/2021-12-06_E2408-Tumor-Normal_exonic_filtered_MUTECT2_calls.txt", sep = ";", header = T) %>% # this is the path for the E2408-Tumor-Normal analysis
  dplyr::distinct(Library, Hugo_Symbol) %>%
  mutate(var = 1) %>%
  reshape2::dcast(Hugo_Symbol ~ Library) %>%
  dplyr::mutate(LY_FL_1095_T1 = NA, LY_FL_1096_T1 = NA, LY_FL_1111_T1 = NA, LY_FL_1118_T1 = NA, LY_FL_1119_T1 = NA, # 13 samples did not return SNV calls from Mutect2 analysis
         LY_FL_1122_T1 = NA, LY_FL_1126_T1 = NA, LY_FL_1129_T1 = NA, LY_FL_1132_T1 = NA, LY_FL_1136_T1 = NA,
         LY_FL_1162_T1 = NA, LY_FL_1167_T1 = NA, LY_FL_1174_T1 = NA)
print("created capseq_E2408")
dplyr::setdiff(genes_common$Gene.Name, capseq_E2408$Hugo_Symbol) # 6 genes within common genes are not mutated within the E2408 Tumor Only cohort
capseq_E2408 <- rbind(capseq_E2408,
                      c("BTG2", rep(NA, ncol(capseq_E2408)-1)),
                      c("CD58", rep(NA, ncol(capseq_E2408)-1)),
                      c("HIST1H2AM", rep(NA, ncol(capseq_E2408)-1)),
                      c("HLA-DMB", rep(NA, ncol(capseq_E2408)-1)),
                      c("MYD88", rep(NA, ncol(capseq_E2408)-1)),
                      c("SGK1", rep(NA, ncol(capseq_E2408)-1))) %>%
  replace(is.na(.), 0)
row.names(capseq_E2408) <- capseq_E2408$Hugo_Symbol
capseq_E2408$Hugo_Symbol <- NULL
capseq_E2408_mat <- as.matrix(capseq_E2408)

### Read in PLOSMED filtered exonic SNV calls
capseq_PLOSMED <- read.table("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/PLOSMED/2021-10-17_PLOSMED_exonic_filtered_MUTECT2_calls.txt", sep = ";", header = T) %>%
  mutate(TIMEPOINT = substr(SAMPLE_ID, 11, 12)) %>%
  filter(TIMEPOINT != "T2") %>%
  dplyr::distinct(SAMPLE_ID, Hugo_Symbol) %>%
  mutate(var = 1) %>%
  reshape2::dcast(Hugo_Symbol ~ SAMPLE_ID) %>%
  dplyr::mutate(LY_FL_679_T1 = NA) # 1 sample did not return SNV calls from Mutect2 analysis
print("created capseq_PLOSMED")
dplyr::setdiff(genes_common$Gene.Name, capseq_PLOSMED$Hugo_Symbol) # 1 gene within common genes are not mutated within the PLOSMED cohort
capseq_PLOSMED <- rbind(capseq_PLOSMED, c("IL4R", rep(NA, ncol(capseq_PLOSMED)-1))) %>% replace(is.na(.), 0)
row.names(capseq_PLOSMED) <- capseq_PLOSMED$Hugo_Symbol
capseq_PLOSMED$Hugo_Symbol <- NULL
capseq_PLOSMED_mat <- as.matrix(capseq_PLOSMED)

#-------------------
#Merge datasets
#-------------------

genes_common_vec <- unlist(genes_common)
capseq_merged_mat <- cbind(capseq_UHN_mat[genes_common_vec,],
                       capseq_E4402_mat[genes_common_vec,],
                       capseq_E2408_mat[genes_common_vec,],
                       capseq_PLOSMED_mat[genes_common_vec,])
print("created capseq_merged_mat")
#print(capseq_merged_mat)

capseq_merged_df <- cbind(capseq_UHN[genes_common_vec,],
                           capseq_E4402[genes_common_vec,],
                           capseq_E2408[genes_common_vec,],
                           capseq_PLOSMED[genes_common_vec,])
print("created capseq_merged_df")
#print(capseq_merged_df)

### Transfer dataframe to new object,
## and remove the extra Hugo_Symbol column in this second dataframe
capseq_merged_no_Hugo <- capseq_merged_df
capseq_merged_no_Hugo$Hugo_Symbol <- NULL #By making Hugo_Symbol null,
                                          #we lose the ability to create a frequency plot downstream later,
                                          #so need to transfer this to a new object

### Examine the merged dataframe
print("Dimensions of the merged dataset including all cohorts")
dim(capseq_merged_df) # 57 x 715
print("Glimpse of the merged dataset")
glimpse(capseq_merged_df)

### Change the class of values in merged dataset from "character" to "num" to input into downstream analyses
capseq_merg_num <- japply(capseq_merged_df, which(sapply(capseq_merged_df, class)=="character"), as.numeric )
#summary(capseq_merg_num) %>% kable() %>% kable_styling()

### Writing out the merged dataframe
capseq_merg_txt = paste0(dir, date, "_gene_vs_sample_SNV_matrix.txt")
capseq_merg_csv = paste0(dir, date, "_gene_vs_sample_SNV_matrix.csv")
write.table(capseq_merg_num, file=capseq_merg_txt, quote=F, row.names=T, sep=";")
write.csv(capseq_merg_num, file = capseq_merg_csv, row.names = FALSE)

### Writing out the merged matrix
capseq_merg_mat_txt = paste0(dir, date, "_gene_vs_sample_SNV_mat_matrix.txt")
capseq_merg_mat_csv = paste0(dir, date, "_gene_vs_sample_SNV_mat_matrix.csv")
write.table(capseq_merged_mat, file=capseq_merg_mat_txt, quote=F, row.names=T, sep=";")
write.csv(capseq_merged_mat, file = capseq_merg_mat_csv, row.names = FALSE)

#----------------------------------------
# Plot mutation frequency in all cohorts
#----------------------------------------

### Combined SNV frequency plot
capseq_merged_freq <- capseq_merged_mat %>%
  reshape2::melt() %>%
  filter(value == "1") %>%
  select(Hugo_Symbol = Var1, Var2) %>%
  dplyr::distinct() %>%
  count(Hugo_Symbol) %>%
  arrange(desc(n)) %>%
  mutate(proportion_mutated = n/ncol(capseq_merged_mat)) %>%
  filter(proportion_mutated > 0.02) %>%
  select(Hugo_Symbol, n, proportion_mutated) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = proportion_mutated)) +
  geom_bar(stat = "identity") +
  ylab("Percentage of samples with mutation") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
        axis.title.x = element_blank())
filename <- paste0(date, "_UHN_E4402_PLOSMED-T-O_E2408-T-O_SNV_combined_freq.png")
ggsave(filename, width = 25, height = 10, units = "cm")
print("Combined SNV frequency plot saved in working directory as:")
print(filename)

### Faceted SNV frequency plot
genes__more_2perc <- capseq_merged_mat %>%
  reshape2::melt() %>%
  filter(value == "1") %>%
  select(Hugo_Symbol = Var1, Var2) %>%
  dplyr::distinct() %>%
  count(Hugo_Symbol) %>%
  arrange(desc(n)) %>%
  mutate(proportion_mutated = n/ncol(capseq_merged_mat)) %>%
  filter(proportion_mutated > 0.02) %>%
  select(Hugo_Symbol, n, proportion_mutated) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  .$Hugo_Symbol

capseq_merged_mat_plot <- capseq_merged_mat %>%
  reshape2::melt() %>%
  filter(value == "1") %>%
  select(Hugo_Symbol = Var1, Var2) %>%
  mutate(cohort = ifelse(Var2 %in% colnames(capseq_UHN), "UHN",
                  ifelse(Var2 %in% colnames(capseq_E4402), "E4402",
                  ifelse(Var2 %in% colnames(capseq_E2408), "E2408",
                  ifelse(Var2 %in% colnames(capseq_PLOSMED), "PLOSMED", Var2))))) %>%
  dplyr::distinct() %>%
  group_by(cohort) %>%
  count(Hugo_Symbol) %>%
  mutate(proportion_mutated = ifelse(cohort == "UHN", n/ncol(capseq_UHN),
                                     ifelse(cohort == "E4402", n/ncol(capseq_E4402),
                                            ifelse(cohort == "E2408", n/ncol(capseq_E2408),
                                                   ifelse(cohort == "PLOSMED", n/ncol(capseq_PLOSMED)))))) %>%
  filter(Hugo_Symbol %in% genes__more_2perc) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, genes__more_2perc)) %>%
  ggplot(aes(x = Hugo_Symbol, y = proportion_mutated)) +
  geom_bar(stat = "identity") +
  facet_grid(cohort ~ .) +
  ylab("Percentage of samples with mutation") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
        axis.title.x = element_blank())
filename <- paste0(date, "_UHN_E4402_PLOSMED-T-O_E2408-T-O_SNV_faceted_freq.png")
ggsave(filename, width = 25, height = 25, units = "cm")
print("Faceted SNV frequency plot saved in working directory as:")
print(filename)
