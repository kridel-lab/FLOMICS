
#---
# This script reads in mutation (capseq) data, creates SNV matrices
#   and plots.
# Authors: Robert Kridel + Victoria Shelton
# Last modified: February 1st, 2023
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
setwd("/your working directory/SNV_Clustering")
dir <- "/your working directory/SNV_Clustering"

#----------------------------------
# Read in data
#----------------------------------

### Read in gene panels
genes_plosmed <- read.csv("panel1.csv") %>%
  filter(PLOS_MED_PANEL == "YES") %>%
  dplyr::select(Gene.Name)
genes_uhn <- read.csv("panel2.csv") %>%
  filter(FLOMICS_PANEL == "YES") %>%
  dplyr::select(Gene.Name)


### Determine the common genes between the two utilized gene panels
genes_common <- intersect(genes_plosmed, genes_uhn)


### Read in filtered exonic and splice-impacted SNV calls (excluding synonymous)
capseq_var <- read.table("SNVcalls.txt", sep = ";", header = TRUE)

#check if there are any missing genes from this SNV dataset
dplyr::setdiff(genes_common$Gene.Name, capseq_var$Hugo_Symbol)

#if gene is targeted, but has no SNVs,
#  ensure values are set to '0' and not 'NA' for this gene
capseq_var <- rbind(capseq_var,
                     c("gene_name", rep(NA, ncol(capseq_var) - 1))) %>%
  replace(is.na(.), 0)
row.names(capseq_var) <- capseq_var$Hugo_Symbol
capseq_var$Hugo_Symbol <- NULL
capseq_var_mat <- as.matrix(capseq_var)



#-------------------
#Merge datasets
#-------------------
genes_common_vec <- unlist(genes_common)
capseq_merged_mat <- cbind(capseq_var_mat[genes_common_vec, ])
capseq_merged_df <- cbind(capseq_var[genes_common_vec, ])


### Transfer dataframe to new object,
## and remove the extra Hugo_Symbol column in this second dataframe
capseq_merged_no_hugo <- capseq_merged_df
capseq_merged_no_hugo$Hugo_Symbol <- NULL
#By making Hugo_Symbol null,
#we lose the ability to create a frequency plot downstream later,
#so need to transfer this to a new object


### Examine the merged dataframe
print("Dimensions of the merged dataset including all cohorts")
dim(capseq_merged_df) # 57 x 713


### Change the class of values in merged dataset
## from "character" to "num" to input into downstream analyses
capseq_merg_num <- japply(capseq_merged_df,
                          which(sapply(capseq_merged_df, class) == "character"),
                          as.numeric)


### Writing out the merged dataframe
capseq_merg_txt <- paste0(dir, date, "_gene_vs_sample_SNV_matrix.txt")
capseq_merg_csv <- paste0(dir, date, "_gene_vs_sample_SNV_matrix.csv")
write.table(capseq_merg_num,
            file = capseq_merg_txt, quote = FALSE, row.names = TRUE, sep = ";")
write.csv(capseq_merg_num, file = capseq_merg_csv, row.names = FALSE)

### Writing out the merged matrix
capseq_merg_mat_txt <- paste0(dir, date, "_gene_vs_sample_SNV_mat_matrix.txt")
capseq_merg_mat_csv <- paste0(dir, date, "_gene_vs_sample_SNV_mat_matrix.csv")
write.table(capseq_merged_mat,
            file = capseq_merg_mat_txt,
                    quote = FALSE, row.names = TRUE, sep = ";")
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
  mutate(proportion_mutated = n / ncol(capseq_merged_mat)) %>%
  filter(proportion_mutated > 0.02) %>%
  select(Hugo_Symbol, n, proportion_mutated) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = proportion_mutated)) +
  geom_bar(stat = "identity") +
  ylab("Percentage of samples with mutation") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                    vjust = 0.5, hjust = 1, face = "italic"),
        axis.title.x = element_blank())
ggsave(paste0(date, "SNV_combined_freq.png"),
         width = 25, height = 10, units = "cm")


### Faceted SNV frequency plot (by cohort)
genes__more_2perc <- capseq_merged_mat %>%
  reshape2::melt() %>%
  filter(value == "1") %>%
  select(Hugo_Symbol = Var1, Var2) %>%
  dplyr::distinct() %>%
  count(Hugo_Symbol) %>%
  arrange(desc(n)) %>%
  mutate(proportion_mutated = n / ncol(capseq_merged_mat)) %>%
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
                  ifelse(Var2 %in% colnames(capseq_PLOSMED), "PLOSMED",
                   Var2))))) %>%
  dplyr::distinct() %>%
  group_by(cohort) %>%
  count(Hugo_Symbol) %>%
  mutate(proportion_mutated = ifelse(cohort == "UHN", n / ncol(capseq_UHN),
                              ifelse(cohort == "E4402", n / ncol(capseq_E4402),
                              ifelse(cohort == "E2408", n / ncol(capseq_E2408),
                              ifelse(cohort == "PLOSMED", n /
                               ncol(capseq_PLOSMED)))))) %>%
  filter(Hugo_Symbol %in% genes__more_2perc) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, genes__more_2perc)) %>%
  ggplot(aes(x = Hugo_Symbol, y = proportion_mutated)) +
  geom_bar(stat = "identity") +
  facet_grid(cohort ~ .) +
  ylab("Percentage of samples with mutation") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                    vjust = 0.5, hjust = 1, face = "italic"),
        axis.title.x = element_blank())
ggsave(paste0(date, "_SNV_faceted_freq.png"),
        width = 25, height = 25, units = "cm")
