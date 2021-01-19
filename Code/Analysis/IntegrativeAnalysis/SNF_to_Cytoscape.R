

library(dplyr)
library(data.table)
library(reshape2)

setwd("~/github/FLOMICS/")

# Read in Cluster labels
cluster_labels <- read.csv("Cluster Labels/InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv") %>%
  select(Source_node = ID, SNFClust)

# set filter
x <- 0.005

# Read in affinity matrices from SNF
dim3_r_m_d <- read.csv(file = "SNF/WCombinedFiltered_6Jan2021.csv") %>%
  melt() %>%
  filter(value != 0.5) %>%
  setnames(old = "X", new = "Source_node") %>%
  setnames(old = "variable", new = "Target_node") %>%
  left_join(cluster_labels) %>%
  mutate(Source_node = substr(Source_node, 1, 9)) %>%
  mutate(Target_node = substr(Target_node, 1, 9)) %>%
  # filter(value > x) %>%
  mutate(BIOPAX_TYPE = "interaction")

write.table(dim3_r_m_d, file = "SNF/WCombinedFiltered_6Jan2021_melted.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#

dim2_m_r <- read.csv(file = "SNF/W2MethylationRNAseq_8Jan2021.csv") %>%
  melt() %>%
  filter(value != 0.5) %>%
  setnames(old = "X", new = "Source_node") %>%
  setnames(old = "variable", new = "Target_node") %>%
  left_join(cluster_labels) %>%
  mutate(Source_node = substr(Source_node, 1, 9)) %>%
  mutate(Target_node = substr(Target_node, 1, 9)) %>%
  # filter(value > x) %>%
  mutate(BIOPAX_TYPE = "interaction")

write.table(dim2_m_r, file = "SNF/W2MethylationRNAseq_8Jan2021_melted.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#

dim2_m_d <- read.csv(file = "SNF/W2MethylationDNAseq_8Jan2021.csv") %>%
  melt() %>%
  filter(value != 0.5) %>%
  setnames(old = "X", new = "Source_node") %>%
  setnames(old = "variable", new = "Target_node") %>%
  left_join(cluster_labels) %>%
  mutate(Source_node = substr(Source_node, 1, 9)) %>%
  mutate(Target_node = substr(Target_node, 1, 9)) %>%
  # filter(value > x) %>%
  mutate(BIOPAX_TYPE = "interaction")

write.table(dim2_m_d, file = "SNF/W2MethylationDNAseq_8Jan2021_melted.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#

dim2_r_d <- read.csv(file = "SNF/W2RNAseqDNAseq_8Jan2021.csv") %>%
  melt() %>%
  filter(value != 0.5) %>%
  setnames(old = "X", new = "Source_node") %>%
  setnames(old = "variable", new = "Target_node") %>%
  left_join(cluster_labels) %>%
  mutate(Source_node = substr(Source_node, 1, 9)) %>%
  mutate(Target_node = substr(Target_node, 1, 9)) %>%
  # filter(value > x) %>%
  mutate(BIOPAX_TYPE = "interaction")

write.table(dim2_r_d, file = "SNF/W2RNAseqDNAseq_8Jan2021_melted.txt", row.names = FALSE, sep = "\t", quote = FALSE)
