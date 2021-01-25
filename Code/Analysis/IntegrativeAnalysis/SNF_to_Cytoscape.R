
#--
# Script to generate network in Cytoscape
#--

library(dplyr)
library(data.table)
library(reshape2)

setwd("~/github/FLOMICS/")

#--
# Generate input for Cytoscape
#--

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

#--
# Steps in Cytoscape
#--

# using Cytoscape 3.8.2
# File -> Import -> Network from File ...
# verify that "Edge Attribute" selected for "value", and "Source Node Attribute" for "SNFClust"
# press OK

# go to "Style" in bar on left side
# change width and height to 30
# change font size to 20
# click on "Fill Color", select Column -> SNFClust, Mapping Type -> Discrete Mapping
#   1: R67, G99, B216 (i.e. color #4363D8)
#   2: R245, G130, B49 (i.e. color #F58231)
# click on "Label Color" -> white
# click on "Shape" -> Ellipse
# click on "Label", select Column -> SNFClust, Mapping Type -> Discrete Mapping
#   1: 1
#   2: 2

# go to "Filter" in bar on left side
# create new filter, "Column Filter", "Choose Column" -> "Edge: value"
# slide to between 0.005 and 0.071 (i.e. max)

# Click "Layout" in Menubar, Edge-Weighted Spring Embedded Layout -> value
                     

