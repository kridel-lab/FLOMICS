#---
# This script reads in the clustering assignment table
## and converts cluster labels to new subtype names
# Authors: Victoria Shelton
# Last modified: April 12th, 2023
#---

#--
# Load in packages
#--
packages <- c("dplyr", "ggplot2", "flexmix", "tidyverse", "gridExtra")
lapply(packages, require, character.only = TRUE)

#--
# Setup working space
#--

date <- Sys.Date()
setwd("~/Kridel Lab/")
dir <- "~/Kridel Lab/"

#----------------------------------
# Read in data
#----------------------------------
### load cluster labels
cluster_assignments <- read.csv("GMM_Cluster_Labels_flexmix_clusters.csv")

#---
# converting labels
#---
new_subtype_labels <- cluster_assignments %>%
  dplyr::mutate(SubtypeAIC = ifelse(ClusterAIC == "C1", "CS",
                             ifelse(ClusterAIC == "C2", "TT",
                             ifelse(ClusterAIC == "C3", "GM",
                             ifelse(ClusterAIC == "C4", "Q",
                             ifelse(ClusterAIC == "C5", "AR", ClusterAIC))))))

new_subtype_labels_2 <- new_subtype_labels %>%
  dplyr::select(-ClusterAIC) %>%
  dplyr::rename(ClusterAIC = SubtypeAIC) %>%
  dplyr::select(PATIENT_ID, SAMPLE_ID, ClusterAIC, cohort)

#write out the new labels
write.csv(new_subtype_labels_2,
 file = paste0(dir, date, "_GMM_Subtype_Labels_flexmix_clusters.csv"),
  row.names = FALSE)
