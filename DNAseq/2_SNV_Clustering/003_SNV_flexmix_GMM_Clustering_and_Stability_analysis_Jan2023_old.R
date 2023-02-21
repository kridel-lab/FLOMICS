#---
# This script reads in mutation (CAPSEQ) data and performs clustering
# Authors: Robert Kridel + Victoria Shelton
# Last modified: February 1st, 2023
#---

#--
# Load in packages
#--
packages <- c("dplyr", "ggplot2", "flexmix", "tidyverse", "gridExtra")
lapply(packages, require, character.only = TRUE)

packages2 <- c("magrittr", "cluster", "cluster.datasets", "cowplot",
               "ggfortify", "dendextend", "factoextra", "FactoMineR",
              "corrplot", "GGally", "ggiraphExtra", "knitr", "kableExtra",
              "taRifx", "ggalluvial", "clValid", "mclust")
lapply(packages2, require, character.only = TRUE)

#--
# Setup working space
#--

date <- Sys.Date()
setwd("/your working directory/SNV_Clustering")
dir <- "/your working directory/SNV_Clustering"

#---
# Read in data
#---

### load SNV dataset
capseq_merg_mat_txt <- paste0("gene_vs_sample_SNV_mat_matrix.txt")
capseq_merged_mat <- read.table(capseq_merg_mat_txt, sep = ";")
print("Number of samples in CAPSEQ Merged dataset")
ncol(capseq_merged_mat)

### load cohort sample IDs
uhn_samples <- read.delim("capseq_UHN_samples.txt",
 header = FALSE, sep = ",")
uhn_samples <- unlist(uhn_samples)

e4402_samples <- read.delim("capseq_E4402_samples.txt",
 header = FALSE, sep = ",")
e4402_samples <- unlist(e4402_samples)

e2408_samples <- read.delim("capseq_E2408_samples.txt",
 header = FALSE, sep = ",")
e2408_samples <- unlist(e2408_samples)

plosmed_samples <- read.delim("capseq_PLOSMED_samples.txt",
 header = FALSE, sep = ",")
plosmed_samples <- unlist(plosmed_samples)

message("Loaded in Cohort sample IDs")

### Read in gene panels to determine comon genes
Genes_plosmed <- read.csv("panel1.csv") %>%
  filter(PLOS_MED_PANEL == "YES") %>%
  dplyr::select(Gene.Name)
genes_uhn <- read.csv("panel2.csv") %>%
  filter(FLOMICS_PANEL == "YES") %>%
  dplyr::select(Gene.Name)

### Determine the common genes between the two utilized gene panels
genes_common <- intersect(genes_plosmed, genes_uhn)

### load AIC and BIC consensus matrices
aic_consensus_mat_txt <- paste0("aic_consensus_mat.txt")
aic_consensus_mat <- read.table(aic_consensus_mat_txt, header = TRUE, sep = ";")
samples <- colnames(aic_consensus_mat)
message("Number of samples in AIC Consensus Matrix")
ncol(aic_consensus_mat)

bic_consensus_mat_txt <- paste0("bic_consensus_mat.txt")
bic_consensus_mat <- read.table(bic_consensus_mat_txt, header = TRUE, sep = ";")
samples <- colnames(bic_consensus_mat)
message("Number of samples in BIC Consensus Matrix")
ncol(bic_consensus_mat)

#--
# Clustering
# Method origin: https://github.com/ecsg-uoy/DLBCLGenomicSubtyping
#--

### loading in functions
source("script-plot-all-mutations-all-cohorts_to_be_loaded.R")
message("loaded in plotting functions")

### set directory again
setwd("/your-working-directory/")

### formating dataset for entry in flexix modeling
muts_all <- t(capseq_merged_mat) %>%
 data.frame() %>%
 mutate(SAMPLE_ID = row.names(.))
ncol(muts_all) #the number of genes
nrow(muts_all) #the nmber of samples
muts_all <- cbind(SAMPLE_ID = muts_all$SAMPLE_ID,
 muts_all[, c(1:(ncol(muts_all) - 1))])
row.names(muts_all) <- NULL
muts_df <- muts_all[, -1]
muts_df <- sapply(muts_df, as.numeric)
message("Dimensions of formated dataset:")
dim(muts_df)

### creating lists to store the cluster stability outputs
tmpar <- matrix(0, 1, 7) #this is a temporary matrix
aic_cluster_stab_arry_5 <- array(tmpar, dim = c(1, 5, 100)) #5 clusters
aic_cluster_stab_arry_6 <- array(tmpar, dim = c(1, 6, 100)) #6 clusters

### creating lists to store the numbers of patients per cluster
aic_cluster_patnum_arry_5 <- array(tmpar, dim = c(1, 5, 100))
aic_cluster_patnum_arry_6 <- array(tmpar, dim = c(1, 6, 100))

### creating lists to store the sum of the -log10(q) scores
aic_cluster_score_arry_5 <- matrix(0, 1, 100)
aic_cluster_score_arry_6 <- matrix(0, 1, 100)

### creating seeds for cluster runs
#list of random numbers for the clustering seed
x5 <- as.integer(runif(n = 100, min = 1, max = 1000))
write.table(x5, file = "5_cluster_clustering_seeds.txt",
 quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ";")

x6 <- as.integer(runif(n = 100, min = 1, max = 1000))
write.table(x6, file = "6_cluster_clustering_seeds.txt",
 quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ";")

### set number of iterations
numtries <- 100

### set overall clustering seed
set.seed(99)

### clustering
RNGversion("4.0.2")
criteria <- "AIC"

### starting iterations for 5 clusters
for(i in 1:numtries){
set.seed(x5[i]) #setting seed for clustering
  muts_all <- t(capseq_merged_mat) %>%
   data.frame() %>%
    mutate(SAMPLE_ID = row.names(.))
  muts_all <- cbind(SAMPLE_ID = muts_all$SAMPLE_ID,
   muts_all[, c(1:(ncol(muts_all) - 1))])
  ex <- initFlexmix(muts_df ~ 1,
                  k = 5,
                  model = FLXMCmvbinary(),
                  control = list(minprior = 0.1),
                  nrep = 1000)
save(ex, file = "clustering_ex_5.RData")

### saving the Number of Components plot
message("finished Clustering")
pdf(file = paste0(date, "_", i, "_components-plot_5-clusters.pdf"))
plot(ex)
dev.off()
message("components plot is saved in working directory as:")
print(paste0(date, "_", i, "_components-plot_5-clusters.pdf"))

### appending cluster assignments
muts_all$ClusterAIC <- as.factor(paste0('C', aic@cluster))
table(muts_all$ClusterAIC)


### generating and saving cluster maps
## and top 10 -log10(q) values for each criteria
genes <- colnames(muts_df)
plt_aic <- heatmap_mutation_extended(muts_all,
 genes, 'ClusterAIC', y_order = 'group', idcol = 'SAMPLE_ID')
grid.arrange(plt_aic)
ggsave(file = paste0(date, "_", i, "_heatmap_aic.png"),
 plt_aic, width = 25, height = 25, units = "cm")

sig_long_out_path <- paste0(date, "_", i, "_sig_long_out_", criteria, ".txt")
sig_long_out_mat <- read.table(sig_long_out_path, sep = ";", header = TRUE) %>%
 dplyr::filter(to_plot == TRUE)
 #the to_plot filter keeps only genes that are considered significant
 # to clustering outcome
message("saved AIC sum of all -log10(q) scores,
 the -log10 of the adjusted p-value (fdr adjusted)")
aic_cluster_score_arry_5[, i] <- sum(sig_long_out_mat$val)
 #summing all the sig scores
sig_long_out_mat_sort <- sig_long_out_mat[with(sig_long_out_mat, order(-val)), ]
 #sorting entries by val in descending order
aic_cluster_score_top_arry_5[, i] <- sig_long_out_mat_sort[1:10, ] %>%
 dplyr::select(val) %>%
  sum() #summing the top 10 entries

### attaching cohort labels to samples
muts_all_sub <- muts_all %>%
  dplyr::mutate(PATIENT_ID = substr(SAMPLE_ID, 1, 9)) %>%
  dplyr::select(PATIENT_ID, SAMPLE_ID, ClusterAIC, ClusterBIC) %>%
  dplyr::mutate(cohort = ifelse(SAMPLE_ID %in% UHNsamples, "UHN",
                         ifelse(SAMPLE_ID %in% E4402samples, "E4402",
                                ifelse(SAMPLE_ID %in% E2408samples, "E2408",
                                       ifelse(SAMPLE_ID %in% PLOSMEDsamples,
                                        "PLOSMED", SAMPLE_ID)))))

write.csv(muts_all_sub,
 file = paste0(date, "_", i, "_5_GMM_Cluster_Labels_flexmix_clusters.csv"),
 row.names = FALSE)
message("saved cluster labels")

#------------------------------------------------------
# Assesing the Stability of Individual Clusters
#------------------------------------------------------

### Table of cluster assignments
aic_clus_assign <- data.frame(matrix(nrow=713)) %>%
 add_column(Samples = colnames(capseq_merged_mat),
  Cluster = ex@models[["5"]]@cluster) %>%
 dplyr::select(-1)

### extracting the samples in each cluster,
## the number of samples and the sample identifiers
aic_clus1 <- aic_clus_assign %>%
 dplyr::filter(Cluster == "1") %>%
  dplyr::select(1)
aic_patnum_1 <- dim(aic_clus1)
aic_clus1_samples <- aic_clus1[, 1]

aic_clus2 <- aic_clus_assign %>%
 dplyr::filter(Cluster == "2") %>%
  dplyr::select(1)
aic_patnum_2 <- dim(aic_clus2)
aic_clus2_samples <- aic_clus2[, 1]

aic_clus3 <- aic_clus_assign %>%
 dplyr::filter(Cluster == "3") %>%
  dplyr::select(1)
aic_patnum_3 <- dim(aic_clus3)
aic_clus3_samples <- aic_clus3[, 1]

aic_clus4 <- aic_clus_assign %>%
 dplyr::filter(Cluster == "4") %>%
  dplyr::select(1)
aic_patnum_4 <- dim(aic_clus4)
aic_clus4_samples <- aic_clus4[, 1]

aic_clus5 <- aic_clus_assign %>%
 dplyr::filter(Cluster == "5") %>%
  dplyr::select(1)
aic_patnum_5 <- dim(aic_clus5)
aic_clus5_samples <- aic_clus5[, 1]

#aic_clus6 <- aic_clus_assign %>% dplyr::filter(Cluster == "6") %>% dplyr::select(1)
#aic_patnum_6 <- dim(aic_clus6) 
#aic_clus6_samples <- aic_clus6[,1]

### storing number of patients
aic_patnum_mat <- cbind(aic_patnum_1, aic_patnum_2, aic_patnum_3, aic_patnum_4,
 aic_patnum_5)#, aic_patnum_6)
aic_cluster_patnum_arry_5[, , i] <- aic_patnum_mat[1, ]
write.csv(aic_patnum_mat,
 file = paste0(date, "_", i, "_5_AIC_Cluster_Patient_Number_Matrix.csv"),
row.names = FALSE)
message("saved aic_patnum_mat")

### AIC CLUS1 MATRIX SUM ###
aic_clus1_sum_mat <- matrix(0, length(aic_clus1_samples),
 length(aic_clus1_samples)) %>%
 set_colnames(aic_clus1_samples) %>%
 set_rownames(aic_clus1_samples) #create matrix to store consensus
for (col_id in aic_clus1_samples){
  for (row_id in aic_clus1_samples){
    aic_clus1_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
  }
}
#sum up the AIC Cluster 1 consensus matrix
aic_clus1_sum <- sum(aic_clus1_sum_mat)
message("made aic_clus1_sum matrix")

### AIC CLUS2 MATRIX SUM ###
aic_clus2_sum_mat <- matrix(0, length(aic_clus2_samples),
 length(aic_clus2_samples)) %>%
 set_colnames(aic_clus2_samples) %>%
 set_rownames(aic_clus2_samples) #create matrix to store consensus
for (col_id in aic_clus2_samples){
  for (row_id in aic_clus2_samples){
    aic_clus2_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
  }
}
#sum up the AIC Cluster 2 consensus matrix
aic_clus2_sum <- sum(aic_clus2_sum_mat)

### AIC CLUS3 MATRIX SUM ###
aic_clus3_sum_mat <- matrix(0, length(aic_clus3_samples),
 length(aic_clus3_samples)) %>%
 set_colnames(aic_clus3_samples) %>%
 set_rownames(aic_clus3_samples) #create matrix to store consensus
for (col_id in aic_clus3_samples){
  for (row_id in aic_clus3_samples){
    aic_clus3_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
  }
}
#sum up the AIC Cluster 3 consensus matrix
aic_clus3_sum <- sum(aic_clus3_sum_mat)

### AIC CLUS4 MATRIX SUM ###
aic_clus4_sum_mat <- matrix(0, length(aic_clus4_samples),
 length(aic_clus4_samples)) %>%
 set_colnames(aic_clus4_samples) %>%
 set_rownames(aic_clus4_samples) #create matrix to store consensus
for (col_id in aic_clus4_samples){
  for (row_id in aic_clus4_samples){
    aic_clus4_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
  }
}
#sum up the AIC Cluster 4 consensus matrix
aic_clus4_sum <- sum(aic_clus4_sum_mat)

### AIC CLUS5 MATRIX SUM ###
aic_clus5_sum_mat <- matrix(0, length(aic_clus5_samples),
 length(aic_clus5_samples)) %>%
 set_colnames(aic_clus5_samples) %>%
 set_rownames(aic_clus5_samples) #create matrix to store consensus
for (col_id in aic_clus5_samples){
  for (row_id in aic_clus5_samples){
    aic_clus5_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
  }
}
#sum up the AIC Cluster 5 consensus matrix
aic_clus5_sum <- sum(aic_clus5_sum_mat)

### AIC CLUS6 MATRIX SUM ###
#aic_clus6_sum_mat <- matrix(0, length(aic_clus6_samples),
# length(aic_clus6_samples)) %>% 
# set_colnames(aic_clus6_samples) %>% 
# set_rownames(aic_clus6_samples) #create matrix to store consensus
#for (col_id in aic_clus6_samples){
#  for (row_id in aic_clus6_samples){
#    aic_clus6_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
#  }
#}
#sum up the AIC Cluster 6 consensus matrix
#aic_clus6_sum <- sum(aic_clus6_sum_mat)


### 5 Cluster Stabilities:
#computing the consensus of AIC Cluster 1
aic_clus1_stab <- ((1/ ((length(aic_clus1_samples)) *
 ((length(aic_clus1_samples)) - 1) / 2)) * (aic_clus1_sum))

#computing the consensus of AIC Cluster 2
aic_clus2_stab <- ((1/ ((length(aic_clus2_samples)) *
 ((length(aic_clus2_samples)) - 1) / 2)) * (aic_clus2_sum))

#computing the consensus of AIC Cluster 3
aic_clus3_stab <- ((1/ ((length(aic_clus3_samples)) *
 ((length(aic_clus3_samples)) - 1) / 2)) * (aic_clus3_sum))

#computing the consensus of AIC Cluster 4
aic_clus4_stab <- ((1/ ((length(aic_clus4_samples)) *
 ((length(aic_clus4_samples)) - 1) / 2)) * (aic_clus4_sum))

#computing the consensus of AIC Cluster 5
aic_clus5_stab <- ((1/ ((length(aic_clus5_samples)) *
 ((length(aic_clus5_samples)) - 1) / 2)) * (aic_clus5_sum))

#computing the consensus of AIC Cluster 6
#aic_clus6_stab <- ((1/ ((length(aic_clus6_samples)) *
# ((length(aic_clus6_samples)) - 1) / 2)) * (aic_clus6_sum))


### A matrix of the cluster stability scores
aic_stab_mat <- cbind(aic_clus1_stab, aic_clus2_stab, aic_clus3_stab,
 aic_clus4_stab, aic_clus5_stab)#, aic_clus6_stab)

# storing the cluster stability scores matrix in the designated array
aic_cluster_stab_arry_5[, , i] <- aic_stab_mat
write.csv(aic_stab_mat,
 file = paste0(date, "_", i, "_5_AIC_Cluster_Stability_Matrix.csv"),
 row.names = FALSE)
message("saved aic_stab_mat")

}

### writing out the 5 cluster stability arrays
write.csv(aic_cluster_stab_arry_5,
 file = paste0(date, "_5_AIC_Cluster_Stability_Array.csv"), row.names = FALSE)
saveRDS(aic_cluster_stab_arry_5, file = "aic_cluster_stab_arry_5.RData")
saveRDS(aic_cluster_patnum_arry_5, file = "aic_cluster_patnum_arry_5.RData")
saveRDS(aic_cluster_score_arry_5, file = "aic_cluster_score_arry_5.RData")
saveRDS(aic_cluster_score_top_arry_5,
 file = "aic_cluster_score_top_arry_5.RData")

#------------------
#------------------

### starting iterations for 6 clusters
for(i in 1:numtries){
set.seed(x6[i]) #setting seed for clustering
  muts_all <- t(capseq_merged_mat) %>%
   data.frame() %>%
    mutate(SAMPLE_ID = row.names(.))
  muts_all <- cbind(SAMPLE_ID = muts_all$SAMPLE_ID,
   muts_all[, c(1:(ncol(muts_all) - 1))])
  ex <- initFlexmix(muts_df ~ 1,
                  k = 6,
                  model = FLXMCmvbinary(),
                  control = list(minprior = 0.1),
                  nrep = 1000)
save(ex, file = "clustering_ex_6.RData")

### saving the Number of Components plot
message("finished Clustering")
pdf(file = paste0(date, "_", i, "_components-plot_6-clusters.pdf"))
plot(ex)
dev.off()
message("components plot is saved in working directory as:")
print(paste0(date, "_", i, "_components-plot_6-clusters.pdf"))

### appending cluster assignments
muts_all$ClusterAIC <- as.factor(paste0('C', aic@cluster))
table(muts_all$ClusterAIC)

### generating and saving cluster maps
## and top 10 -log10(q) values for each criteria
genes <- colnames(muts_df)
plt_aic <- heatmap_mutation_extended(muts_all,
 genes, 'ClusterAIC', y_order = 'group', idcol = 'SAMPLE_ID')
grid.arrange(plt_aic)
ggsave(file = paste0(date, "_", i, "_heatmap_aic.png"),
 plt_aic, width = 25, height = 25, units = "cm")

sig_long_out_path <- paste0(date, "_", i, "_sig_long_out_", criteria, ".txt")
sig_long_out_mat <- read.table(sig_long_out_path, sep = ";", header = TRUE) %>%
 dplyr::filter(to_plot == TRUE)
 #the to_plot filter keeps only genes that are considered significant
 # to clustering outcome
message("saved AIC sum of all -log10(q) scores,
 the -log10 of the adjusted p-value (fdr adjusted)")
aic_cluster_score_arry_6[, i] <- sum(sig_long_out_mat$val)
 #summing all the sig scores
sig_long_out_mat_sort <- sig_long_out_mat[with(sig_long_out_mat, order(-val)), ]
 #sorting entries by val in descending order
aic_cluster_score_top_arry_6[, i] <- sig_long_out_mat_sort[1:10, ] %>%
 dplyr::select(val) %>%
  sum() #summing the top 10 entries


### attaching cohort labels to samples
muts_all_sub <- muts_all %>%
  dplyr::mutate(PATIENT_ID = substr(SAMPLE_ID, 1, 9)) %>%
  dplyr::select(PATIENT_ID, SAMPLE_ID, ClusterAIC, ClusterBIC) %>%
  dplyr::mutate(cohort = ifelse(SAMPLE_ID %in% UHNsamples, "UHN",
                         ifelse(SAMPLE_ID %in% E4402samples, "E4402",
                                ifelse(SAMPLE_ID %in% E2408samples, "E2408",
                                       ifelse(SAMPLE_ID %in% PLOSMEDsamples,
                                        "PLOSMED", SAMPLE_ID)))))

write.csv(muts_all_sub,
 file = paste0(date, "_", i, "_6_GMM_Cluster_Labels_flexmix_clusters.csv"),
 row.names = FALSE)
message("saved cluster labels")

#------------------------------------------------------
# Assesing the Stability of Individual Clusters
#------------------------------------------------------

### Table of cluster assignments
aic_clus_assign <- data.frame(matrix(nrow = 713)) %>%
 add_column(Samples = colnames(capseq_merged_mat),
  Cluster = ex@models[["6"]]@cluster) %>%
 dplyr::select(-1)

### extracting the samples in each cluster,
## the number of samples and the sample identifiers
aic_clus1 <- aic_clus_assign %>%
 dplyr::filter(Cluster == "1") %>%
  dplyr::select(1)
aic_patnum_1 <- dim(aic_clus1)
aic_clus1_samples <- aic_clus1[, 1]

aic_clus2 <- aic_clus_assign %>%
 dplyr::filter(Cluster == "2") %>%
  dplyr::select(1)
aic_patnum_2 <- dim(aic_clus2)
aic_clus2_samples <- aic_clus2[, 1]

aic_clus3 <- aic_clus_assign %>%
 dplyr::filter(Cluster == "3") %>%
  dplyr::select(1)
aic_patnum_3 <- dim(aic_clus3)
aic_clus3_samples <- aic_clus3[, 1]

aic_clus4 <- aic_clus_assign %>%
 dplyr::filter(Cluster == "4") %>%
  dplyr::select(1)
aic_patnum_4 <- dim(aic_clus4)
aic_clus4_samples <- aic_clus4[, 1]

aic_clus5 <- aic_clus_assign %>%
 dplyr::filter(Cluster == "5") %>%
  dplyr::select(1)
aic_patnum_5 <- dim(aic_clus5)
aic_clus5_samples <- aic_clus5[, 1]

aic_clus6 <- aic_clus_assign %>%
 dplyr::filter(Cluster == "6") %>%
  dplyr::select(1)
aic_patnum_6 <- dim(aic_clus6)
aic_clus6_samples <- aic_clus6[, 1]

### storing number of patients
aic_patnum_mat <- cbind(aic_patnum_1, aic_patnum_2, aic_patnum_3, aic_patnum_4,
 aic_patnum_5, aic_patnum_6)
aic_cluster_patnum_arry_5[, , i] <- aic_patnum_mat[1, ]
write.csv(aic_patnum_mat,
 file = paste0(date, "_", i, "_6_AIC_Cluster_Patient_Number_Matrix.csv"),
row.names = FALSE)
message("saved aic_patnum_mat")

### AIC CLUS1 MATRIX SUM ###
aic_clus1_sum_mat <- matrix(0, length(aic_clus1_samples),
 length(aic_clus1_samples)) %>%
 set_colnames(aic_clus1_samples) %>%
 set_rownames(aic_clus1_samples) #create matrix to store consensus
for (col_id in aic_clus1_samples){
  for (row_id in aic_clus1_samples){
    aic_clus1_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
  }
}
#sum up the AIC Cluster 1 consensus matrix
aic_clus1_sum <- sum(aic_clus1_sum_mat)
message("made aic_clus1_sum matrix")

### AIC CLUS2 MATRIX SUM ###
aic_clus2_sum_mat <- matrix(0, length(aic_clus2_samples),
 length(aic_clus2_samples)) %>%
 set_colnames(aic_clus2_samples) %>%
 set_rownames(aic_clus2_samples) #create matrix to store consensus
for (col_id in aic_clus2_samples){
  for (row_id in aic_clus2_samples){
    aic_clus2_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
  }
}
#sum up the AIC Cluster 2 consensus matrix
aic_clus2_sum <- sum(aic_clus2_sum_mat)

### AIC CLUS3 MATRIX SUM ###
aic_clus3_sum_mat <- matrix(0, length(aic_clus3_samples),
 length(aic_clus3_samples)) %>%
 set_colnames(aic_clus3_samples) %>%
 set_rownames(aic_clus3_samples) #create matrix to store consensus
for (col_id in aic_clus3_samples){
  for (row_id in aic_clus3_samples){
    aic_clus3_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
  }
}
#sum up the AIC Cluster 3 consensus matrix
aic_clus3_sum <- sum(aic_clus3_sum_mat)

### AIC CLUS4 MATRIX SUM ###
aic_clus4_sum_mat <- matrix(0, length(aic_clus4_samples),
 length(aic_clus4_samples)) %>%
 set_colnames(aic_clus4_samples) %>%
 set_rownames(aic_clus4_samples) #create matrix to store consensus
for (col_id in aic_clus4_samples){
  for (row_id in aic_clus4_samples){
    aic_clus4_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
  }
}
#sum up the AIC Cluster 4 consensus matrix
aic_clus4_sum <- sum(aic_clus4_sum_mat)

### AIC CLUS5 MATRIX SUM ###
aic_clus5_sum_mat <- matrix(0, length(aic_clus5_samples),
 length(aic_clus5_samples)) %>%
 set_colnames(aic_clus5_samples) %>%
 set_rownames(aic_clus5_samples) #create matrix to store consensus
for (col_id in aic_clus5_samples){
  for (row_id in aic_clus5_samples){
    aic_clus5_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
  }
}
#sum up the AIC Cluster 5 consensus matrix
aic_clus5_sum <- sum(aic_clus5_sum_mat)

### AIC CLUS6 MATRIX SUM ###
aic_clus6_sum_mat <- matrix(0, length(aic_clus6_samples),
 length(aic_clus6_samples)) %>% 
 set_colnames(aic_clus6_samples) %>% 
 set_rownames(aic_clus6_samples) #create matrix to store consensus
for (col_id in aic_clus6_samples){
  for (row_id in aic_clus6_samples){
    aic_clus6_sum_mat[row_id, col_id] <- aic_consensus_mat[row_id, col_id]
  }
}
#sum up the AIC Cluster 6 consensus matrix
aic_clus6_sum <- sum(aic_clus6_sum_mat)


### 6 Cluster Stabilities:
#computing the consensus of AIC Cluster 1
aic_clus1_stab <- ((1/ ((length(aic_clus1_samples)) *
 ((length(aic_clus1_samples)) - 1) / 2)) * (aic_clus1_sum))

#computing the consensus of AIC Cluster 2
aic_clus2_stab <- ((1 / ((length(aic_clus2_samples)) *
 ((length(aic_clus2_samples)) - 1) / 2)) * (aic_clus2_sum))

#computing the consensus of AIC Cluster 3
aic_clus3_stab <- ((1/ ((length(aic_clus3_samples)) *
 ((length(aic_clus3_samples)) - 1) / 2)) * (aic_clus3_sum))

#computing the consensus of AIC Cluster 4
aic_clus4_stab <- ((1/ ((length(aic_clus4_samples)) *
 ((length(aic_clus4_samples)) - 1) / 2)) * (aic_clus4_sum))

#computing the consensus of AIC Cluster 5
aic_clus5_stab <- ((1/ ((length(aic_clus5_samples)) *
 ((length(aic_clus5_samples)) - 1) / 2)) * (aic_clus5_sum))

#computing the consensus of AIC Cluster 6
aic_clus6_stab <- ((1/ ((length(aic_clus6_samples)) *
 ((length(aic_clus6_samples)) - 1) / 2)) * (aic_clus6_sum))


### A matrix of the cluster stability scores
aic_stab_mat <- cbind(aic_clus1_stab, aic_clus2_stab, aic_clus3_stab,
 aic_clus4_stab, aic_clus5_stab, aic_clus6_stab)

# storing the cluster stability scores matrix in the designated array
aic_cluster_stab_arry_6[, , i] <- aic_stab_mat
write.csv(aic_stab_mat,
 file = paste0(date, "_", i, "_6_AIC_Cluster_Stability_Matrix.csv"),
 row.names = FALSE)
message("saved aic_stab_mat")

}

### writing out the 6 cluster stability arrays
write.csv(aic_cluster_stab_arry_6,
 file = paste0(date, "_6_AIC_Cluster_Stability_Array.csv"), row.names = FALSE)
saveRDS(aic_cluster_stab_arry_6, file = "aic_cluster_stab_arry_6.RData")
saveRDS(aic_cluster_patnum_arry_6, file = "aic_cluster_patnum_arry_6.RData")
saveRDS(aic_cluster_score_arry_6, file = "aic_cluster_score_arry_6.RData")
saveRDS(aic_cluster_score_top_arry_6,
 file = "aic_cluster_score_top_arry_6.RData")



#--
# Score sheets
#--

### 5 CLUSTERS ###
aic_score_mat_5 <- matrix(nrow = 100, ncol = 10)
colnames(aic_score_mat_5) <- c("weighted_mean", "mean", "median", "max", "min",
 "sum", "max_index", "min_index", "sum of -log10(q) scores",
  "sum of top 10 gene -log10(q) scores")
#computing the weighted averages
aic_weighted_arry_5 <- aic_cluster_stab_arry_5 * aic_cluster_patnum_arry_5
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
aic_weighted_arry_5[is.nan(aic_weighted_arry_5)] <- 0
aic_score_mat_5[, 1] <- apply(aic_weighted_arry_5, 3, sum) / 713
#the means
aic_score_mat_5[, 2] <- apply(aic_cluster_stab_arry_5, 3, mean, na.rm = TRUE)
#the medians
aic_score_mat_5[, 3] <- apply(aic_cluster_stab_arry_5, 3, median, na.rm = TRUE)
#the max stability
aic_score_mat_5[, 4] <- apply(aic_cluster_stab_arry_5, 3, max, na.rm = TRUE)
#the min stability
aic_score_mat_5[, 5] <- apply(aic_cluster_stab_arry_5, 3, min, na.rm = TRUE)
#the sum of all stabilities
aic_score_mat_5[, 6] <- apply(aic_cluster_stab_arry_5, 3, sum, na.rm = TRUE)
#the cluster with the max stability
aic_score_mat_5[, 7] <- apply(aic_cluster_stab_arry_5, 3, which.max)
#the cluster with the min stability
aic_score_mat_5[, 8] <- apply(aic_cluster_stab_arry_5, 3, which.min)
#the sum of -log10(q) scores
aic_score_mat_5[, 9] <- aic_cluster_score_arry_5
# the sum of top 10 -log10(q) scores
aic_score_mat_5[, 10] <- aic_cluster_score_top_arry_5

aic_score_mat_5_2 <- aic_score_mat_5 %>%
  as.data.frame() %>%
  dplyr::mutate(ranked_weighted_mean = rank(aic_score_mat_5[, 1]),
                ranked_unweighted_mean = rank(aic_score_mat_5[, 2]),
                ranked_median = rank(aic_score_mat_5[, 3]),
                ranked_max = rank(aic_score_mat_5[, 4]),
                ranked_min = rank(aic_score_mat_5[, 5]),
                ranked_sum = rank(aic_score_mat_5[, 6]),
                ranked_sigscore = rank(aic_score_mat_5[, 9]),
                ranked_topsigscore = rank(aic_score_mat_5[, 10]))

aic_score_mat_5_3 <- aic_score_mat_5_2 %>%
  dplyr::mutate(Stabilityscore =
   rowSums(aic_score_mat_5_2[, c(11, 13:15, 17)]))

write.csv(aic_score_mat_4, file = "aic_cluster_stability_score_sheet.csv", row.names = T)



### Matrix to store the score sheet of the clustering models
aic_score_mat_5 <- matrix(nrow = 20, ncol = 9)
colnames(aic_score_mat_5) <- c("weighted_mean", "mean", "median", "max", "min", "sum", "max_index", "min_index", "sum of all gene -log10(q) scores")

### Computing weighted average cluster stablility 
aic_weighted_arry_5 <- aic_cluster_stab_arry_5*aic_cluster_patnum_arry_5 # multiplying cluster stability and number of samples
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
aic_weighted_arry_5[is.nan(aic_weighted_arry_5)] <- 0
aic_score_mat_5[,1] <- apply(aic_weighted_arry_5, 3, sum)/713 # computing the weighted average

### Computing other possibile stability summary metrics to use for scoring confidence in the clustering model
aic_score_mat_5[,2] <-apply(aic_cluster_stab_arry_5, 3, mean, na.rm = T) # the means
aic_score_mat_5[,3] <-apply(aic_cluster_stab_arry_5, 3, median, na.rm = T) # the medians
aic_score_mat_5[,7] <-apply(aic_cluster_stab_arry_5, 3, which.max) #t he cluster with the max stability
aic_score_mat_5[,4] <-apply(aic_cluster_stab_arry_5, 3, max, na.rm = T) # the max stability
aic_score_mat_5[,8] <-apply(aic_cluster_stab_arry_5, 3, which.min) #t he cluster with the min stability
aic_score_mat_5[,5] <-apply(aic_cluster_stab_arry_5, 3, min, na.rm = T) # the min stability
aic_score_mat_5[,6] <-apply(aic_cluster_stab_arry_5, 3, sum, na.rm = T) # sum of stabilities
aic_score_mat_5[,9] <-aic_cluster_score_arry_5 # sums of top 10 -log10(q) scores

### Learning max values in each metric
which.max(aic_score_mat_5[,1])
which.max(aic_score_mat_5[,2])
which.max(aic_score_mat_5[,3])
which.max(aic_score_mat_5[,4])
which.max(aic_score_mat_5[,5])
which.max(aic_score_mat_5[,6])
which.max(aic_score_mat_5[,9])

### Ranking the clustering models based on their metric scores, column by column (not accumulative yet)
aic_score_mat_5_2 <- aic_score_mat_5 %>%
  as.data.frame() %>%
  dplyr::mutate(ranked_weighted_mean = rank(aic_score_mat_5[,1]),
                ranked_unweighted_mean = rank(aic_score_mat_5[,2]),
                ranked_median = rank(aic_score_mat_5[,3]),
                ranked_max = rank(aic_score_mat_5[,4]),
                ranked_min = rank(aic_score_mat_5[,5]), 
                ranked_sum = rank(aic_score_mat_5[,6]),
                ranked_sigscore = rank(aic_score_mat_5[,9]))

### Summing the rows of rank scores for different scoring regimes
aic_score_mat_5_3 <- aic_score_mat_5_2 %>%
  dplyr::mutate(weighted_mean_scores = rowSums(aic_score_mat_5_2[,c(10,12:16)]),
                unweighted_mean_scores = rowSums(aic_score_mat_5_2[,c(11:16)]),
                weighted_mean_no_sum_scores = rowSums(aic_score_mat_5_2[,c(10, 12:14, 16)]))

### Learning which model is the top pick of the chosen scoring regime
which.max(aic_score_mat_5_3$weighted_mean_no_sum_scores) 

### Writing out the final score sheet summarization of the analysis
aic_score_mat_5_4 <- aic_score_mat_5_3[,c(1,3:5,9:10,12:14,16)] %>%
  dplyr::rename('scores (out of 80)' = weighted_mean_no_sum_scores)

write.table(aic_score_mat_4, file="2022-05_30_aic_cluster_stability_score_sheet.txt", quote=F, row.names=T, col.names=T, sep=";")
write.csv(aic_score_mat_4, file = "2022-05-30_aic_cluster_stability_score_sheet.csv", row.names = T)


