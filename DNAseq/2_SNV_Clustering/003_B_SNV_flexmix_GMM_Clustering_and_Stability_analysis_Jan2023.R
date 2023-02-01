#---
# This script reads in mutation (CAPSEQ) data and performs clustering
# Authors: Robert Kridel + Victoria Shelton
# Last modified: Jan 9th, 2023
#---

#--
# Load in packages
#--
packages <- c("dplyr", "ggplot2", "flexmix", "tidyverse", "gridExtra")
lapply(packages, require, character.only = TRUE)

packages2 <- c("magrittr", "cluster", "cluster.datasets", "cowplot",
               "ggfortify", "dendextend", "factoextra", "FactoMineR", "corrplot",
               "GGally", "ggiraphExtra", "knitr", "kableExtra", "taRifx", "ggalluvial",
               "clValid", "mclust")
lapply(packages2, require, character.only = TRUE)

#--
# Setup working space
#--

date <- Sys.Date()
#setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/") #need to create these directories on the cluster
#dir <- "/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/"
setwd("~/Desktop/Kridel Lab/ConsensusClustering-Dec2022/")
dir <- "~/Desktop/Kridel Lab/ConsensusClustering-Dec2022/"

#---
# Read in data
#---

### load SNV dataset
#capseq_merg_mat_txt = paste0("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/May2022/2022-05-15_gene_vs_sample_E2408_TO_SNV_mat_matrix.txt")
capseq_merg_mat_txt = paste0("~/Desktop/Kridel Lab/ConsensusClustering/2022-05-15_gene_vs_sample_E2408_TO_SNV_mat_matrix.txt")
capseq_merged_mat <- read.table(capseq_merg_mat_txt, sep = ";")
print("Number of samples in CAPSEQ Merged dataset")
ncol(capseq_merged_mat) # n=715

## ADDITION: 2022-12-20
### remove duplicate samples from dataset
capseq_merged_mat <- capseq_merged_mat %>%
  dplyr::select(-c(LY_FL_1156_T1, LY_FL_1135_T1))
print("Number of samples in CAPSEQ Merged dataset")
ncol(capseq_merged_mat) # n=713

### load cohort sample IDs
#UHNsamples <- read.delim("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/May2022/capseq_UHN_samples.txt", header = F, sep = ",")
UHNsamples <- read.delim("~/Desktop/Kridel Lab/ConsensusClustering/capseq_UHN_samples.txt", header = F, sep = ",")
UHNsamples <- unlist(UHNsamples)

#E4402samples <- read.delim("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/May2022/capseq_E4402_samples.txt", header = F, sep = ",")
E4402samples <- read.delim("~/Desktop/Kridel Lab/ConsensusClustering/capseq_E4402_samples.txt", header = F, sep = ",")
E4402samples <- unlist(E4402samples)

#E2408samples <- read.delim("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/May2022/capseq_E2408_samples.txt", header = F, sep = ",")
E2408samples <- read.delim("~/Desktop/Kridel Lab/ConsensusClustering/capseq_E2408_samples.txt", header = F, sep = ",")
E2408samples <- unlist(E2408samples)

#PLOSMEDsamples <- read.delim("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/May2022/capseq_PLOSMED_samples.txt", header = F, sep = ",")
PLOSMEDsamples <- read.delim("~/Desktop/Kridel Lab/ConsensusClustering/capseq_PLOSMED_samples.txt", header = F, sep = ",")
PLOSMEDsamples <- unlist(PLOSMEDsamples)

message("Loaded in Cohort sample IDs")

### Read in gene panels to determine comon genes
#genes_PLOSMED <- read.csv("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/PROBES/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
genes_PLOSMED <- read.csv("~/Desktop/Kridel Lab/sample annotation/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
  dplyr::filter(PLOS_MED_PANEL == "YES") %>% dplyr::select(Gene.Name) # n=86 ?..Gene.Name
message("collected genes_PLOSMED")
print(genes_PLOSMED)

#genes_UHN <- read.csv("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/PROBES/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
genes_UHN <- read.csv("~/Desktop/Kridel Lab/sample annotation/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
  dplyr::filter(FLOMICS_PANEL == "YES") %>% dplyr::select(Gene.Name) # n=71 ?..Gene.Name
message("collected genes_UHN")
print(genes_UHN)

### Determine the common genes between the two utilized gene panels
genes_common <- intersect(genes_PLOSMED, genes_UHN) # n=57
message("collected common genes between genes_PLOSMED and genes_UHN")
print(genes_common)

### load AIC and BIC consensus matrices
#aic_consensus_mat_txt = paste0("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/May2022/2022-05-22_aic_consensus_mat.txt")
aic_consensus_mat_txt = paste0("~/Desktop/Kridel Lab/ConsensusClustering-Dec2022/2022-12-23/2022-12-23_aic_consensus_mat.txt")
aic_consensus_mat <- read.table(aic_consensus_mat_txt, header = T, sep = ";")
samples <- colnames(aic_consensus_mat)
message("Number of samples in AIC Consensus Matrix")
ncol(aic_consensus_mat) # n=713

#bic_consensus_mat_txt = paste0("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/May2022/2022-05-22_bic_consensus_mat.txt")
bic_consensus_mat_txt = paste0("~/Desktop/Kridel Lab/ConsensusClustering-Dec2022/2022-12-23/2022-12-23_bic_consensus_mat.txt")
bic_consensus_mat <- read.table(bic_consensus_mat_txt, header = T, sep = ";")
samples <- colnames(bic_consensus_mat)
message("Number of samples in BIC Consensus Matrix")
ncol(bic_consensus_mat) # n=713

#--
# Clustering
# method: https://github.com/ecsg-uoy/DLBCLGenomicSubtyping
#--

### loading in functions
#source("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/VS-script-plot-all-mutations-all-cohorts_to_be_loaded.R")
source("~/Desktop/Kridel Lab/filtering_muts/VS-script-plot-all-mutations-all-cohorts_to_be_loaded_3.R")
message("loaded in plotting functions")

### set directory again
#setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/May2022/cluster_stability")
setwd("~/Desktop/Kridel Lab/ConsensusClustering-Dec2022/ClusterStability")

### formating dataset for entry in flexix modeling
muts_all <- t(capseq_merged_mat) %>% data.frame() %>% mutate(SAMPLE_ID = row.names(.))
ncol(muts_all) # 58
nrow(muts_all) #713
muts_all <- cbind(SAMPLE_ID = muts_all$SAMPLE_ID, muts_all[,c(1:(ncol(muts_all)-1))])
row.names(muts_all) <- NULL
muts_df <- muts_all[, -1]
muts_df <- sapply(muts_df, as.numeric)
message("Dimensions of formated dataset:")
dim(muts_df)
#2022-05-30 = 715 57
#2022-12-20 = 713 57, removed duplicates, -c(LY_FL_1156_T1, LY_FL_1135_T1)

### creating lists to store the cluster stability outputs
tmpar <- matrix(0, 1, 7) #this is a temporary matrix
aic_cluster_stab_arry_5 <- array(tmpar, dim = c(1, 5, 20))
bic_cluster_stab_arry_5 <- array(tmpar, dim = c(1, 2, 20))
icl_cluster_stab_arry_5 <- array(tmpar, dim = c(1, 2, 20))

aic_cluster_stab_arry_6 <- array(tmpar, dim = c(1, 6, 20))
bic_cluster_stab_arry_6 <- array(tmpar, dim = c(1, 2, 20))
icl_cluster_stab_arry_6 <- array(tmpar, dim = c(1, 2, 20))

### creating lists to store the numbers of patients per cluster 
aic_cluster_patnum_arry_5 <- array(tmpar, dim = c(1, 5, 20))
bic_cluster_patnum_arry_5 <- array(tmpar, dim = c(1, 2, 20))
icl_cluster_patnum_arry_5 <- array(tmpar, dim = c(1, 2, 20))

aic_cluster_patnum_arry_6 <- array(tmpar, dim = c(1, 6, 20))
bic_cluster_patnum_arry_6 <- array(tmpar, dim = c(1, 2, 20))
icl_cluster_patnum_arry_6 <- array(tmpar, dim = c(1, 2, 20))

### creating lists to store the sum of the -log10(q) scores
aic_cluster_score_arry_5 <- matrix(0, 1, 20)
bic_cluster_score_arry_5 <- matrix(0, 1, 20)
icl_cluster_score_arry_5 <- matrix(0, 1, 20)

aic_cluster_score_arry_6 <- matrix(0, 1, 20)
bic_cluster_score_arry_6 <- matrix(0, 1, 20)
icl_cluster_score_arry_6 <- matrix(0, 1, 20)

### creating seeds for cluster runs
x3 <- as.integer(runif(n = 20, min = 1, max = 1000)) #list of random numbers for the clustering seed
write.table(x3, file="clustering_seeds.txt", quote=F, row.names=F, col.names=F, sep=";")

#new 5 cluster seeds = 671, 504, 890, 459,  93, 497, 25, 319, 888, 729, 106, 572, 312, 80, 310, 338, 273, 302, 984, 138
x3 <- c(671, 504, 890, 459,  93, 497, 25, 319, 888, 729, 106, 572, 312, 80, 310, 338, 273, 302, 984, 138)

#new 6 cluster seeds = 288, 781, 713, 278, 408, 691, 698, 148, 215, 743, 262, 201, 260, 512, 300, 511, 99, 540, 495, 980
x3 <- c(288, 781, 713, 278, 408, 691, 698, 148, 215, 743, 262, 201, 260, 512, 300, 511, 99, 540, 495, 980)


### set number of iterations
numtries = 20

### set overall clustering seed
set.seed(99)

### clustering
RNGversion("4.0.2")

### starting iterations
for(i in 1:numtries){
set.seed(x3[i]) #setting seed for clustering
  set.seed(671)
  muts_all <- t(capseq_merged_mat) %>% data.frame() %>% mutate(SAMPLE_ID = row.names(.))
  muts_all <- cbind(SAMPLE_ID = muts_all$SAMPLE_ID, muts_all[,c(1:(ncol(muts_all)-1))])
  ex <- initFlexmix(muts_df ~ 1,
                  k = c(2,5),
                  model = FLXMCmvbinary(),
                  control = list(minprior = 0.1),
                  nrep = 1000)
save(ex, file = "clustering_ex_December2022.RData")

### saving the Number of Components plot
message("finished Clustering")
pdf(file = paste0(date, "_", i, "_all_mutations_all_cohorts_E2408-TO_GMM_num-components-plot_5-clusters.pdf"))
plot(ex)
dev.off()
message("GMM Clustering Number of Components plot is saved in working directory as:")
print(paste0(date, "_", i, "_all_mutations_all_cohorts_E2408-TO_GMM_num-components-plot.pdf"))

aic <- getModel(ex, which = "AIC")
bic <- getModel(ex, which = "BIC")
icl <- getModel(ex, which = "ICL")
message("retrieved AIC, BIC & ICL models")
save(aic, file = "clustering_aic_December2022.RData")

### appending cluster assignments
muts_all$ClusterAIC <- as.factor(paste0('C', aic@cluster))
table(muts_all$ClusterAIC)
#best outcome of seed 65
#C1  C2  C3  C4  C5 
#96 215  78 144 182
## NEW ADDITION 2022-12-20
#C1  C2  C3  C4  C5 
#89 268 141  91 124

muts_all$ClusterBIC <- as.factor(paste0('C', bic@cluster))
table(muts_all$ClusterBIC)
#C1  C2 
#176 539 
## NEW ADDITION 2022-12-20
#C1  C2 
#538 175

muts_all$ClusterICL <- as.factor(paste0('C', icl@cluster))
table(muts_all$ClusterICL)
#C1  C2 
#176 539
## NEW ADDITION 2022-12-20
#C1  C2 
#538 175

### generating and saving cluster maps and top 10 -log10(q) values for each criteria
genes <- colnames(muts_df)

# AIC
criteria <- "AIC"
plt_aic <- heatmap_mutation_extended(muts_all, genes, 'ClusterAIC', y_order = 'group', idcol = 'SAMPLE_ID')
grid.arrange(plt_aic)
ggsave(file = paste0(date, "_", i, "_heatmap_aic.png"), plt_aic, width = 25, height = 25, units = "cm")

sig_long_out_path = paste0("~/Desktop/Kridel Lab/ConsensusClustering-Dec2022/ClusterStability/", date, "_", i, "_sig_long_out_", criteria, ".txt")
sig_long_out_mat <- read.table(sig_long_out_path, sep = ";", header = T) %>% dplyr::filter(to_plot == TRUE) #the to_plot filter keeps only genes that are considered significant to clustering outcome (usually filter in loaded functions, but this is a safe guard)
message("saved AIC sum of all -log10(q) scores, the -log10 of the adjusted p-value (fdr adjusted)")
sig_long_out_mat_sort <- sig_long_out_mat[with(sig_long_out_mat,order(-val)),] #sorting entries by val in descending order
aic_cluster_score_arry_5[,i] <- sig_long_out_mat_sort[1:10,] %>% dplyr::select(val) %>% sum() #summing the top 10 entries

# BIC
criteria <- "BIC"
plt_bic <- heatmap_mutation_extended(muts_all, genes, 'ClusterBIC', y_order = 'group', idcol = 'SAMPLE_ID')
grid.arrange(plt_bic)
ggsave(file = paste0(date, "_", i, "_heatmap_bic.png"), plt_bic, width = 25, height = 25, units = "cm")

sig_long_out_path = paste0("~/Desktop/Kridel Lab/ConsensusClustering-Dec2022/ClusterStability/", date, "_", i, "_sig_long_out_", criteria, ".txt")
sig_long_out_mat <- read.table(sig_long_out_path, sep = ";", header = T) %>% dplyr::filter(to_plot == TRUE) #the to_plot filter keeps only genes that are considered significant to clustering outcome (usually filter in loaded functions, but this is a safe guard)
message("saved BIC sum of all -log10(q) scores, the -log10 of the adjusted p-value (fdr adjusted)")
sig_long_out_mat_sort <- sig_long_out_mat[with(sig_long_out_mat,order(-val)),] #sorting entries by val in descending order
bic_cluster_score_arry_5[,i] <- sig_long_out_mat_sort[1:10,] %>% dplyr::select(val) %>% sum() #summing the top 10 entries

# ICL
criteria <- "ICL"
plt_icl <- heatmap_mutation_extended(muts_all, genes, 'ClusterICL', y_order = 'group', idcol = 'SAMPLE_ID')
grid.arrange(plt_icl)
ggsave(file = paste0(date, "_", i, "_heatmap_icl.png"), plt_icl, width = 25, height = 25, units = "cm")
message("saved all model heatmaps")

sig_long_out_path = paste0("~/Desktop/Kridel Lab/ConsensusClustering-Dec2022/ClusterStability/", date, "_", i, "_sig_long_out_", criteria, ".txt")
sig_long_out_mat <- read.table(sig_long_out_path, sep = ";", header = T) %>% dplyr::filter(to_plot == TRUE) #the to_plot filter keeps only genes that are considered significant to clustering outcome (usually filter in loaded functions, but this is a safe guard)
message("saved ICl sum of all -log10(q) scores, the -log10 of the adjusted p-value (fdr adjusted)")
sig_long_out_mat_sort <- sig_long_out_mat[with(sig_long_out_mat,order(-val)),] #sorting entries by val in descending order
icl_cluster_score_arry_5[,i] <- sig_long_out_mat_sort[1:10,] %>% dplyr::select(val) %>% sum() #summing the top 10 entries
message("saved -log10(q) scores, summed top 10")

### attaching cohort labels to samples
muts_all_sub <- muts_all %>%
  dplyr::mutate(PATIENT_ID = substr(SAMPLE_ID, 1, 9)) %>%
  dplyr::select(PATIENT_ID, SAMPLE_ID, ClusterAIC, ClusterBIC) %>%
  dplyr::mutate(cohort = ifelse(SAMPLE_ID %in% UHNsamples, "UHN",
                         ifelse(SAMPLE_ID %in% E4402samples, "E4402",
                                ifelse(SAMPLE_ID %in% E2408samples, "E2408",
                                       ifelse(SAMPLE_ID %in% PLOSMEDsamples, "PLOSMED", SAMPLE_ID)))))

write.csv(muts_all_sub, file = paste0(date, "_", i, "_GMM_Cluster_Labels_flexmix_clusters.csv"), row.names = FALSE)
message("saved cluster labels")

### Table of how many samples are in each cluster, from each cohort
# AIC
table(muts_all_sub$ClusterAIC, muts_all_sub$cohort)
# BIC
table(muts_all_sub$ClusterBIC, muts_all_sub$cohort)

### Table displaying the proportion of samples from each cluster and from each cohort
# AIC
round(prop.table(table(muts_all_sub$ClusterAIC, muts_all_sub$cohort), margin = 2),2)*100
# BIC
round(prop.table(table(muts_all_sub$ClusterBIC, muts_all_sub$cohort), margin = 2),2)*100

### Looking the chisquare statistic
# AIC
chisq.test(muts_all_sub$ClusterAIC, muts_all_sub$cohort)
# BIC 
chisq.test(muts_all_sub$ClusterBIC, muts_all_sub$cohort)

### A "AIC Cluster label" x "BIC Cluster label" table
table(muts_all_sub$ClusterAIC, muts_all_sub$ClusterBIC)

#------------------------------------------------------
# Assesing the Stability of Individual Clusters
#------------------------------------------------------

### AIC Cluster Assignment ###

### Table of cluster assignments
aic_clus_assign <- data.frame(matrix(nrow=713)) %>% add_column(Samples = colnames(capseq_merged_mat), Cluster = ex@models[["5"]]@cluster) %>% dplyr::select(-1)

### extracting the samples in each cluster, the number of samples and the sample identifiers
aic_clus1 <- aic_clus_assign %>% dplyr::filter(Cluster == "1") %>% dplyr::select(1)
aic_patnum_1 <- dim(aic_clus1) 
aic_clus1_samples <- aic_clus1[,1]

aic_clus2 <- aic_clus_assign %>% dplyr::filter(Cluster == "2") %>% dplyr::select(1)
aic_patnum_2 <- dim(aic_clus2) 
aic_clus2_samples <- aic_clus2[,1]

aic_clus3 <- aic_clus_assign %>% dplyr::filter(Cluster == "3") %>% dplyr::select(1)
aic_patnum_3 <- dim(aic_clus3) 
aic_clus3_samples <- aic_clus3[,1]

aic_clus4 <- aic_clus_assign %>% dplyr::filter(Cluster == "4") %>% dplyr::select(1)
aic_patnum_4 <- dim(aic_clus4) 
aic_clus4_samples <- aic_clus4[,1]

aic_clus5 <- aic_clus_assign %>% dplyr::filter(Cluster == "5") %>% dplyr::select(1)
aic_patnum_5 <- dim(aic_clus5) 
aic_clus5_samples <- aic_clus5[,1]

#aic_clus6 <- aic_clus_assign %>% dplyr::filter(Cluster == "6") %>% dplyr::select(1)
#aic_patnum_6 <- dim(aic_clus6) 
#aic_clus6_samples <- aic_clus6[,1]

#aic_clus7 <- aic_clus_assign %>% dplyr::filter(Cluster == "7") %>% dplyr::select(1)
#aic_patnum_7 <- dim(aic_clus7) 
#aic_clus7_samples <- aic_clus7[,1]

### storing number of patients
aic_patnum_mat <- cbind(aic_patnum_1, aic_patnum_2, aic_patnum_3, aic_patnum_4, aic_patnum_5)#, aic_patnum_6, aic_patnum_7)
aic_cluster_patnum_arry_5[,,i] <- aic_patnum_mat[1,]
write.csv(aic_patnum_mat, file = paste0(date, "_", i, "_AIC_Cluster_Patient_Number_Matrix.csv"), row.names = FALSE)
message("saved aic_patnum_mat")

### AIC CLUS1 MATRIX SUM ###
aic_clus1_sum_mat <- matrix(0, length(aic_clus1_samples), length(aic_clus1_samples)) %>% set_colnames(aic_clus1_samples) %>% set_rownames(aic_clus1_samples) #create matrix to store consensus

for (col_id in aic_clus1_samples){
  #print("column name:")
  #print(col_id)
  for (row_id in aic_clus1_samples){
    #print("row name:")
    #print(row_id)
    aic_clus1_sum_mat[row_id,col_id] <- aic_consensus_mat[row_id,col_id]
  }
}

aic_clus1_sum <- sum(aic_clus1_sum_mat) #sum up the AIC Cluster 1 consensus matrix
message("made aic_clus1_sum matrix")

### AIC CLUS2 MATRIX SUM ###
aic_clus2_sum_mat <- matrix(0, length(aic_clus2_samples), length(aic_clus2_samples)) %>% set_colnames(aic_clus2_samples) %>% set_rownames(aic_clus2_samples) #create matrix to store consensus

for (col_id in aic_clus2_samples){

  for (row_id in aic_clus2_samples){

    aic_clus2_sum_mat[row_id,col_id] <- aic_consensus_mat[row_id,col_id]
  }
}

aic_clus2_sum <- sum(aic_clus2_sum_mat) #sum up the AIC Cluster 2 consensus matrix

### AIC CLUS3 MATRIX SUM ###
aic_clus3_sum_mat <- matrix(0, length(aic_clus3_samples), length(aic_clus3_samples)) %>% set_colnames(aic_clus3_samples) %>% set_rownames(aic_clus3_samples) #create matrix to store consensus

for (col_id in aic_clus3_samples){

  for (row_id in aic_clus3_samples){

    aic_clus3_sum_mat[row_id,col_id] <- aic_consensus_mat[row_id,col_id]
  }
}

aic_clus3_sum <- sum(aic_clus3_sum_mat) #sum up the AIC Cluster 3 consensus matrix

### AIC CLUS4 MATRIX SUM ###
aic_clus4_sum_mat <- matrix(0, length(aic_clus4_samples), length(aic_clus4_samples)) %>% set_colnames(aic_clus4_samples) %>% set_rownames(aic_clus4_samples) #create matrix to store consensus

for (col_id in aic_clus4_samples){

  for (row_id in aic_clus4_samples){

    aic_clus4_sum_mat[row_id,col_id] <- aic_consensus_mat[row_id,col_id]
  }
}

aic_clus4_sum <- sum(aic_clus4_sum_mat) #sum up the AIC Cluster 4 consensus matrix

### AIC CLUS5 MATRIX SUM ###
aic_clus5_sum_mat <- matrix(0, length(aic_clus5_samples), length(aic_clus5_samples)) %>% set_colnames(aic_clus5_samples) %>% set_rownames(aic_clus5_samples) #create matrix to store consensus

for (col_id in aic_clus5_samples){

  for (row_id in aic_clus5_samples){

    aic_clus5_sum_mat[row_id,col_id] <- aic_consensus_mat[row_id,col_id]
  }
}

aic_clus5_sum <- sum(aic_clus5_sum_mat) #sum up the AIC Cluster 5 consensus matrix

### AIC CLUS6 MATRIX SUM ###
#aic_clus6_sum_mat <- matrix(0, length(aic_clus6_samples), length(aic_clus6_samples)) %>% set_colnames(aic_clus6_samples) %>% set_rownames(aic_clus6_samples) #create matrix to store consensus

#for (col_id in aic_clus6_samples){

#  for (row_id in aic_clus6_samples){

#    aic_clus6_sum_mat[row_id,col_id] <- aic_consensus_mat[row_id,col_id]
#  }
#}

#aic_clus6_sum <- sum(aic_clus6_sum_mat) #sum up the AIC Cluster 6 consensus matrix

### AIC CLUS7 MATRIX SUM ###
#aic_clus7_sum_mat <- matrix(0, length(aic_clus7_samples), length(aic_clus7_samples)) %>% set_colnames(aic_clus7_samples) %>% set_rownames(aic_clus7_samples) #create matrix to store consensus

#for (col_id in aic_clus7_samples){

#  for (row_id in aic_clus7_samples){

#    aic_clus7_sum_mat[row_id,col_id] <- aic_consensus_mat[row_id,col_id]
#  }
#}

#aic_clus7_sum <- sum(aic_clus7_sum_mat) #sum up the AIC Cluster 7 consensus matrix

### AIC Cluster Stabilities:
aic_clus1_stab <- ((1/ ((length(aic_clus1_samples)) * ((length(aic_clus1_samples)) -1) /2)) * (aic_clus1_sum)) #computing the consensus of AIC Cluster 1
message("AIC CLUSTER 1 STABILITY:")
print(aic_clus1_stab)

aic_clus2_stab <- ((1/ ((length(aic_clus2_samples)) * ((length(aic_clus2_samples)) -1) /2)) * (aic_clus2_sum)) #computing the consensus of AIC Cluster 2
message("AIC CLUSTER 2 STABILITY:")
print(aic_clus2_stab)

aic_clus3_stab <- ((1/ ((length(aic_clus3_samples)) * ((length(aic_clus3_samples)) -1) /2)) * (aic_clus3_sum)) #computing the consensus of AIC Cluster 3
message("AIC CLUSTER 3 STABILITY:")
print(aic_clus3_stab)

aic_clus4_stab <- ((1/ ((length(aic_clus4_samples)) * ((length(aic_clus4_samples)) -1) /2)) * (aic_clus4_sum)) #computing the consensus of AIC Cluster 4
message("AIC CLUSTER 4 STABILITY:")
print(aic_clus4_stab)

aic_clus5_stab <- ((1/ ((length(aic_clus5_samples)) * ((length(aic_clus5_samples)) -1) /2)) * (aic_clus5_sum)) #computing the consensus of AIC Cluster 5
message("AIC CLUSTER 5 STABILITY:")
print(aic_clus5_stab)

#aic_clus6_stab <- ((1/ ((length(aic_clus6_samples)) * ((length(aic_clus6_samples)) -1) /2)) * (aic_clus6_sum)) #computing the consensus of AIC Cluster 6
#print("AIC CLUSTER 6 STABILITY:")
#print(aic_clus6_stab)

#aic_clus7_stab <- ((1/ ((length(aic_clus7_samples)) * ((length(aic_clus7_samples)) -1) /2)) * (aic_clus7_sum)) #computing the consensus of AIC Cluster 7
#print("AIC CLUSTER 7 STABILITY:")
#print(aic_clus7_stab)

### A matrix of the cluster stability scores
aic_stab_mat <- cbind(aic_clus1_stab, aic_clus2_stab, aic_clus3_stab, aic_clus4_stab, aic_clus5_stab)#, aic_clus6_stab, aic_clus7_stab)

aic_cluster_stab_arry_5[,,i] <- aic_stab_mat # storing the cluster stability scores matrix in the designated array
write.csv(aic_stab_mat, file = paste0(date, "_", i, "_AIC_Cluster_Stability_Matrix.csv"), row.names = FALSE)
message("saved aic_stab_mat")


### BIC Cluster Assignment ###
bic_clus_assign <- data.frame(matrix(nrow=713)) %>% add_column(Samples = samples, Cluster = ex@models[["2"]]@cluster) %>% dplyr::select(-1)

bic_clus1 <- bic_clus_assign %>% dplyr::filter(Cluster == "1") %>% dplyr::select(1)
bic_patnum_1 <- dim(bic_clus1) 
bic_clus1_samples <- bic_clus1[,1]

bic_clus2 <- bic_clus_assign %>% dplyr::filter(Cluster == "2") %>% dplyr::select(1)
bic_patnum_2 <- dim(bic_clus2) 
bic_clus2_samples <- bic_clus2[,1]

### storing number of patients
bic_patnum_mat <- cbind(bic_patnum_1, bic_patnum_2)
bic_cluster_patnum_arry_5[,,i] <- bic_patnum_mat[1,]
write.csv(bic_patnum_mat, file = paste0(date, "_", i, "_BIC_Cluster_Patient_Number_Matrix.csv"), row.names = FALSE)
message("saved bic_patnum_mat")

### BIC CLUS1 MATRIX SUM ###
bic_clus1_sum_mat <- matrix(0, length(bic_clus1_samples), length(bic_clus1_samples)) %>% set_colnames(bic_clus1_samples) %>% set_rownames(bic_clus1_samples) # create matrix to store consensus

for (col_id in bic_clus1_samples){

  for (row_id in bic_clus1_samples){

    bic_clus1_sum_mat[row_id,col_id] <- bic_consensus_mat[row_id,col_id]
  }
}

bic_clus1_sum <- sum(bic_clus1_sum_mat) # sum up the BIC Cluster 1 consensus matrix


### BIC CLUS2 MATRIX SUM ###
bic_clus2_sum_mat <- matrix(0, length(bic_clus2_samples), length(bic_clus2_samples)) %>% set_colnames(bic_clus2_samples) %>% set_rownames(bic_clus2_samples) # create matrix to store consensus

for (col_id in bic_clus2_samples){

  for (row_id in bic_clus2_samples){

    bic_clus2_sum_mat[row_id,col_id] <- bic_consensus_mat[row_id,col_id]
  }
}

bic_clus2_sum <- sum(bic_clus2_sum_mat) # sum up the BIC Cluster 2 consensus matrix

### BIC Cluster Stabilities:
bic_clus1_stab <- ((1/ ((length(bic_clus1_samples)) * ((length(bic_clus1_samples)) -1) /2)) * (bic_clus1_sum)) # computing the consensus of BIC Cluster 1
message("BIC CLUSTER 1 STABILITY:")
print(bic_clus1_stab)

bic_clus2_stab <- ((1/ ((length(bic_clus2_samples)) * ((length(bic_clus2_samples)) -1) /2)) * (bic_clus2_sum)) # computing the consensus of BIC Cluster 2
message("BIC CLUSTER 2 STABILITY:")
print(bic_clus2_stab)

bic_stab_mat <- cbind(bic_clus1_stab, bic_clus2_stab)
bic_cluster_stab_arry_5[,,i] <- bic_stab_mat
write.csv(bic_stab_mat, file = paste0(date, "_", i,  "_BIC_Cluster_Stability_Matrix.csv"), row.names = FALSE)

message("saved bic_stab_mat")

###### -----------------------

### ICL Cluster Assignment ###
# icl_clus_assign <- data.frame(matrix(nrow=713)) %>% add_column(Samples = samples, Cluster = ex@models[["2"]]@cluster) %>% dplyr::select(-1)
# 
# icl_clus1 <- icl_clus_assign %>% dplyr::filter(Cluster == "1") %>% dplyr::select(1)
# icl_patnum_1 <- dim(icl_clus1) 
# icl_clus1_samples <- icl_clus1[,1]
# 
# icl_clus2 <- icl_clus_assign %>% dplyr::filter(Cluster == "2") %>% dplyr::select(1)
# icl_patnum_2 <- dim(icl_clus2) 
# icl_clus2_samples <- icl_clus2[,1]
# 
# #storing number of patients
# icl_patnum_mat <- cbind(icl_patnum_1, icl_patnum_2)
# icl_cluster_patnum_arry_5[,,i] <- icl_patnum_mat[1,]
# write.csv(icl_patnum_mat, file = paste0(date, "_", i, "_ICL_Cluster_Patient_Number_Matrix.csv"), row.names = FALSE)
# print("saved icl_patnum_mat")
# 
# 
# ### ICL CLUS1 MATRIX SUM ###
# icl_clus1_sum_mat <- matrix(0, length(icl_clus1_samples), length(icl_clus1_samples)) %>% set_colnames(icl_clus1_samples) %>% set_rownames(icl_clus1_samples)#create matrix to store consensus
# 
# for (col_id in icl_clus1_samples){

#   for (row_id in icl_clus1_samples){

#     icl_clus1_sum_mat[row_id,col_id] <- icl_consensus_mat[row_id,col_id]
#   }
# }
# 
# icl_clus1_sum <- sum(icl_clus1_sum_mat) #sum up the ICL Cluster 1 consensus matrix
# 
# 
# ### ICL CLUS2 MATRIX SUM ###
# icl_clus2_sum_mat <- matrix(0, length(icl_clus2_samples), length(icl_clus2_samples)) %>% set_colnames(icl_clus2_samples) %>% set_rownames(icl_clus2_samples)#create matrix to store consensus
# 
# for (col_id in icl_clus2_samples){

#   for (row_id in icl_clus2_samples){

#     icl_clus2_sum_mat[row_id,col_id] <- icl_consensus_mat[row_id,col_id]
#   }
# }
# 
# icl_clus2_sum <- sum(icl_clus2_sum_mat) #sum up the ICL Cluster 2 consensus matrix
# 
# ### ICL Cluster Stabilities:
# icl_clus1_stab <- ((1/ ((length(icl_clus1_samples)) * ((length(icl_clus1_samples)) -1) /2)) * (icl_clus1_sum)) #computing the consensus of ICL Cluster 1
# print("ICL CLUSTER 1 STABILITY:")
# print(icl_clus1_stab)
# 
# icl_clus2_stab <- ((1/ ((length(icl_clus2_samples)) * ((length(icl_clus2_samples)) -1) /2)) * (icl_clus2_sum)) #computing the consensus of ICL Cluster 2
# print("ICL CLUSTER 2 STABILITY:")
# print(icl_clus2_stab)
# 
# icl_stab_mat <- cbind(icl_clus1_stab, icl_clus2_stab)
# icl_cluster_stab_arry_5[,,i] <- icl_stab_mat
# write.csv(icl_stab_mat, file = paste0(date, "_", i,  "_ICL_Cluster_Stability_Matrix.csv"), row.names = FALSE)
# 
# print("saved icl_stab_mat")

###### -----------------------

}

### writing out the stability arrays
setwd("~/Desktop/Kridel Lab/ConsensusClustering-Dec2022/ClusterStability")
#write.csv(icl_cluster_stab_arry, file = paste0(date, "_ICL_Cluster_Stability_Array.csv"), row.names = FALSE)
write.csv(bic_cluster_stab_arry, file = paste0(date, "_BIC_Cluster_Stability_Array.csv"), row.names = FALSE)
write.csv(aic_cluster_stab_arry, file = paste0(date, "_AIC_Cluster_Stability_Array.csv"), row.names = FALSE)

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


