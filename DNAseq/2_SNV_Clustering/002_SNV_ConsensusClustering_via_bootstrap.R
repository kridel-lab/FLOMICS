#---
# This script reads in mutation (CAPSEQ) data and performs clustering
# (updated to not include duplicated samples)
# Author: Victoria Shelton
# Last modified: February 1st, 2023
#---

#--
# Load in packages
#--
packages <- c("dplyr", "ggplot2", "flexmix", "tidyverse", "gridExtra")
lapply(packages, require, character.only = TRUE)

packages2 <- c("magrittr", "cluster", "cluster.datasets", "cowplot", "boot",
               "ggfortify", "dendextend", "factoextra", "FactoMineR", "corrplot",
               "GGally", "ggiraphExtra", "knitr", "kableExtra", "taRifx", "tidymodels",
               "gplots")
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
capseq_merged_mat <- read.table("2022-05-15_gene_vs_sample_E2408_TO_SNV_mat_matrix.txt", header = T, sep = ";")
print(capseq_merged_mat)

### remove the sample duplictaes
capseq_merged_mat <- capseq_merged_mat %>%
  dplyr::select(-c(LY_FL_1156_T1, LY_FL_1135_T1))
print("Number of samples in CAPSEQ Merged dataset")
ncol(capseq_merged_mat) # n=713

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/Dec2022")
dir <- "/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_VariantAnalysis/SNV_Clustering/UHN_E4402_PLOSMED-T-O_E2408-T-O/Dec2022"

### Extract samples
samples <- colnames(capseq_merged_mat)

### Change the class of values in merged dataset from "character" to "num" to input into downstream analyses
capseq_merg_num <- japply(capseq_merged_mat, which(sapply(capseq_merged_mat, class)=="character"), as.numeric )
summary(capseq_merg_num) %>% kable() %>% kable_styling()
print("made capseq_merg_num")

muts_all <- t(capseq_merged_mat) %>% data.frame() %>% mutate(SAMPLE_ID = row.names(.))
muts_all <- cbind(SAMPLE_ID = muts_all["SAMPLE_ID"], muts_all[,c(1:(ncol(muts_all)-1))])
row.names(muts_all) <- NULL
muts_df <- muts_all[, -1]
muts_df <- sapply(muts_df, as.numeric)
print("made muts_df")
dim(muts_df)


#--
# Consensus clustering following the method used in: https://doi.org/10.1182/blood.2019003535
#--

tmpar <- matrix(0, 713, 713) #this is a temporary matrix
ar_inclu <- array(tmpar, dim = c(713, 713, 1)) #creating an array to store the inclusion matrices.
ar_aic_conn <- array(tmpar, dim = c(713, 713, 1)) #creating an array to store the connectivity matrices.
ar_bic_conn <- array(tmpar, dim = c(713, 713, 1)) #creating an array to store the connectivity matrices.
clus_num_aic <- c() #creating a vector to store the number of aic clusters per bootstrap run
clus_num_bic <- c() #creating a vector to store the number of bic clusters per bootstrap run
print("made all matrices to be filled by bootstrap")

set.seed(99)
x2 <- as.integer(runif(n = 1000, min = 1, max = 1000)) #list of random numbers for the clustering seed
write.table(x2, file="clustering_seeds.txt", quote=F, row.names=F, col.names=F, sep=";")
numtries = 1000
for(i in 1:numtries){
  comment1 <- paste0("Bootstrap Run ", i)
  print(comment1)
  # bootstrap resample from the dataset: resample should be 713 samples, with replacement, so duplicates are expected
   muts_boot <- rsample::bootstraps(muts_df,
                            times = 1,
                            apparent = F)
  resamp <- c(muts_boot[[1]][[1]][["in_id"]]) #sample indices
  resamp_df <- capseq_merged_mat[,resamp] #collecting the SNV info for the bootstrapped samples
  ss <- samples[resamp]

  resamp_df_t <- t(resamp_df) %>% data.frame() %>% mutate(SAMPLE_ID = row.names(.))
  resamp_df_t <- cbind(SAMPLE_ID = resamp_df_t["SAMPLE_ID"], resamp_df_t[,c(1:(ncol(resamp_df_t)-1))])
  row.names(resamp_df_t) <- NULL
  resamp_df_t_noid <-resamp_df_t[, -1]
  resamp_df_t_noid <- sapply(resamp_df_t_noid, as.numeric)
  com <- paste0("Dimensions of Bootstrap ", i, ":")
  print(com)
  print(dim(resamp_df_t_noid))

  # create an inclusion matrix, where 1 means sample was in bootstrap resample,
  # and 0 means sample was not in bootstrap resample
  ### 1) make an empty matrix (713 x 713)
  inclu_mat <- matrix(, nrow = 713, ncol = 713)
  colnames(inclu_mat) <- samples
  rownames(inclu_mat) <- samples

  ### 2) Populate the matrix with "1"s for every pair of samples found in the boostrap resample
  for (col_id in ss){
    #print("column name:")
    #print(col_id)
    for (row_id in ss){
      #print("row name:")
      #print(row_id)
      inclu_mat[row_id,col_id] = 1
    }
    }
  inclu_mat[is.na(inclu_mat)] <- 0 #convert all 'NA' values to '0'
  comment2 <- paste0("Inclusion matrix ", i, " filled")
  print(comment2)
  #i = 1
  #Appending the inclusion matrix to the inclusion array
  ar_inclu <- array(c(ar_inclu, inclu_mat), dim = c(713, 713, (i + 1)))
  #ar_inclu[2,,2]

  ### 3) Fit the mixture model as before using the same list of variables on the bootstrapped sample
  set.seed(x2[i]) #setting seed for clustering
  boot_ex <- initFlexmix(resamp_df_t_noid ~ 1,
                    k = 2:10,
                    model = FLXMCmvbinary(),
                    control = list(minprior = 0.1),
                    nrep = 1000)
  aic <- getModel(boot_ex, which = "AIC")
  clus_num_aic <- c(clus_num_aic, aic@k0) #to return the number of clusters
  bic <- getModel(boot_ex, which = "BIC")
  clus_num_bic <- c(clus_num_bic, bic@k0) #to return the number of clusters
  icl <- getModel(boot_ex, which = "ICL")
  print("retrieved AIC, BIC & ICL models")

  ##### AIC #####
  ### 4) Extract cluster assignment probabilities for AIC
  aic_prob <- aic@posterior[["scaled"]]
  #rownames(aic_prob) <- samples
  rownames(aic_prob) <- ss
  colnames(aic_prob) <- seq(1,  ncol(aic_prob), 1)
  print("Extracted cluster assignment probabilities for AIC")

  ### Populating the interim connectivity matrices
  aic_add_array <- array(tmpar, dim = c(713, 713, 1)) #creating an array to store 'product' connectivity matrices
  for (col_num2 in seq(1, ncol(aic_prob), 1)){
    #print(col_num2)
    tmp_mat <- matrix(, nrow = 713, ncol = 713)
    colnames(tmp_mat) <- samples
    rownames(tmp_mat) <- samples
    for (col_id in ss){
      #print("column name:")
      #print(col_id)
      for (row_id in ss){
        #print("row name:")
        #print(row_id)
        #col_id = "LY_FL_001_T1"
        #row_id = "LY_FL_001_T1"
        clus_prob <- ((aic_prob[row_id, col_num2])*(aic_prob[col_id, col_num2]))
        tmp_mat[row_id,col_id] = clus_prob
        }

      }
    tmp_mat[is.na(tmp_mat)] <- 0
    assign(  paste0("aic_conn_mat_interim_", col_num2), tmp_mat )
    aic_add_array <- array(c(aic_add_array, tmp_mat), dim = c(713, 713, (col_num2 + 1)))
    }

  print("Populated the AIC interim connectivity matrices")

  sum_aic_interim_array <- rowSums(aic_add_array, dims = 2) #summing the 'product' matrices
  print("Summed the interim connectivity matrices for AIC")

  comment3 <- paste0("AIC Connectivity matrix ", i, " filled")
  print(comment3)

  #i = 1
  #Appending the inclusion matrix to the inclusion array
  ar_aic_conn <- array(c(ar_aic_conn, sum_aic_interim_array), dim = c(713, 713, (i + 1)))
  print("Appended the AIC inclusion matrix to the AIC inclusion array")

  ##### BIC #####
  ### 5) Extract cluster assignment probabilities for BIC
  bic_prob <- bic@posterior[["scaled"]]
  rownames(bic_prob) <- ss
  colnames(bic_prob) <- seq(1,  ncol(bic_prob), 1)
  print("Extracted cluster assignment probabilities for BIC")

  ### Populating the interim connectivity matrices
  bic_add_array <- array(tmpar, dim = c(713, 713, 1)) #creating an array to store 'product' connectivity matrices
  for (col_num2 in seq(1, ncol(bic_prob), 1)){
    #print(col_num2)
    tmp_mat <- matrix(, nrow = 713, ncol = 713)
    colnames(tmp_mat) <- samples
    rownames(tmp_mat) <- samples
    for (col_id in ss){
      #print("column name:")
      #print(col_id)
      for (row_id in ss){
        #print("row name:")
        #print(row_id)
        #col_id = "LY_FL_001_T1"
        #row_id = "LY_FL_001_T1"
        clus_prob <- ((bic_prob[row_id, col_num2])*(bic_prob[col_id, col_num2]))
        tmp_mat[row_id,col_id] = clus_prob
      }

    }
    tmp_mat[is.na(tmp_mat)] <- 0
    assign(  paste0("bic_conn_mat_interim_", col_num2), tmp_mat )
    bic_add_array <- array(c(bic_add_array, tmp_mat), dim = c(713, 713, (col_num2 + 1)))

  }
  print("Populated the BIC interim connectivity matrices")

  sum_bic_interim_array <- rowSums(bic_add_array, dims = 2) #summing the 'product' matrices
  print("Summed the interim connectivity matrices for BIC")

  comment3 <- paste0("BIC Connectivity matrix ", i, " filled")
  print(comment3)

  #i = 1
  #Appending the inclusion matrix to the inclusion array
  ar_bic_conn <- array(c(ar_bic_conn, sum_bic_interim_array), dim = c(713, 713, (i + 1)))
  print("Appended the BIC inclusion matrix to the BIC inclusion array")


}
# This works correctly and well :)) yay!
print("Completed bootstrapping")

### Summing the inclusion matrices
I_mat <- rowSums(ar_inclu, dims = 2) %>% set_colnames(samples) %>% set_rownames(samples) #summing the inclusion matrices
### Summing the AIC connectivity 'addition' matrices (the number of matrices should match the number of bootstrap samples)
aic_C_mat <- rowSums(ar_aic_conn, dims = 2) %>% set_colnames(samples) %>% set_rownames(samples)#summing the connectivity 'product' matrices
### Summing the BIC connectivity 'addition' matrices (the number of matrices should match the number of bootstrap samples)
bic_C_mat <- rowSums(ar_bic_conn, dims = 2) %>% set_colnames(samples) %>% set_rownames(samples)#summing the connectivity 'product' matrices
print("Completed summing the AIC and BIC connectivity 'addition' matrices")

### Generate AIC concensus matrix
aic_consensus_mat <- ((aic_C_mat)/(I_mat))
aic_consensus_mat[is.na(aic_consensus_mat)] <- 0 #'NA's become '0's
print("Created AIC concensus matrix")
a_mat_save = paste0(date, "_aic_consensus_mat.txt")
write.table(aic_consensus_mat , file=a_mat_save , quote=F, row.names=T, sep=";")
print("Saved AIC concensus matrix under this filename: ")
print(a_mat_save)

### Generate BIC concensus matrix
bic_consensus_mat <- ((bic_C_mat)/(I_mat))
bic_consensus_mat[is.na(bic_consensus_mat)] <- 0 #'NA's become '0's
print("Created BIC concensus matrix")
b_mat_save = paste0(date, "_bic_consensus_mat.txt")
write.table(bic_consensus_mat , file=b_mat_save , quote=F, row.names=T, sep=";")
print("Saved BIC concensus matrix under this filename: ")

### Let's try visualizing the consensus matrices
pdf(file = paste0(date, "_aic_consensuscluster.pdf"))
hm <- heatmap(aic_consensus_mat, symm = T)
dev.off()
print("saved AIC Consensus Cluster heatmap")
#-

pdf(file = paste0(date, "_bic_consensuscluster.pdf"))
hm <- heatmap(bic_consensus_mat, symm = T)
dev.off()
print("saved BIC Consensus Cluster heatmap")

#--
#  Looking at the proportion of the number of clusters identified by AIC and BIC over bootstrapping.
#--
### AIC # of clusters
aic_clus_num_freq <- as.data.frame(table(clus_num_aic))

aic_clus_num_total <- sum(aic_clus_num_freq["Freq"])
aic_clus_num_freq["Proportion"] <- format(round((((aic_clus_num_freq["Freq"])/(aic_clus_num_total))*100), 2), nsmall = 2)
aic_clus_num_freq_pro <- aic_clus_num_freq %>% dplyr::select(clus_num_aic, Proportion)
aic_clus_num_freq_pro <- sapply(aic_clus_num_freq_pro, as.numeric)

pdf(file = paste0(date, "_aic_cluster_barplot.pdf"))
ggplot(aic_clus_num_freq_pro, aes(x = clus_num_aic, y = Proportion)) + geom_bar(stat = "identity") + geom_text(aes(label=Proportion), vjust=1.6, color="white", size=3.5)
dev.off()
print("saved AIC cluster proportion plot")


### BIC # of clusters
bic_clus_num_freq <- as.data.frame(table(clus_num_bic))
bic_clus_num_total <- sum(bic_clus_num_freq["Freq"])
bic_clus_num_freq["Proportion"] <- format(round((((bic_clus_num_freq["Freq"])/(bic_clus_num_total))*100), 2), nsmall = 2)
bic_clus_num_freq_pro <- bic_clus_num_freq %>% dplyr::select(clus_num_bic, Proportion)
bic_clus_num_freq_pro <- as.data.frame(sapply(bic_clus_num_freq_pro, as.numeric))

pdf(file = paste0(date, "_bic_cluster_barplot.pdf"))
ggplot(bic_clus_num_freq_pro, aes(x = clus_num_bic, y = Proportion)) + geom_bar(stat = "identity") + geom_text(aes(label=Proportion), vjust=1.6, color="white", size=3.5)
dev.off()
print("saved BIC cluster proportion plot")
