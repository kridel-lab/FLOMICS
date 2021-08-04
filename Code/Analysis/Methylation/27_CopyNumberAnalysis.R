# Updated 3 Aug 2021
# Date: 12 September 2019
# Function: Find the association between copy number data and cluster labels using a chi-squared
#           contingency table and return results with adjusted p values for multiple testing.
# Author: Anjali Silva

# Input:
# CopyNumberData: A dataframe of size N x d, where N provides copy number sites (peak) in rows 
#                 and d provides the samples.
# ClusterLabels: A matrix with N rows and one column. The N rownames should correspond 
#                to sample identity. The column should provide the cluster label for each sample. 

# Output:
# SignificantResults: Significant results from Chi-squared contingency table tests from associating 
#                     cluster labels with each copy number site. Not adjusted for multiple testing.
# SignificantResults_Adjusted: Significant results from Chi-squared contingency table tests adjusted
#                              for multiple testing. 

CopyNumberAnalysis27 <- function(CopyNumberData, 
                                 ClusterLabels) {
  
  # Match copy number data with cluster order
  Match_Copynumber_ByCluster <- match(rownames(ClusterLabels), colnames(CopyNumberData))
  CopyNumberData_Ordered <- CopyNumberData[, Match_Copynumber_ByCluster[- 
                            which(is.na(Match_Copynumber_ByCluster) == TRUE)]]
  Match_Cluster_ByCopyNumber <- match(colnames(CopyNumberData), rownames(ClusterLabels))
  ClusterLabels_Ordered <- ClusterLabels[Match_Cluster_ByCopyNumber[- 
                           which(is.na(Match_Cluster_ByCopyNumber) == TRUE)], ]
  
  # For a data frame with combined cluster and copy number data
  Combined_Dataframe <- data.frame(ClusterLabels_Ordered, t(CopyNumberData_Ordered))
  cat("\n Analyzing data for", nrow(Combined_Dataframe)," samples across", 
      nrow(CopyNumberData_Ordered),"sites. \n")
  names(Combined_Dataframe) <- c("Cluster", rownames(CopyNumberData_Ordered))

  # Get association between each position and cluster labels  
  Chisq.test_CopyNumberData <- lapply(c(2:60), 
           function(i) stats::chisq.test(table(Combined_Dataframe$Cluster, Combined_Dataframe[[i]])))
  Chisq.test_CopyNumberData_Sig_PAdjusted <- 
    Chisq.test_CopyNumberData_Sig <- 
    matrix(unlist(lapply(c(2:60), function(w) 
           if(Chisq.test_CopyNumberData[[w - 1]]$p.value < 0.05) 
           { c(colnames(Combined_Dataframe)[w], 
             Chisq.test_CopyNumberData[[w - 1]]$p.value)})), 
             ncol = 2, byrow = TRUE)
  
  # Adjusting P value for multiple testing
  Chisq.test_CopyNumberData_Sig_PAdjusted[, 2] <- 
    stats::p.adjust(p = as.numeric(Chisq.test_CopyNumberData_Sig[, 2]), 
                                   method = "bonferroni", 
                                   n = length(as.numeric(Chisq.test_CopyNumberData_Sig[, 2])))

  # Printing tables
  Comparing_Tables <- lapply(1:length(Chisq.test_CopyNumberData_Sig_PAdjusted[, 1]), 
                             function(u) table(Combined_Dataframe$Cluster, Combined_Dataframe[, 
                             which(colnames(Combined_Dataframe) == 
                             Chisq.test_CopyNumberData_Sig_PAdjusted[u, 1])])) 
  names(Comparing_Tables) <- c(Chisq.test_CopyNumberData_Sig_PAdjusted[, 1])
  
  
  RESULTS <- list(SignificantResults = Chisq.test_CopyNumberData_Sig,
                  SignificantResults_Adjusted = Chisq.test_CopyNumberData_Sig_PAdjusted,
                  CrossTabulation_Tables = Comparing_Tables)
  
  class(RESULTS) <- "CopyNumberAnalysis"
  return(RESULTS)
}
# [END]
