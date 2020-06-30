# Updated 9 April 2019
# Function: Summarize descriptive statistics of all datasets (Beta and Mvalues)
# Author: Anjali Silva

# Input:
# BetaMatrixProbes: A matrix of beta values for probes x patients, with probes in rows and patients as columns
# BetaMatrixGenes: A matrix of beta values for genes x patients, with genes in rows and patients as columns
# ClinicalCategoryToVisualize: Specify column of ClinicalFile to use. Options "STAGE", "SEX", "SITE_BIOPSY", "TYPE_BIOPSY", "INSTITUTION", "COO", "TYPE", "TRANSLOC_14_18"
# MvalueMatrixProbes: A matrix of M-values for probes x patients, with probes in rows and patients as columns
# MvalueMatrixGenes A matrix of M-values for genes x patients, with genes in rows and patients as columns

# Output: 
# BetaMatrixProbesSummary: list with summary statistics for Beta Matrix of probes including overall mean, patient means, overall median, patient medians, range, quantile, IQR, variance, standard deviation, median absolute deviation
# BetaMatrixGenesSummary: list with summary statistics for Beta Matrix of genes including overall mean, patient means, overall median, patient medians, range, quantile, IQR, variance, standard deviation, median absolute deviation
# MvalueMatrixProbesSummary: list with summary statistics for Mvalue Matrix of probes including overall mean, patient means, overall median, patient medians, range, quantile, IQR, variance, standard deviation, median absolute deviation
# MvalueMatrixGenesSummary: list with summary statistics for Mvalue Matrix of probes including overall mean, patient means, overall median, patient medians, range, quantile, IQR, variance, standard deviation, median absolute deviation
# 6_SummaryStatistics_BoxPlot.p* 

SummaryStatistics <- function(BetaMatrixProbes = NA, 
                              BetaMatrixGenes = NA, 
                              ClinicalCategoryToVisualize = NA, 
                              ClinicalFile = NA, 
                              MvalueMatrixProbes = NA, 
                              MvalueMatrixGenes = NA, 
                              ProduceImages = "Yes", 
                              PNGorPDF = "png") {
  
  # Loading needed packages
  # RegularPckgs=c("pastecs","ggpubr")
  library(pastecs)
  library(ggpubr)
  library(BiocGenerics)
  
  # Define a function for finding mode
  # https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  # Measure of central tendency (mean, median, # mode)
  # Measure of variablity (range, quantile, interquartile range, variance and standard deviation, median absolute deviation)
  if(typeof(BetaMatrixProbes) == "double") {
    overall_mean_BetaMatrixProbes <- mean(as.matrix(BetaMatrixProbes))
    patient_means_BetaMatrixProbes <- colMeans(as.matrix(BetaMatrixProbes))
    overall_median_BetaMatrixProbes <- median(as.matrix(BetaMatrixProbes))
    patient_median_BetaMatrixProbes <- apply(as.matrix(BetaMatrixProbes), 2, median)
    # Mode commented out as it takes a long time
    # overall_mode_BetaMatrixProbes <- Mode(round(BetaMatrixProbes,2))
    # patient_mode_BetaMatrixProbes <- apply(BetaMatrixProbes,2,Mode)
    range_BetaMatrixProbes <- range(as.matrix(BetaMatrixProbes))
    quantile_BetaMatrixProbes <- quantile(as.matrix(BetaMatrixProbes))
    IQR_BetaMatrixProbes <- IQR(as.matrix(BetaMatrixProbes))
    var_BetaMatrixProbes <- var(as.matrix(BetaMatrixProbes))
    sd_BetaMatrixProbes <- sd(as.matrix(BetaMatrixProbes))
    # median absolute deviation
    mad_BetaMatrixProbes <- mad(as.matrix(BetaMatrixProbes)) 
    
    if(ProduceImages == "Yes") {
      # Checking spread of data 
      # Boxplot organized by Type
      # This boxplot is specific to type: FL, DLBCL, RLN
      TypeLymphomVector <- c(which(substr(colnames(BetaMatrixProbes), 4, 5 ) == "FL"), 
                             which(substr(colnames(BetaMatrixProbes), 4, 5 ) == "DL"),
                             which(substr(colnames(BetaMatrixProbes), 4, 5 ) == "RL"))
      ColVector <- c(rep(1, length(which(substr(colnames(BetaMatrixProbes), 4, 5 ) == "FL"))), 
                     rep(2, length(which(substr(colnames(BetaMatrixProbes), 4, 5 ) == "DL"))), 
                     rep(3, length(which(substr(colnames(BetaMatrixProbes), 4, 5 ) == "RL"))))
      grDevices::png(paste0(pathNow, "/img/6_SummaryStatistics_BoxPlot.", PNGorPDF))
      boxplot(BetaMatrixProbes[, TypeLymphomVector], las = 2, col = ColVector + 1, 
              main = "Boxplot showing spread of SWAN normalized Beta values")
      grDevices::dev.off()
      
      # QQ plots is used to check whether the data is normally distributed
      # if (!require(car)) install.packages("car")
      # png(paste0(pathNow,"/img/6_SummaryStatistics_QQPlot.",PNGorPDF))
      # qqPlot(BetaMatrixProbes)
      # dev.off()
    }
  }
  
  if(typeof(as.matrix(BetaMatrixGenes)) == "double") {
    overall_mean_BetaMatrixGenes <- mean(as.matrix(BetaMatrixGenes))
    patient_means_BetaMatrixGenes <- colMeans(as.matrix(BetaMatrixGenes))
    overall_median_BetaMatrixGenes <- median(as.matrix(BetaMatrixGenes))
    patient_median_BetaMatrixGenes <- apply(as.matrix(BetaMatrixGenes), 2, median)
    
    range_BetaMatrixGenes <- range(as.matrix(BetaMatrixGenes))
    quantile_BetaMatrixGenes <- quantile(as.matrix(BetaMatrixGenes))
    IQR_BetaMatrixGenes <- IQR(as.matrix(BetaMatrixGenes))
    var_BetaMatrixGenes <- var(as.matrix(BetaMatrixGenes))
    sd_BetaMatrixGenes <- sd(as.matrix(BetaMatrixGenes))
    # median absolute deviation
    mad_BetaMatrixGenes <- mad(as.matrix(BetaMatrixGenes)) 
  } else {
    overall_mean_BetaMatrixGenes <- patient_means_BetaMatrixGenes <- 
    overall_median_BetaMatrixGenes <- patient_median_BetaMatrixGenes <- 
    range_BetaMatrixGenes <- quantile_BetaMatrixGenes <- 
    IQR_BetaMatrixGenes <- var_BetaMatrixGenes <- 
    sd_BetaMatrixGenes <- mad_BetaMatrixGenes <- 
    "BetaMatrixGenes Not Provided"
  }
  
  if(typeof(as.matrix(MvalueMatrixProbes)) == "double") {
    overall_mean_MvalueMatrixProbes <- mean(as.matrix(MvalueMatrixProbes))
    patient_means_MvalueMatrixProbes <- colMeans(as.matrix(MvalueMatrixProbes))
    overall_median_MvalueMatrixProbes <- median(as.matrix(MvalueMatrixProbes))
    patient_median_MvalueMatrixProbes <- apply(as.matrix(MvalueMatrixProbes),2,median)
    
    range_MvalueMatrixProbes <- range(as.matrix(MvalueMatrixProbes))
    quantile_MvalueMatrixProbes <- quantile(as.matrix(MvalueMatrixProbes))
    IQR_MvalueMatrixProbes <- IQR(as.matrix(MvalueMatrixProbes))
    var_MvalueMatrixProbes <- var(as.matrix(MvalueMatrixProbes))
    sd_MvalueMatrixProbes <- sd(as.matrix(MvalueMatrixProbes))
    # median absolute deviation
    mad_MvalueMatrixProbes <- mad(as.matrix(MvalueMatrixProbes)) 
  }
  
  if(typeof(as.matrix(MvalueMatrixGenes)) == "double") {
    overall_mean_MvalueMatrixGenes <- mean(as.matrix(MvalueMatrixGenes))
    patient_meansMvalueMatrixGenes <- colMeans(as.matrix(MvalueMatrixGenes))
    overall_median_MvalueMatrixGenes <- median(as.matrix(MvalueMatrixGenes))
    patient_median_MvalueMatrixGenes <- apply(as.matrix(MvalueMatrixGenes),2,median)
    
    range_MvalueMatrixGenes <- range(as.matrix(MvalueMatrixGenes))
    quantile_MvalueMatrixGenes <- quantile(as.matrix(MvalueMatrixGenes))
    IQR_MvalueMatrixGenes <- IQR(as.matrix(MvalueMatrixGenes))
    var_MvalueMatrixGenes <- var(as.matrix(MvalueMatrixGenes))
    sd_MvalueMatrixGenes <- sd(as.matrix(MvalueMatrixGenes))
    # median absolute deviation
    mad_MvalueMatrixGenes <- mad(as.matrix(MvalueMatrixGenes)) 
  } else {
    overall_mean_MvalueMatrixGenes <- patient_meansMvalueMatrixGenes <- 
    overall_median_MvalueMatrixGenes <- patient_median_MvalueMatrixGenes <-
    range_MvalueMatrixGenes <- quantile_MvalueMatrixGenes <- 
    IQR_MvalueMatrixGenes <- var_MvalueMatrixGenes <- 
    sd_MvalueMatrixGenes <- mad_MvalueMatrixGenes <- 
    "MvalueMatrixGenes Not Provided"
  }
  
  if(! is.na(ClinicalCategoryToVisualize)) {
    
    if(ClinicalCategoryToVisualize == "STAGE") {
      
      mean_advanced_beta <- mean(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                             ClinicalCategoryToVisualize)] == "ADVANCED")])
      median_advanced_beta <- median(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                                 ClinicalCategoryToVisualize)] == "ADVANCED")])
      range_advanced_beta <- range(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                               ClinicalCategoryToVisualize)] == "ADVANCED")])
      
      mean_advanced_mvalue <- mean(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                                 ClinicalCategoryToVisualize)] == "ADVANCED")])
      median_advanced_mvalue <- median(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                                     ClinicalCategoryToVisualize)] == "ADVANCED")])
      range_advanced_mvalue <- range(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                                   ClinicalCategoryToVisualize)] == "ADVANCED")])
      
      mean_limited_beta <- mean(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                            ClinicalCategoryToVisualize)] == "LIMITED")])
      median_limited_beta <- median(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                                ClinicalCategoryToVisualize)] == "LIMITED")])
      range_limited_beta <- range(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                              ClinicalCategoryToVisualize)] == "LIMITED")])
      
      mean_limited_mvalue <- mean(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                                ClinicalCategoryToVisualize)] == "LIMITED")])
      medianlimited_mvalue <- median(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                                   ClinicalCategoryToVisualize)] == "LIMITED")])
      range_limited_mvalue <- range(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                                  ClinicalCategoryToVisualize)] == "LIMITED")])
      
      RESULTS <- list(BetaMatrixProbesSummary = list(OverallMean = overall_mean_BetaMatrixProbes,
                                                     PatientMeans = patient_means_BetaMatrixProbes, 
                                                     OverallMedian = overall_median_BetaMatrixProbes, 
                                                     PatientMedian = patient_median_BetaMatrixProbes, 
                                                     Range = range_BetaMatrixProbes, 
                                                     Quantile = quantile_BetaMatrixProbes, 
                                                     IQR = IQR_BetaMatrixProbes, 
                                                     Variance = var_BetaMatrixProbes, 
                                                     StandardDeviation = sd_BetaMatrixProbes,
                                                     MedianAbsoluteDeviation = mad_BetaMatrixProbes),
                      BetaMatrixGenesSummary = list(OverallMean = overall_mean_BetaMatrixGenes,
                                                    PatientMeans = patient_means_BetaMatrixGenes, 
                                                    OverallMedian = overall_median_BetaMatrixGenes, 
                                                    PatientMedian = patient_median_BetaMatrixGenes, 
                                                    Range = range_BetaMatrixGenes, 
                                                    Quantile = quantile_BetaMatrixGenes, 
                                                    IQR = IQR_BetaMatrixGenes, 
                                                    Variance = var_BetaMatrixGenes, 
                                                    StandardDeviation = sd_BetaMatrixGenes, 
                                                    MedianAbsoluteDeviation = mad_BetaMatrixGenes),
                      MvalueMatrixProbesSummary = list(OverallMean = overall_mean_MvalueMatrixProbes,
                                                       PatientMeans = patient_means_MvalueMatrixProbes, 
                                                       OverallMedian = overall_median_MvalueMatrixProbes, 
                                                       PatientMedian = patient_median_MvalueMatrixProbes, 
                                                       Range = range_MvalueMatrixProbes, 
                                                       Quantile = quantile_MvalueMatrixProbes, 
                                                       IQR = IQR_MvalueMatrixProbes, 
                                                       Variance = var_MvalueMatrixProbes, 
                                                       StandardDeviation = sd_MvalueMatrixProbes, 
                                                       MedianAbsoluteDeviation = mad_MvalueMatrixProbes),
                      MvalueMatrixGenesSummary =  list(OverallMean = overall_mean_MvalueMatrixGenes,
                                                       PatientMeans = patient_meansMvalueMatrixGenes,
                                                       OverallMedian = overall_median_MvalueMatrixGenes, 
                                                       PatientMedian = patient_median_MvalueMatrixGenes, 
                                                       Range = range_MvalueMatrixGenes, 
                                                       Quantile = quantile_MvalueMatrixGenes, 
                                                       IQR = IQR_MvalueMatrixGenes, 
                                                       Variance = var_MvalueMatrixGenes, 
                                                       StandardDeviation = sd_MvalueMatrixGenes, 
                                                       MedianAbsoluteDeviation = mad_MvalueMatrixGenes),
                      ClinicalCategorySpecificSummary = list(mean_advanced_beta, 
                                                             median_advanced_beta, 
                                                             range_advanced_beta, 
                                                             mean_advanced_mvalue, 
                                                             median_advanced_mvalue, 
                                                             range_advanced_mvalue,
                                                             mean_limited_beta, 
                                                             median_limited_beta, 
                                                             range_limited_beta, 
                                                             mean_limited_mvalue, 
                                                             medianlimited_mvalue, 
                                                             range_limited_mvalue))
      
    } else if(ClinicalCategoryToVisualize == "TYPE") {
      
      mean_fl_beta <- mean(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                       ClinicalCategoryToVisualize)] == "FL")])
      median_fl_beta <- median(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                           ClinicalCategoryToVisualize)] == "FL")])
      range_fl_beta <- range(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                         ClinicalCategoryToVisualize)] == "FL")])
      
      mean_fl_mavlue <- mean(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                           ClinicalCategoryToVisualize)] == "FL")])
      median_fl_mavlue <- median(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                               ClinicalCategoryToVisualize)] == "FL")])
      range_fl_mavlue <- range(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                             ClinicalCategoryToVisualize)] == "FL")])
      
      mean_DLBCL_beta <- mean(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                          ClinicalCategoryToVisualize)] == "DLBCL")])
      median_DLBCL_beta <- median(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                              ClinicalCategoryToVisualize)] == "DLBCL")])
      range_DLBCL_beta <- range(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                            ClinicalCategoryToVisualize)] == "DLBCL")])
      
      mean_DLBCL_mavlue <- mean(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                              ClinicalCategoryToVisualize)] == "DLBCL")])
      median_DLBCL_mavlue <- median(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                                  ClinicalCategoryToVisualize)] == "DLBCL")])
      range_DLBCL_mavlue <- range(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                                ClinicalCategoryToVisualize)] == "DLBCL")])
      
      mean_RLN_beta <- mean(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                        ClinicalCategoryToVisualize)] == "RLN")])
      median_RLN_beta <- median(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                            ClinicalCategoryToVisualize)] == "RLN")])
      range_RLN_beta <- range(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                          ClinicalCategoryToVisualize)] == "RLN")])
      
      mean_RLN_mvalue <- mean(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                            ClinicalCategoryToVisualize)] == "RLN")])
      median_RLN_mvalue <- median(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                                ClinicalCategoryToVisualize)] == "RLN")])
      range_RLN_mvalue <- range(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                                                              ClinicalCategoryToVisualize)] == "RLN")])
      
      RESULTS <- list(BetaMatrixProbesSummary = list(OverallMean = overall_mean_BetaMatrixProbes,
                                                     PatientMeans = patient_means_BetaMatrixProbes, 
                                                     OverallMedian = overall_median_BetaMatrixProbes, 
                                                     PatientMedian = patient_median_BetaMatrixProbes, 
                                                     Range = range_BetaMatrixProbes, 
                                                     Quantile = quantile_BetaMatrixProbes, 
                                                     IQR = IQR_BetaMatrixProbes, 
                                                     Variance = var_BetaMatrixProbes, 
                                                     StandardDeviation = sd_BetaMatrixProbes, 
                                                     MedianAbsoluteDeviation = mad_BetaMatrixProbes),
                      BetaMatrixGenesSummary = list(OverallMean = overall_mean_BetaMatrixGenes,
                                                    PatientMeans = patient_means_BetaMatrixGenes, 
                                                    OverallMedian = overall_median_BetaMatrixGenes, 
                                                    PatientMedian = patient_median_BetaMatrixGenes, 
                                                    Range = range_BetaMatrixGenes, 
                                                    Quantile = quantile_BetaMatrixGenes, 
                                                    IQR = IQR_BetaMatrixGenes, 
                                                    Variance = var_BetaMatrixGenes, 
                                                    StandardDeviation = sd_BetaMatrixGenes, 
                                                    MedianAbsoluteDeviation = mad_BetaMatrixGenes),
                      MvalueMatrixProbesSummary = list(OverallMean = overall_mean_MvalueMatrixProbes,
                                                       PatientMeans = patient_means_MvalueMatrixProbes, 
                                                       OverallMedian = overall_median_MvalueMatrixProbes, 
                                                       PatientMedian = patient_median_MvalueMatrixProbes, 
                                                       Range = range_MvalueMatrixProbes, 
                                                       Quantile = quantile_MvalueMatrixProbes, 
                                                       IQR = IQR_MvalueMatrixProbes, 
                                                       Variance = var_MvalueMatrixProbes, 
                                                       StandardDeviation = sd_MvalueMatrixProbes, 
                                                       MedianAbsoluteDeviation = mad_MvalueMatrixProbes),
                      MvalueMatrixGenesSummary =  list(OverallMean = overall_mean_MvalueMatrixGenes,
                                                       PatientMeans = patient_meansMvalueMatrixGenes, 
                                                       OverallMedian = overall_median_MvalueMatrixGenes, 
                                                       PatientMedian = patient_median_MvalueMatrixGenes, 
                                                       Range = range_MvalueMatrixGenes, 
                                                       Quantile = quantile_MvalueMatrixGenes, 
                                                       IQR = IQR_MvalueMatrixGenes, 
                                                       Variance = var_MvalueMatrixGenes, 
                                                       StandardDeviation = sd_MvalueMatrixGenes, 
                                                       MedianAbsoluteDeviation = mad_MvalueMatrixGenes),
                      ClinicalCategorySpecificSummary = list(mean_fl_beta, 
                                                             median_fl_beta, 
                                                             range_fl_beta, 
                                                             mean_fl_mavlue, 
                                                             median_fl_mavlue, 
                                                             range_fl_mavlue, 
                                                             mean_DLBCL_beta, 
                                                             median_DLBCL_beta, 
                                                             range_DLBCL_beta, 
                                                             mean_DLBCL_mavlue, 
                                                             median_DLBCL_mavlue, 
                                                             range_DLBCL_mavlue, 
                                                             mean_RLN_beta, 
                                                             median_RLN_beta, 
                                                             range_RLN_beta, 
                                                             mean_RLN_mvalue, 
                                                             median_RLN_mvalue, 
                                                             range_RLN_mvalue))
    }
  } else {
    RESULTS <- list(BetaMatrixProbesSummary = list(OverallMean = overall_mean_BetaMatrixProbes,
                                                   PatientMeans = patient_means_BetaMatrixProbes, 
                                                   OverallMedian = overall_median_BetaMatrixProbes, 
                                                   PatientMedian = patient_median_BetaMatrixProbes, 
                                                   Range = range_BetaMatrixProbes, 
                                                   Quantile = quantile_BetaMatrixProbes, 
                                                   IQR = IQR_BetaMatrixProbes, 
                                                   Variance = var_BetaMatrixProbes, 
                                                   StandardDeviation = sd_BetaMatrixProbes, 
                                                   MedianAbsoluteDeviation = mad_BetaMatrixProbes),
                    BetaMatrixGenesSummary = list(OverallMean = overall_mean_BetaMatrixGenes,
                                                  PatientMeans = patient_means_BetaMatrixGenes, 
                                                  OverallMedian = overall_median_BetaMatrixGenes, 
                                                  PatientMedian = patient_median_BetaMatrixGenes, 
                                                  Range = range_BetaMatrixGenes, 
                                                  Quantile = quantile_BetaMatrixGenes, 
                                                  IQR = IQR_BetaMatrixGenes, 
                                                  Variance = var_BetaMatrixGenes, 
                                                  StandardDeviation = sd_BetaMatrixGenes, 
                                                  MedianAbsoluteDeviation = mad_BetaMatrixGenes),
                    MvalueMatrixProbesSummary = list(OverallMean = overall_mean_MvalueMatrixProbes,
                                                     PatientMeans = patient_means_MvalueMatrixProbes, 
                                                     OverallMedian = overall_median_MvalueMatrixProbes, 
                                                     PatientMedian = patient_median_MvalueMatrixProbes, 
                                                     Range =range_MvalueMatrixProbes, 
                                                     Quantile =quantile_MvalueMatrixProbes, 
                                                     IQR = IQR_MvalueMatrixProbes, 
                                                     Variance = var_MvalueMatrixProbes, 
                                                     StandardDeviation = sd_MvalueMatrixProbes,
                                                     MedianAbsoluteDeviation = mad_MvalueMatrixProbes),
                    MvalueMatrixGenesSummary =  list(OverallMean = overall_mean_MvalueMatrixGenes,
                                                     PatientMeans = patient_meansMvalueMatrixGenes, 
                                                     OverallMedian = overall_median_MvalueMatrixGenes, 
                                                     PatientMedian = patient_median_MvalueMatrixGenes, 
                                                     Range = range_MvalueMatrixGenes, 
                                                     Quantile = quantile_MvalueMatrixGenes, 
                                                     IQR = IQR_MvalueMatrixGenes, 
                                                     Variance = var_MvalueMatrixGenes, 
                                                     StandardDeviation = sd_MvalueMatrixGenes, 
                                                     MedianAbsoluteDeviation = mad_MvalueMatrixGenes))
  }
  
  class(RESULTS) <- "SummaryStatistics_ASilva"
  return(RESULTS) 
}



