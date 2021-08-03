# Updated 3 Aug 2021
# Updated 9 April 2019
# Function: Summarize descriptive statistics of all datasets (Beta and Mvalues).
# Author: Anjali Silva

# Input:
# BetaMatrixProbes: A matrix of beta values for probes x patients, with probes in rows 
#                   and patients as columns.
# BetaMatrixGenes: A matrix of beta values for genes x patients, with genes in rows and patients as columns.
# ClinicalCategoryToVisualize: Specify column of ClinicalFile to use. Options "STAGE", "SEX", "SITE_BIOPSY", 
#                              "TYPE_BIOPSY", "INSTITUTION", "COO", "TYPE", "TRANSLOC_14_18".
# MvalueMatrixProbes: A matrix of M-values for probes x patients, with probes in rows and patients as columns.
# MvalueMatrixGenes A matrix of M-values for genes x patients, with genes in rows and patients as columns.

# Output: 
# BetaMatrixProbesSummary: List with summary statistics for Beta Matrix of probes including overall mean, 
#                          patient means, overall median, patient medians, range, quantile, IQR, variance,
#                          standard deviation, median absolute deviation.
# BetaMatrixGenesSummary: List with summary statistics for Beta Matrix of genes including overall mean, 
#                         patient means, overall median, patient medians, range, quantile, IQR, variance,
#                         standard deviation, median absolute deviation.
# MvalueMatrixProbesSummary: List with summary statistics for Mvalue Matrix of probes including overall
#                            mean, patient means, overall median, patient medians, range, quantile, IQR,
#                            variance, standard deviation, median absolute deviation
# MvalueMatrixGenesSummary: List with summary statistics for Mvalue Matrix of probes including overall mean, 
#                           patient means, overall median, patient medians, range, quantile, IQR, variance,
#                           standard deviation, median absolute deviation.

# Visuals saved to img folder
# 6_SummaryStatistics_BoxPlot.p* 

SummaryStatistics6 <- function(BetaMatrixProbes = NA, 
                              BetaMatrixGenes = NA, 
                              ClinicalCategoryToVisualize = NA, 
                              ClinicalFile = NA, 
                              MvalueMatrixProbes = NA, 
                              MvalueMatrixGenes = NA, 
                              ProduceImages = "Yes", 
                              PNGorPDF = "png") {
  
  # Loading needed packages
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
    overallMeanBetaMatrixProbes <- mean(as.matrix(BetaMatrixProbes))
    patientMeansBetaMatrixProbes <- colMeans(as.matrix(BetaMatrixProbes))
    overallMedianBetaMatrixProbes <- median(as.matrix(BetaMatrixProbes))
    patientMedianBetaMatrixProbes <- apply(as.matrix(BetaMatrixProbes), 2, median)
    # Mode commented out as it takes a long time
    # overall_mode_BetaMatrixProbes <- Mode(round(BetaMatrixProbes,2))
    # patient_mode_BetaMatrixProbes <- apply(BetaMatrixProbes,2,Mode)
    rangeBetaMatrixProbes <- range(as.matrix(BetaMatrixProbes))
    quantileBetaMatrixProbes <- quantile(as.matrix(BetaMatrixProbes))
    IQRBetaMatrixProbes <- IQR(as.matrix(BetaMatrixProbes))
    varBetaMatrixProbes <- var(as.matrix(BetaMatrixProbes))
    sdBetaMatrixProbes <- sd(as.matrix(BetaMatrixProbes))
    # median absolute deviation
    madBetaMatrixProbes <- mad(as.matrix(BetaMatrixProbes)) 
    
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
    overallMeanBetaMatrixGenes <- mean(as.matrix(BetaMatrixGenes))
    patientMeansBetaMatrixGenes <- colMeans(as.matrix(BetaMatrixGenes))
    overallMedianBetaMatrixGenes <- median(as.matrix(BetaMatrixGenes))
    patientMedianBetaMatrixGenes <- apply(as.matrix(BetaMatrixGenes), 2, median)
    
    rangeBetaMatrixGenes <- range(as.matrix(BetaMatrixGenes))
    quantileBetaMatrixGenes <- quantile(as.matrix(BetaMatrixGenes))
    IQRBetaMatrixGenes <- IQR(as.matrix(BetaMatrixGenes))
    varBetaMatrixGenes <- var(as.matrix(BetaMatrixGenes))
    sdBetaMatrixGenes <- sd(as.matrix(BetaMatrixGenes))
    # median absolute deviation
    madBetaMatrixGenes <- mad(as.matrix(BetaMatrixGenes)) 
  } else {
    overallMeanBetaMatrixGenes <- patientMeansBetaMatrixGenes <- 
    overallMedianBetaMatrixGenes <- patientMedianBetaMatrixGenes <- 
    rangeBetaMatrixGenes <- quantileBetaMatrixGenes <- 
    IQRBetaMatrixGenes <- varBetaMatrixGenes <- 
    sdBetaMatrixGenes <- madBetaMatrixGenes <- 
    "BetaMatrixGenes Not Provided"
  }
  
  if(typeof(as.matrix(MvalueMatrixProbes)) == "double") {
    overallMeanMvalueMatrixProbes <- mean(as.matrix(MvalueMatrixProbes))
    patientMeansMvalueMatrixProbes <- colMeans(as.matrix(MvalueMatrixProbes))
    overallMedianMvalueMatrixProbes <- median(as.matrix(MvalueMatrixProbes))
    patientMedianMvalueMatrixProbes <- apply(as.matrix(MvalueMatrixProbes), 2, median)
    
    rangeMvalueMatrixProbes <- range(as.matrix(MvalueMatrixProbes))
    quantileMvalueMatrixProbes <- quantile(as.matrix(MvalueMatrixProbes))
    IQRMvalueMatrixProbes <- IQR(as.matrix(MvalueMatrixProbes))
    varMvalueMatrixProbes <- var(as.matrix(MvalueMatrixProbes))
    sdMvalueMatrixProbes <- sd(as.matrix(MvalueMatrixProbes))
    # median absolute deviation
    madMvalueMatrixProbes <- mad(as.matrix(MvalueMatrixProbes)) 
  }
  
  if(typeof(as.matrix(MvalueMatrixGenes)) == "double") {
    overallMeanMvalueMatrixGenes <- mean(as.matrix(MvalueMatrixGenes))
    patientMeansMvalueMatrixGenes <- colMeans(as.matrix(MvalueMatrixGenes))
    overallMedianMvalueMatrixGenes <- median(as.matrix(MvalueMatrixGenes))
    patientMedianMvalueMatrixGenes <- apply(as.matrix(MvalueMatrixGenes), 2, median)
    
    rangeMvalueMatrixGenes <- range(as.matrix(MvalueMatrixGenes))
    quantileMvalueMatrixGenes <- quantile(as.matrix(MvalueMatrixGenes))
    IQRMvalueMatrixGenes <- IQR(as.matrix(MvalueMatrixGenes))
    varMvalueMatrixGenes <- var(as.matrix(MvalueMatrixGenes))
    sdMvalueMatrixGenes <- sd(as.matrix(MvalueMatrixGenes))
    # median absolute deviation
    madMvalueMatrixGenes <- mad(as.matrix(MvalueMatrixGenes)) 
  } else {
    overallMeanMvalueMatrixGenes <- patientMeansMvalueMatrixGenes <- 
    overallMedianMvalueMatrixGenes <- patientMedianMvalueMatrixGenes <-
    rangeMvalueMatrixGenes <- quantileMvalueMatrixGenes <- 
    IQRMvalueMatrixGenes <- varMvalueMatrixGenes <- 
    sdMvalueMatrixGenes <- madMvalueMatrixGenes <- 
    "MvalueMatrixGenes Not Provided"
  }
  
  if(! is.na(ClinicalCategoryToVisualize)) {
    
    if(ClinicalCategoryToVisualize == "STAGE") {
      
      meanAdvancedBeta <- mean(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                            ClinicalCategoryToVisualize)] == "ADVANCED")])
      medianAdvancedBeta <- median(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                            ClinicalCategoryToVisualize)] == "ADVANCED")])
      rangeAdvancedBeta <- range(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                          ClinicalCategoryToVisualize)] == "ADVANCED")])
      
      meanAdvancedMvalue <- mean(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                              ClinicalCategoryToVisualize)] == "ADVANCED")])
      medianAdvancedMvalue <- median(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                                ClinicalCategoryToVisualize)] == "ADVANCED")])
      rangeAdvancedMvalue <- range(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                               ClinicalCategoryToVisualize)] == "ADVANCED")])
      
      meanLimitedBeta <- mean(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                           ClinicalCategoryToVisualize)] == "LIMITED")])
      medianLimitedBeta <- median(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                             ClinicalCategoryToVisualize)] == "LIMITED")])
      rangeLimitedBeta <- range(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                            ClinicalCategoryToVisualize)] == "LIMITED")])
      
      meanLimitedMvalue <- mean(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                             ClinicalCategoryToVisualize)] == "LIMITED")])
      medianlimitedMvalue <- median(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                              ClinicalCategoryToVisualize)] == "LIMITED")])
      rangeLimitedMvalue <- range(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                              ClinicalCategoryToVisualize)] == "LIMITED")])
      
      RESULTS <- list(BetaMatrixProbesSummary = list(OverallMean = overallMeanBetaMatrixProbes,
                                                     PatientMeans = patientMeansBetaMatrixProbes, 
                                                     OverallMedian = overallMedianBetaMatrixProbes, 
                                                     PatientMedian = patientMedianBetaMatrixProbes, 
                                                     Range = rangeBetaMatrixProbes, 
                                                     Quantile = quantileBetaMatrixProbes, 
                                                     IQR = IQRBetaMatrixProbes, 
                                                     Variance = varBetaMatrixProbes, 
                                                     StandardDeviation = sdBetaMatrixProbes,
                                                     MedianAbsoluteDeviation = madBetaMatrixProbes),
                      BetaMatrixGenesSummary = list(OverallMean = overallMeanBetaMatrixGenes,
                                                    PatientMeans = patientMeansBetaMatrixGenes, 
                                                    OverallMedian = overallMedianBetaMatrixGenes, 
                                                    PatientMedian = patientMedianBetaMatrixGenes, 
                                                    Range = rangeBetaMatrixGenes, 
                                                    Quantile = quantileBetaMatrixGenes, 
                                                    IQR = IQRBetaMatrixGenes, 
                                                    Variance = varBetaMatrixGenes, 
                                                    StandardDeviation = sdBetaMatrixGenes, 
                                                    MedianAbsoluteDeviation = madBetaMatrixGenes),
                      MvalueMatrixProbesSummary = list(OverallMean = overallMeanMvalueMatrixProbes,
                                                       PatientMeans = patientMeansMvalueMatrixProbes, 
                                                       OverallMedian = overallMedianMvalueMatrixProbes, 
                                                       PatientMedian = patientMedianMvalueMatrixProbes, 
                                                       Range = rangeMvalueMatrixProbes, 
                                                       Quantile = quantileMvalueMatrixProbes, 
                                                       IQR = IQRMvalueMatrixProbes, 
                                                       Variance = varMvalueMatrixProbes, 
                                                       StandardDeviation = sdMvalueMatrixProbes, 
                                                       MedianAbsoluteDeviation = madMvalueMatrixProbes),
                      MvalueMatrixGenesSummary =  list(OverallMean = overallMeanMvalueMatrixGenes,
                                                       PatientMeans = patientMeansMvalueMatrixGenes,
                                                       OverallMedian = overallMedianMvalueMatrixGenes, 
                                                       PatientMedian = patientMedianMvalueMatrixGenes, 
                                                       Range = rangeMvalueMatrixGenes, 
                                                       Quantile = quantileMvalueMatrixGenes, 
                                                       IQR = IQRMvalueMatrixGenes, 
                                                       Variance = varMvalueMatrixGenes, 
                                                       StandardDeviation = sdMvalueMatrixGenes, 
                                                       MedianAbsoluteDeviation = madMvalueMatrixGenes),
                      ClinicalCategorySpecificSummary = list(meanAdvancedBeta, 
                                                             medianAdvancedBeta, 
                                                             rangeAdvancedBeta, 
                                                             meanAdvancedMvalue, 
                                                             medianAdvancedMvalue, 
                                                             rangeAdvancedMvalue,
                                                             meanLimitedBeta, 
                                                             medianLimitedBeta, 
                                                             rangeLimitedBeta, 
                                                             meanLimitedMvalue, 
                                                             medianlimitedMvalue, 
                                                             rangeLimitedMvalue))
      
    } else if(ClinicalCategoryToVisualize == "TYPE") {
      
      meanFLbeta <- mean(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                    ClinicalCategoryToVisualize)] == "FL")])
      medianFLbeta <- median(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                        ClinicalCategoryToVisualize)] == "FL")])
      rangeFLbeta <- range(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                       ClinicalCategoryToVisualize)] == "FL")])
      
      meanFLmavlue <- mean(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                        ClinicalCategoryToVisualize)] == "FL")])
      medianFLmavlue <- median(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                          ClinicalCategoryToVisualize)] == "FL")])
      rangeFLmavlue <- range(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                         ClinicalCategoryToVisualize)] == "FL")])
      
      meanDLBCLbeta <- mean(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                         ClinicalCategoryToVisualize)] == "DLBCL")])
      medianDLBCLbeta <- median(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                           ClinicalCategoryToVisualize)] == "DLBCL")])
      rangeDLBCLbeta <- range(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                          ClinicalCategoryToVisualize)] == "DLBCL")])
      
      meanDLBCLmvalue <- mean(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                           ClinicalCategoryToVisualize)] == "DLBCL")])
      medianDLBCLmvalue <- median(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                             ClinicalCategoryToVisualize)] == "DLBCL")])
      rangeDLBCLmvalue <- range(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                            ClinicalCategoryToVisualize)] == "DLBCL")])
      
      meanRLNbeta <- mean(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                       ClinicalCategoryToVisualize)] == "RLN")])
      medianRLNbeta <- median(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                         ClinicalCategoryToVisualize)] == "RLN")])
      rangeRLNbeta <- range(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                        ClinicalCategoryToVisualize)] == "RLN")])
      
      meanRLNmvalue <- mean(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                         ClinicalCategoryToVisualize)] == "RLN")])
      medianRLNmvalue <- median(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                           ClinicalCategoryToVisualize)] == "RLN")])
      rangeRLNmvalue <- range(MvalueMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == 
                          ClinicalCategoryToVisualize)] == "RLN")])
      
      RESULTS <- list(BetaMatrixProbesSummary = list(OverallMean = overallMeanBetaMatrixProbes,
                                                     PatientMeans = patientMeansBetaMatrixProbes, 
                                                     OverallMedian = overallMedianBetaMatrixProbes, 
                                                     PatientMedian = patientMedianBetaMatrixProbes, 
                                                     Range = rangeBetaMatrixProbes, 
                                                     Quantile = quantileBetaMatrixProbes, 
                                                     IQR = IQRBetaMatrixProbes, 
                                                     Variance = varBetaMatrixProbes, 
                                                     StandardDeviation = sdBetaMatrixProbes, 
                                                     MedianAbsoluteDeviation = madBetaMatrixProbes),
                      BetaMatrixGenesSummary = list(OverallMean = overallMeanBetaMatrixGenes,
                                                    PatientMeans = patientMeansBetaMatrixGenes, 
                                                    OverallMedian = overallMedianBetaMatrixGenes, 
                                                    PatientMedian = patientMedianBetaMatrixGenes, 
                                                    Range = rangeBetaMatrixGenes, 
                                                    Quantile = quantileBetaMatrixGenes, 
                                                    IQR = IQRBetaMatrixGenes, 
                                                    Variance = varBetaMatrixGenes, 
                                                    StandardDeviation = sdBetaMatrixGenes, 
                                                    MedianAbsoluteDeviation = madBetaMatrixGenes),
                      MvalueMatrixProbesSummary = list(OverallMean = overallMeanMvalueMatrixProbes,
                                                       PatientMeans = patientMeansMvalueMatrixProbes, 
                                                       OverallMedian = overallMedianMvalueMatrixProbes, 
                                                       PatientMedian = patientMedianMvalueMatrixProbes, 
                                                       Range = rangeMvalueMatrixProbes, 
                                                       Quantile = quantileMvalueMatrixProbes, 
                                                       IQR = IQRMvalueMatrixProbes, 
                                                       Variance = varMvalueMatrixProbes, 
                                                       StandardDeviation = sdMvalueMatrixProbes, 
                                                       MedianAbsoluteDeviation = madMvalueMatrixProbes),
                      MvalueMatrixGenesSummary =  list(OverallMean = overallMeanMvalueMatrixGenes,
                                                       PatientMeans = patientMeansMvalueMatrixGenes, 
                                                       OverallMedian = overallMedianMvalueMatrixGenes, 
                                                       PatientMedian = patientMedianMvalueMatrixGenes, 
                                                       Range = rangeMvalueMatrixGenes, 
                                                       Quantile = quantileMvalueMatrixGenes, 
                                                       IQR = IQRMvalueMatrixGenes, 
                                                       Variance = varMvalueMatrixGenes, 
                                                       StandardDeviation = sdMvalueMatrixGenes, 
                                                       MedianAbsoluteDeviation = madMvalueMatrixGenes),
                      ClinicalCategorySpecificSummary = list(meanFLbeta, 
                                                             medianFLbeta, 
                                                             rangeFLbeta, 
                                                             meanFLmavlue, 
                                                             medianFLmavlue, 
                                                             rangeFLmavlue, 
                                                             meanDLBCLbeta, 
                                                             medianDLBCLbeta, 
                                                             rangeDLBCLbeta, 
                                                             meanDLBCLmvalue, 
                                                             medianDLBCLmvalue, 
                                                             rangeDLBCLmvalue, 
                                                             meanRLNbeta, 
                                                             medianRLNbeta, 
                                                             rangeRLNbeta, 
                                                             meanRLNmvalue, 
                                                             medianRLNmvalue, 
                                                             rangeRLNmvalue))
    }
  } else {
    RESULTS <- list(BetaMatrixProbesSummary = list(OverallMean = overallMeanBetaMatrixProbes,
                                                   PatientMeans = patientMeansBetaMatrixProbes, 
                                                   OverallMedian = overallMedianBetaMatrixProbes, 
                                                   PatientMedian = patientMedianBetaMatrixProbes, 
                                                   Range = rangeBetaMatrixProbes, 
                                                   Quantile = quantileBetaMatrixProbes, 
                                                   IQR = IQRBetaMatrixProbes, 
                                                   Variance = varBetaMatrixProbes, 
                                                   StandardDeviation = sdBetaMatrixProbes, 
                                                   MedianAbsoluteDeviation = madBetaMatrixProbes),
                    BetaMatrixGenesSummary = list(OverallMean = overallMeanBetaMatrixGenes,
                                                  PatientMeans = patientMeansBetaMatrixGenes, 
                                                  OverallMedian = overallMedianBetaMatrixGenes, 
                                                  PatientMedian = patientMedianBetaMatrixGenes, 
                                                  Range = rangeBetaMatrixGenes, 
                                                  Quantile = quantileBetaMatrixGenes, 
                                                  IQR = IQRBetaMatrixGenes, 
                                                  Variance = varBetaMatrixGenes, 
                                                  StandardDeviation = sdBetaMatrixGenes, 
                                                  MedianAbsoluteDeviation = madBetaMatrixGenes),
                    MvalueMatrixProbesSummary = list(OverallMean = overallMeanMvalueMatrixProbes,
                                                     PatientMeans = patientMeansMvalueMatrixProbes, 
                                                     OverallMedian = overallMedianMvalueMatrixProbes, 
                                                     PatientMedian = patientMedianMvalueMatrixProbes, 
                                                     Range =rangeMvalueMatrixProbes, 
                                                     Quantile =quantileMvalueMatrixProbes, 
                                                     IQR = IQRMvalueMatrixProbes, 
                                                     Variance = varMvalueMatrixProbes, 
                                                     StandardDeviation = sdMvalueMatrixProbes,
                                                     MedianAbsoluteDeviation = madMvalueMatrixProbes),
                    MvalueMatrixGenesSummary =  list(OverallMean = overallMeanMvalueMatrixGenes,
                                                     PatientMeans = patientMeansMvalueMatrixGenes, 
                                                     OverallMedian = overallMedianMvalueMatrixGenes, 
                                                     PatientMedian = patientMedianMvalueMatrixGenes, 
                                                     Range = rangeMvalueMatrixGenes, 
                                                     Quantile = quantileMvalueMatrixGenes, 
                                                     IQR = IQRMvalueMatrixGenes, 
                                                     Variance = varMvalueMatrixGenes, 
                                                     StandardDeviation = sdMvalueMatrixGenes, 
                                                     MedianAbsoluteDeviation = madMvalueMatrixGenes))
  }
  
  class(RESULTS) <- "SummaryStatistics_ASilva"
  return(RESULTS) 
}
# [END]


