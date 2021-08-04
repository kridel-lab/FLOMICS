# Updated 3 Aug 2021
# Updated 3 May 2019
# Function: Fit a multinomial model for feature selection using penalized logistic regression
#           via lasso regression. Features (probes) are then selected for BetaMatrix and 
#           MvalueMatrix. Currently developed only for "TYPE".
# Ref: http://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/
# Author: Anjali Silva

# Input:
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients 
#             as columns.
# MvalueMatrix: A matrix of probes x patients, with probes in rows and patients as columns.
# FigureGenerate: Produce images or not, options = "Yes" or "No"; default "Yes". 
# PNGorPDF: Output format of the image, options = "png" or "pdf"; default "png".

# Output:
# BetaMatrix_Lambda.min: Beta matrix with probes corresponding to non zero regression coefficients
#                        using lambda.min as the best lambda.
# MvalueMatrix_Lambda.min: Mvalue matrix with probes corresponding to non zero regression coefficients
#                          using lambda.min as the best lambda.
# Probes: Probe names corresponding to non zero regression coefficients using lambda.min as the best lambda.

GlmnetFeatureSelection20 <- function(BetaMatrix,
                                    MvalueMatrix, 
                                    FigureGenerate = "Yes", 
                                    PNGorPDF = "png") {
  
 # Loading needed packages
 # LoadCheckPkg(RegularPckgs = c("tidyverse", "caret", "glmnet", 
 # "plotmo", "pROC", "BootValidation"))
  library(tidyverse)
  library(caret)
  library(glmnet)
  library(plotmo)
  library(pROC)
  # library(BootValidation)
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  TypeLymphomVector <- c(which(substr(colnames(BetaMatrix), 4, 5 ) == "FL"), 
                         which(substr(colnames(BetaMatrix), 4, 5 ) == "DL"), 
                         which(substr(colnames(BetaMatrix), 4, 5 ) == "RL"))
  ColVector <- c(rep(1, length(which(substr(colnames(BetaMatrix), 4, 5 ) == "FL"))), 
                 rep(2, length(which(substr(colnames(BetaMatrix), 4, 5 ) == "DL"))), 
                 rep(3, length(which(substr(colnames(BetaMatrix), 4, 5 ) == "RL"))))
  
  # Find the optimal value of lambda that minimizes the cross-validation error
  set.seed(123) 
  cv.lasso <- glmnet::cv.glmnet(x = t(x = BetaMatrix[, TypeLymphomVector]), 
                                y = ColVector, 
                                alpha = 1, 
                                family = "multinomial")
  # cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "multinomial", type.multinomial = "grouped", 
  #                       nfolds = 10, weights = weights)
  
  if (PNGorPDF == "png") {
    grDevices::png(paste0(pathNow, "/img/20_GlmnetFeatureSelection.", PNGorPDF))
  }
  if (PNGorPDF == "pdf") { 
    grDevices::pdf(paste0(pathNow, "/img/20_GlmnetFeatureSelection.", PNGorPDF))
  }
  plot(cv.lasso)
  grDevices::dev.off()
  
  # cv.lasso$lambda.min : this value will give the most accurate model
  # coef(cv.lasso, cv.lasso$lambda.min) was used
  # cv.lasso$lambda.1se: value of lambda that gives the simplest model but also lies within one standard error of the optimal value of lambda
  #   coef(cv.lasso, cv.lasso$lambda.1se)[[1]] was not used

  # Determine the regression coefficients using lambda.min as the best lambda
  coefs_lambda.min <- list(FL = rownames(coef(cv.lasso, cv.lasso$lambda.min)[[1]])[which(
                                coef(cv.lasso, cv.lasso$lambda.min)[[1]] != 0)[- 1]],
                           DL = rownames(coef(cv.lasso, cv.lasso$lambda.min)[[2]])[
                                which(coef(cv.lasso, cv.lasso$lambda.min)[[2]] != 0)[- 1]],
                           RL = rownames(coef(cv.lasso, cv.lasso$lambda.min)[[3]])[
                                which(coef(cv.lasso, cv.lasso$lambda.min)[[3]] != 0)[- 1]]) 
  
  # Getting the corresponding probe names from BetaMatrix
  BetaMatrix_coefs_lambda.min <- BetaMatrix[sapply(1:length(unique(unlist(coefs_lambda.min))), 
                                            function(x) 
                                            which(rownames(BetaMatrix) == unique(unlist(coefs_lambda.min))[x])), ]
  MvalueMatrix_coefs_lambda.min <- MvalueMatrix[sapply(1:length(unique(unlist(coefs_lambda.min))), 
                                                function(x) 
                                                which(rownames(BetaMatrix) == unique(unlist(coefs_lambda.min))[x])), ] 
  
  # Standardize coefficients
  # See https://stats.stackexchange.com/questions/14853/variable-importance-from-glmnet/211396#211396
  # std_coefs <- coefs.df
  # row.names(std_coefs) <- std_coefs$genes
  # std_coefs <- std_coefs[-1,2:4]
  # sds <- apply(x, 2, sd)
  # std_coefs <- std_coefs * sds
  # std_coefs$CNS.coefs <- as.numeric(std_coefs$CNS.coefs)
  # std_coefs$SYST.coefs <- as.numeric(std_coefs$SYST.coefs)
  # std_coefs$NO.coefs <- as.numeric(std_coefs$NO.coefs)

  RESULTS <- list(BetaMatrix_Lambda.min = BetaMatrix_coefs_lambda.min,
                  MvalueMatrix_Lambda.min = MvalueMatrix_coefs_lambda.min,
                  Probes = coefs_lambda.min)
  
  class(RESULTS) <- "GlmnetFeatureSelection_ASilva"
  return(RESULTS)
}
# [END]
