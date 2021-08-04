# Updated 3 Aug 2021
# 23 October 2019
# Function: Detecting cell-type specific differential methylation. This function is under construction.
# Ref: # https://bioconductor.org/packages/release/bioc/vignettes/TOAST/inst/doc/TOAST.html
# Author: Anjali Silva


# Input:
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns.
#            The ncol(BetaMatrix) should match ncol(RNAseqCountMatrix), ncol(MvalueMatrix), nrow(ClinicalFile)
#            and nrow(QCMatrix).
# ClinicalFile: A data frame of clinical categories, that has size patients x clinical categories. The
#               nrow(ClinicalFile) should match ncol(RNAseqCountMatrix), ncol(MvalueMatrix), ncol(BetaMatrix),
#               and nrow(QCMatrix).
# NumCellTypes: An integer specifying the nuymber of cell types. 


# Output: 
# This function is under construction.

ApplyingTOAST29 <- function(BetaMatrix,
                          ClinicalFile,
                          NumCellTypes) {
  
  library(TOAST)
  # Reference-free deconvolution using RefFreeEWAS
  library(RefFreeEWAS)
  
  # Select the top 1000 most variant features
  refinx <- TOAST::findRefinx(BetaMatrix, nmarker = 1000, sortBy = "var")
  Y <- BetaMatrix[refinx,]
  # Do reference-free deconvolution on the RA dataset.
  # Y = Mu Omega^T
  outT <- RefFreeEWAS::RefFreeCellMix(Y, mu0 = RefFreeEWAS::RefFreeCellMixInitialize(Y, K = NumCellTypes))
  # names(outT) # "Mu" "Omega" "incrementalChangeSummary"
  estProp_RF <- outT$Omega
  # dim(estProp_RF) # 170   6
  
  CellTypeMap6 <- mclust::map(estProp_RF)
  # table(CellTypeMap6)
  #  1  2  3  4  5  6 
  # 19 15 38 13 83  2
  # Likely, better to have 5 cell types?
  
  
  # Do reference-free deconvolution using improved-RF implemented with RefFreeCellMix. 
  NumCellTypes = 6
  set.seed(1234)
  # Improve reference-free deconvolution using cross-cell type differential analysis
  outRF1_6 <- TOAST::csDeconv(BetaMatrix, NumCellTypes, TotalIter = 30, bound_negative = TRUE) 
  # Comparing udpated RF estimations versus RB results.
  estProp_RF_improved_6 <- TOAST::assignCellType(input = outRF1_6$estProp,
                                               reference = estProp_RF) 
  mean(diag(cor(estProp_RF_improved_6, estProp_RF))) # 0.5310188
  CellTypeMap_improved_6 <- mclust::map(estProp_RF_improved_6)
  table(CellTypeMap_improved_6)
  
  
  ##### ##### ##### ##### ##### ###
  # Try with 5 cell types
  NumCellTypes <- 5
  # Y = Mu Omega^T
  outT5 <- RefFreeEWAS::RefFreeCellMix(Y, mu0 = RefFreeEWAS::RefFreeCellMixInitialize(Y, K = NumCellTypes))
  names(outT5) # "Mu" "Omega" "incrementalChangeSummary"
  estProp_RF5 <- outT5$Omega
  dim(estProp_RF5) # 170   5
  
  CellTypeMap5 <- mclust::map(estProp_RF5)
  table(CellTypeMap5)
  #  1  2  3  4  5 
  #  20 14 39 13 84 
  
  
  # Do reference-free deconvolution using improved-RF implemented with RefFreeCellMix. 
  NumCellTypes = 5
  set.seed(1234)
  # Improve reference-free deconvolution using cross-cell type differential analysis
  outRF1_5 <- TOAST::csDeconv(BetaMatrix, NumCellTypes, TotalIter = 60, bound_negative = TRUE) 
  # Comparing udpated RF estimations versus RB results.
  estProp_RF_improved_5 <- TOAST::assignCellType(input = outRF1_5$estProp,
                                                 reference = estProp_RF5) 
  mean(diag(cor(estProp_RF_improved_5, estProp_RF5))) # 0.6901539
  CellTypeMap_improved_5 <- mclust::map(estProp_RF_improved_5)
  table(CellTypeMap_improved_5)
  # CellTypeMap_improved_5
  # 1  2  3  4  5 
  # 35 55  9 13 58
  
  colnames(estProp_RF_improved_5) <- c("CellType1", "CellType2", 
                                       "CellType3", "CellType4", "CellType5")
  
  boxplot(estProp_RF_improved_5)
  par(mfrow = c(2,3))
  hist(estProp_RF_improved_5[,1], main ="Cell Type 1")
  hist(estProp_RF_improved_5[,2], main ="Cell Type 2")
  hist(estProp_RF_improved_5[,3], main ="Cell Type 3")
  hist(estProp_RF_improved_5[,4], main ="Cell Type 4")
  hist(estProp_RF_improved_5[,5], main ="Cell Type 5")
  range(estProp_RF_improved_5) # -4.581426e-17  8.576891e-01
  

  # Detect cell type-specific and cross-cell type differential signals
  TYPE_design <- data.frame(disease = factor(ClinicalFile$TYPE, 
                                             levels = c("DLBCL", "FL", "RLN"), 
                                             labels = c("0", "1", "2")))
  
  # Generate design matrix from input phenotypes and proportions.
  TYPE_Design_out <- TOAST::makeDesign(design = TYPE_design, Prop = estProp_RF_improved_5)
  # Fit model with proportions and phenotypes
  TYPE_fitted_model <- TOAST::fitModel(Design_out = TYPE_Design_out, Y = BetaMatrix)
  # Error in solve.default(t(W) %*% W) : 
  #   system is computationally singular: reciprocal condition number = 1.02816e-36
  # Does not work
  
  STAGE_original <- factor(ClinicalFile$STAGE, 
                           levels = c("ADVANCED", "LIMITED"),
                           labels = c("1", "0"))
  STAGE <- STAGE_original[! is.na(STAGE_original)]
  STAGE_design <- data.frame(stage = STAGE)
  ## columns of proportion matrix should have names
  
  # Make model design using the design (phenotype) data frame and proportion matrix.
  STAGE_estProp_RF_improved_5_design_out <- TOAST::makeDesign(STAGE_design, 
                                                              estProp_RF_improved_5[! is.na(STAGE_original),])
  names(STAGE_estProp_RF_improved_5_design_out)
  # "design_matrix"  "Prop"           "design"         "all_coefs"      "all_cell_types" "formula"   
  STAGE_estProp_RF_improved_5_design_out$all_cell_types # "CellType1" "CellType2" "CellType3" "CellType4" "CellType5"
  STAGE_estProp_RF_improved_5_design_out$Prop
  
  # Fit linear models for raw data and the design generated from Design_out().
  STAGE_estProp_RF_improved_5_fitted_model <- TOAST::fitModel(STAGE_estProp_RF_improved_5_design_out, 
                                                              BetaMatrix[, ! is.na(STAGE_original)])
  # print all the cell type names
  names(STAGE_estProp_RF_improved_5_fitted_model)
  # "Design_out"     "N"              "coefs"          "coefs_var"      "Y"              "Ypred"         
  # "resi"           "all_coefs"      "all_cell_types" "MSE"            "model_names"   
  STAGE_estProp_RF_improved_5_fitted_model$all_cell_types # The names of all cell types
  STAGE_estProp_RF_improved_5_fitted_model$coefs # Estimated coefficients (beta) in the model
  STAGE_estProp_RF_improved_5_fitted_model$Ypred # Predicted Y from the fitted model
  STAGE_estProp_RF_improved_5_fitted_model$all_coefs # The names of all phenotypes
  STAGE_estProp_RF_improved_5_fitted_model$MSE # Estimated mean squared error
  STAGE_estProp_RF_improved_5_fitted_model$model_names # 	The names of all terms in the fitted model.
  
  # Detecting cell type-specific differential signals
  # For example, testing stage effect in CellType1
  # Testing differential signals for specified phenotype and cell type(s).
  STAGE_estProp_RF_improved_5_res_table_CellType1 <- TOAST::csTest(STAGE_estProp_RF_improved_5_fitted_model, 
                                                                   coef = "stage", 
                                                                   cell_type = "CellType1")
  head(STAGE_estProp_RF_improved_5_res_table_CellType1, 3)
  range(STAGE_estProp_RF_improved_5_res_table_CellType1$fdr) # 0.05850281 0.99999989
  length(STAGE_estProp_RF_improved_5_res_table_CellType1$p_value < 0.05) # 559,110
  
  
  # CellType2
  STAGE_estProp_RF_improved_5_res_table_CellType2 <- TOAST::csTest(STAGE_estProp_RF_improved_5_fitted_model, 
                                                                   coef = "stage", 
                                                                   cell_type = "CellType2")
  head(STAGE_estProp_RF_improved_5_res_table_CellType2, 3)
  range(STAGE_estProp_RF_improved_5_res_table_CellType2$fdr) # 0.08787088 0.99999851
  length(STAGE_estProp_RF_improved_5_res_table_CellType2$p_value < 0.05) # 559110
  
  # CellType3
  STAGE_estProp_RF_improved_5_res_table_CellType3 <- TOAST::csTest(STAGE_estProp_RF_improved_5_fitted_model, 
                                                                   coef = "stage", 
                                                                   cell_type = "CellType3")
  head(STAGE_estProp_RF_improved_5_res_table_CellType3, 3)
  range(STAGE_estProp_RF_improved_5_res_table_CellType3$fdr) # 1.143896e-05 9.999974e-01
  length(STAGE_estProp_RF_improved_5_res_table_CellType3$p_value < 0.05) # 559,110
  
  # CellType4
  STAGE_estProp_RF_improved_5_res_table_CellType4 <- TOAST::csTest(STAGE_estProp_RF_improved_5_fitted_model, 
                                                                   coef = "stage", 
                                                                   cell_type = "CellType4")
  head(STAGE_estProp_RF_improved_5_res_table_CellType4, 3)
  range(STAGE_estProp_RF_improved_5_res_table_CellType4$fdr) # 1.619254e-07 9.999974e-01
  length(STAGE_estProp_RF_improved_5_res_table_CellType4$p_value < 0.05) # 559110
  
  # CellType5
  STAGE_estProp_RF_improved_5_res_table_CellType5 <- TOAST::csTest(STAGE_estProp_RF_improved_5_fitted_model, 
                                                                   coef = "stage", 
                                                                   cell_type = "CellType5")
  head(STAGE_estProp_RF_improved_5_res_table_CellType5, 3)
  range(STAGE_estProp_RF_improved_5_res_table_CellType5$fdr) # 0.0005062402 0.9999978604
  length(STAGE_estProp_RF_improved_5_res_table_CellType5$p_value < 0.05) # 559110
  
  
  
  #--
  # Compare estimateCellCounts with purity estimates from methylation arrays
  #--  
  
  # read purity estimates from TGL
  purity <- read.delim(file = "purity_Kridel_ifpf-DLBC.txt")
  purity$SampleID
  dim(purity) # 177   2
  length(unique(purity$SampleID)) # 172
  purity_names <- c(substr(purity$SampleID[1:78], 5, 13), 
                    substr(purity$SampleID[79:92], 1, 9),
                    substr(purity$SampleID[93:107], 1, 10),
                    substr(purity$SampleID[108:177], 1, 9))
  purity[, 1] <- purity_names
  length(which(substr(colnames(Output_Remove_2_BetaMatrix), 4, 5) == "FL"))
  Output_Remove_2_BetaMatrix_editnames <- Output_Remove_2_BetaMatrix
  
  colnames(Output_Remove_2_BetaMatrix_editnames)[which
    (substr(colnames(Output_Remove_2_BetaMatrix_editnames), 4, 5) == "FL")] <-
    substr(colnames(Output_Remove_2_BetaMatrix_editnames)[which
    (substr(colnames(Output_Remove_2_BetaMatrix_editnames), 4, 5) == "FL")], 1, 9)
  
  
  match_beta_bypurity <- match(colnames(Output_Remove_2_BetaMatrix_editnames), purity[, 1])
  purity_edited <- purity[match_beta_bypurity[!is.na(match_beta_bypurity)], ]
  dim(purity_edited) # 170   2
  length(unique(purity_edited$SampleID)) # 170
  purity_edited
  
  estProp_RF_improved_5_edited <- estProp_RF_improved_5 
  rownames(estProp_RF_improved_5_edited)[which
    (substr(rownames(estProp_RF_improved_5_edited), 4, 5) == "FL")] <-
    substr(rownames(estProp_RF_improved_5_edited)[which
    (substr(rownames(estProp_RF_improved_5_edited), 4, 5) == "FL")], 1, 9)
  
  match_purity_byProp <- match(rownames(estProp_RF_improved_5_edited), purity_edited[,1])
  purity_edited <- purity_edited[match_purity_byProp, ]
  
  # combine data 
  merged_purity_estProp <- cbind(purity_edited, 
                                 estProp_RF_improved_5_edited, 
                                 Output_Remove_2_ClinicalFile$TYPE, 
                                 Output_Remove_2_ClinicalFile$STAGE, 
                                 as.character(Output_Remove_2_ClinicalFile$TRANSLOC_14_18) )
  dim(merged_purity_estProp) # 170  10
  colnames(merged_purity_estProp)[8] <- "TYPE"
  colnames(merged_purity_estProp)[9] <- "STAGE"
  colnames(merged_purity_estProp)[10] <-"TRANSLOC1418"
  head(merged_purity_estProp)

  # plotting only CellType1
  ggplot2::ggplot(merged_purity_estProp, 
                  aes(x = purity, y = CellType1)) +
                  geom_point() +
                  xlim(0,1) +
                  ylim(0,1) +
                  stat_smooth(method = "lm") +
                  theme_bw() +
                  labs(y = "proportion")
  
  # plotting only CellType2
  ggplot2::ggplot(merged_purity_estProp, 
                  aes(x = purity, y = CellType2)) +
                  geom_point() +
                  xlim(0,1) +
                  ylim(0,1) +
                  stat_smooth(method = "lm") +
                  theme_bw() +
                  labs(y = "proportion")
  
  # plotting only CellType3
  ggplot2::ggplot(merged_purity_estProp, 
                  aes(x = purity, y = CellType3)) +
                  geom_point() +
                  xlim(0,1) +
                  ylim(0,1) +
                  stat_smooth(method = "lm") +
                  theme_bw() +
                  labs(y = "proportion")
  
  # plotting all Cell Types
  ggplot2::ggplot(as.data.frame(merged_purity_estProp), 
                  aes(x = purity, y = value, color = variable)) +
                  geom_point(aes(y = CellType1, col = "CellType1")) + 
                  geom_point(aes(y = CellType2, col = "CellType2")) +
                  geom_point(aes(y = CellType3, col = "CellType3")) +
                  geom_point(aes(y = CellType4, col = "CellType4")) +
                  geom_point(aes(y = CellType5, col = "CellType5")) +
                  xlim(0,1) +
                  ylim(0,1) +
                  #stat_smooth(method = "lm") +
                  theme_bw() +
                  labs(y = "proportion") +
                  theme_classic()
  
  
  stats::cor.test(merged_purity_estProp$Purity, 
                  merged_purity_estProp$estPropCellType1,
                  method = c("spearman"))
  # p-value = 0.001925
  #  rho 0.2362496 

  #--
  # Compare immune cell subsets by STAGE
  #--
  
  # Melting data
  merged_purity_estProp_melted <- reshape2::melt(merged_purity_estProp[, -2])
  head(merged_purity_estProp_melted)
  # plotting by STAGE
  ggplot2::ggplot(merged_purity_estProp_melted, aes(x = variable, y = value, color = STAGE)) +
                  geom_boxplot() +
                  stat_compare_means(aes(group = STAGE)) +
                  labs(y = "proportion") +
                  theme_bw() +
                  theme_classic() 

  # plotting by TYPE
  ggplot2::ggplot(merged_purity_estProp_melted, aes(x = variable, y = value, color = TYPE)) +
                  geom_boxplot() +
                  stat_compare_means(aes(group = TYPE)) +
                  labs(y = "proportion") +
                  theme_bw() +
                  theme_classic() 
  
  ggplot2::ggplot(merged_purity_estProp_melted, 
                  aes(x = variable, y = value, color = TRANSLOC1418)) +
                  geom_boxplot() +
                  stat_compare_means(aes(group = TRANSLOC1418)) +
                  labs(y = "proportion") +
                  theme_bw() +
                  theme_classic()
  
  
  return(NULL)
  
}
# [END]  
  
