# Updated 12 Feb 2019
# Function: Summarize Beta values (mainly used for visualization) by genes not probes
# Author: Anjali Silva

# Input:
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns
# ClinicalFile:  File of patients and corresponding clinical category; matrix of size patients x clinical categories

# Output: 
# SummarizedBetaMatrix: A matrix of size genes x patients, with summarized probe data

SummarizationBeta <- function(BetaMatrix, 
                              ClinicalFile) {
  
  # Loading needed packages
  # RegularPckgs=c("minfi","IlluminaHumanMethylationEPICanno.ilm10b2.hg19","WGCNA")
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  library(WGCNA)
  
  # getting the annotation file (also saved in Methylation data folder)
  ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

  # rownames(BetaMatrix) # to get the probenames in the Mvalues dataset

  # Matching rownames in methylation file with corresponding probe location in annotation file
  # "Name defines the column with probe names in annotation file
  match_ids_methylation_beta<-match(rownames(BetaMatrix), ann[, which(colnames(ann) == "Name")])
  
  # Using matched rownames, get the gene names from annotation file for corresponding probe names in methylation file
  # colnames(ann) # UCSC_RefGene_Name define genes, which is column 22

  # some annotations are missing. This was confirmed using excel sheet too. This is due to the way annotation file is.
  # UCSC_RefGene_Name defines the name of the gene in annotation file 
  missing_annotation_methylation_beta <- which(ann[match_ids_methylation_beta, 
                                                   which(colnames(ann) == "UCSC_RefGene_Name")] == "")

  # getting annotations with missing gene rows removed
  ann_missing_removed_beta <- ann[match_ids_methylation_beta[- missing_annotation_methylation_beta], ]
  # gettting probes, such that probes with missing annotation removed
  beta_lymph_missing_removed <- BetaMatrix[-missing_annotation_methylation_beta, ]

  # only taking first feature
  ann_missing_removed_truncated = sub ("\\;.*", "" , ann_missing_removed_beta[, which(colnames(ann_missing_removed_beta) == "UCSC_RefGene_Name")]) 
  
  # CollapseRows allows to select one representative row per group
  beta_missingremoved_collapsed_truncated <- collapseRows(beta_lymph_missing_removed,
                                                          rowGroup = ann_missing_removed_truncated,
                                                          rowID = rownames(beta_lymph_missing_removed),
                                                          method = "Average")
  
  beta_missingremoved_collapsed2_truncated <- data.frame(beta_missingremoved_collapsed_truncated$group2row, 
                                                         beta_missingremoved_collapsed_truncated$datETcollapsed)
  
  colnames(beta_missingremoved_collapsed2_truncated)[1] <- "NAME"
  
  beta_missingremoved_collapsed2_truncated$selectedRowID <- NULL
  
  beta_missingremoved_collapsed2_truncated <- beta_missingremoved_collapsed2_truncated[, order(names(beta_missingremoved_collapsed2_truncated))]
  
  beta_missingremoved_collapsed2_truncated <- beta_missingremoved_collapsed2_truncated[, colnames(beta_missingremoved_collapsed2_truncated) != "NAME"]
  
  
  RESULTS <- list(SummarizedBetaMatrix = beta_missingremoved_collapsed2_truncated)
  
  class(RESULTS) <- "SummarizationBeta_ASilva"
  return(RESULTS) 
}
