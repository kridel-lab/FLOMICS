# Updated 3 Aug 2021
# Updated 12 Feb 2019
# Function: Summarize Beta values (mainly used for visualization) by genes not probes.
# Author: Anjali Silva

# Input:
# BetaMatrix: A matrix of beta values for probes x patients, with probes 
#             in rows and patients as columns.
# ClinicalFile:  File of patients and corresponding clinical category; 
#                matrix of size patients x clinical categories.

# Output: 
# SummarizedBetaMatrix: A matrix of size genes x patients, with summarized probe data.

SummarizationBeta5 <- function(BetaMatrix, 
                               ClinicalFile) {
  
  # Loading needed packages
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  library(WGCNA)
  
  # getting the annotation file (also saved in Methylation data folder)
  ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

  # rownames(BetaMatrix) # to get the probenames in the Mvalues dataset

  # Matching rownames in methylation file with corresponding probe location in annotation file
  # "Name defines the column with probe names in annotation file
  matchIdsMethylationBeta <- match(rownames(BetaMatrix), ann[, which(colnames(ann) == "Name")])
  
  # Using matched rownames, get the gene names from annotation file for corresponding probe names in methylation file
  # colnames(ann) # UCSC_RefGene_Name define genes, which is column 22

  # some annotations are missing. This was confirmed using excel sheet too. This is due to the way annotation file is.
  # UCSC_RefGene_Name defines the name of the gene in annotation file 
  missingAnnotationMethylationBeta <- which(ann[matchIdsMethylationBeta, 
                                      which(colnames(ann) == "UCSC_RefGene_Name")] == "")

  # getting annotations with missing gene rows removed
  annMissingRemovedBeta <- ann[matchIdsMethylationBeta[- missingAnnotationMethylationBeta], ]
  # gettting probes, such that probes with missing annotation removed
  betaLymphMissingRemoved <- BetaMatrix[- missingAnnotationMethylationBeta, ]

  # only taking first feature
  annMissingRemovedTruncated <- sub ("\\;.*", "" , annMissingRemovedBeta[ , 
                                which(colnames(annMissingRemovedBeta) == "UCSC_RefGene_Name")]) 
  
  # CollapseRows allows to select one representative row per group
  betaMissingremovedCollapsedTruncated <- WGCNA::collapseRows(betaLymphMissingRemoved,
                                                              rowGroup = annMissingRemovedTruncated,
                                                              rowID = rownames(betaLymphMissingRemoved),
                                                              method = "Average")
  
  betaMissingremovedCollapsed2Truncated <- data.frame(betaMissingremovedCollapsedTruncated$group2row, 
                                                      betaMissingremovedCollapsedTruncated$datETcollapsed)
  
  colnames(betaMissingremovedCollapsed2Truncated)[1] <- "NAME"
  
  betaMissingremovedCollapsed2Truncated$selectedRowID <- NULL
  
  betaMissingremovedCollapsed2Truncated <- betaMissingremovedCollapsed2Truncated[, 
                                           order(names(betaMissingremovedCollapsed2Truncated))]
  
  betaMissingremovedCollapsed2Truncated <- betaMissingremovedCollapsed2Truncated[, 
                                           colnames(betaMissingremovedCollapsed2Truncated) != "NAME"]
  
  
  RESULTS <- list(SummarizedBetaMatrix = betaMissingremovedCollapsed2Truncated)
  
  class(RESULTS) <- "SummarizationBeta_ASilva"
  return(RESULTS) 
}
# [END]
