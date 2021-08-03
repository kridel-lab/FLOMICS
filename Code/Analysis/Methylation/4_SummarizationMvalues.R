# Updated 3 Aug 2021
# Updated 12 Feb 2019
# Function: Summarize M values (mainly used for statistical analysis) by genes rather than probes. 
# Author: Anjali Silva

# Input:
# MvalueMatrix: A matrix of probes x patients, with probes in rows and patients as columns.
# ClinicalFile: File of patients and corresponding clinical category; matrix of size 
#               patients x clinical categories.

# Output: 
# SummarizedMvalueMatrix: A matrix of size genes x patients, with summarized probe data for M values.

SummarizationMvalues4 <- function(MvalueMatrix, 
                                  ClinicalFile) {
  
  # Loading needed packages
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  library(WGCNA)
  
  # Summarize probes to genes
  ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  
  # Matching rownames in methylation file with corresponding location in annotation file
  matchIdsMethylation <- match(rownames(MvalueMatrix), ann$Name)
  
  # Using matched rownames, get the gene names from annotation file for corresponding probe names in methylation file
  # colnames(ann) # UCSC_RefGene_Name define genes, which is column 22
  #length(ann[matchIdsMethylation,22]) # 569817
  
  # some annotations are missing. This was confirmed using excel sheet too. This is due to the way annotation file is.
  missingAnnotationMethylation <- which(ann[matchIdsMethylation, 22] == "")
  # length(missingAnnotationMethylation) # 160094
  
  # getting annotations with missing gene rows removed
  annMissingRemoved <- ann[matchIdsMethylation[- missingAnnotationMethylation], ]
  # dim(annMissingRemoved)[1] # 46; this makes sense as 569817 - 160094 = 409723
  # gettting probes, such that probes with missing annotation removed
  MvalueMatrixMissingRemoved <- MvalueMatrix[- missingAnnotationMethylation, ]
  # dim(MvalueMatrixMissingRemoved) #  409723     173
  # match(rownames(MvalueMatrixMissingRemoved), annMissingRemoved[,4]) # in order
  
  # only taking first feature
  annMissingRemovedTruncated <- sub ("\\;.*", "" , annMissingRemoved[ , 22] ) 
  # length(annMissingRemovedTruncated) # 409723
  
  # CollapseRows allows to select one representative row per group
  MvaluesMissingremovedCollapsedTruncated <- WGCNA::collapseRows(MvalueMatrixMissingRemoved,
                                                             rowGroup = annMissingRemovedTruncated,
                                                             rowID = rownames(MvalueMatrixMissingRemoved),
                                                             method = "Average")
  MvaluesMissingremovedCollapsed2Truncated <- data.frame(MvaluesMissingremovedCollapsedTruncated$group2row, 
                                                            MvaluesMissingremovedCollapsedTruncated$datETcollapsed)
  
  colnames(MvaluesMissingremovedCollapsed2Truncated)[1] <- "NAME"
  
  MvaluesMissingremovedCollapsed2Truncated$selectedRowID <- NULL
  
  MvaluesMissingremovedCollapsed2Truncated <- MvaluesMissingremovedCollapsed2Truncated[ , 
                                                 order(names(MvaluesMissingremovedCollapsed2Truncated))]
  
  MvaluesMissingremovedCollapsed2Truncated <- MvaluesMissingremovedCollapsed2Truncated[ , 
                                                 colnames(MvaluesMissingremovedCollapsed2Truncated) != "NAME"]
  
  # dim(MvaluesMissingremovedCollapsed2Truncated) #  26589    30
  
  # tail(MvaluesMissingremovedCollapsed2Truncated)
  
  RESULTS <- list(SummarizedMvalueMatrix = MvaluesMissingremovedCollapsed2Truncated)
  
  class(RESULTS) <- "SummarizationMvalues_ASilva"
  return(RESULTS) 
}
# [END]
