# Updated 12 Feb 2019
# Function: Summarize M values (mainly used for statistical analysis) by genes rather than probes. 
# Author: Anjali Silva

# Input:
# MvalueMatrix: A matrix of probes x patients, with probes in rows and patients as columns.
# ClinicalFile: File of patients and corresponding clinical category; matrix of size 
#               patients x clinical categories.

# Output: 
# SummarizedMvalueMatrix: A matrix of size genes x patients, with summarized probe data for M values.

SummarizationMvalues <- function(MvalueMatrix, 
                                 ClinicalFile) {
  
  # Loading needed packages
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  library(WGCNA)
  
  # Summarize probes to genes
  ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  
  # Matching rownames in methylation file with corresponding location in annotation file
  match_ids_methylation <- match(rownames(MvalueMatrix), ann$Name)
  
  # Using matched rownames, get the gene names from annotation file for corresponding probe names in methylation file
  # colnames(ann) # UCSC_RefGene_Name define genes, which is column 22
  #length(ann[match_ids_methylation,22]) # 569817
  
  # some annotations are missing. This was confirmed using excel sheet too. This is due to the way annotation file is.
  missing_annotation_methylation <- which(ann[match_ids_methylation, 22] == "")
  # length(missing_annotation_methylation) # 160094
  
  # getting annotations with missing gene rows removed
  ann_missing_removed <- ann[match_ids_methylation[- missing_annotation_methylation], ]
  # dim(ann_missing_removed)[1] # 46; this makes sense as 569817 - 160094 = 409723
  # gettting probes, such that probes with missing annotation removed
  MvalueMatrix_missing_removed <- MvalueMatrix[- missing_annotation_methylation, ]
  # dim(MvalueMatrix_missing_removed) #  409723     173
  # match(rownames(MvalueMatrix_missing_removed), ann_missing_removed[,4]) # in order
  
  # only taking first feature
  ann_missing_removed_truncated <- sub ("\\;.*", "" , ann_missing_removed[ , 22] ) 
  # length(ann_missing_removed_truncated) # 409723
  
  # CollapseRows allows to select one representative row per group
  Mvalues_missingremoved_collapsed_truncated <- WGCNA::collapseRows(MvalueMatrix_missing_removed,
                                                             rowGroup = ann_missing_removed_truncated,
                                                             rowID = rownames(MvalueMatrix_missing_removed),
                                                             method = "Average")
  Mvalues_missingremoved_collapsed2_truncated <- data.frame(Mvalues_missingremoved_collapsed_truncated$group2row, 
                                                            Mvalues_missingremoved_collapsed_truncated$datETcollapsed)
  
  colnames(Mvalues_missingremoved_collapsed2_truncated)[1] <- "NAME"
  
  Mvalues_missingremoved_collapsed2_truncated$selectedRowID <- NULL
  
  Mvalues_missingremoved_collapsed2_truncated <- Mvalues_missingremoved_collapsed2_truncated[ , 
                                                 order(names(Mvalues_missingremoved_collapsed2_truncated))]
  
  Mvalues_missingremoved_collapsed2_truncated <- Mvalues_missingremoved_collapsed2_truncated[ , 
                                                 colnames(Mvalues_missingremoved_collapsed2_truncated) != "NAME"]
  
  # dim(Mvalues_missingremoved_collapsed2_truncated) #  26589    30
  
  # tail(Mvalues_missingremoved_collapsed2_truncated)
  
  RESULTS <- list(SummarizedMvalueMatrix = Mvalues_missingremoved_collapsed2_truncated)
  
  class(RESULTS) <- "SummarizationMvalues_ASilva"
  return(RESULTS) 

  
}



