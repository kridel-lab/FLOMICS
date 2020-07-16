# Updated 27 April 2019
# Function: Obtain probes corresponding to CategoryToVisualize for Beta and Mvalue matrices.
# Author: Anjali Silva

# Input:
# CategoryToVisualize: Specifies which column needs to be visualized from AnnotationFile.
#                      Need exact spelling and caps/noCaps for CategoryToVisualize, e.g.: CategoryToVisualize="DMR"
#                      Options "Relation_to_Island", "chr", "Regulatory_Feature_Group", "X450k_Enhancer", 
#                      "Phantom4_Enhancers", "Phantom5_Enhancers".
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns.
# MvalueMatrix: A matrix of M values for probes x patients, with probes in rows and patients as columns.
# AnnotationFile: Matrix of annotations for all the probes found in BetaMatrix. It is of size probes x annotations.
# ClinicalFile:  File of patients and corresponding clinical category; matrix of size patients x clinical categories
# ProvideWithinCategories: Should indicate "Yes" or "No", as to whether categories within each category need to be investigated.

# Output:
# Matrix: A list containing BetaMatrix_CategoryToVisualize, MvalueMatrix_CategoryToVisualize, Category
#         BetaMatrix_CategoryToVisualize: Beta Matrix with probes corresponding to CategoryToVisualize 
#         if ProvideWithinCategories is set to "No". Otherwise, return Beta Matrix with probes 
#         corresponding to different categories within ProvideWithinCategories. 
#         MvalueMatrix_CategoryToVisualize Category: Return input CategoryToVisualize if 
#         ProvideWithinCategories is set to "No". Otherwise, return the different categories within
#         ProvideWithinCategories.  
# ProvideWithinCategories: Return if the input was "Yes" or "No".


IsolateEntries <- function(CategoryToVisualize, 
                           BetaMatrix, 
                           MvalueMatrix,
                           AnnotationFile, 
                           ClinicalFile, 
                           ProvideWithinCategories = "Yes") {
  
  # Loading needed packages
  # LoadCheckPkg(RegularPckgs=c("stringi","data.table","ggpubr"))
  # "stringi" used to remove extra white space, if present, e.g. "EN " to "EN"
  library(stringi)
  library(data.table)
  library(ggpubr)
  
  # Remove empty entries in the clinical file
  ClinicalFile[ClinicalFile == ""] <- NA
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  # Matching rownames in methylation file with corresponding location in annotation file
  match_ids_methylation_beta <- match(rownames(BetaMatrix), AnnotationFile[, 1])
  # match(rownames(BetaMatrix), rownames(MvalueMatrix))
  
  # Bringing annotation file in the same order as BetaMatrix
  AnnotationFile2 <- AnnotationFile[match_ids_methylation_beta, ] 
  
  # Identify column number matching the CategoryToVisualize
  columnNumber <- which(colnames(AnnotationFile2) == CategoryToVisualize)
  
  # Only take in the non empty entries of the selected column
  non_empty_entries <- which(AnnotationFile2[, columnNumber] != "")

  Matrix <- list(BetaMatrix_CategoryToVisualize = BetaMatrix[non_empty_entries, ], 
                 MvalueMatrix_CategoryToVisualize = MvalueMatrix[non_empty_entries, ],
                 Category = CategoryToVisualize)
  
  if (ProvideWithinCategories == "Yes") {
    # Defining function for plotting categories within
    WithinCategories_Yes <- function(AnnotationCategory, BetaMatrix, non_empty_entries, ClinicalFile) {
      
      # AnnotationCategory is the category within the annotation file that is being analyzed
      # E.g. table(AnnotationFile2[non_empty_entries,columnNumber]) gives CDMR   DMR  RDMR
      # Then set AnnotationCategory="CDMR"
      
      # rows corresponding to the current AnnotationCategory
      correspondingRows <- which(AnnotationFile2[non_empty_entries, columnNumber] == AnnotationCategory)
      
      # Initialize average_cat3
      # average_cat3<-0
      
      # Looking at all probes belonging to ClinicalCategory (e.g. if ClinicalCategory = STAGE, then look at both ADVANCED LIMITED )
      # ClinicalCategory options "SEX","STAGE","SITE_BIOPSY","TYPE_BIOPSY","COO","TRANSLOC_14_18","INSTITUTION","EPIC_QC","TYPE"
      
      # vectorNames<-unique(na.omit(ClinicalFile[,which(colnames(ClinicalFile)==ClinicalCategory)]))
      # average_cat1 <- as.matrix(BetaMatrix[correspondingRows,which(stri_trim(ClinicalFile[c(which(!is.na(ClinicalFile[,which(colnames(ClinicalFile)==ClinicalCategory)])==TRUE)),which(colnames(ClinicalFile)==ClinicalCategory)])==vectorNames[1])])
      # average_cat2 <- as.matrix(BetaMatrix[correspondingRows,which(stri_trim(ClinicalFile[c(which(!is.na(ClinicalFile[,which(colnames(ClinicalFile)==ClinicalCategory)])==TRUE)),which(colnames(ClinicalFile)==ClinicalCategory)])==vectorNames[2])])
      # if(length(vectorNames)==3){
      #   average_cat3 <- as.matrix(BetaMatrix[correspondingRows,which(stri_trim(ClinicalFile[c(which(!is.na(ClinicalFile[,which(colnames(ClinicalFile)==ClinicalCategory)])==TRUE)),which(colnames(ClinicalFile)==ClinicalCategory)])==vectorNames[3])])
      # }else if(length(vectorNames)==4){
      #    average_cat3 <- as.matrix(BetaMatrix[correspondingRows,which(stri_trim(ClinicalFile[c(which(!is.na(ClinicalFile[,which(colnames(ClinicalFile)==ClinicalCategory)])==TRUE)),which(colnames(ClinicalFile)==ClinicalCategory)])==vectorNames[3])])
      # average_cat4 <- as.matrix(BetaMatrix[correspondingRows,which(stri_trim(ClinicalFile[c(which(!is.na(ClinicalFile[,which(colnames(ClinicalFile)==ClinicalCategory)])==TRUE)),which(colnames(ClinicalFile)==ClinicalCategory)])==vectorNames[4])])
      # }else if(length(vectorNames)>=5){
      #   stop("More than 4 categories under ClinicalCategory. Code needs to be updated")
      # }
      
      RESULTS <- list(BetaMatrix_AnnotationCategory = BetaMatrix[correspondingRows, ],
                      MvalueMatrix_AnnotationCategory = MvalueMatrix[correspondingRows, ],
                      Category = AnnotationCategory)
      class(RESULTS) <- paste(AnnotationCategory)
      return(RESULTS)
    }
    
    CategoriesToPlot <- unique(AnnotationFile2[non_empty_entries, columnNumber])
    Matrix <- lapply(CategoriesToPlot,
                     WithinCategories_Yes, 
                     BetaMatrix = BetaMatrix,
                     non_empty_entries = non_empty_entries,
                     ClinicalFile = ClinicalFile)
  }
  
  
  RESULTS <- list(Matrix = Matrix,
                  ProvideWithinCategories = ProvideWithinCategories)
  
  class(RESULTS) <- "IsolateEntries_ASilva"
  return(RESULTS)
}
