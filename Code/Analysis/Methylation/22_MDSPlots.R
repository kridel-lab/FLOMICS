# Updated 17 May 2019
# Function: Generate MDS plots for a specified ClinicalCategoryToVisualize. 
# Author: Anjali Silva

# Input:
# ClinicalCategoryToVisualize: Specify column of ClinicalFile to use. Options "STAGE", "SEX", 
#                              "SITE_BIOPSY", "TYPE_BIOPSY", "INSTITUTION", "COO", "TYPE", "TRANSLOC_14_18".
#                               Default: "TYPE".
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns
# AnnotationFile: Matrix of annotations for all the probes found in BetaMatrix. It is of size probes x annotations.
# ClinicalFile: File with patient sample names along rows, and categories on columns. 
# SampleSheet: Illumina methylation sample sheet as obtained from 1_QCRemoveSamples().
# FigureGenerate: Should figures be generated, default = "No". Options: "Yes", "No".
# PNGorPDF: Output format of the image, options = "png" or "pdf".

# Output: None, plots are saved

# Visuals saved to img folder
# 22_MvalueMDS_", ClinicalCategoryToVisualize,*

MDSplots <- function(ClinicalCategoryToVisualize = "TYPE", 
                     BetaMatrix, 
                     AnnotationFile, 
                     ClinicalFile, 
                     SampleSheet, 
                     FigureGenerate = "Yes", 
                     PNGorPDF = "png") {

  # Loading needed packages
  # RegularPckgs=c("minfi","missMethyl","limma")
  library(limma)
  library(missMethyl)
  library(minfi)
  
  # Remove empty entries in the clinical file
  ClinicalFile[ClinicalFile == ""] <- NA
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  # Matching rownames in methylation file with corresponding location in annotation file
  match_ids_methylation_beta <- match(rownames(BetaMatrix), AnnotationFile[ , 1])
  
  # Bringing annotation file in the same order as BetaMatrix
  AnnotationFile2 <- AnnotationFile[match_ids_methylation_beta, ] 
  
  # Identify column number matching the AnnotationCategoryToVisualize
  columnNumber <- which(colnames(AnnotationFile2) == AnnotationCategoryToVisualize)
  
  # Only take in the non empty entries of the selected column
  non_empty_entries <- which(AnnotationFile2[ , columnNumber] != "")
  
  
  # Obtain the unique clinical categories 
  vectorNames <- unique(na.omit(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]))
  
  # all probes and samples only corresponding to those specified in ClinicalCategoryToVisualize
  average_cat11 <- as.matrix(BetaMatrix[ , which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == TRUE)),
                                                                       which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == vectorNames[1])])
  average_cat22 <- as.matrix(BetaMatrix[ , which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == TRUE)),
                                                                       which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == vectorNames[2])])
  if (length(vectorNames) > 2) {
    average_cat33 <- as.matrix(BetaMatrix[ , which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[ , which(colnames(ClinicalFile)==ClinicalCategoryToVisualize)]) == TRUE)),
                                                                        which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == vectorNames[3])])
  }
  if (length(vectorNames) > 3) {
    average_cat44 <- as.matrix(BetaMatrix[correspondingRows, which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == TRUE)), 
                                                                                         which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == vectorNames[4])])
  }
  

  if(FigureGenerate == "Yes") {
    if(PNGorPDF == "pdf") {   
      pdf(paste0(pathNow, "/img/22_MvalueMDS_", ClinicalCategoryToVisualize,".pdf"), 
          width = 60, height = 60, pointsize = 50)
    }else{
      png(paste0(pathNow, "/img/22_MvalueMDS_", ClinicalCategoryToVisualize,".png"))
    }
    par(mfrow = c(2, length(vectorNames)))
    
    limma::plotMDS(average_cat11, labels = substr(colnames(average_cat11), 4, 10), col = 1)
    limma::plotMDS(average_cat22, labels = substr(colnames(average_cat22), 4, 10), col = 2)
    if (length(vectorNames) > 2) {
      limma::plotMDS(average_cat33, labels = substr(colnames(average_cat33), 4, 10), col = 3)
    }
    if (length(vectorNames) > 3) {
      limma::plotMDS(average_cat44, labels = substr(colnames(average_cat44), 4, 10), col = 4)
    }
    limma::plotMDS(BetaMatrix, 
                   labels = substr(colnames(BetaMatrix), 4, 10), 
                   col = as.integer(factor(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)])))
    graphics::legend("topright", legend = vectorNames, pch = 16, cex = 0.8, 
                     col = unique(as.integer(factor(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]))))
    grDevices::dev.off()
  }
  
  return(NULL)
}


