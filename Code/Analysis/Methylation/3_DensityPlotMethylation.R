# Updated 3 Aug 2021
# Updated 5 Feb 2020
# Updated 11 Feb 2019
# Function: Density plot based on annotation category, clusters or clinical category. 
# Author: Anjali Silva


# Input:
# AnnotationCategoryToVisualize: Specifies which column needs to be visualized from AnnotationFile.
#                      Need exact spelling and caps/noCaps for AnnotationCategoryToVisualize, 
#                      e.g.: AnnotationCategoryToVisualize="DMR"
#                      Options "UCSC_RefGene_Group", Relation_to_Island", "chr", "Regulatory_Feature_Group", "X450k_Enhancer", 
#                      "Phantom4_Enhancers", "Phantom5_Enhancers". By default set to "None", in which case 
#                      all probes are visualized with respect to selected Clinical Category. If set to "None",
#                      then argument PlotWithinAnnotationCategories should be set to "No".
# ClusterLabels: A vector of cluster labels that is the same length as the number of samples (dimensionality) 
#                of Beta matrix. If provided, ClinicalCategoryToVisualize will be ignored. Default to "NA". 
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns
# AnnotationFile: Matrix of annotations for all the probes found in BetaMatrix. It is of size probes x annotations.
# ClinicalFile:  File of patients and corresponding clinical category; matrix of size patients x clinical categories
# ClinicalCategoryToVisualize: Specify column of ClinicalFile to use. Options "STAGE", "SEX", "SITE_BIOPSY", 
#                              "TYPE_BIOPSY", "INSTITUTION", "COO", "TYPE", "TRANSLOC_14_18".
# PlotWithinAnnotationCategories: Should indicate "Yes" or "No", as to whether categories within each category 
#                                 to visualize need to be investigated. By default set to "No".
# FigureGenerate: Should figures be generated, default = "No". Options: "Yes", "No".
# ImageName: A character vector indicating a unique for saving the plots. 
# PNGorPDF: Output format of the image, options = "png" or "pdf".

# Output: None, plots are saved

MethylationDensityPlot3 <- function(AnnotationCategoryToVisualize = "None", 
                                   ClusterLabels = NA, 
                                   PlotWithinAnnotationCategories = "No", 
                                   ClinicalCategoryToVisualize = "TYPE", 
                                   BetaMatrix, 
                                   AnnotationFile = NA, 
                                   ClinicalFile = NA, 
                                   SampleSheet = NA, 
                                   FigureGenerate = "Yes", 
                                   ImageName = NA, 
                                   PNGorPDF = "png") {
  
  # Loading needed packages
  library(minfi)
  library(missMethyl)
  library(limma)
  library(stringi)
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  if(typeof(AnnotationFile) == "list") {
    geneRegion <- sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Group)
    geneRegion[which(geneRegion == "")] <- "NA"
    AnnotationFile$UCSC_RefGene_Group <- geneRegion
  }
  
  # Check points
  if (AnnotationCategoryToVisualize == "None" && PlotWithinAnnotationCategories != "No") {
    stop(paste("Conflicting argument entries for AnnotationCategoryToVisualize:", 
               AnnotationCategoryToVisualize, "and PlotWithinAnnotationCategories:", 
               PlotWithinAnnotationCategories,". \n Fix and rerun."));
  }
  #if (AnnotationCategoryToVisualize != "None" && PlotWithinAnnotationCategories == "No") {
  #  stop(paste("Conflicting argument entries for AnnotationCategoryToVisualize:", AnnotationCategoryToVisualize, "and PlotWithinAnnotationCategories:", PlotWithinAnnotationCategories,". \n Fix and rerun."));
  #}
  
  # if ClusterLabels provided
  if (sum(is.na(ClusterLabels)) == 0) {
    
    # all probes and samples only corresponding to those specified in Clusters
    averageCat11 <- as.matrix(BetaMatrix[ , which(ClusterLabels == 1)])
    averageCat22 <- as.matrix(BetaMatrix[ , which(ClusterLabels == 2)])
    if (max(unique(ClusterLabels)) > 2) {
      averageCat33 <- as.matrix(BetaMatrix[ , which(ClusterLabels == 3)])
    }
    if (max(unique(ClusterLabels)) > 3) {
      averageCat44 <- as.matrix(BetaMatrix[ , which(ClusterLabels == 4)])
    }
    if (max(unique(ClusterLabels)) > 4) {
      average_cat55 <- as.matrix(BetaMatrix[ , which(ClusterLabels == 5)])
    }
    if (max(unique(ClusterLabels)) > 5) {
      average_cat66 <- as.matrix(BetaMatrix[ , which(ClusterLabels == 6)])
    }
    if (max(unique(ClusterLabels)) > 6) {
      stop("\nCode only support upto 6 clusters. Please alter code");
    }
    
    if (FigureGenerate == "Yes") {
      # Plotting density distributions
      if (PNGorPDF == "pdf") {
        pdf(paste0(pathNow,"/img/3_DensityBeta_ClusterLabels_", ImageName, "_only.pdf"), 
            width = 60, height = 60, pointsize = 50)
      } else {
        png(paste0(pathNow,"/img/3_DensityBeta_ClusterLabels_", ImageName, "only.png"))
      }
      
      if (max(unique(ClusterLabels)) <= 2) {
        par(mfrow = c(1, 2))        
      } else if (max(unique(ClusterLabels)) <= 4) {
        par(mfrow = c(2, 2))     
      } else {
        par(mfrow = c(2, 3)) 
      }
      
      minfi::densityPlot(dat = as.matrix(averageCat11), main = paste0("Density of Cluster 1"))
      minfi::densityPlot(dat = as.matrix(averageCat22), main = paste0("Density of Cluster 2"))
      if (max(unique(ClusterLabels)) > 2) {
        minfi::densityPlot(dat = as.matrix(averageCat33), main = paste0("Density of Cluster 3"))
      }
      if (max(unique(ClusterLabels)) > 3) {
        minfi::densityPlot(dat=as.matrix(averageCat44), main = paste0("Density of Cluster 4"))      
      }
      if (max(unique(ClusterLabels)) > 4) {
        minfi::densityPlot(dat = as.matrix(average_cat55), main = paste0("Density of Cluster 5"))      
      }
      if (max(unique(ClusterLabels)) > 5 ) {
        minfi::densityPlot(dat = as.matrix(average_cat66), main = paste0("Density of Cluster 6"))      
      }
      dev.off()
    }
    
    
    # Remove empty entries in the clinical file
    ClinicalFile[ClinicalFile == ""] <- NA
    
    # Matching rownames in methylation file with corresponding location in annotation file
    match_ids_methylation_beta <- match(rownames(BetaMatrix), AnnotationFile$V1)
    
    # Bringing annotation file in the same order as BetaMatrix
    AnnotationFile2 <- data.frame(AnnotationFile[match_ids_methylation_beta, ]) 
    
    # Identify column number matching the AnnotationCategoryToVisualize
    if(AnnotationCategoryToVisualize != "None") {
      columnNumber = which(colnames(AnnotationFile2) == AnnotationCategoryToVisualize)
      
      # Only take in the non empty entries of the selected column
      nonEmptyEntries <- which(AnnotationFile2[ , columnNumber] != "")
      
      # Defining functions to visualize
      eachplotYes <- function(ClinicalCategory, BetaMatrix, nonEmptyEntries, ClinicalFile, AnnotationCategory) {
        # AnnotationCategory is the category within the annotation file that is being analyzed
        # E.g. table(AnnotationFile2[nonEmptyEntries,columnNumber]) gives CDMR   DMR  RDMR
        # Then set AnnotationCategory="CDMR"
        
        # columns corresponding to the current cluster
        correspondingColumns <- which(ClusterLabels == ClinicalCategory)
        
        # Looking at averge for all probes belonging to each Annotation Category
        # AnnotationCategoryToVisualize = "Relation_to_Island"
        # Island  N_Shore OpenSea N_Shelf S_Shelf S_Shore
        vectorNames <- unique(AnnotationFile2[nonEmptyEntries, columnNumber])
        
        # all probes corresponding to AnnotationCategory and samples corresponding to those specified in vectorNames
        averageCat1 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[1]), 
                                             correspondingColumns])
        averageCat2 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[2]), 
                                             correspondingColumns])
        if (length(vectorNames) > 2) {
          averageCat3 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[3]), 
                                               correspondingColumns])
        }
        if (length(vectorNames) > 3) {
          averageCat4 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[4]), 
                                               correspondingColumns])
        }
        if (length(vectorNames) > 4) {
          averageCat5 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[5]), 
                                               correspondingColumns])
        }
        if (length(vectorNames) > 5) {
          averageCat6 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[6]), 
                                               correspondingColumns])
        }
        if (length(vectorNames) > 6) {
          averageCat7 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[7]), 
                                               correspondingColumns])
        }
        
        if(FigureGenerate == "Yes") {
          # Plotting density distributions
          if(PNGorPDF == "pdf") {
            pdf(paste0(pathNow, "/img/3_DensityBeta_Cluster", ClinicalCategory,"_", AnnotationCategoryToVisualize,".pdf"), 
                width = 60, height = 60, pointsize = 50)
          }else {
            png(paste0(pathNow, "/img/3_DensityBeta_Cluster", ClinicalCategory,"_", AnnotationCategoryToVisualize,".png"))
          }
          
          par(mfrow = c(2, trunc(length(vectorNames) / 2, 1)))
          
          minfi::densityPlot(dat = as.matrix(averageCat1), main = paste0("Density of ", vectorNames[1]))
          minfi::densityPlot(dat = as.matrix(averageCat2), main = paste0("Density of ", vectorNames[2]))
          if (length(vectorNames) > 2) {
            minfi::densityPlot(dat = as.matrix(averageCat3), main = paste0("Density of ", vectorNames[3]))
          }
          if (length(vectorNames) > 3) {
            minfi::densityPlot(dat = as.matrix(averageCat4), main = paste0("Density of ", vectorNames[4]))
          }
          if (length(vectorNames) > 4) {
            minfi::densityPlot(dat = as.matrix(averageCat5), main = paste0("Density of ", vectorNames[5]))
          }
          if (length(vectorNames) > 5) {
            minfi::densityPlot(dat = as.matrix(averageCat6), main = paste0("Density of ", vectorNames[6]))
          }
          if (length(vectorNames) > 6) {
            minfi::densityPlot(dat = as.matrix(averageCat7), main = paste0("Density of ", vectorNames[7]))
          }
          dev.off()
        }
      }
      
      if (FigureGenerate == "Yes") {
        # CategoriesToPlot<-unique(AnnotationFile2[nonEmptyEntries,columnNumber])
        CategoriesToPlot <- unique(na.omit(ClusterLabels))
        plots <- lapply(CategoriesToPlot, 
                        eachplotYes, 
                        BetaMatrix = BetaMatrix,
                        non_empty_entries = nonEmptyEntries,
                        ClinicalFile = ClinicalFile,
                        AnnotationCategory = AnnotationCategoryToVisualize)
      }
    }
  } else {
      # Remove empty entries in the clinical file
      ClinicalFile[ClinicalFile == ""] <- NA
      
      # Matching rownames in methylation file with corresponding location in annotation file
      match_ids_methylation_beta <- match(rownames(BetaMatrix), AnnotationFile$V1)
      
      # Bringing annotation file in the same order as BetaMatrix
      AnnotationFile2 <- data.frame(AnnotationFile[match_ids_methylation_beta, ]) 
      
      # Identify column number matching the AnnotationCategoryToVisualize
      columnNumber <- which(colnames(AnnotationFile2) == AnnotationCategoryToVisualize)
      
      # Only take in the non empty entries of the selected column
      nonEmptyEntries <- which(AnnotationFile2[ , columnNumber] != "")
      
      # Defining functions to visualize
      eachplotYes <- function(ClinicalCategory, 
                               BetaMatrix, 
                               non_empty_entries,
                               ClinicalFile, 
                               AnnotationCategory, 
                               columnNumber) {
        # AnnotationCategory is the category within the annotation file that is being analyzed
        # E.g. table(AnnotationFile2[nonEmptyEntries,columnNumber]) gives CDMR   DMR  RDMR
        # Then set AnnotationCategory="CDMR"
        
        # columns corresponding to the current ClinicalCategory
        correspondingColumns <- which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == TRUE)), 
                                                             which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == ClinicalCategory)
        
        
        # Looking at averge for all probes belonging to each Annotation Category
        # AnnotationCategoryToVisualize = "Relation_to_Island"
        # Island  N_Shore OpenSea N_Shelf S_Shelf S_Shore
        vectorNames <- unique(AnnotationFile2[nonEmptyEntries, columnNumber])
        
        # all probes corresponding to AnnotationCategory and samples corresponding to those specified in vectorNames
        averageCat1 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[1]), 
                                             correspondingColumns])
        averageCat2 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[2]), 
                                             correspondingColumns])
        if (length(vectorNames) > 2) {
          averageCat3 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[3]), 
                                               correspondingColumns])
        }
        if (length(vectorNames) > 3) {
          averageCat4 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[4]), 
                                               correspondingColumns])
        }
        if (length(vectorNames) > 4) {
          averageCat5 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[5]), 
                                               correspondingColumns])
        }
        if (length(vectorNames) > 5) {
          averageCat6 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[6]), 
                                               correspondingColumns])
        }
        if (length(vectorNames) > 6) {
          averageCat7 <- as.matrix(BetaMatrix[which(AnnotationFile2[nonEmptyEntries, columnNumber] == vectorNames[7]), 
                                               correspondingColumns])
        }
        
        if(FigureGenerate == "Yes") {
          # Plotting density distributions
          if(PNGorPDF == "pdf") {
            pdf(paste0(pathNow, "/img/3_DensityBeta_", ClinicalCategory,"_", AnnotationCategoryToVisualize,".pdf"), 
                width = 60, height = 60, pointsize = 50)
          }else {
            png(paste0(pathNow, "/img/3_DensityBeta_", ClinicalCategory,"_", AnnotationCategoryToVisualize,".png"))
          }
          
          par(mfrow = c(2, trunc(length(vectorNames) / 2, 1)))
          
          minfi::densityPlot(dat = as.matrix(averageCat1), main = paste0("Density of ", vectorNames[1]))
          minfi::densityPlot(dat = as.matrix(averageCat2), main = paste0("Density of ", vectorNames[2]))
          if (length(vectorNames) > 2) {
            minfi::densityPlot(dat = as.matrix(averageCat3), main = paste0("Density of ", vectorNames[3]))
          }
          if (length(vectorNames) > 3) {
            minfi::densityPlot(dat = as.matrix(averageCat4), main = paste0("Density of ", vectorNames[4]))
          }
          if (length(vectorNames) > 4) {
            minfi::densityPlot(dat = as.matrix(averageCat5), main = paste0("Density of ", vectorNames[5]))
          }
          if (length(vectorNames) > 5) {
            minfi::densityPlot(dat = as.matrix(averageCat6), main = paste0("Density of ", vectorNames[6]))
          }
          if (length(vectorNames) > 6) {
            minfi::densityPlot(dat = as.matrix(averageCat7), main = paste0("Density of ", vectorNames[7]))
          }
          dev.off()
        }
      }
      
      eachplotNo <- function(AnnotationCategoryToVisualize, 
                              BetaMatrix, 
                              non_empty_entries, 
                              ClinicalFile, 
                              ClinicalCategoryToVisualize) {
        
        # Obtain the unique clinical categories 
        vectorNames <- unique(na.omit(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]))
        
        # all probes and samples only corresponding to those specified in ClinicalCategoryToVisualize
        averageCat11 <- as.matrix(BetaMatrix[ , which(stri_trim(ClinicalFile[c(which(!is.na(ClinicalFile[,which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == TRUE)), 
                                                                              which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == vectorNames[1])])
        averageCat22 <- as.matrix(BetaMatrix[ , which(stri_trim(ClinicalFile[c(which(!is.na(ClinicalFile[,which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == TRUE)), 
                                                                              which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == vectorNames[2])])
        if (length(vectorNames) > 2) {
          averageCat33 <- as.matrix(BetaMatrix[ , which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == TRUE)), 
                                                                                which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == vectorNames[3])])
        }
        if (length(vectorNames) > 3) {
          averageCat44 <- as.matrix(BetaMatrix[correspondingRows, which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == TRUE)), 
                                                                                                which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]) == vectorNames[4])])
        }
        
        if (FigureGenerate == "Yes") {
          # Plotting density distributions
          if (PNGorPDF == "pdf") {
            pdf(paste0(pathNow,"/img/3_DensityBeta_",ClinicalCategoryToVisualize,"_only.pdf"), width = 60, height = 60, pointsize = 50)
          } else {
            png(paste0(pathNow,"/img/3_DensityBeta_",ClinicalCategoryToVisualize,"_only.png"))
          }
          
          if (length(vectorNames) == 1) {
            par(mfrow = c(1, 1))        
          } else if (length(vectorNames) == 2) {
            par(mfrow = c(1, 2))        
          } else if(length(vectorNames) == 3) {
            par(mfrow = c(3, 1))        
          }else {
            par(mfrow = c(2, 2)) 
          }
          minfi::densityPlot(dat = as.matrix(averageCat11), main = paste0("Density of ", vectorNames[1]))
          minfi::densityPlot(dat = as.matrix(averageCat22), main = paste0("Density of ", vectorNames[2]))
          if (length(vectorNames) > 2) {
            minfi::densityPlot(dat = as.matrix(averageCat33), 
                               main = paste0("Density of ", vectorNames[3]))
          }
          if (length(vectorNames) > 3) {
            minfi::densityPlot(dat = as.matrix(averageCat44), main = paste0("Density of ", vectorNames[4]))
          }
          dev.off()
        }
      }
      
      if (PlotWithinAnnotationCategories == "Yes" && FigureGenerate == "Yes") {
        # CategoriesToPlot<-unique(AnnotationFile2[nonEmptyEntries,columnNumber])
        CategoriesToPlot <- unique(na.omit(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategoryToVisualize)]))
        plots <- lapply(CategoriesToPlot, eachplotYes, BetaMatrix = BetaMatrix,
                        non_empty_entries = nonEmptyEntries,
                        ClinicalFile = ClinicalFile,
                        AnnotationCategory = AnnotationCategoryToVisualize,
                        columnNumber = columnNumber)
      } else if (PlotWithinAnnotationCategories == "No" && FigureGenerate == "Yes"){
        plots <- eachplotNo(AnnotationCategoryToVisualize = AnnotationCategoryToVisualize,
                             BetaMatrix = BetaMatrix,
                             non_empty_entries = nonEmptyEntries,
                             ClinicalFile = ClinicalFile,
                             ClinicalCategory = ClinicalCategoryToVisualize)
      }
  }
  
  return(NULL)
}
# [END]
