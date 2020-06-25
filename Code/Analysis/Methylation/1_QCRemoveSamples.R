# Updated 5 Feb 2020
# Updated 8 Nov 2019
# Updated 11 Feb 2019
# Function: Remove samples based on clinical annotations EPIC_DATA, EPIC_QC, and classify 
#           based on TIME_POINT. Returns updated ClinicalFile and sheet.
# Author: Anjali Silva

# Input:
# ClinicalFile: A data frame with observations (e.g., patient IDs) along rows, and categories
#              (e.g., SEX SITE_BIOPSY, etc.) on columns. 
# Path: Provide the path where Illumina methylation sample sheet and IDAT files are located.
#       Default it set to "/Users/anjalisilva/Desktop/UHN/FLOMICS/Git/Methylation/Pipeline/extdata".

# Output: 
# ClinicalFile_updSamples_Ordered_T1: Clinical file updated by removing samples, keeping only T1 samples passing QC
# ClinicalFile_updSamples_Ordered_T2:  Clinical file updated by removing samples, keeping both T1 and T2 samples passing QC
# Sheet_updSamples_Ordered_T1: metharray sheet updated by removing samples, keeping only T1 samples passing QC
# Sheet_updSamples_Ordered_T2: metharray sheet updated by removing samples, keeping both T1 and T2 samples passing QC


QCRemoveSamples <- function(ClinicalFile, 
                            Path = "/Users/anjalisilva/Desktop/UHN/FLOMICS/GitHub/Methylation/Pipeline/extdata") {

  library(data.table)
  library(dplyr)
  library(minfi)
  # Check points
  # cat("\n Dimension of ClinicalFile",dim(ClinicalFile))
  
  ########################################## #
  # Remove samples based on clinical annotations
  ########################################## #
  
  # Check if ClinicalFile has EPIC data
  epic_data <- which(ClinicalFile[, which(colnames(ClinicalFile) == "EPIC_DATA")] == "TRUE") 
  
  # Further analysis is only performed if EPIC data is present in Clinical File
  if (length(epic_data) > 0) {
    
    # ClinicalFile with QC passed T1 timepoints only
    ClinicalFile_RmdSamples_T1 <- ClinicalFile %>%
      filter(EPIC_QC != "Bad" & EPIC_INCLUDE == "YES" & TIME_POINT == "T1")
    
    # ClinicalFile with QC passed T1 and T2 timepoints 
    ClinicalFile_RmdSamples_T2 <- ClinicalFile %>%
      filter(EPIC_QC != "Bad" & EPIC_INCLUDE == "YES")
    
    # Check column for EPIC_QC data
    quality <- which(ClinicalFile[, which(colnames(ClinicalFile) == "EPIC_QC")] == "Bad")
    
    # Check column for INCLUDE data  
    include <- which(ClinicalFile[, which(colnames(ClinicalFile) == "EPIC_INCLUDE")] == "NO")
    
  } else {
    # If no EPIC data present
    stop("\n No EPIC_DATA present in ClinicalFile provided. \n");
  }
  
  # If samples to be removed present based on EPIC_QC, EPIC_INCLUDE from Clinical File
  message("\nBased on 'EPIC_DATA = TRUE', EPIC_QC == BAD', 'EPIC_INCLUDE == NO', \n the following ", 
          length(unique(c(quality, include))) ," sample(s) is (are) removed from Clinical File: 
          \n", as.character(ClinicalFile$SAMPLE_ID[unique(c(quality, include))]), ".\n")

  
  ########################################## #
  # Defining the path to reading the an Illumina methylation sample sheet
  ########################################## #
  
  # sheet <- read.metharray.sheet(base =file.path(SampleSheet_Path, "extdata"))
  sheet <- minfi::read.metharray.sheet(base = Path)
  # class(RGset) # "RGChannelSet"
  # RGChannelSet: raw data from the IDAT files; this data is organized at the probe 
  # (not CpG locus) level. This data has two channels: Red and Green.
  
  # Adjuste sample sheet sample names
  sheet$Sample_Name <- gsub("KRI_", "", sheet$Sample_Name)
  sheet$Sample_Name[which(substr(sheet$Sample_Name, 10, 12) == "")] <-
    paste0(sheet$Sample_Name[which(substr(sheet$Sample_Name, 10, 12) == "")], "_T1")
  
  # Saving original sample sheet
  sheet_original <- sheet
  # Manually remove "LY_FL_159_T1" and keep only "LY_FL_159_T1_rep"
  sheet <- sheet[- which(sheet$Sample_Name == "LY_FL_159_T1"), ]
  sheet$Sample_Name[which(sheet$Sample_Name == "LY_FL_159_T1_rep")] <- "LY_FL_159_T1"
  
  # Order sample sheet sample names by ClinicalFile_RmdSamples_T1
  OrderSamples_T1 <- match(ClinicalFile_RmdSamples_T1$SAMPLE_ID, sheet$Sample_Name)
  
  # sheet$Sample_Name[OrderSamples[!is.na(OrderSamples)]]
  Sheet_updSamples_Ordered_T1 <- sheet[OrderSamples_T1[! is.na(OrderSamples_T1)], ]
  
  # Order sample sheet sample names by ClinicalFile_RmdSamples_T1
  OrderSamples_T2 <- match(ClinicalFile_RmdSamples_T2$SAMPLE_ID, sheet$Sample_Name)
  # sheet$Sample_Name[OrderSamples[!is.na(OrderSamples)]]
  Sheet_updSamples_Ordered_T2 <- sheet[OrderSamples_T2[! is.na(OrderSamples_T2)], ]
  
  # Generate a clinical file with all samples
  ClinicalFile_AllSamples <- ClinicalFile[match(sheet$Sample_Name, ClinicalFile$SAMPLE_ID), ]

  RESULTS <- list(ClinicalFile_updSamples_170_Ordered_T1 = ClinicalFile_RmdSamples_T1,
                  ClinicalFile_updSamples_170_Ordered_T2 = ClinicalFile_RmdSamples_T2,
                  ClinicalFile_updSamples_176 = ClinicalFile_AllSamples,
                  Sheet_updSamples_170_Ordered_T1 = Sheet_updSamples_Ordered_T1,
                  Sheet_updSamples_170_Ordered_T2 = Sheet_updSamples_Ordered_T2,
                  Sheet_AllSamples_177 = sheet_original)
  
  class(RESULTS) <- "SamplesToRemove_ASilva"
  return(RESULTS)
}

