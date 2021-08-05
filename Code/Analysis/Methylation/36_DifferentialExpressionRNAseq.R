# Updated 3 Aug 2021
# 3 March 2020
# Function: Perform differential expression analysis given RNAseq count matrix and contrast. 
# Ref: https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# Author: Anjali Silva and Robert Kridel

# Input:
# RNAseqCountMatrix: A matrix of RNAseq counts of ENSEMBL IDs x samples that has been NOT 
#                    been normalized but features filtered based on users choice. 
# ContrastColumnName: Type of differential analysis to be done. ContrastColumnName options
#                     include "STAGE", "SEX", "SITE_BIOPSY", "TYPE_BIOPSY", "INSTITUTION",
#                     "COO", "TYPE", "TRANSLOC_14_18", "EPIC_QC", "CLUSTER". Default: "TYPE".
# AnnotationFileEnsemblToGenes: A data frame containing Ensemble ID in first column and gene  
#                               name in the next column, including type which lists if a 
#                               protein-coding gene or not.  
# ClinicalFile: File with patient sample names, and categories. 
# ClusterLabels: If "ContrastColumnName" == "CLUSTER", vector of integers with length equal 
#                to ncol(RNAseqCountMatrix) should be provided. Default is NA. 
# LogFCcutoff: An integer indicating the log fold change cut off to be used. All entries with a 
#              log fold change cut off of greater than LogFCcutoff is selected for further 
#              analysis. Default: 2. 
# FDRcutoff: A numeric value indicating the FDR cut off. All entries with a FDR value less than
#            FDRcutoff will be selected. Default: 0.05. 
# ProduceImages: Produce images or not, options = "Yes" or "No"
# ProduceHeatMap: A character vector specifying "Yes" or "No", as to if heatmap should be produced.
#                 Only available if ClusterLabels are provided. efault is "No". 
# PNGorPDF: Output format of the image, options = "png" or "pdf"


# Output:
# TopTagsEachContrast: A list of length(contrasts), with all tags for each contrast obtained using
#                      edgeR::topTags with Benjamini & Hochberg adjustment method for p-values for 
#                      multiple testing. IDs sorted by p-values, smallest at top. Contains logFC, 
#                      logCPM, F, PValue, FDR, ensemblGeneId for each ID. 
# AllDEIDsTable: A list of length(contrasts), with table for each contrast that lists only the 
#                differentially expressed IDs that are SIGNIFICANT (with a FDR of 0.05 and logFC > 2). 
#                Contains logFC, logCPM, F, PValue, FDR, ensemblGeneId for each ID. IDs 
#                sorted by p-values, smallest at top.
# EachContrastPosSignificant: A list of length(contrasts), with table for each contrast that lists only the 
#                            significant tags (FDR of 0.05 and logFC > 2) with (+) fold change. IDs 
#                            sorted by p-values, smallest at top.
# EachContrastNegSignificant: A list of length(contrasts), with table for each contrast that lists only the 
#                            significant tags (FDR of 0.05 and logFC > 2) with (-) fold change. IDs 
#                            sorted by p-values, smallest at top.

# AllDEIDs: A list containing unique ENSEMBL IDs from all different contrasts. E.g., if contrast is 
#           TYPE, differentially expressed ENSEMBL IDs from FL-DLBCL, FL-RLN, RLN-DLBCL will be
#           combined together in this list. 
# RawCountsRNAseqMatrixDEgenes: A matrix of size differentially expressed ENSEMBL IDs x samples
#                               containing the raw counts. 
# NormalizedRNAseqMatrixDEgenes: A matrix of size differentially expressed ENSEMBL IDs x samples
#                               containing the normalized counts. 


# Images
# 36_logFCsVSaverageCountSize.p*
# 36_EnhancedVolcano_Contrast_", colnames(contrasts)[i], ".p*
# 36_VolcanoPlot_GGplot_Contrast_", colnames(contrasts)[i], ".p
# 36_Heatmap.p*


DifferentialExpressionRNAseq <- function(RNAseqCountMatrix, 
                                         ContrastColumnName = "TYPE", 
                                         AnnotationFileEnsemblToGenes = NA,
                                         ClinicalFile, 
                                         ClusterLabels = NA, 
                                         LogFCcutoff = 2, 
                                         FDRcutoff = 0.05, 
                                         ProduceImages = "Yes", 
                                         ProduceHeatMap = "No",
                                         PNGorPDF = "png") {
  
  # Loading needed packages
  library(dplyr)
  library(edgeR)
  library(EnhancedVolcano)
  library(ggpubr)
  library(grDevices)
  library(KEGGprofile)
  library(limma)
  library(RColorBrewer)
  library(stringi)
  library(pheatmap)
  
  pathNow <- getwd() # getting the path
  
  # Remove empty entries in the clinical file
  ClinicalFile[ClinicalFile == ""] <- NA
  
  # check user input
  if(ncol(RNAseqCountMatrix) != nrow(ClinicalFile)) {
    stop("The number of samples in RNAseqCountMatrix should match those in ClinicalFile.")
  }
  
  # if cluster labels are provided check the length
  if((all(is.na(ClusterLabels)) == FALSE)  && (ncol(RNAseqCountMatrix) != length(ClusterLabels))) {
    stop("The number of samples in RNAseqCountMatrix should match length of ClusterLabels.")
  }
  
  # AnnotationFileEnsemblToGenes is needed to select protein-coding genes
  if(all(is.na(AnnotationFileEnsemblToGenes)) == TRUE) {
    stop("\n AnnotationFileEnsemblToGenes must be provided since only analyzing protein-coding genes.")
  }
  
  
  # Indicate cut-offs for logFC and FDR
  cat("\n LogFCcutoff is set to > ", LogFCcutoff, " and FDRcutoff is set to < ", FDRcutoff, ".")

  
  
  # Keep only protein coding ENSEMBLIDs for now

  AnnotationFileENSEMBLIDsProteinCoding <- AnnotationFileEnsemblToGenes %>% 
    dplyr::filter(type == "protein_coding")
  matchRNAseqAnnotation <- match(AnnotationFileENSEMBLIDsProteinCoding$gene, 
                                 rownames(RNAseqCountMatrix))
  RNAseqCountMatrix <- RNAseqCountMatrix[matchRNAseqAnnotation[! is.na(matchRNAseqAnnotation)] , ]                               
  cat("\n Only analyze protein coding ENSEMBLIDs. Therefore, ", 
      nrow(RNAseqCountMatrix), "ENSEMBL IDs analyzed.")
  
  
  #  create a DGEList object:
  exprsRNAseq <- edgeR::DGEList(counts = RNAseqCountMatrix)
  exprsRNAseqNormFacs <- edgeR::calcNormFactors(exprsRNAseq)
  
  
  # Set contrasts 
  if(ContrastColumnName == "STAGE") {
    # check if NA values are present, if so remove those patients 
    statusContrastName <- factor(stringi::stri_trim( ClinicalFile[which(ClinicalFile[ ,
                                     which(colnames(ClinicalFile) == "TYPE")] == "FL"), 
                                     which(colnames(ClinicalFile) == "STAGE")]), 
                                     levels = c("ADVANCED", "LIMITED"))
    exprsRNAseqNormFacs <- exprsRNAseqNormFacs[, which(ClinicalFile[ , 
                                     which(colnames(ClinicalFile) == "TYPE")] == "FL")]
    designRNAseq <- stats::model.matrix (~0 + statusContrastName)
    colnames(designRNAseq) <- levels(statusContrastName)
    contrasts <- limma::makeContrasts("ADVANCED-LIMITED", levels = designRNAseq)
  } else if(ContrastColumnName == "SEX") {
    if(length(which(is.na(ClinicalFile[ , which(colnames(ClinicalFile) == "SEX")]) == TRUE)) > 0) {
      statusContrastName <- factor(stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[, 
                                                which(colnames(ClinicalFile) == "SEX")]) == TRUE)), 
                                                          which(colnames(ClinicalFile) == "SEX")]),  
                                   levels = c("M", "F"), labels = c("Male", "Female"))
      exprsRNAseqNormFacs <- exprsRNAseqNormFacs[, - c(which(is.na(ClinicalFile[, 
                                                which(colnames(ClinicalFile) == "SEX")]) == TRUE))]
      
    } else {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[ , 
                                   which(colnames(ClinicalFile) == "SEX")]),  levels = c("M", "F"), 
                                   labels = c("Male", "Female"))
    }
    designRNAseq <- stats::model.matrix(~0 + statusContrastName)
    colnames(designRNAseq) <- levels(statusContrastName)
    contrasts <- limma::makeContrasts("Female-Male", levels = designRNAseq)
  } else if(ContrastColumnName == "SITE_BIOPSY") {
    if(length(which(is.na(ClinicalFile[ , which(colnames(ClinicalFile) == "SITE_BIOPSY")]) == TRUE)) > 0) {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[ , 
                                            which(colnames(ClinicalFile) == "SITE_BIOPSY")]) == TRUE)), 
                                            which(colnames(ClinicalFile) == "SITE_BIOPSY")]),  
                                            levels = c("EN", "LN"))
      exprsRNAseqNormFacs <- exprsRNAseqNormFacs[, - c(which(is.na(ClinicalFile[ , 
                                            which(colnames(ClinicalFile) == "SITE_BIOPSY")]) == TRUE))]
    }else {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[, 
                                            which(colnames(ClinicalFile) == "SITE_BIOPSY")] == TRUE))), 
                                            which(colnames(ClinicalFile) == "SITE_BIOPSY")]), 
                                            levels = c("EN", "LN"))
    }
    designRNAseq <- stats::model.matrix(~0 + statusContrastName)
    colnames(designRNAseq) <- levels(statusContrastName)
    contrasts <- limma::makeContrasts("LN-EN", levels = designRNAseq)
  } else if(ContrastColumnName == "TYPE_BIOPSY") {
    if(length(which(is.na(ClinicalFile[, which(colnames(ClinicalFile) == "TYPE_BIOPSY")]) == TRUE)) > 0) {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[, 
                                            which(colnames(ClinicalFile) == "TYPE_BIOPSY")]) == TRUE)), 
                                            which(colnames(ClinicalFile) == "TYPE_BIOPSY")]), 
                                            levels = c("CORE", "TISSUE"))
      exprsRNAseqNormFacs <- exprsRNAseqNormFacs[, - c(which(is.na(ClinicalFile[ , 
                                            which(colnames(ClinicalFile) == "TYPE_BIOPSY")]) == TRUE))]
      
    } else {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[, 
                                                      which(colnames(ClinicalFile) == "TYPE_BIOPSY")] == TRUE))), 
                                                      which(colnames(ClinicalFile) == "TYPE_BIOPSY")]), 
                                                      levels = c("CORE", "TISSUE"))
    }
    designRNAseq <- stats::model.matrix(~0 + statusContrastName)
    colnames(designRNAseq) <- levels(statusContrastName)
    contrasts <- limma::makeContrasts("TISSUE-CORE", levels = designRNAseq)
  } else if(ContrastColumnName == "INSTITUTION") {
    # https://support.bioconductor.org/p/44216/
    if(length(which(is.na(ClinicalFile[,which(colnames(ClinicalFile) == "INSTITUTION")]) == TRUE)) > 0) {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[ , 
                                                      which(colnames(ClinicalFile) == "INSTITUTION")]) == TRUE)), 
                                                      which(colnames(ClinicalFile) == "INSTITUTION")]), 
                                                      levels = c("BCCA", "UHN","JGH"))
      exprsRNAseqNormFacs <- exprsRNAseqNormFacs[, - c(which(is.na(ClinicalFile[ , 
                                                      which(colnames(ClinicalFile) == "INSTITUTION")]) == TRUE))]
    } else {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[ , 
                                                    which(colnames(ClinicalFile) == "INSTITUTION")]),
                                                    levels = c("BCCA", "UHN","JGH"))
    }
    designRNAseq <- stats::model.matrix(~0 + statusContrastName)
    colnames(designRNAseq) <- levels(statusContrastName)
    contrasts <- limma::makeContrasts("UHN-BCCA", "JGH-BCCA", "JGH-UHN", levels = designRNAseq)
  } else if(ContrastColumnName == "EPIC_QC") {
    if(length(which(is.na(ClinicalFile[,which(colnames(ClinicalFile) == "EPIC_QC")]) == TRUE)) > 0) {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[, 
                                                      which(colnames(ClinicalFile) == "EPIC_QC")]) == TRUE)), 
                                                      which(colnames(ClinicalFile) == "EPIC_QC")]), 
                                                      levels = c("Fair", "Good", "Poor"))
      exprsRNAseqNormFacs <- exprsRNAseqNormFacs[, - c(which(is.na(ClinicalFile[ , 
                                                 which(colnames(ClinicalFile) == "EPIC_QC")]) == TRUE))]
      
    } else {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[ , which(colnames(ClinicalFile) == "EPIC_QC")]), 
                                   levels = c("Fair", "Good", "Poor"))
    }
    designRNAseq <- stats::model.matrix(~0 + statusContrastName)
    colnames(designRNAseq) <- levels(statusContrastName)
    contrasts <- limma::makeContrasts("Good-Fair", "Poor-Good", "Poor-Fair", levels = designRNAseq)
  } else if(ContrastColumnName == "COO") {
    # https://support.bioconductor.org/p/44216/ # contrasts for multiple comparisons
    if(length(which(is.na(ClinicalFile[, which(colnames(ClinicalFile) == "COO")]) == TRUE)) > 0) {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[, 
                                                      which(colnames(ClinicalFile) == "COO")]) == TRUE)), 
                                                      which(colnames(ClinicalFile) == "COO")]), 
                                                      levels = c("GCB", "ABC"))
      exprsRNAseqNormFacs <- exprsRNAseqNormFacs[, - c(which(is.na(ClinicalFile[ , 
                                                 which(colnames(ClinicalFile) == "COO")]) == TRUE))]
     } else {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[which(ClinicalFile[ , 
                                                      which(colnames(ClinicalFile) == "TYPE")] == "DLBCL"), 
                                                      which(colnames(ClinicalFile) == "COO")]),
                                                      levels = c("GCB", "ABC"))
    }
    designRNAseq <- stats::model.matrix(~0 + statusContrastName)
    colnames(designRNAseq) <- levels(statusContrastName)
    contrasts <- limma::makeContrasts("ABC-GCB", levels = designRNAseq)
  } else if(ContrastColumnName == "TYPE") {
    if(length(which(is.na(ClinicalFile[, which(colnames(ClinicalFile) == "TYPE")]) == TRUE)) > 0) {
      if (length(unique(ClinicalFile[, which(colnames(ClinicalFile) == "TYPE")])) == 2) {
        statusContrastName <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[, 
                                                        which(colnames(ClinicalFile) == "TYPE")]) == TRUE)), 
                                                        which(colnames(ClinicalFile) == "TYPE")]), 
                                                        levels = c("DLBCL", "FL","RLN"))
        exprsRNAseqNormFacs <- exprsRNAseqNormFacs[, - c(which(is.na(ClinicalFile[ , 
                                                   which(colnames(ClinicalFile) == "TYPE")]) == TRUE))]
      } else {
        statusContrastName <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[ , 
                                                        which(colnames(ClinicalFile) == "TYPE")]) == TRUE)), 
                                                        which(colnames(ClinicalFile) == "TYPE")]), 
                                                        levels = c("DLBCL", "FL", "RLN"))}
    } else {
      if (length(unique(ClinicalFile[,which(colnames(ClinicalFile) == "TYPE")])) == 2) {
        statusContrastName <- factor(stringi::stri_trim(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE")]))
      } else {
        statusContrastName <- factor(stringi::stri_trim(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE")]), 
                                     levels = c("DLBCL", "FL", "RLN"))
      }
    }
    designRNAseq <- stats::model.matrix(~0 + statusContrastName)
    colnames(designRNAseq) <- levels(statusContrastName)
    if (length(unique(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE")])) == 2) {
      contrasts <- limma::makeContrasts(paste0(unique(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE")])[1], "-", 
                                        unique(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE")])[2]), 
                                        levels = designRNAseq)
    } else {
      contrasts <- limma::makeContrasts("FL-RLN", "DLBCL-RLN", "DLBCL-FL", levels = designRNAseq)
    }
  } else if (ContrastColumnName == "TRANSLOC_14_18") {
    if (length(which(is.na(ClinicalFile[ , which(colnames(ClinicalFile) == "TRANSLOC_14_18")]) == TRUE)) > 0) {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[ , 
                                                      which(colnames(ClinicalFile) == "TRANSLOC_14_18")]) == TRUE)),
                                                      which(colnames(ClinicalFile) == "TRANSLOC_14_18")]),  
                                                      levels = c("0", "1"), 
                                                      labels = c("NoTranslocation", "Translocation"))
      exprsRNAseqNormFacs <- exprsRNAseqNormFacs[, - c(which(is.na(ClinicalFile[ , 
                                                which(colnames(ClinicalFile) == "TRANSLOC_14_18")]) == TRUE))]
    } else {
      statusContrastName <- factor(stringi::stri_trim(ClinicalFile[which(! is.na(ClinicalFile[ , 
                                                      which(colnames(ClinicalFile) == "TRANSLOC_14_18")]) == TRUE),
                                                      which(colnames(ClinicalFile) == "TRANSLOC_14_18")]), 
                                                      levels = c("0", "1"), 
                                                      labels = c("NoTranslocation", "Translocation"))
    }
    designRNAseq <- stats::model.matrix(~0 + statusContrastName)
    colnames(designRNAseq) <- levels(statusContrastName)
    contrasts <- limma::makeContrasts("Translocation-NoTranslocation", levels = designRNAseq)
  } else if(ContrastColumnName == "CLUSTER") {
    if (sum(is.na(ClusterLabels)) > 0) {
      stop("ContrastColumnName is set to CLUSTER, but ClusterLabels not provided. Provide the ClusterLabels");}
    if (length(unique(ClusterLabels)) == 4) {
      statusContrastName <- factor(ClusterLabels, levels = c("1", "2", "3", "4"), 
                                   labels = c("one", "two", "three", "four"))
    } else if(length(unique(ClusterLabels)) == 3) {
      statusContrastName <- factor(ClusterLabels, levels = c("1", "2", "3"), 
                                   labels = c("one", "two", "three"))
    } else if(length(unique(ClusterLabels)) == 2) {
      statusContrastName <- factor(ClusterLabels, levels = c("1", "2"), 
                                   labels = c("one", "two"))
    } else{
      stop("ClusterLabels have more than 4 clusters. Currently only upto 4 clusters supported");}

    designRNAseq <- stats::model.matrix(~0 + statusContrastName)
    colnames(designRNAseq) <- levels(statusContrastName)
    if (length(unique(ClusterLabels)) == 4) {
      contrasts <- limma::makeContrasts("one-two", "one-three", "one-four", "two-three", 
                                        "two-four", "three-four", 
                                        levels = designRNAseq)
    } else if (length(unique(ClusterLabels)) == 3) {
      contrasts <- limma::makeContrasts("one-two","one-three","two-three", levels = designRNAseq)
    } else if(length(unique(ClusterLabels)) == 2) {
      contrasts <- limma::makeContrasts("one-two", levels = designRNAseq)
    }
    
  } else {
    stop("Not a valid ContrastColumnName. Check doumentation for valid ContrastColumnNames.")
  }
  
  
  
  
  # Estimate dispersion.
  # Estimate the genewise dispersion estimates over all genes, allowing for a possible 
  # abundance trend.
  # Allowing gene-specific dispersion is necessary in order that differential
  # expression is not driven by outliers. Therefore the tagwise dispersions are
  # strongly recommended in model fitting and testing for differential expression.
  # This method estimates common dispersion, trended dispersions and tagwise dispersions
  # in one run and is recommended.
  dispEstimates <- edgeR::estimateDisp(y = exprsRNAseqNormFacs, 
                                       design = designRNAseq, 
                                       robust = TRUE)
  
  # The dispersion estimates can be viewed in a BCV plot
  # Plot the genewise biological coefficient of variation (BCV) against gene abundance 
  # (in log2 counts per million).
  # edgeR::plotBCV(dispEstimates)
  
  # quasi-likelihood (QL) dispersions can be estimated and visualized 
  dispEstimateQL <- edgeR::glmQLFit(y = dispEstimates, 
                                    design = designRNAseq, 
                                    robust = TRUE)
  # plotQLDisp(dispEstimateQL)
  
  # Fit generalized linear model
  # Such a model is an extension of classical linear models to non-normally distributed 
  # response data.
  fitGenLModel <- edgeR::glmQLFit(y = dispEstimates, 
                                  design = designRNAseq)
  
  # Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests
  
  
  
  
  # Differential gene expression analysis
  if(ProduceImages == "Yes") {
    if (PNGorPDF == "pdf") {
      pdf(paste0(pathNow,"/img/36_logFCsVSaverageCountSize.pdf"),
          width = 60, height = 60, pointsize = 50)
    } else {
      png(paste0(pathNow,"/img/36_logFCsVSaverageCountSize.png"))
    }
    
    if (ncol(contrasts) == 1) {
      par(mfrow = c(1, 1))
    } else if (ncol(contrasts) == 2) {
      par(mfrow = c(1, 2))
    } else if (ncol(contrasts) <= 4) {
      par(mfrow = c(2, 2))
    } else if (ncol(contrasts) <= 6) {
      par(mfrow = c(2, 3))
    }
  }
  
  
  # Fit model and plot log CPM against logFC
  topTagsPvalue <- filterFCandFDR <- IDsavelogFCandFDRup <- IDsavelogFCandFDRdown <-
    IDsavelogFCandFDRupGeneList <- IDsavelogFCandFDRdownGeneList <- list()
  for (i in 1:ncol(contrasts)) {
    # Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
    # Conduct genewise statistical tests for a given coefficient or contrast.
    fitGenLModQLTest <- edgeR::glmQLFTest(glmfit = fitGenLModel, 
                                          contrast = contrasts[, i])
    
    if(ProduceImages == "Yes") {
      # Plot all the logFCs against average count size, highlighting the DE genes
      limma::plotMD(fitGenLModQLTest, main = paste0("Contrast: ", colnames(contrasts)[i]))
      graphics::abline(h = c(-1, 1), col = "blue")
    }
    
    # Extracts the most differentially expressed genes ranked by p-value 
    # Issue present with duplicated issues. Not sure why.
    topTagsPvalue[[i]] <- data.frame(edgeR::topTags(fitGenLModQLTest, 
                                                    n = Inf, 
                                                    adjust.method = "BH", 
                                                    sort.by = "PValue")) %>% 
                                                    mutate(ensemblGeneId = row.names(.))
    
    
    # AllDifferentilalENSEMBLids[67]
    # Filter based on prespecified logFC and FDR cutoff values. 
    # Note for logFC, the absolute value is taken.
    filterFCandFDR[[i]] <- topTagsPvalue[[i]] %>% dplyr::filter(abs(logFC) > LogFCcutoff & FDR < FDRcutoff) 
    cat("\n Contrast:", colnames(contrasts)[i], "has total of", nrow(filterFCandFDR[[i]]), 
        "ENSEMBL IDs that passed logFC > ", LogFCcutoff, "and FDR < ", FDRcutoff, ".")
    
    # Seperate filtered values based on + or - logFC
    IDsavelogFCandFDRup[[i]] <- filterFCandFDR[[i]] %>% dplyr::filter(logFC > 0) 
    IDsavelogFCandFDRupGeneList[[i]] <- IDsavelogFCandFDRup[[i]]$ensemblGeneId
    cat("\n Contrast:", colnames(contrasts)[i], "has", nrow(IDsavelogFCandFDRup[[i]]), 
        "ensembl IDs (+ve)logFC and significant.")
    IDsavelogFCandFDRdown[[i]] <- filterFCandFDR[[i]] %>% dplyr::filter(logFC < 0) 
    IDsavelogFCandFDRdownGeneList[[i]] <- IDsavelogFCandFDRdown[[i]]$ensemblGeneId
    cat("\n Contrast:", colnames(contrasts)[i], "has", nrow(IDsavelogFCandFDRdown[[i]]), 
        "ensembl IDs (-ve)logFC and significant.")
  }
  
  if(ProduceImages == "Yes") {
    grDevices::dev.off()
  }
  
  # name list with contrast names
  names(topTagsPvalue) <- names(filterFCandFDR) <- names(IDsavelogFCandFDRup) <- 
    names(IDsavelogFCandFDRdown) <- colnames(contrasts)
  
  if(ProduceImages == "Yes") {
    # Enhanced volcano plot
    for (i in 1:ncol(contrasts)) {
      par(mfrow = c(1, 1))
      
      # https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
      EnhancedVolcano::EnhancedVolcano(toptable = topTagsPvalue[[i]],
                                       lab = topTagsPvalue[[i]]$ensemblGeneId,
                                       x = 'logFC',
                                       y = 'PValue',
                                       pCutoff = c(0.01, 0.05),
                                       FCcutoff = 2.0,
                                       title = paste0("VolcanoPlot Contrast:", colnames(contrasts)[i]),
                                       xlim = c(-15, 15),
                                       gridlines.major = FALSE,
                                       gridlines.minor = FALSE)
      
      if (PNGorPDF == "pdf") {
        ggplot2::ggsave(paste0(pathNow, "/img/36_EnhancedVolcano_Contrast_", colnames(contrasts)[i], ".pdf"))
      } else {
        ggplot2::ggsave(paste0(pathNow, "/img/36_EnhancedVolcano_Contrast_", colnames(contrasts)[i], ".png"))
      }
    } 
    
    # GGplot2 volcano plot
    for (i in 1:ncol(contrasts)) {
      grDevices::dev.off()
      par(mfrow = c(1, 1))
      VolcanoPlot <- ggplot2::ggplot(topTagsPvalue[[i]]) +
                                    geom_point(aes(x = logFC, y= -log10(PValue))) +
                                    ggtitle(paste0("Contrast: ", colnames(contrasts)[i])) +
                                    labs(y = "-log10(adjusted P.Value)", x = "logFC") +
                                    geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
                                               col = "#4d9221", size = 1) +
                                    geom_hline(yintercept = -log10(0.001), linetype = "dashed", 
                                               col = "#762a83", size = 1) +
                                    theme_bw()  +
                                    theme(text = element_text(size=15), 
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(), 
                                          panel.background = element_rect(colour = "black", size = 1.5), 
                                          axis.title =  element_text(face = "bold"),
                                          aspect.ratio = 1, legend.position = "none", 
                                          axis.text.x = element_text(face = "bold"), 
                                          axis.text.y = element_text(face = "bold")) 
      
      ggplot2::ggsave(paste0(pathNow, "/img/36_VolcanoPlot_GGplot_Contrast_", colnames(contrasts)[i], ".png"))
    }
  }
  
  
  
  
  
  EachContrastAllDifferentilalENSEMBLids <- 
    EachContrastRNAseqCountMatrixAllDEgenes <- 
    EachContrastData2RNAseqFilterededgeRNorm <- list()
  # plot heatmpa - may take a long time depending on number of IDs; only done if ClusterLabels provided
  if(ProduceHeatMap == "Yes" && (all(is.na(ClusterLabels)) == FALSE)) {
    for (i in 1:ncol(contrasts)) {
      
      
      UniqueDifferentilalUpENSEMBLids <- unique(unlist(IDsavelogFCandFDRupGeneList[[i]])) 
      UniqueDifferentilalDownENSEMBLids <- unique(unlist(IDsavelogFCandFDRdownGeneList[[i]])) 
      EachContrastAllDifferentilalENSEMBLids[[i]] <- AllDifferentilalENSEMBLids <- 
        unique(c(UniqueDifferentilalUpENSEMBLids, UniqueDifferentilalDownENSEMBLids)) 
      
      # Issue with duplicated issues beyond second decimal
      # substr(AllDifferentilalENSEMBLids[which(substr(AllDifferentilalENSEMBLids,18, 18) == ".")], 1, 17)
      
      # get count matrix for DE ids
      EachContrastRNAseqCountMatrixAllDEgenes[[i]] <- RNAseqCountMatrixAllDEgenes <- 
        RNAseqCountMatrix[match(AllDifferentilalENSEMBLids, rownames(RNAseqCountMatrix)), ] 
      

    # order patients by cluster 
    orderColumns <- c(which(ClusterLabels == 1),
                      which(ClusterLabels == 2),
                      which(ClusterLabels == 3),
                      which(ClusterLabels == 4),
                      which(ClusterLabels == 5))
    
    annotation_col = data.frame(
      Cluster = factor(ClusterLabels[orderColumns]),
      Disease = factor(ClinicalFile$TYPE[orderColumns]),
      Stage = factor(ClinicalFile$STAGE[orderColumns]),
      Sex = factor(ClinicalFile$SEX[orderColumns]),
      Translocation = factor(ClinicalFile$TRANSLOC_14_18[orderColumns]),
      TypeBiopsy =  ClinicalFile$TYPE_BIOPSY[orderColumns],
      SiteBiopsy =  ClinicalFile$SITE_BIOPSY[orderColumns])
    rownames(annotation_col) = colnames(RNAseqCountMatrixAllDEgenes[, orderColumns])
    
    ann_colors = list(
      Cluster = c("1" = "#4363d8", "2" = "#f58231"), 
      Disease = c(DLBCL = "#d6604d", FL = "#66bd63", RLN = "#4575b4"), #option 5
      Stage = c(ADVANCED = "#762a83", LIMITED = "#c2a5cf"),
      Sex = c("F" = "#b35806", "M" = "#fdb863"),
      Translocation = c("0" = "#e0e0e0","1" = "#878787"),
      TypeBiopsy = c("TISSUE" = "#a6dba0", "CORE" = "#878787"),
      SiteBiopsy = c("LN" = "#f1b6da" , "EN" = "#c51b7d"))
    
    mypalette <- RColorBrewer::brewer.pal(11, "RdYlBu")
    morecols <- grDevices::colorRampPalette(mypalette)
    
    grDevices::dev.off()
    cat("\n Printing HeatmapRawCounts for contrast", colnames(contrasts)[i])
    par(mfrow = c(1, 1))
    if (PNGorPDF == "pdf") {
      pdf(paste0(pathNow, "/img/36_HeatmapRawCounts", nrow(RNAseqCountMatrixAllDEgenes), 
                 "DEgenes_contrast", colnames(contrasts)[i], ".pdf"),
          width = 60, height = 60, pointsize = 50)
    } else {
      png(paste0(pathNow, "/img/36_HeatmapRawCounts", nrow(RNAseqCountMatrixAllDEgenes), 
                 "DEgenes_contrast", colnames(contrasts)[i], ".png"))
    }
    
    # Plot the heatmap
    pheatmap::pheatmap((log(RNAseqCountMatrixAllDEgenes + 0.01))[, orderColumns],
                       show_colnames = T, 
                       show_rownames = F,
                       fontface = "italic", 
                       legend = T,
                       annotation_colors = ann_colors, 
                       border_color = "black", 
                       cluster_cols = FALSE,
                       cluster_rows = TRUE,
                       annotation_col = annotation_col, 
                       color =  rev(morecols(50)))
    
    grDevices::dev.off()
    
    

    
    
    # calculate normalized values for DE IDs
    exprsRNAseq <- edgeR::DGEList(counts = RNAseqCountMatrixAllDEgenes, 
                                  group = factor(ClusterLabels))
    exprsRNAseq$counts <- exprsRNAseq$counts
    Data2RNAseqFilterededgeRNormFactors <- edgeR::calcNormFactors(exprsRNAseq)
    EachContrastData2RNAseqFilterededgeRNorm [[i]] <- Data2RNAseqFilterededgeRNorm <- 
                                                        edgeR::cpm(Data2RNAseqFilterededgeRNormFactors, 
                                                        normalized.lib.sizes = TRUE, log = TRUE)
    # dim(Data2RNAseqFilterededgeRNorm) # IDs x samples
    par(mfrow = c(1, 1))
    cat("\n Printing HeatmapNormalizedCounts for contrast", colnames(contrasts)[i])
    if (PNGorPDF == "pdf") {
      pdf(paste0(pathNow, "/img/36_HeatmapNormalizedCounts", nrow(RNAseqCountMatrixAllDEgenes), 
                 "DEgenes_contrast", colnames(contrasts)[i], ".pdf"),
          width = 60, height = 60, pointsize = 50)
    } else {
      png(paste0(pathNow, "/img/36_HeatmapNormalizedCounts", nrow(RNAseqCountMatrixAllDEgenes),  
                 "DEgenes_contrast", colnames(contrasts)[i], ".png"))
    }
    # Plot the heatmap
    pheatmap::pheatmap(Data2RNAseqFilterededgeRNorm[, orderColumns],
                       #scale = "column",
                       show_colnames = T, 
                       show_rownames = F,
                       fontface = "italic", 
                       legend = T,
                       annotation_colors = ann_colors, 
                       border_color = "black", 
                       cluster_cols = FALSE,
                       cluster_rows = TRUE,
                       annotation_col = annotation_col, 
                       color =  rev(morecols(50)))
    
    grDevices::dev.off()
    
    
    }

  }
  
  # Do pathway analysis if AnnotationFileEnsemblToGenes provided 
  #if (! all(is.na(AnnotationFileEnsemblToGenes))) {
    # AnnotationFileEnsemblToGenes <- as.data.frame(AnnotationFileEnsemblToGenes)
    # MatchedGeneNames <- AnnotationFileEnsemblToGenes[match(sub("\\..*", "", AllDifferentilalENSEMBLids), 
    #                                         AnnotationFileEnsemblToGenes$gene), ]$name
    # cbind(AllDifferentilalENSEMBLids, AllDifferentilalENSEMBLids)
    
    # gsub("\\.(?=[^.]*\\.)", "", "ENSG00000197889.5", perl=TRUE)
    # gsub("\\.(?=[^.]*\\.)", "", "ENSG00000205755.6.1", perl=TRUE)
    # gsub("\\.(?=[.$]*\\.)", "", "ENSG00000205755.6.1", perl=TRUE)
    
  #}
  
  notRun <- function() {
    go <- goana(fitGenLModel, species = "Hs")
  }
  
  
  RESULTS <- list(TopTagsEachContrast = topTagsPvalue,
                  AllDEIDsTable = filterFCandFDR,
                  EachContrastPosSignificant = IDsavelogFCandFDRup,
                  EachContrastNegSignificant = IDsavelogFCandFDRdown,
                  AllDEIDs = EachContrastAllDifferentilalENSEMBLids,
                  RawCountsRNAseqMatrixDEgenes = EachContrastRNAseqCountMatrixAllDEgenes, 
                  NormalizedRNAseqMatrixDEgenes = EachContrastData2RNAseqFilterededgeRNorm)
  class(RESULTS) <- "DifferentialExpressionRNAseq_ASilva"
  return(RESULTS)
}
# [END]
