# Updated 3 Aug 2021
# Updated 28 July 2020
# Updated 12 April 2019
# Function: Huet23GeneModel
# Author: Anjali Silva

# Input:
# BetaMatrixSummarized: A matrix of beta values for genes (ONLY) x patients, with genes in
#                       rows and patients as columns.
# MvalueMatrixSummarized: A matrix of M values for genes (ONLY) x patients, with genes in 
#                         rows and patients as columns.
# RNAseqCountMatrixNormalized: 
# ClinicalFile: File with patient sample names, and categories. 
# FigureGenerate: Produce images or not, options = "Yes" or "No"; default "Yes".
# PNGorPDF: Output format of the image, options = "png" or "pdf"; default "png".


# Output:
# Probe23MatrixArrangedByType: matrix of Beta values of probes corresponding to 23 genes x samples. 
# Gene23MatrixArrangedByType: matrix of Beta values of 23 genes x samples. 

# Visuals saved to img folder
# 16_Huet23GeneModel_NonSummarized_Beta_ClusterSamples*
# 16_Huet23GeneModel_NonSummarized_Beta_NoClustering*
# 16_Huet23GeneModel_NonSummarized_Beta_FLSamplesOnly_Clustering*

Huet23GeneModel16 <- function(BetaMatrixSummarized = NA, 
                              BetaMatrixNotSummarized = NA, 
                              MvalueMatrixSummarized = NA, 
                              MethylationAnnotationFile = NA,
                              RNAseqCountMatrixNormalized,
                              RNAseqAnnotationFile,
                              ClinicalFile, 
                              FigureGenerate = "Yes", 
                              PNGorPDF = "png") {

  # Loading needed packages
  # LoadCheckPkg(RegularPckgs=c("gplots","ggplot2","pheatmap","RColorBrewer","IlluminaHumanMethylationEPICanno.ilm10b2.hg19"))
  library(gplots)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(minfi)
  
  # Getting the path of files
  pathNow <- getwd()
  
  # Based on Huet et al., 2018 paper 
  goodGenes <- c("SHISA8", "ABCB1", "METRNL", "VCL", "ALDH2", "RGS10")
  badGenes <- c("GADD45A", "FOXO1", "DCAF12", "CXCR4", "AFF3", "ORAI2", 
           "E2F5", "KIAA0040", "TAGAP", "TCF4", "SEMA4B", "FCRL2", "PRDM15", 
           "EML6", "USP44", "RASSF6", "VPREB1")
  
  notRun <- function() {
    # Getting location of genes in the summarized Beta values matrix that correspond to 23 genes
    GoodGenesLocation <- unlist(sapply(c(1:length(goodGenes)), function(i) 
      which(rownames(BetaMatrixSummarized) == goodGenes[i])))
    BadGeneLocation <- unlist(sapply(c(1:length(badGenes)), function(i) 
      which(rownames(BetaMatrixSummarized) == badGenes[i])))
    
    # Visualizations using summarized beta values
    if(FigureGenerate == "Yes") {
      
      set.seed(1) # altering set.seed() did not alter the results
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      
      if (PNGorPDF == "png") {
        grDevices::png(paste0(pathNow, "/img/16_Huet23GeneModel_Summarized_Beta_Values.", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/16_Huet23GeneModel_Summarized_Beta_Values.", PNGorPDF))
      }
      gplots::heatmap.2(as.matrix(BetaMatrixSummarized[c(BadGeneLocation,GoodGenesLocation), ]), 
                        dendrogram = "column",
                        trace = "none",
                        scale = "none",
                        Rowv = F, col = rev(redgreen(75)), 
                        key = T, density.info = "density", 
                        RowSideColors = c(rep("maroon", length(BadGeneLocation)), 
                                          rep("grey", length(GoodGenesLocation))))
      graphics::legend(y = .8, x = -0.5, xpd = TRUE,     
                       legend = c("Bad Genes", "Good Genes"),
                       col = c(unique(c(rep("maroon", length(BadGeneLocation)),
                                        rep("grey", length(GoodGenesLocation))))), 
                       lty= 1, lwd = 5, cex = .7)
      grDevices::dev.off()
    }
  
    # Getting location of probes in non-summarized Beta values corresponding to 23 genes
    
    # Getting the annotation file (also saved in Methylation data folder)
    ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
    
    # Matching rownames in methylation file with corresponding probe location in annotation file
    match_ids_methylation_beta <- match(rownames(BetaMatrix), ann[, which(colnames(ann) == "Name")])
    
    # Some annotations are missing. This was confirmed using excel sheet too. This is due to the 
    # way annotation file is.
    # column 22 is UCSC_RefGene_Name
    missing_annotation_methylation_beta <- which(ann[match_ids_methylation_beta,which(colnames(ann) == 
                                                                          "UCSC_RefGene_Name")] == "")
    
    # Getting annotations with missing gene rows removed
    ann_missing_removed_beta <- ann[match_ids_methylation_beta[-missing_annotation_methylation_beta],]
    
    # Gettting probes, such that probes with missing annotation removed
    beta_lymph_missing_removed <- BetaMatrix[-missing_annotation_methylation_beta,]
    
    # Only taking first feature
    ann_missing_removed_truncated <- sub ("\\;.*", "" , ann_missing_removed_beta[, 
                                          which(colnames(ann_missing_removed_beta) == "UCSC_RefGene_Name")] ) 
   
    GoodGenesLocation_Unsummarized <- sapply(c(1:length(goodGenes)), 
                                             function(i)  which(ann_missing_removed_truncated == goodGenes[i]))
    BadGeneLocation_Unsummarized <- sapply(c(1:length(badGenes)), 
                                           function(i) which(ann_missing_removed_truncated == badGenes[i]))
    
    # Visualizations using non-summarized beta values
    if(FigureGenerate == "Yes") {
      
      if (PNGorPDF == "png") {
        # Clustering samples only then view heatmap
        grDevices::png(paste0(pathNow, "/img/16_Huet23GeneModel_NonSummarized_Beta_ClusterSamples.", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/16_Huet23GeneModel_NonSummarized_Beta_ClusterSamples.", PNGorPDF))
      }
      gplots::heatmap.2(as.matrix(beta_lymph_missing_removed[c(unlist(BadGeneLocation_Unsummarized),
                                                               unlist(GoodGenesLocation_Unsummarized)), ]), 
                        dendrogram = 'column', Rowv = FALSE, trace = 'none', col = rev(redgreen(75)),
                        key = T, density.info="density", 
                        RowSideColors = col_vector[c(rep(1, length(unlist(BadGeneLocation_Unsummarized))),
                                                     rep(2, length(unlist(GoodGenesLocation_Unsummarized)))) + 4], 
                        ColSideColors = col_vector[c(ClinicalFile$TYPE) + 1])
      graphics::legend(y=.8, x=-0.5, xpd=TRUE,     
             legend = c("Bad Genes", "Good Genes", "FL samples", "DLBCL samples", "RLN samples"),
             col = c(unique(col_vector[c(rep(1, length(unlist(BadGeneLocation_Unsummarized))),
                                         rep(2, length(unlist(GoodGenesLocation_Unsummarized)))) + 4]), 
                     unique(col_vector[c(ClinicalFile$TYPE) + 1])), 
             lty = 1, lwd = 5, cex =.7)
      grDevices::dev.off()
      
      # Without clustering, look at heatmap using only by sample type 
      TypeLymphomVector <- c(which(substr(colnames(BetaMatrixNotSummarized), 4, 5) == "FL"), 
                             which(substr(colnames(BetaMatrixNotSummarized), 4, 5) == "DL"),
                             which(substr(colnames(BetaMatrixNotSummarized), 4, 5) == "RL"))
      if (PNGorPDF == "png") {
        grDevices::png(paste0(pathNow,"/img/16_Huet23GeneModel_NonSummarized_Beta_NoClustering.", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {  
        grDevices::pdf(paste0(pathNow,"/img/16_Huet23GeneModel_NonSummarized_Beta_NoClustering.", PNGorPDF))
      }
      gplots::heatmap.2(as.matrix(beta_lymph_missing_removed[c(unlist(BadGeneLocation_Unsummarized),
                                                               unlist(GoodGenesLocation_Unsummarized)),
                                                             TypeLymphomVector]),
                        dendrogram='column', 
                        Rowv=FALSE, 
                        Colv=FALSE, 
                        trace='none', 
                        col = rev(redgreen(75)), 
                        key = T, 
                        density.info="density", 
                        RowSideColors = col_vector[c(rep(1, length(unlist(BadGeneLocation_Unsummarized))),
                                                     rep(2, length(unlist(GoodGenesLocation_Unsummarized)))) + 4], 
                        ColSideColors = col_vector[c(ClinicalFile$TYPE)[TypeLymphomVector] + 1] )
      graphics::legend(y = 1.2, x = .25, xpd = TRUE,     
             legend = c("Bad Genes", "Good Genes", "FL samples", "DLBCL samples", "RLN samples"),
             col = c(unique(col_vector[c(rep(1, length(unlist(BadGeneLocation_Unsummarized))),
                                         rep(2, length(unlist(GoodGenesLocation_Unsummarized)))) + 4]), 
                     unique(col_vector[c(ClinicalFile$TYPE)[TypeLymphomVector] + 1])), 
             lty = 1, lwd = 5, cex=.7)
      grDevices::dev.off()
      
      # Without clustering, look at heatmap using only by sample type FL
      TypeLymphomVector_FL <- c(which(substr(colnames(BetaMatrixNotSummarized), 4, 5) == "FL"))
      if (PNGorPDF == "png") {
        png(paste0(pathNow, "/img/16_Huet23GeneModel_NonSummarized_Beta_FLSamplesOnly_Clustering.", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        pdf(paste0(pathNow, "/img/16_Huet23GeneModel_NonSummarized_Beta_FLSamplesOnly_Clustering.", PNGorPDF))
      }
      # pheatmap(as.matrix(beta_lymph_missing_removed[c(unlist(BadGeneLocation_Unsummarized),
      # unlist(GoodGenesLocation_Unsummarized)),TypeLymphomVector_FL]), cluster_rows = FALSE,  
      # scale ="none", show_colnames = T, fontface="italic", legend = T, border_color = "black", 
      # color =  rev(redgreen(1000)) )
      gplots::heatmap.2(as.matrix(beta_lymph_missing_removed[c(unlist(BadGeneLocation_Unsummarized),
                                                               unlist(GoodGenesLocation_Unsummarized)),
                                                             TypeLymphomVector_FL]), 
                        dendrogram = 'column', 
                        Rowv = FALSE,trace = 'none', col = rev(redgreen(75)), key = T, 
                        density.info = "density", 
                        RowSideColors = col_vector[c(rep(1, length(unlist(BadGeneLocation_Unsummarized))),
                                                     rep(2, length(unlist(GoodGenesLocation_Unsummarized)))) + 4], 
                        ColSideColors = col_vector[c(ClinicalFile$TYPE)[TypeLymphomVector_FL] + 1] )
      graphics::legend(y = 1.2, x = .25, xpd = TRUE,      
             legend = c("Bad Genes", "Good Genes", "FL samples"),
             col = c(unique(col_vector[c(rep(1, length(unlist(BadGeneLocation_Unsummarized))),
                                         rep(2, length(unlist(GoodGenesLocation_Unsummarized)))) + 4]),
                     unique(col_vector[c(ClinicalFile$TYPE)[TypeLymphomVector_FL] + 1])), 
             lty= 1, lwd = 5, cex = .7)
      grDevices::dev.off()
    }
  }
  
  
  # RNAseq analysis
  # Annotate RNA matrix with gene names
  RNAseqCountMatrixNormalizedName <- data.frame(RNAseqCountMatrixNormalized,
                                                name = ENSMBLid[match(rownames(RNAseqCountMatrixNormalized), 
                                                                      RNAseqAnnotationFile$gene), ])
  # retain the same column names as 
  colnames(RNAseqCountMatrixNormalizedName)[(ncol(RNAseqCountMatrixNormalized) + 1):ncol(RNAseqCountMatrixNormalizedName)] <- 
    colnames(RNAseqAnnotationFile)
  
  # Create matrix 23 genes x patients
  
  goodGenesMatched <- match(goodGenes, RNAseqCountMatrixNormalizedName$name)
  badGenesMatched <- match(badGenes, RNAseqCountMatrixNormalizedName$name)
  
  # giving all the bad genes (1:15) a coefficient of 1
  # giving all good genes (16:20) a coefficient of -1
  
  scoreSamples <- sapply(c(1:ncol(RNAseqCountMatrixNormalizedName[goodGenesMatched, 
                               1:ncol(RNAseqCountMatrixNormalized)])), 
                               function(sample) 
                               (sum((1 * RNAseqCountMatrixNormalizedName[badGenesMatched, sample]), na.rm = TRUE) + 
                                sum((- 1 * RNAseqCountMatrixNormalizedName[goodGenesMatched, sample]), na.rm = TRUE)))
  
  names(scoreSamples) <- colnames(RNAseqCountMatrixNormalized)

  
  
  if(all(is.na(BetaMatrixSummarized)) != TRUE) {
  
  RESULTS <- list(Probe23MatrixArrangedByType = as.matrix(beta_lymph_missing_removed[c(unlist(BadGeneLocation_Unsummarized),
                                                                                       unlist(GoodGenesLocation_Unsummarized)),
                                                                                     TypeLymphomVector]),
                  Gene23MatrixArrangedByType = as.matrix(BetaMatrixSummarized[c(BadGeneLocation,
                                                                                GoodGenesLocation),
                                                                              TypeLymphomVector]))
  } else {
    RESULTS <- list(scores = scoreSamples,
                    goodGenesMatched = RNAseqCountMatrixNormalizedName$name[goodGenesMatched],
                    badGenesMatched = RNAseqCountMatrixNormalizedName$name[badGenesMatched],
                    RNAseqMatrix = RNAseqCountMatrixNormalizedName[c(goodGenesMatched, badGenesMatched), ])
  }
  
  
  class(RESULTS) <- "Huet23GeneModel_ASilva"
  return(RESULTS)
}
# [END]
