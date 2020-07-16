# 11 March 2020
# # Function: Produce violin plots of RNAseq counts of all samples against TYPE and 
#             STAGE, including mean, median and variance. 
# Author: Anjali Silva

# Input
# RNAseqCountMatrix: A matrix of RNAseq counts that has been NOT been normalized but 
#                    features filtered based on users choice. 
# ClinicalFile: File with patient sample names, and categories. 
# CategoryToVisualize: Category to visualize from the Clinical File;
#                      Options: "SAMPLE_ID", SAMPLE_ID_TGL", "INCLUDE", "TIME_POINT", 
#                     "RES_ID","LY_FL_ID", "OTHER_ID", "SEX","SITE_BIOPSY", "SITE_EXTRANODAL", 
#                     "TYPE_BIOPSY", "TRANSLOC_14_18", "INSTITUTION", "TYPE", "STAGE", 
#                     "COO", "EPIC_QC" or "Cluster". Default: "TYPE".
# ClusterLabels: A vector of integers with length equal to ncol(RNAseqCountMatrix). Default is NA. 
# ProduceImages: Produce images or not, options = "Yes" or "No".
# PNGorPDF: Output format of the image, options = "png" or "pdf".

# Output
# 37_ViolinPlot_MeanRNAseq_", CategoryToVisualize.p*
# 37_ViolinPlot_MedianRNAseq_", CategoryToVisualize.p*
# 37_ViolinPlot_VarianceRNAseq", CategoryToVisualize.p*


ViolinPlotsRNAseq <- function(RNAseqCountMatrix, 
                              ClinicalFile, 
                              CategoryToVisualize = "TYPE", 
                              ClusterLabels = NA,
                              ProduceImages = "Yes",
                              PNGorPDF = "png") {
    
  # Loading needed packages
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(grDevices)
  
  # checking user input
  if ((CategoryToVisualize == "CLUSTER") && (all(is.na(ClusterLabels)) == "TRUE")) {
    stop("\n CategoryToVisualize is set to CLUSTER, but no ClusterLabels provided.")
  }
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  
  if (CategoryToVisualize == "CLUSTER") {
    # Creating a matrix with values for CategoryToVisualize and RowMeans of probes using beta values
    clinicalCategoryPlotValues <- data.frame(Cluster = factor(ClusterLabels),
                                             MeanRNAseq = colMeans(RNAseqCountMatrix),
                                             MedianRNAseq = colMeans(RNAseqCountMatrix),
                                             VarianceRNAseq = colVars(RNAseqCountMatrix))
    
    # Define comparisons
    if(length(unique(ClusterLabels)) == 2) {
      ComparisonOptions <- list(c("1", "2"))
    } else if(length(unique(ClusterLabels)) == 3) {
      ComparisonOptions <- list(c("1", "2"),
                                       c("1", "3"), 
                                       c("2", "3"))
    } else if(length(unique(ClusterLabels)) == 4) {
      ComparisonOptions <- list(c("1", "2"), 
                                       c("1", "3"), 
                                       c("1", "4"), 
                                       c("2", "3"), 
                                       c("2", "4"), 
                                       c("3", "4"))
    } else if(length(unique(ClusterLabels)) == 5) {
      ComparisonOptions <- list(c("1", "2"), 
                                       c("1", "3"), 
                                       c("1", "4"), 
                                       c("1", "5"), 
                                       c("2", "3"), 
                                       c("2", "4"), 
                                       c("2", "5"), 
                                       c("3", "4"), 
                                       c("3", "5"))
    }
    
    # Define colours
    coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                        '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080')
    
    ViolinPlot <- ggplot2::ggplot(clinicalCategoryPlotValues, 
                                  aes(x = Cluster, y = MeanRNAseq, fill = Cluster)) +
                                  geom_violin(trim = FALSE) + 
                                  geom_boxplot(width = 0.1) +
                                  scale_y_continuous("Average Expression Using RNAseq Expected Counts") +
                                  scale_fill_manual(values = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                  scale_x_discrete(name = paste(CategoryToVisualize)) +
                                  theme_bw() + theme(panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank()) +
                                  theme(aspect.ratio = 1) # +
    # geom_violin(trim=FALSE) # to get rid of boxplot inside
    ggplot2::ggsave(paste0(pathNow, "/img/37_ViolinPlot_MeanRNAseq_", CategoryToVisualize, ".", PNGorPDF))
    
    p6 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Cluster", 
                           y = "MeanRNAseq", fill = "Cluster",
                           palette = coloursBarPlot[sort(unique(ClusterLabels))],
                           ylab = " Mean RNAseq Expected Count Value") + 
                           stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                           stat_compare_means(paired = FALSE)
    ggplot2::ggsave(paste0(pathNow, "/img/37_ViolinPlot_MeanRNAseq_", CategoryToVisualize, ".", PNGorPDF))
    
    p7 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Cluster", 
                           y = "MedianRNAseq", fill = "Cluster",
                           palette = coloursBarPlot[sort(unique(ClusterLabels))],
                           ylab = " Median RNAseq Expected Count Value") + 
                           stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                           stat_compare_means(paired = FALSE)
    ggplot2::ggsave(paste0(pathNow, "/img/37_ViolinPlot_MedianRNAseq_", CategoryToVisualize, ".", PNGorPDF))
    
    p8 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Cluster", 
                           y = "VarianceRNAseq", fill = "Cluster",
                           palette = coloursBarPlot[sort(unique(ClusterLabels))],
                           ylab = " Variance of RNAseq Expected Count Value") + 
                           stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                           stat_compare_means(paired = FALSE)
    ggplot2::ggsave(paste0(pathNow, "/img/37_ViolinPlot_VarianceRNAseq", CategoryToVisualize, ".", PNGorPDF))
    
    
    
  } else {
    # Creating a matrix with values for CategoryToVisualize and RowMeans of probes using beta values
    clinicalCategoryPlotValues <- data.frame(Category = ClinicalFile[, which(colnames(ClinicalFile) == CategoryToVisualize)],
                                             MeanRNAseq = colMeans(RNAseqCountMatrix),
                                             MedianRNAseq = colMeans(RNAseqCountMatrix),
                                             VarianceRNAseq = colVars(RNAseqCountMatrix))
    
    # Remove NA values, if present
    if(length(which(is.na(clinicalCategoryPlotValues) == TRUE)) > 0) {
      clinicalCategoryPlotValues <- clinicalCategoryPlotValues[- which(is.na(clinicalCategoryPlotValues) == TRUE), ]
    }
    
    # Change violin plot colors by groups
    if(CategoryToVisualize == "TYPE") {
      ColourChoice <- c("red", "green", "blue")
    } else if(CategoryToVisualize == "STAGE") {
      ColourChoice <- c("chocolate4", "orange")
      
    }
    
    ViolinPlot <- ggplot2::ggplot(clinicalCategoryPlotValues, 
                                  aes(x = Category, y = MeanRNAseq, fill = Category)) +
                                  geom_violin(trim = FALSE) + 
                                  geom_boxplot(width = 0.1) +
                                  scale_y_continuous("Average Expression Using RNAseq Expected Counts") +
                                  scale_fill_manual(values = ColourChoice) +
                                  scale_x_discrete(name = paste(CategoryToVisualize)) +
                                  theme_bw() + theme(panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank()) +
                                  theme(aspect.ratio = 1) # +
                                # geom_violin(trim=FALSE) # to get rid of boxplot inside
    ggplot2::ggsave(paste0(pathNow, "/img/37_ViolinPlot_MeanRNAseq_", CategoryToVisualize, ".", PNGorPDF))
    
    #ViolinPlot2 <- ggplot(clinicalCategoryPlotValues, 
    #                      aes(x = Category, y = MeanBeta, fill = Category)) +
    #                      geom_violin(trim = FALSE) + 
    #                      geom_jitter(width = 0.1) +
    #                      scale_y_continuous("Average methylation using beta values") +
    #                      scale_x_discrete(name = paste(CategoryToVisualize)) +
    #                      theme_bw() + theme(panel.grid.major = element_blank(), 
    #                                         panel.grid.minor = element_blank())+
    #                      theme(aspect.ratio = 1) # +
    # geom_violin(trim=FALSE) # to get rid of boxplot inside
    # ggsave(paste0(pathNow,"/img/37_ViolinPlot_WithJitter_MeanBetaValue_", CategoryToVisualize, ".", PNGorPDF))
    
    # Setting the number of comparisons for ggpubr
    if(length(unique(clinicalCategoryPlotValues$Category)) == 1) {
      ComparisonOptions <- list(names(table(clinicalCategoryPlotValues$Category))[1])
    } else if(length(unique(clinicalCategoryPlotValues$Category)) == 2) {
      ComparisonOptions <- list(names(table(clinicalCategoryPlotValues$Category))[1:2])
    } else if(length(unique(clinicalCategoryPlotValues$Category)) == 3) {
      ComparisonOptions <- list(names(table(clinicalCategoryPlotValues$Category))[1:2],
                                names(table(clinicalCategoryPlotValues$Category))[2:3],
                                names(table(clinicalCategoryPlotValues$Category))[c(1, 3)])
    } else if(length(unique(clinicalCategoryPlotValues$Category)) == 4) {
      ComparisonOptions <- list(names(table(clinicalCategoryPlotValues$Category))[1:2],
                                names(table(clinicalCategoryPlotValues$Category))[c(1, 3)], 
                                names(table(clinicalCategoryPlotValues$Category))[c(1, 4)], 
                                names(table(clinicalCategoryPlotValues$Category))[c(2, 3)], 
                                names(table(clinicalCategoryPlotValues$Category))[c(2, 4)], 
                                names(table(clinicalCategoryPlotValues$Category))[c(3, 4)])
    }
    
    if(length(unique(clinicalCategoryPlotValues$Category)) == 1) {
      # if only one value present
      p6 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Category", 
                           y = "MeanRNAseq", fill = "Category",
                           palette = ColourChoice,
                           ylab = " Mean RNAseq Expected Count Value") 
      ggplot2::ggsave(paste0(pathNow, "/img/37_ViolinPlot_MeanRNAseq_", CategoryToVisualize, ".", PNGorPDF))
      
      p7 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Category", 
                              y = "MedianRNAseq", fill = "Category",
                              palette = ColourChoice,
                              ylab = " Median RNAseq Expected Count Value") 
      ggplot2::ggsave(paste0(pathNow, "/img/37_ViolinPlot_MedianRNAseq_", CategoryToVisualize, ".", PNGorPDF))
      
      p8 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Category", 
                              y = "VarianceRNAseq", fill = "Category",
                              palette = ColourChoice,
                              ylab = " Variance of RNAseq Expected Count Value") 
      ggplot2::ggsave(paste0(pathNow, "/img/37_ViolinPlot_VarianceRNAseq", CategoryToVisualize, ".", PNGorPDF))
    } else {
      p6 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Category", 
                              y = "MeanRNAseq", fill = "Category",
                              palette = ColourChoice,
                              ylab = " Mean RNAseq Expected Count Value") + 
                              stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                              stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/37_ViolinPlot_MeanRNAseq_", CategoryToVisualize, ".", PNGorPDF))
      
      p7 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Category", 
                              y = "MedianRNAseq", fill = "Category",
                              palette = ColourChoice,
                              ylab = " Median RNAseq Expected Count Value") + 
                              stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                              stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/37_ViolinPlot_MedianRNAseq_", CategoryToVisualize, ".", PNGorPDF))
      
      p8 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Category", 
                              y = "VarianceRNAseq", fill = "Category",
                              palette = ColourChoice,
                              ylab = " Variance of RNAseq Expected Count Value") + 
                              stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                              stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/37_ViolinPlot_VarianceRNAseq", CategoryToVisualize, ".", PNGorPDF))
    }
  }

  
  RESULTS <- list(DataTable = clinicalCategoryPlotValues)
  
  class(RESULTS) <- "ViolinPlotsRNAseq_ASilva"
  return(RESULTS)
}
  
