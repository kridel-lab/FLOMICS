# Updated 17 June 2020
# Function: Create box plots with siginificance shown, for beta and M values for a given
#           clinical category or cluster. 
# Author: Anjali Silva

# Input:
# BetaMatrix: Matrix of beta values for probes x patients, with probes in rows and patients as columns.
# MvalueMatrix: Matrix of M values of size probes x patients. 
# ClinicalFile: File of patients and corresponding clinical category; 
#               Matrix of size patients x clinical categories.
# CategoryToVisualize: Category to visualize from the Clinical File;
#                      Options: "SAMPLE_ID", SAMPLE_ID_TGL", "INCLUDE", "TIME_POINT", "RES_ID","LY_FL_ID",
#                               "OTHER_ID", "SEX","SITE_BIOPSY", "SITE_EXTRANODAL", "TYPE_BIOPSY", 
#                               "TRANSLOC_14_18", "INSTITUTION", "TYPE", "STAGE", "COO", "EPIC_QC"
#                     Else  "CLUSTER" maybe provided, in which case ClusterLabels argument should be provided. 
# ClusterLabels: A vector of integers of length ncol(BetaMatrix), indicating cluster membership.
#                Default NA.   
# TumorPurity: A numeric vector of length equalling number of samples, indicating the purity level of each sample. 
# PlotClustersWithinCategories: If no needed, leave as "NA". If plot clusters within categories is required, "Yes".
#                               Default = "NA". If "Yes" need to provide ClusterLabels argument.
# PNGorPDF: Output format of the image, options = "png" or "pdf".

# Output: 

# Visuals saved to img folder
# 10_BoxPlot_MeanBetaValue_", CategoryToVisualize, ".".p*

BoxPlotsMethylation <- function(BetaMatrix, 
                                MvalueMatrix, 
                                ClinicalFile, 
                                CategoryToVisualize = "TYPE", 
                                ClusterLabels = NA,
                                TumorPurity,
                                PlotClustersWithinCategories = NA, 
                                ProduceImages = "Yes",
                                PNGorPDF = "png") {
    
  # Loading needed packages
  # LoadCheckPkg(RegularPckgs=c("tidyverse","ggplot2","ggpubr"))
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  
  # checking user input
  if (CategoryToVisualize == "CLUSTER" && all(is.na(ClusterLabels)) == "TRUE") {
    stop("\n CategoryToVisualize is set to CLUSTER, but no ClusterLabels provided.")
  }
  
  if ((all(is.na(PlotClustersWithinCategories)) == FALSE) && all(is.na(ClusterLabels)) == "TRUE") {
    stop("\n PlotClustersWithinCategories is set to 'Yes', but no ClusterLabels provided.")
  }
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  TumorPurity <- TumorPurity[match(colnames(BetaMatrix), rownames(TumorPurity)), ]
  
  
  if (CategoryToVisualize == "CLUSTER") {
    # Creating a matrix with values for CategoryToVisualize and RowMeans of probes using beta values
    clinicalCategoryPlotValues <- data.frame(Cluster = factor(ClusterLabels),
                                             TumorPurity = TumorPurity,
                                             MeanBeta = colMeans(BetaMatrix),
                                             MedianBeta = colMedians(BetaMatrix),
                                             VarianceBeta = colVars(BetaMatrix),
                                             MeanMvalue = colMeans(MvalueMatrix))
    
    # Define comparisons
    if(length(unique(ClusterLabels)) == 2) {
      ComparisonOptions <- list(c("1","2"))
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
    
    BoxPlot <- ggplot2::ggplot(clinicalCategoryPlotValues, 
                               aes(x = Cluster, y =  MeanBeta, fill = Cluster)) +
                               geom_boxplot(width = 0.1) +
                               labs(y = "Mean Beta Value", 
                                    x = paste(CategoryToVisualize),
                                    fill = paste(CategoryToVisualize)) +
                               scale_fill_manual(values = coloursBarPlot[sort(unique(ClusterLabels))]) +
                               theme_bw() + 
                               theme(panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank()) +
                               theme(aspect.ratio = 1) # +
    # geom_violin(trim=FALSE) # to get rid of boxplot inside
    ggsave(paste0(pathNow, "/img/10_BoxPlot_MeanBeta_", CategoryToVisualize, ".", PNGorPDF))
    
    p6 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Cluster", 
                            y = "MeanBeta", fill = "Cluster",
                            palette = coloursBarPlot[sort(unique(ClusterLabels))],
                            ylab = " Mean Beta Value") + 
                            stat_compare_means(comparisons = ComparisonOptions) + 
                            # Add pairwise comparisons p-value
                            stat_compare_means(paired = FALSE)
    ggsave(paste0(pathNow, "/img/10_BoxPlot_MeanBeta_", CategoryToVisualize, ".", PNGorPDF))
    
    p7 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Cluster", 
                            y = "MedianBeta", fill = "Cluster",
                            palette = coloursBarPlot[sort(unique(ClusterLabels))],
                            ylab = " Median Beta Value") + 
                            stat_compare_means(comparisons = ComparisonOptions) + 
                            # Add pairwise comparisons p-value
                            stat_compare_means(paired = FALSE)
    ggsave(paste0(pathNow, "/img/10_BoxPlot_MedianBeta_", CategoryToVisualize, ".", PNGorPDF))
    
    p8 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Cluster", 
                            y = "VarianceBeta", fill = "Cluster",
                            palette = coloursBarPlot[sort(unique(ClusterLabels))],
                            ylab = " Variance of Beta Value") + 
                            stat_compare_means(comparisons = ComparisonOptions) + 
                            # Add pairwise comparisons p-value
                            stat_compare_means(paired = FALSE)
    ggsave(paste0(pathNow, "/img/10_BoxPlot_VarianceBeta", CategoryToVisualize, ".", PNGorPDF))
    
    p9 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Cluster", 
                            y = "MeanMvalue", fill = "Cluster",
                            palette = coloursBarPlot[sort(unique(ClusterLabels))],
                            ylab = " Mean of MValue") + 
                            stat_compare_means(comparisons = ComparisonOptions) + 
                            # Add pairwise comparisons p-value
                            stat_compare_means(paired = FALSE)
    ggsave(paste0(pathNow, "/img/10_BoxPlot_MeanMvalue", CategoryToVisualize, ".", PNGorPDF))
    
    p10 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Cluster", 
                            y = "TumorPurity", fill = "Cluster",
                            palette = coloursBarPlot[sort(unique(ClusterLabels))],
                            ylab = "Tumor Purity") + 
                            stat_compare_means(comparisons = ComparisonOptions) + 
                            # Add pairwise comparisons p-value
                            stat_compare_means(paired = FALSE)
    ggsave(paste0(pathNow, "/img/10_BoxPlot_Purity", CategoryToVisualize, ".", PNGorPDF))
    
  } else {
  
    
    
    # If no need to plot clusters within categories
    if(is.na(PlotClustersWithinCategories) == TRUE) {
      # Creating a matrix with values for CategoryToVisualize and RowMeans of probes using beta values
      clinicalCategoryPlotValues <- data.frame(Category = factor(ClinicalFile[, which(colnames(ClinicalFile) == CategoryToVisualize)]),
                                               TumorPurity = TumorPurity,
                                               MeanBeta = colMeans(BetaMatrix),
                                               MedianBeta = colMedians(BetaMatrix),
                                               VarianceBeta = colVars(BetaMatrix),
                                               MeanMvalue = colMeans(MvalueMatrix))
                                              
      # Remove NA values, if present
      if(length(which(is.na(clinicalCategoryPlotValues) == TRUE)) > 0) {
        clinicalCategoryPlotValues <- clinicalCategoryPlotValues[- which(is.na(clinicalCategoryPlotValues$Category) == TRUE), ]
      }
      
      # Change violin plot colors by groups
      if(CategoryToVisualize == "TYPE") {
        ColourChoice <- c("#d6604d", "#66bd63", "#4575b4")
      } else if(CategoryToVisualize == "STAGE") {
        ColourChoice <- c("#762a83","#c2a5cf")
      } else if(CategoryToVisualize == "SEX") {
        ColourChoice <- c("#b35806","#fdb863")
      } else if(CategoryToVisualize == "SITE_BIOPSY") {
        ColourChoice <- c("#f1b6da","#c51b7d")
      } else if(CategoryToVisualize == "TYPE_BIOPSY") {
        ColourChoice <- c("#a6dba0", "#878787")
      } else if(CategoryToVisualize == "TRANSLOC_14_18") {
        ColourChoice <- c("#e0e0e0","#878787")
      } else {
        ColourChoice <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                            '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                            '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                            '#000075', '#808080')
      }
      
      
      
      # #                      Options: "SAMPLE_ID", SAMPLE_ID_TGL", "INCLUDE", "TIME_POINT", "RES_ID","LY_FL_ID",
      #                               "OTHER_ID", "SEX","SITE_BIOPSY", "SITE_EXTRANODAL", "TYPE_BIOPSY", "TRANSLOC_14_18", 
      #                               "INSTITUTION", "TYPE", "STAGE", "COO", "EPIC_QC"
      #                     Else  "CLUSTER" maybe provided, in which case ClusterLabels argument should be provided. 
      
      
      BoxPlot <- ggplot2::ggplot(clinicalCategoryPlotValues, 
                                 aes(x = Category, y = MeanBeta, fill = Category)) +
                                 geom_boxplot() +
                                 labs(y = "Mean Beta Value", 
                                      fill = paste(CategoryToVisualize), 
                                      x = paste(CategoryToVisualize)) +
                                 scale_fill_manual(values = ColourChoice) +
                                 theme_bw() + 
                                 theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) +
                                 theme(aspect.ratio = 1) # +
                        # geom_violin(trim=FALSE) # to get rid of boxplot inside
      ggplot2::ggsave(paste0(pathNow, "/img/10_BoxPlot_MeanBetaValue_", CategoryToVisualize, ".", PNGorPDF))
      
      
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
        
      if (length(unique(clinicalCategoryPlotValues$Category)) == 1) {
        p6 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, 
                                x = "Category", 
                                y = "MeanBeta", fill = "Category",
                                palette = ColourChoice,
                                xlab = paste(CategoryToVisualize),
                                add = "boxplot", ylab = " Mean Beta Value") 
          
        ggplot2::ggsave(paste0(pathNow, "/img/10_BoxPlot_MeanBetaValue_", CategoryToVisualize, ".", PNGorPDF))
        
        p7 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Category", 
                               y = "MedianBeta", fill = "Category",
                               palette = ColourChoice,
                               xlab = paste(CategoryToVisualize),
                               add = "boxplot", ylab = " Median Beta Value")
        
        ggplot2::ggsave(paste0(pathNow, "/img/10_BoxPlot_MedianBetaValue_", CategoryToVisualize, ".", PNGorPDF))
        
        p8 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Category", 
                               y = "VarianceBeta", fill = "Category",
                               palette = ColourChoice,
                               xlab = paste(CategoryToVisualize),
                               add = "boxplot", ylab = " Variance of Beta Values") 
          
        ggplot2::ggsave(paste0(pathNow, "/img/10_BoxPlot_VarianceOfBetaValue_", CategoryToVisualize, ".", PNGorPDF))
        
        p9 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, 
                                x = "Category", y = "MeanMvalue", fill = "Category",
                                palette = ColourChoice,
                                xlab = paste(CategoryToVisualize),
                                add = "boxplot", ylab="Mean Mvalue") 
        
        ggplot2::ggsave(paste0(pathNow, "/img/10_BoxPlot_MeanMValue_", CategoryToVisualize, ".", PNGorPDF))
        
      } else { 
        p6 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Category", 
                                y = "MeanBeta", fill = "Category",
                                palette = ColourChoice,
                                xlab = paste(CategoryToVisualize),
                                add = "boxplot", ylab = " Mean Beta Value") + 
                                stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                                stat_compare_means(paired = FALSE)
        
        ggplot2::ggsave(paste0(pathNow, "/img/10_BoxPlot_MeanBetaValue_", CategoryToVisualize, ".", PNGorPDF))
        
        p7 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Category", 
                                y = "MedianBeta", fill = "Category",
                                palette = ColourChoice,
                                xlab = paste(CategoryToVisualize),
                                add = "boxplot", ylab = " Median Beta Value") + 
                                stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                                stat_compare_means(paired = FALSE)
                              
        ggsave(paste0(pathNow, "/img/10_BoxPlot_MedianBetaValue_", CategoryToVisualize, ".", PNGorPDF))
        
        p8 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, x = "Category", 
                                y = "VarianceBeta", fill = "Category",
                                palette = ColourChoice,
                                xlab = paste(CategoryToVisualize),
                                add = "boxplot", ylab = " Variance of Beta Values") + 
                                stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                                stat_compare_means(paired = FALSE)
        
        ggplot2::ggsave(paste0(pathNow, "/img/10_BetaPlot_VarianceOfBetaValue_", CategoryToVisualize, ".", PNGorPDF))
        
        p9 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, 
                                x = "Category", y = "MeanMvalue", fill = "Category",
                                palette = ColourChoice,
                                xlab = paste(CategoryToVisualize),
                                add = "boxplot", ylab = "Mean Mvalue") + 
                                stat_compare_means(comparisons = ComparisonOptions)+ # Add pairwise comparisons p-value
                                stat_compare_means(paired = FALSE)
        
        ggplot2::ggsave(paste0(pathNow, "/img/10_BetaPlot_MeanMValue_", CategoryToVisualize, ".", PNGorPDF))
        
        p10 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, 
                                x = "Category", y = "TumorPurity", fill = "Category",
                                palette = ColourChoice,
                                xlab = paste(CategoryToVisualize),
                                add = "boxplot", ylab = "Tumor Purity") + 
                                stat_compare_means(comparisons = ComparisonOptions)+ # Add pairwise comparisons p-value
                                stat_compare_means(paired = FALSE)
        
        ggplot2::ggsave(paste0(pathNow, "/img/10_BetaPlot_Purity_", CategoryToVisualize, ".", PNGorPDF))
      }
      # If need to plot clusters within categories
    } else if(is.na(PlotClustersWithinCategories) == FALSE) { 
      # Creating a matrix with values for CategoryToVisualize and RowMeans of probes using beta values
      clinicalCategoryPlotValues <- data.frame(Cluster = factor(ClusterLabels),
                                               Category = ClinicalFile[ , which(colnames(ClinicalFile) == CategoryToVisualize)],
                                               TumorPurity = TumorPurity,
                                               MeanBeta = colMeans(BetaMatrix),
                                               MedianBeta = colMedians(BetaMatrix),
                                               VarianceBeta = colVars(BetaMatrix),
                                               MeanMvalue = colMeans(MvalueMatrix))
      
      # Remove NA values, if present
      if(length(which(is.na(clinicalCategoryPlotValues) == TRUE)) > 0) {
        clinicalCategoryPlotValues <- clinicalCategoryPlotValues[- which(is.na(clinicalCategoryPlotValues$Category) == TRUE), ]
      }
      
      # Define colours
      coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                          '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                          '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                          '#000075', '#808080')
      
      
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
      
        p6 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, 
                                x = "Category", 
                                y = "MeanBeta", 
                                fill = "Cluster",
                                add = "boxplot", 
                                palette = coloursBarPlot[sort(unique(ClusterLabels))],
                                ylab = " Mean Beta Value",
                                xlab = paste(CategoryToVisualize)) + 
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(paired = FALSE)
        
        ggplot2::ggsave(paste0(pathNow, "/img/10_BoxPlot_MeanBetaValue_", CategoryToVisualize, ".", PNGorPDF))
        
        p7 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, 
                                x = "Category", 
                                y = "MedianBeta", 
                                fill = "Cluster",
                                palette = coloursBarPlot[sort(unique(ClusterLabels))],
                                add = "boxplot", 
                                ylab = " Median Beta Value",
                                xlab = paste(CategoryToVisualize)) + 
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(paired = FALSE)
        
        ggsave(paste0(pathNow, "/img/10_BoxPlot_MedianBetaValue_", CategoryToVisualize, ".", PNGorPDF))
        
        p8 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, 
                                x = "Category", 
                                y = "VarianceBeta", 
                                fill = "Cluster",
                                palette = coloursBarPlot[sort(unique(ClusterLabels))],
                                add = "boxplot", 
                                ylab = " Variance of Beta Values",
                                xlab = paste(CategoryToVisualize)) + 
                                stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                                stat_compare_means(paired = FALSE)
                              
        ggplot2::ggsave(paste0(pathNow, "/img/10_BoxPlot_VarianceOfBetaValue_", CategoryToVisualize, ".", PNGorPDF))
        
        p9 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, 
                                x = "Category", 
                                y = "MeanMvalue", 
                                fill = "Cluster",
                                palette = coloursBarPlot[sort(unique(ClusterLabels))],
                                add = "boxplot", 
                                ylab ="Mean Mvalue", 
                                xlab = paste(CategoryToVisualize)) + 
                                stat_compare_means(comparisons = ComparisonOptions)+ # Add pairwise comparisons p-value
                                stat_compare_means(paired = FALSE)
        
        ggplot2::ggsave(paste0(pathNow, "/img/10_BoxPlot_MeanMValue_", CategoryToVisualize, ".", PNGorPDF))
        
        p10 <- ggpubr::ggboxplot(clinicalCategoryPlotValues, 
                                x = "Category", 
                                y = "TumorPurity", 
                                fill = "Cluster",
                                palette = coloursBarPlot[sort(unique(ClusterLabels))],
                                add = "boxplot", 
                                ylab ="Tumor Purity", 
                                xlab = paste(CategoryToVisualize)) + 
                                stat_compare_means(comparisons = ComparisonOptions)+ # Add pairwise comparisons p-value
                                stat_compare_means(paired = FALSE)
        
        ggplot2::ggsave(paste0(pathNow, "/img/10_Purity_", CategoryToVisualize, ".", PNGorPDF))
      
      
    }
    
  }  
  
  RESULTS <- list(DataTable = clinicalCategoryPlotValues)
  
  class(RESULTS) <- "BoxPlotsMehtylation_ASilva"
  return(RESULTS)
}
  



# Updated 3 April 2019
# Function: Generate Boxplots Based on Type (DLBCL, FL, RLN) ONLY given GeneOrProbeName
# Author: Anjali Silva

# Input:
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns.
# ClinicalFile: File with patient sample names, and categories. 
# CategoryToVisualize: Category to visualize from the Clinical File;
#                      Options: "SAMPLE_ID", SAMPLE_ID_TGL", "INCLUDE", "TIME_POINT", "RES_ID","LY_FL_ID",           
#                      "OTHER_ID", "SEX","SITE_BIOPSY" SITE_EXTRANODAL", "TYPE_BIOPSY", "TRANSLOC_14_18", 
#                      "INSTITUTION", "TYPE", "STAGE", "COO", "EPIC_QC" 
# GeneOrProbeName: Name of the gene or probe to be plotted (e.g. "CTSS");
#                  Note if a gene name is provided, need to provide the summarized beta matrix for BetaMatrix. 
# PNGorPDF: Output format of the image, options = "png" or "pdf"

# Output:
# DataFrame: 

# Visuals saved to img folder
# Boxplot 



BoxPlotMethylationGeneProbe <- function(BetaMatrix, 
                                        ClinicalFile, 
                                        CategoryToVisualize, 
                                        GeneOrProbeName = "none", 
                                        PNGorPDF) {
  
  
  # Loading needed packages
  # LoadCheckPkg(RegularPckgs=c("tidyr","ggplot2","plotly"))
  library(tidyr)
  library(ggplot2)
  library(plotly)
  
  # Getting the beta values
  if(CategoryToVisualize == "TYPE") {
    
    if (GeneOrProbeName == "none") {
      testdata2 <- data.frame(names = c (rep("FL", length(which(ClinicalFile$TYPE == "FL"))), 
                                         rep("DLBCL", length(which(ClinicalFile$TYPE == "DLBCL"))), 
                                         rep("RLN", length(which(ClinicalFile$TYPE == "RLN")))),
                              values = c(colMeans(BetaMatrix[ , which(ClinicalFile$TYPE == "FL")]),
                                         colMeans(BetaMatrix[ , which(ClinicalFile$TYPE == "DLBCL")]),
                                         colMeans(BetaMatrix[ , which(ClinicalFile$TYPE == "RLN")])))
    } else{
      testdata2 <- data.frame(names = c (rep("FL", length(which(ClinicalFile$TYPE == "FL"))), 
                                         rep("DLBCL", length(which(ClinicalFile$TYPE == "DLBCL"))), 
                                         rep("RLN", length(which(ClinicalFile$TYPE == "RLN")))),
                              values = c(unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName), 
                                                           which(ClinicalFile$TYPE == "FL")]),
                                         unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName), 
                                                           which(ClinicalFile$TYPE == "DLBCL")]),
                                         unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName),
                                                           which(ClinicalFile$TYPE == "RLN")])))
    }
    
    
    p2 <- ggplot(testdata2, aes(x = names, y = values, color = factor(names))) +
                 geom_boxplot(size = 2, outlier.size = 3) +
                 scale_y_continuous(name = "Methylation Beta Value") +
                 geom_jitter(aes(color = names), size = 3, alpha = 1, width = 0.2)+
                 labs(x = paste0(CategoryToVisualize)) +
                 theme_bw() + 
                 theme(text = element_text(size = 20), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.text.x = element_text(face = "bold"),
                       axis.text.y = element_text(face = "bold")) +
                 scale_color_manual(values = c("#d6604d", "#66bd63", "#4575b4"),
                                    name = paste0(CategoryToVisualize)) +
                 theme(aspect.ratio = 1, 
                      legend.position = "right", 
                      panel.background = element_rect(colour = "black", size=1.5),  
                      axis.title =  element_text(face = "bold"))
    
    pathNow <- getwd()
    ggsave(paste0(pathNow, "/img/12_BoxPlot=", GeneOrProbeName, "_", 
                  CategoryToVisualize, ".", PNGorPDF))
    
  } else if(CategoryToVisualize == "STAGE") {
    
    if (GeneOrProbeName == "none") {
      testdata2 <- data.frame(names = c(rep("ADVANCED", length(which(ClinicalFile$STAGE == "ADVANCED"))),
                                        rep("LIMITED", length(which(ClinicalFile$STAGE == "LIMITED")))),
                              values = c(unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName), 
                                                           which(ClinicalFile$STAGE == "ADVANCED")]),
                                         unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName), 
                                                           which(ClinicalFile$STAGE == "LIMITED")])))
      
    } else {
      testdata2 <- data.frame(names = c(rep("ADVANCED", length(which(ClinicalFile$STAGE == "ADVANCED"))),
                                      rep("LIMITED",length(which(ClinicalFile$STAGE == "LIMITED")))),
                              values = c(unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName), 
                                                         which(ClinicalFile$STAGE == "ADVANCED")]),
                                       unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName), 
                                                         which(ClinicalFile$STAGE =="LIMITED")])))
      
    }
    
    p2 <- ggplot(testdata2, aes(x = names, y = values, color = factor(names))) +
                 geom_boxplot(size = 2, outlier.size = 3) +
                 scale_y_continuous(name = "Methylation Beta Value") +
                 geom_jitter(aes(color = names), 
                             size = 3, alpha = 1, width = 0.2) +
                 labs(x = paste0(CategoryToVisualize)) +
                 theme_bw() + theme(text = element_text(size = 20), 
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), 
                                    axis.text.x = element_text(face="bold"), 
                                    axis.text.y = element_text(face="bold")) +
                 scale_color_manual(values = c("#fddbc7", "#b2182b", "red"),
                                    name = paste0(CategoryToVisualize) )+
                 theme(aspect.ratio = 1, legend.position = "right", 
                       panel.background = element_rect(colour = "black", size=1.5),  
                       axis.title =  element_text(face = "bold"))
    
    pathNow <- getwd()
    ggsave(paste0(pathNow, "/img/12_BoxPlot=", GeneOrProbeName, "_", CategoryToVisualize, ".", PNGorPDF))
    
  } else if(CategoryToVisualize == "SEX") {
    testdata2 <- data.frame(names = c(rep("F",length(which(ClinicalFile$SEX == "F"))), 
                                      rep("M",length(which(ClinicalFile$SEX == "M")))),
                            values = c(unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName), 
                                                         which(ClinicalFile$SEX == "F")]),
                                     unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName), 
                                                       which(ClinicalFile$SEX == "M")])))
    p2 <- ggplot(testdata2, aes(x = names, y = values, color = factor(names))) +
                 geom_boxplot(size = 2, outlier.size = 3) +
                 scale_y_continuous(name = "Methylation Beta Value") +
                 geom_jitter(aes(color = names), size = 3, alpha = 1, width = 0.2)+
                 labs(x = paste0(CategoryToVisualize))+
                 theme_bw() + theme(text = element_text(size = 20), 
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), 
                                    axis.text.x = element_text(face = "bold"), 
                                    axis.text.y = element_text(face = "bold") )+
                 scale_color_manual(values = c("#fddbc7", "#b2182b","red"),
                                    name = paste0(CategoryToVisualize) )+
                 theme(aspect.ratio = 1, legend.position = "right", 
                       panel.background = element_rect(colour = "black", size=1.5), 
                       axis.title =  element_text(face = "bold"))
              
    pathNow <- getwd()
    ggsave(paste0(pathNow, "/img/12_BoxPlot=", GeneOrProbeName, "_", CategoryToVisualize, ".", PNGorPDF))
  } else if(CategoryToVisualize == "TYPE_BIOPSY") {
    testdata2 <- data.frame(names = c(rep("CORE", length(which(ClinicalFile$TYPE_BIOPSY == "CORE"))), 
                                      rep("TISSUE",length(which(ClinicalFile$TYPE_BIOPSY == "TISSUE")))),
                            values = c(unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName), 
                                                         which(ClinicalFile$TYPE_BIOPSY == "CORE")]),
                                   unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName),
                                                     which(ClinicalFile$TYPE_BIOPSY == "TISSUE")])))
    p2 <- ggplot(testdata2, aes(x = names, y = values, color = factor(names))) +
                 geom_boxplot(size = 2, outlier.size = 3) +
                 scale_y_continuous(name = "Methylation Beta Value") +
                 geom_jitter(aes(color = names), size = 3, alpha = 1, width = 0.2) +
                 labs(x = paste0(CategoryToVisualize)) +
                 theme_bw() + theme(text = element_text(size = 20), 
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), 
                                    axis.text.x = element_text(face = "bold"), 
                                    axis.text.y = element_text(face="bold")) +
                 scale_color_manual(values = c("#fddbc7", "#b2182b", "red"), name = paste0(CategoryToVisualize)) +
                 theme(aspect.ratio=1, legend.position="right",
                       panel.background = element_rect(colour = "black", size = 1.5),  
                       axis.title =  element_text(face = "bold"))
    
    pathNow <- getwd()
    ggsave(paste0(pathNow, "/img/12_BoxPlot=", GeneOrProbeName, "_", CategoryToVisualize, ".", PNGorPDF))
  } else if(CategoryToVisualize == "COO") {
    testdata2 <- data.frame(names = c(rep("ABC", length(which(ClinicalFile$COO == "ABC"))), 
                                      rep("GCB", length(which(ClinicalFile$COO == "GCB")))),
                            values = c(unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName), 
                                                         which(ClinicalFile$COO == "ABC")]),
                                   unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName),
                                                     which(ClinicalFile$COO == "GCB")])))
    p2 <- ggplot(testdata2, aes(x = names, y = values, color = factor(names))) +
                 geom_boxplot(size = 2, outlier.size = 3) +
                 scale_y_continuous(name = "Methylation Beta Value") +
                 geom_jitter(aes(color = names), size = 3, alpha = 1, width = 0.2) +
                 labs(x = paste0(CategoryToVisualize)) +
                 theme_bw() + theme(text = element_text(size = 20), 
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), 
                                    axis.text.x = element_text(face = "bold"), 
                                    axis.text.y = element_text(face = "bold")) +
                 scale_color_manual(values = c("#fddbc7", "#b2182b", "red"),
                                    name = paste0(CategoryToVisualize)) +
                 theme(aspect.ratio = 1, legend.position = "right", 
                       panel.background = element_rect(colour = "black", size = 1.5),  
                       axis.title = element_text(face = "bold"))
    
    pathNow <- getwd()
    ggsave(paste0(pathNow, "/img/12_BoxPlot=", GeneOrProbeName, "_", CategoryToVisualize, ".", PNGorPDF))
  } else if(CategoryToVisualize == "INSTITUTION") {
    testdata2 <- data.frame(names = c(rep("BCCA",length(which(ClinicalFile$INSTITUTION == "BCCA"))), 
                                    rep("JGH",length(which(ClinicalFile$INSTITUTION == "JGH"))), 
                                    rep("UHN",length(which(ClinicalFile$INSTITUTION == "UHN")))),
                          values = c(unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName),
                                                       which(ClinicalFile$INSTITUTION == "BCCA")]),
                                   unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName),
                                                     which(ClinicalFile$INSTITUTION == "JGH")]),
                                   unlist(BetaMatrix[which(rownames(BetaMatrix) == GeneOrProbeName),
                                                     which(ClinicalFile$INSTITUTION == "UHN")])))
    p2 <- ggplot(testdata2, aes(x = names, y = values, color = factor(names))) +
                 geom_boxplot(size = 2, outlier.size = 3) +
                 scale_y_continuous(name = "Methylation Beta Value") +
                 geom_jitter(aes(color = names), size = 3, alpha = 1, width = 0.2)+
                 labs(x = paste0(CategoryToVisualize)) +
                 theme_bw() + theme(text = element_text(size = 20), 
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), 
                                    axis.text.x = element_text(face = "bold"), 
                                    axis.text.y = element_text(face = "bold")) +
                 scale_color_manual(values = c("#fddbc7", "#b2182b", "red"),
                                    name = paste0(CategoryToVisualize)) +
                 theme(aspect.ratio = 1, legend.position = "right", 
                       panel.background = element_rect(colour = "black", size = 1.5),  
                       axis.title =  element_text(face = "bold"))
    
    pathNow<-getwd()
    ggsave(paste0(pathNow, "/img/12_BoxPlot=", GeneOrProbeName, "_", CategoryToVisualize, ".", PNGorPDF))
  } else {
    cat("\n Not available for this category")
    testdata2 <- NA
  }
  
  RESULTS <- list(DataFrame = testdata2)
  
  class(RESULTS) <- "BoxPlotMethylationGeneProbe_ASilva"
  return(RESULTS)
}






