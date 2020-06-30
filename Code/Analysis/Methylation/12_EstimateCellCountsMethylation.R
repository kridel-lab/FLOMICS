# Created 25 June 2020
# Function: Create box plots with siginificance shown, for fraction of tumor represented by cell 
# type as estimated from FlowSorted.Blood.EPIC::estimateCellCounts2() function. 
# Author: Anjali Silva 

# Input:
# CellTypesFile: A matrix of size patients x cell types. 
# RGChannelSet: Only need to provide if CellTypesFile argument is not provided. RGChannelSet should
#               an object of class "RGChannelSet" or "RGChannelSetExtended" returned via
#               read.metharray.exp {minfi} function. This can be obatained from running function 2, 
#               LoadingMethylData() function in pipeline. The RGChannelSet will be used to generate
#               cell types using FlowSorted.Blood.EPIC::estimateCellCounts2() function. Default NA.
# ClinicalFile: File of patients and corresponding clinical category; 
#               Matrix of size patients x clinical categories.
# ClusterLabels: A vector of integers of length ncol(BetaMatrix), indicating cluster membership.
#                Default NA.   
# FigureGenerate: Should figures be generated, default = "Yes". Options: "Yes", "No".
# PNGorPDF: Output format of the image, options = "png" or "pdf".

# Output: 

# Visuals saved to img folder
# 12_BoxPlot_MeanBetaValue_", CategoryToVisualize, ".".p*

estimateCellCountsMethylation <- function(CellTypesFile = "NA",
                                          RGChannelSet = "NA",
                                          ClinicalFile, 
                                          ClusterLabels = NA,
                                          FigureGenerate = "Yes",
                                          PNGorPDF = "png") {
  

  # Loading needed packages
  # LoadCheckPkg(RegularPckgs=c("celltypes450")
  library(celltypes450)
  library(ggplot2)
  library(ggpubr)
  
  # check user input
  if((all(is.na(CellTypesFile)) == TRUE) && (all(is.na(RGChannelSet)) == TRUE)) {
    stop("\n Need to either provide CellTypesFile or RGChannelSet")
  }

 if((all(is.na(CellTypesFile)) == FALSE) && (nrow(CellTypesFile) != nrow(ClinicalFile))) {
   stop("\n nrow(CellTypesFile) does not equal nrow(ClinicalFile). The number of samples in 
        both files should be the same.")
 }
  
  # RGChannelSet is provided then run estimateCellCounts2()
  if (all(is.na(RGChannelSet)) == FALSE) {
  countsEPIC <- FlowSorted.Blood.EPIC::estimateCellCounts2(
                                  rgSet = RGChannelSet,
                                  compositeCellType = "Blood",
                                  processMethod = "preprocessNoob",
                                  probeSelect = "IDOL",
                                  cellTypes = c("CD8T", "CD4T", "NK", "Bcell","Mono", "Neu"),
                                  referencePlatform = "IlluminaHumanMethylationEPIC",
                                  referenceset = NULL,
                                  IDOLOptimizedCpGs = IDOLOptimizedCpGs,
                                  returnAll = TRUE,
                                  meanPlot = TRUE,
                                  verbose = FALSE)
  CellTypesFile <- countsEPIC
  }
  
  # Getting the path to save files, which should contain a folder called "img"
  pathNow <- getwd()
  
  # Create table from Clinical File
  Table <- CellTypesFile %>%
           data.frame() %>% mutate(SAMPLE_ID = row.names(.)) %>%
           left_join(ClinicalFile_T1[, c("SAMPLE_ID", "STAGE", "TYPE", "SEX", "SITE_BIOPSY", 
                                  "TYPE_BIOPSY", "INSTITUTION", "EPIC_QC")])
  Table <- data.frame(Table, "TRANSLOC_14_18" = factor(ClinicalFile$TRANSLOC_14_18))
  
  
  # If cluster labels are provided
  if(all(is.na(ClusterLabels)) == FALSE) {
    TableCluster <- data.frame(Table, "CLUSTER" = factor(ClusterLabels))
    
    
  # plotting by category
  for (category in c("STAGE", "TYPE", "SEX", "SITE_BIOPSY", 
                 "TYPE_BIOPSY", "INSTITUTION", "TRANSLOC_14_18", "EPIC_QC",
                 "CLUSTER")) {
    cat("\n", category)
    
    # pick colours for category
    if(category == "TYPE") {
      ColourChoice <- c("#d6604d", "#66bd63", "#4575b4")
    } else if(category == "STAGE") {
      ColourChoice <- c("#762a83","#c2a5cf")
    } else if(category == "SEX") {
      ColourChoice <- c("#b35806","#fdb863")
    } else if(category == "SITE_BIOPSY") {
      ColourChoice <- c("#f1b6da","#c51b7d")
    } else if(category == "TYPE_BIOPSY") {
      ColourChoice <- c("#a6dba0", "#878787")
    } else if(category == "TRANSLOC_14_18") {
      ColourChoice <- c("#e0e0e0","#878787")
    } else {
      ColourChoice <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                        '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080')
    }
    
    # if NA values are present for a given category, remove those
    if(length(which(is.na(TableCluster[, which(colnames(TableCluster) == category)]) == TRUE)) > 0) {
      # remove NA value entries
      TableMod <- TableCluster[- which(is.na(TableCluster[, which(colnames(TableCluster) == category)]) == TRUE), ]
      TableMod[, which(colnames(TableMod) == "TRANSLOC_14_18")] <- factor(TableMod$TRANSLOC_14_18)
      
    } else {
      TableMod <- TableCluster
      TableMod[, which(colnames(TableMod) == "TRANSLOC_14_18")] <- factor(TableMod$TRANSLOC_14_18)
      
    }
      
      # identifying comparison option for plot
      CategoriesVector <- TableMod[, which(colnames(TableMod) == category)]
      if(length(unique(CategoriesVector)) == 1) {
        ComparisonOptions <- list(names(table(CategoriesVector))[1])    
      } else if(length(unique(CategoriesVector)) == 2) {
        ComparisonOptions <- list(names(table(CategoriesVector))[1:2])
      } else if(length(unique(CategoriesVector)) == 3) {
        ComparisonOptions <- list(names(table(CategoriesVector))[1:2],
                                  names(table(CategoriesVector))[2:3],
                                  names(table(CategoriesVector))[c(1, 3)])
      } else if(length(unique(clinicalCategoryPlotValues$Category)) == 4) {
        ComparisonOptions <- list(names(table(CategoriesVector))[1:2],
                                  names(table(CategoriesVector))[c(1, 3)], 
                                  names(table(CategoriesVector))[c(1, 4)], 
                                  names(table(CategoriesVector))[c(2, 3)], 
                                  names(table(CategoriesVector))[c(2, 4)], 
                                  names(table(CategoriesVector))[c(3, 4)])
      }
      
      Table2 <- TableMod %>%
        reshape2::melt(id.vars = c("SAMPLE_ID", "STAGE", "TYPE", "SEX", "SITE_BIOPSY", 
                                   "TYPE_BIOPSY", "INSTITUTION", "TRANSLOC_14_18", "EPIC_QC",
                                   "CLUSTER")) 
      
      if (category == "STAGE") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                              add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                              xlab = "Cell type",
                              font.label = list(size = 20, color = "black"), 
                              palette = ColourChoice) +  
                              ggtitle("Cell type vs. fraction of tumor") +
                              stat_compare_means(comparisons = ComparisonOptions) + 
                              # Add pairwise comparisons p-value
                              stat_compare_means(aes(group = STAGE))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "TYPE") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = TYPE))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "SEX") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = SEX))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "SITE_BIOPSY") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = SITE_BIOPSY))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "TYPE_BIOPSY") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = TYPE_BIOPSY))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "INSTITUTION") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = INSTITUTION))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "TRANSLOC_14_18") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = "TRANSLOC_14_18",
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = TRANSLOC_14_18))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "EPIC_QC") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = EPIC_QC))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      
      if (category == "CLUSTER") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = CLUSTER))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }

      
    } 
  } else { # if cluster labels are not provided
    
    TableCluster <- data.frame(Table)
    
    # plotting by category
    for (category in c("STAGE", "TYPE", "SEX", "SITE_BIOPSY", 
                       "TYPE_BIOPSY", "INSTITUTION", "TRANSLOC_14_18", 
                       "EPIC_QC")) {
      cat("\n", category)
      
      # pick colours for category
      if(category == "TYPE") {
        ColourChoice <- c("#d6604d", "#66bd63", "#4575b4")
      } else if(category == "STAGE") {
        ColourChoice <- c("#762a83","#c2a5cf")
      } else if(category == "SEX") {
        ColourChoice <- c("#b35806","#fdb863")
      } else if(category == "SITE_BIOPSY") {
        ColourChoice <- c("#f1b6da","#c51b7d")
      } else if(category == "TYPE_BIOPSY") {
        ColourChoice <- c("#a6dba0", "#878787")
      } else if(category == "TRANSLOC_14_18") {
        ColourChoice <- c("#e0e0e0","#878787")
      } else {
        ColourChoice <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                          '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                          '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                          '#000075', '#808080')
      }
      
      # if NA values are present for a given category, remove those
      if(length(which(is.na(TableCluster[, which(colnames(TableCluster) == category)]) == TRUE)) > 0) {
        # remove NA value entries
        TableMod <- TableCluster[- which(is.na(TableCluster[, which(colnames(TableCluster) == category)]) == TRUE), ]
      } else {
        TableMod <- TableCluster
      }
      
      # identifying comparison option for plot
      CategoriesVector <- TableMod[, which(colnames(TableMod) == category)]
        if(length(unique(CategoriesVector)) == 1) {
          ComparisonOptions <- list(names(table(CategoriesVector))[1])    
        } else if(length(unique(CategoriesVector)) == 2) {
          ComparisonOptions <- list(names(table(CategoriesVector))[1:2])
        } else if(length(unique(CategoriesVector)) == 3) {
          ComparisonOptions <- list(names(table(CategoriesVector))[1:2],
                                    names(table(CategoriesVector))[2:3],
                                    names(table(CategoriesVector))[c(1, 3)])
        } else if(length(unique(clinicalCategoryPlotValues$Category)) == 4) {
          ComparisonOptions <- list(names(table(CategoriesVector))[1:2],
                                    names(table(CategoriesVector))[c(1, 3)], 
                                    names(table(CategoriesVector))[c(1, 4)], 
                                    names(table(CategoriesVector))[c(2, 3)], 
                                    names(table(CategoriesVector))[c(2, 4)], 
                                    names(table(CategoriesVector))[c(3, 4)])
        }
        
      Table2 <- TableMod %>%
          reshape2::melt(id.vars = c("SAMPLE_ID", "STAGE", "TYPE", "SEX", "SITE_BIOPSY", 
                                     "TYPE_BIOPSY", "INSTITUTION", "TRANSLOC_14_18", "EPIC_QC")) 
        
      if (category == "STAGE") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = STAGE))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "TYPE") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = TYPE))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "SEX") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = SEX))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "SITE_BIOPSY") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = SITE_BIOPSY))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "TYPE_BIOPSY") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = TYPE_BIOPSY))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "INSTITUTION") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = INSTITUTION))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "TRANSLOC_14_18") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = TRANSLOC_14_18))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
      if (category == "EPIC_QC") {
        p1 <- ggpubr::ggboxplot(Table2, x = "variable", y = "value", fill = category,
                                add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                                xlab = "Cell type",
                                font.label = list(size = 20, color = "black"), 
                                palette = ColourChoice) +  
                                ggtitle("Cell type vs. fraction of tumor") +
                                stat_compare_means(comparisons = ComparisonOptions) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = EPIC_QC))
        ggsave(paste0(pathNow, "/img/12_BoxPlot_FractionOfTumor_", category, ".", PNGorPDF))
      }
      
        
        
      } 
    }

  RESULTS <- list(DataTable = Table2)
  
  class(RESULTS) <- "BoxPlotsMehtylation_ASilva"
  return(RESULTS)
}
      
    
  