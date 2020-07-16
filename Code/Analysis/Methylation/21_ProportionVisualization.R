# Updated 3 April 2019
# Function: Visualize proportion of data falling within a category or ClusterLabels. 
#           - If PlotWithinCategories is set, but no ClusterLabels are provided, 
#             PlotWithinCategories is set to "No", and ClinicalCategory = "STAGE",  
#             then values for PlotWithinCategories for STAGE can be visualized. 
#           - If PlotWithinCategories is set, but no ClusterLabels are provided and
#             PlotWithinCategories is set to "Yes", then values for PlotWithinCategories
#             can be visualized. E.g., if CategoryToVisualize = "Relation_to_Island"
#             and PlotWithinCategories = "Yes", and ClinicalCategory = "STAGE", 
#             then values for STAGE in terms of "Island", "Shore", etc may be visualized.
#           - If PlotWithinCategories is NA, but only ClusterLabels are provided and
#             PlotWithinCategories is set to "No", then values for only the clusters 
#             can be visualized, in terms of entire Beta Matrix (unless filtered). 
#           - If PlotWithinCategories is set, and ClusterLabels are provided and
#             PlotWithinCategories is set to "Yes", then values for only the clusters 
#             can be visualized, in terms of "Island", "Shore", etc.
# Author: Anjali Silva


# Input:
# CategoryToVisualize: Specifies which column needs to be visualized from AnnotationFile.
#                      Need exact spelling and caps/noCaps for CategoryToVisualize, e.g.: 
#                      CategoryToVisualize = "DMR", "UCSC_RefGene_Group",Options "Relation_to_Island", 
#                      "chr", "Regulatory_Feature_Group", "X450k_Enhancer", "Phantom4_Enhancers", 
#                      "Phantom5_Enhancers". Default NA. Must be provided, unless ClusterLabels 
#                       are given. 
# ClusterLabels: If proportions are being visualized for clusters, provide a vector of integers 
#                indicating cluster membership.
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns.
# AnnotationFile: Matrix of annotations for all the probes found in BetaMatrix. It is of size probes x annotations.
# ClinicalFile: File of patients and corresponding clinical category; matrix of size patients x clinical categories.
# ClinicalCategory: Specify column of ClinicalFile to use. Options "STAGE", "SEX", "SITE_BIOPSY", "TYPE_BIOPSY", 
#                   "INSTITUTION", "COO", "TYPE", "TRANSLOC_14_18". CLUSTER is also possible if provided. 
# PlotWithinCategories: Should indicate "Yes" or "No", as to whether categories within each category to 
#                       visualize (CategoryToVisualize) will to be investigated. Deafult is "No".
#                       E.g., if CategoryToVisualize = "Relation_to_Island", and if PlotWithinCategories = "Yes",
#                       then plots will be generated seperately for "Island", "Shore", etc. If 
#                       PlotWithinCategories = "No", then plots will be generated for all probes falling within
#                       "Relation_to_Island". 
# PNGorPDF: Output format of the image, options = "png" or "pdf"



# Output: 
# Proportion0.35and0.6: A data frame of proportion of data values that fell between 0.35 and 0.6 
#                      for each person. 
# TestResults: Output from comparing means of groups present. If only two groups to compare, then
#              t-test will be performed. Otherwise, One-Way ANOVA test is performed. 


# Visuals saved to img folder
# 21_Proportion_Average_*",ClinicalCategory.p*
# 21_Proportion_Btw0.2and0.8_ggplotPoint_",ClinicalCategory,.p*
# 21_Proportion_Btw0.2and0.8_Boxplot_",ClinicalCategory,.p*
# 21_Proportion_Btw0.2and0.8_ScatterPlot_",ClinicalCategory,.p*
# 21_ggboxplot_MeanBetaValueAcrossProbes_ANOVA_,ClinicalCategory,.p*
# 21_ggboxplot_MedianBetaValueAcrossProbes_ANOVA_,ClinicalCategory,.p*


ProportionVisualization <- function(CategoryToVisualize = NA, 
                                    ClusterLabels = NA, 
                                    BetaMatrix, 
                                    AnnotationFile = NA, 
                                    ClinicalFile = NA, 
                                    ClinicalCategory = NA, 
                                    PlotWithinCategories = "No", 
                                    FigureGenerate = "No", 
                                    PNGorPDF = "png", 
                                    ImageName = paste0("21_ProportionVisualization_", date())) {
  
  # Loading needed packages
  # LoadCheckPkg(RegularPckgs=c("ggplot2","stringi","data.table","ggpubr"))
  # "stringi" used to remove extra white space, if present, e.g. "EN " to "EN"
  library(ggplot2)
  library(stringi)
  library(data.table)
  library(ggpubr)
  library(dplyr)
  
  # Remove empty entries in the clinical file
  ClinicalFile[ClinicalFile == ""] <- NA
  
  # check user input
  if((all(is.na(CategoryToVisualize)) == TRUE)  && (all(is.na(ClusterLabels)) == TRUE)) {
    stop("CategoryToVisualize must be provided, unless ClusterLabels are given.")
  }
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  

  
  
  # If cluster labels are provided 
  if (sum(is.na(ClusterLabels)) == 0) {
    
    cat("\n Cluster labels are provided.")
    
    
    # Define colours for clusters
    coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                        '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080')
    
    # All probes and samples only corresponding to those specified in ClinicalCategoryToVisualize
    Cluster1Betas <- rowMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
    Cluster1BetasMean <- rowMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
    
    Cluster2Betas <- rowMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
    Cluster2BetasMean <- rowMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
    
    
    # Setting the number of comparisons for ggpubr plots
    # https://www.r-bloggers.com/add-p-values-and-significance-levels-to-ggplots/
    ComparisonOptionsCluster <- list(c("one", "two"))
    
    if (max(unique(ClusterLabels)) == 3) {
      Cluster3Betas <- rowMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      Cluster3BetasMean <- rowMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      
      
      # Setting the number of comparisons for ggpubr plots
      ComparisonOptionsCluster <- list(c("one", "two"),
                                       c("one", "three"),
                                       c("two", "three"))
    } else if (max(unique(ClusterLabels)) == 4) {
      Cluster3Betas <- rowMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      Cluster3BetasMean <- rowMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      
      Cluster4Betas <- rowMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 4)]))
      Cluster4BetasMean <- rowMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 4)]))
      
      
      # Setting the number of comparisons for ggpubr plots
      ComparisonOptionsCluster <- list(c("one", "two"),
                                       c("one", "three"),
                                       c("one", "four"),
                                       c("two", "three"),
                                       c("two", "four"),
                                       c("three", "four"))
    } else if (max(unique(ClusterLabels)) == 5) {
      Cluster3Betas <- rowMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      Cluster4Betas <- rowMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 4)]))
      Cluster5Betas <- rowMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 5)]))
      
      Cluster3BetasMean <- rowMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      Cluster4BetasMean <- rowMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 4)]))
      Cluster5BetasMean <- rowMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 5)]))
      
      
      # Setting the number of comparisons for ggpubr plots
      ComparisonOptionsCluster <- list(c("one", "two"),
                                       c("one", "three"),
                                       c("one", "four"),
                                       c("one", "five"),
                                       c("two", "three"),
                                       c("two", "four"),
                                       c("two", "five"),
                                       c("three", "four"),
                                       c("three", "five"),
                                       c("four", "five"))
    } else if (max(unique(ClusterLabels)) > 5) {
      stop("\nCode only supports upto 5 clusters. Please alter code");
    }
    
    
    
    # Plotting variability across clusters
    if (FigureGenerate == "Yes" && max(unique(ClusterLabels)) == 2) {
      
      # All probes and samples only corresponding to those specified in ClinicalCategoryToVisualize
      Cluster1Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
      Cluster2Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
      
      Cluster1Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
      Cluster2Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
      
      Cluster1Betas_Mean <- colMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
      Cluster2Betas_Mean <- colMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
      
      
      Cluster1Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 1)]
      Cluster2Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 2)]
      
      Cluster1Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 1)]
      Cluster2Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 2)]
      
      # Create data frame
      data_patient_var <- data.frame(Cluster = c(rep("one", length(Cluster1Betas_Var)),
                                                 rep("two", length(Cluster2Betas_Var))),
                                     variance = c(Cluster1Betas_Var, Cluster2Betas_Var),
                                     medians = c(Cluster1Betas_Median, Cluster2Betas_Median),
                                     means = c(Cluster1Betas_Mean, Cluster2Betas_Mean),
                                     institute = c(Cluster1Institute, Cluster2Institute),
                                     samples = c(Cluster1Samples, Cluster2Samples))
      
      
      # Plot median beta value for patients in each cluster       
      p0_mean <- ggplot2::ggplot(data_patient_var, aes(x = Cluster, 
                                                       y = means, 
                                                       color = factor(Cluster))) +
                                 # geom_point(size=3)
                                 scale_y_continuous(name = "Mean Beta values") +
                                 geom_jitter(aes(color = Cluster), 
                                            size = 3, 
                                            alpha = 1, 
                                            width = 0.2) +
                                 labs(x = "Cluster") +
                                 # annotate("text", x = 1.5, y = max(data_patient_var$variance) + 0.01, size = 7) +
                                 scale_color_manual(values = coloursBarPlot[1:2], 
                                                   name = "Cluster" ) +
                                 theme_bw() + 
                                 theme(axis.title =  element_text(face = "bold"), 
                                       aspect.ratio = 1, 
                                       legend.position = "right", 
                                       text = element_text(size = 20), 
                                       panel.grid.major = element_blank(), 
                                       panel.grid.minor = element_blank(), 
                                       panel.background = element_rect(colour = "black", size=1.5), 
                                       axis.text.x = element_text(face = "bold"), 
                                       axis.text.y = element_text(face = "bold"))
      ggplot2::ggsave(paste0(pathNow,"/img/21_ggplotPoint_MeanBetaValue_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      # Plot median beta value for patients in each cluster 
      p0_median <- ggplot2::ggplot(data_patient_var, aes(x = Cluster, 
                                                         y = medians, 
                                                         color = factor(Cluster))) +
                                   # geom_point(size=3)
                                   scale_y_continuous(name = "Median Beta values") +
                                   geom_jitter(aes(color = Cluster), 
                                               size = 3, 
                                               alpha = 1, 
                                               width = 0.2) +
                                   labs(x = "Cluster") +
                                   # annotate("text", x = 1.5, y = max(data_patient_var$variance) + 0.01, size = 7) +
                                   scale_color_manual(values = coloursBarPlot[1:2], 
                                                      name = "Cluster" ) +
                                   theme_bw() + 
                                   theme(axis.title =  element_text(face = "bold"), 
                                         aspect.ratio = 1, 
                                         legend.position = "right", 
                                         text = element_text(size = 20), 
                                         panel.grid.major = element_blank(), 
                                         panel.grid.minor = element_blank(), 
                                         panel.background = element_rect(colour = "black", size=1.5), 
                                         axis.text.x = element_text(face = "bold"), 
                                         axis.text.y = element_text(face = "bold"))
      ggplot2::ggsave(paste0(pathNow,"/img/21_ggplotPoint_MedianBetaValue_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      # Plotting median beta value for patients in each cluster, seperated by institute
      # Setting the number of comparisons for ggpubr
      if(length(unique(data_patient_var$institute)) == 2) {
        ComparisonOptionsInstitute <- list(names(table(data_patient_var$institute))[1:2])
      }else if(length(unique(data_patient_var$institute)) == 3) {
        ComparisonOptionsInstitute <- list(names(table(data_patient_var$institute))[1:2], 
                                           names(table(data_patient_var$institute))[2:3],
                                           names(table(data_patient_var$institute))[c(1, 3)])
      }else if(length(unique(data_patient_var$institute)) == 4) {
        ComparisonOptionsInstitute <- list( names(table(data_patient_var$institute))[1:2], 
                                            names(table(data_patient_var$institute))[c(1, 3)], 
                                            names(table(data_patient_var$institute))[c(1, 4)], 
                                            names(table(data_patient_var$institute))[c(2, 3)], 
                                            names(table(data_patient_var$institute))[c(2, 4)], 
                                            names(table(data_patient_var$institute))[c(3, 4)])
      }
      
      p0_boxplot_median_institute <- ggpubr::ggboxplot(data_patient_var, 
                                                       x = "institute", 
                                                       y = "medians", 
                                                       fill = "Cluster",
                                                       add = "boxplot", 
                                                       xlab = "institute", 
                                                       ylab = "Median Beta Values",
                                                       palette = coloursBarPlot[1:2]) + 
                                                       stat_compare_means(comparisons = ComparisonOptionsInstitute) + 
                                                       # Add pairwise comparisons p-value
                                                       stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow,"/img/21_ggboxplot_MedianBetaValueAcrossClusters_Institute_", 
                    ImageName, ".", PNGorPDF))
      
      
      # Plot variance of beta value for patients in each cluster 
      p0_variance <- ggplot2::ggplot(data_patient_var, aes(x = Cluster, 
                                                  y = variance, 
                                                  color=factor(Cluster))) +
                                     # geom_point(size=3)+
                                     scale_y_continuous(name = "Varaince of Beta values") +
                                     geom_jitter(aes(color = Cluster), size = 3, alpha = 1, width = 0.2) +
                                     labs(x = "Cluster") +
                                     # annotate("text", x = 1.5, y = max(data_patient_var$variance) + 0.01, size = 7) +
                                     theme_bw() + 
                                     theme(text = element_text(size = 20), 
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  axis.text.x = element_text(face="bold"), 
                                                  axis.text.y = element_text(face="bold") ) +
                                     scale_color_manual(values = coloursBarPlot[1:2], 
                                                       name = "Cluster" ) +
                                     theme(aspect.ratio = 1,
                                          legend.position = "right", 
                                          panel.background = element_rect(colour = "black", size=1.5),
                                          axis.title =  element_text(face = "bold"))
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggplotPoint_Variance_Cluster_", 
                             ImageName, ".", PNGorPDF))
      
      # Setting the number of comparisons for ggpubr
      p0_violin <- ggpubr::ggboxplot(data_patient_var, 
                                     x = "Cluster", 
                                     y = "variance", 
                                     fill = "Cluster",
                                     add = "boxplot", 
                                     xlab = "Cluster", 
                                     ylab = "Variance",
                                     palette = coloursBarPlot[1:2]) + 
                                     stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                                     # Add pairwise comparisons p-value
                                     stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggboxplot_Variance_Violin_Cluster_", 
                             ImageName, ".", PNGorPDF))
      
      p0_scatter <- ggpubr::ggscatter(data_patient_var, 
                                      x = "Cluster", 
                                      y = "variance",
                                      color = "Cluster", 
                                      palette = coloursBarPlot[1:2]) + 
                                      stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                                      # Add pairwise comparisons p-value
                                      stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggplotPoint_Variance_Scatter_Cluster_", 
                             ImageName, ".", PNGorPDF))
      
    }
    
    if (FigureGenerate == "Yes" && max(unique(ClusterLabels)) == 3) {
      # All probes and samples only corresponding to those specified in ClinicalCategoryToVisualize
      Cluster1Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
      Cluster2Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
      Cluster3Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      
      Cluster1Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
      Cluster2Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
      Cluster3Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      
      Cluster1Betas_Mean <- colMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
      Cluster2Betas_Mean <- colMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
      Cluster3Betas_Mean <- colMeans(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      
      
      Cluster1Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 1)]
      Cluster2Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 2)]
      Cluster3Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 3)]
      
      Cluster1Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 1)]
      Cluster2Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 2)]
      Cluster3Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 3)]
      
      # Create data frame
      data_patient_var <- data.frame(Cluster = c(rep("one", length(Cluster1Betas_Var)),
                                                 rep("two", length(Cluster2Betas_Var)),
                                                 rep("three", length(Cluster3Betas_Var))),
                                     variance = c(Cluster1Betas_Var, 
                                                  Cluster2Betas_Var, 
                                                  Cluster3Betas_Var),
                                     medians = c(Cluster1Betas_Median, 
                                                 Cluster2Betas_Median, 
                                                 Cluster3Betas_Median),
                                     means = c(Cluster1Betas_Mean, 
                                               Cluster2Betas_Mean, 
                                               Cluster3Betas_Mean),
                                     institute = c(Cluster1Institute, 
                                                   Cluster2Institute, 
                                                   Cluster3Institute),
                                     samples = c(Cluster1Samples, 
                                                 Cluster2Samples, 
                                                 Cluster3Samples))
      data_patient_var$Cluster <- factor(data_patient_var$Cluster, 
                                         levels = c("one", "two", "three"), 
                                         labels = c("one", "two", "three"))

      # Plot mean beta value for patients in each cluster 
      p0_mean <- ggplot2::ggplot(data_patient_var, 
                                   aes(x = Cluster, 
                                       y = means, color = factor(Cluster))) +
                                # geom_point(size=3)
                                scale_y_continuous(name = "Mean Beta values") +
                                geom_jitter(aes(color = Cluster),
                                            size = 3, alpha = 1, width = 0.2) +
                                labs(x = "Cluster") +
                                # annotate("text", x = 1.5, y = max(data_patient_var$variance) + 0.01, size = 7) +
                                scale_color_manual(values = coloursBarPlot[1:3], 
                                                   name = "Cluster" ) +
                                theme_bw() + 
                                theme(axis.title =  element_text(face = "bold"), 
                                      aspect.ratio = 1, 
                                      legend.position = "right", 
                                      text = element_text(size = 20), 
                                      panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), 
                                      panel.background = element_rect(colour = "black", size=1.5), 
                                      axis.text.x = element_text(face = "bold"), 
                                      axis.text.y = element_text(face = "bold"))
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggplotPoint_MeanBetaValue_Cluster_", 
                             ImageName, ".", PNGorPDF))
      
            
      # Plot median beta value for patients in each cluster 
      p0_median <- ggplot2::ggplot(data_patient_var, 
                                   aes(x = Cluster, 
                                       y = medians, color=factor(Cluster))) +
                                  # geom_point(size=3)
                                  scale_y_continuous(name = "Median Beta values") +
                                  geom_jitter(aes(color = Cluster),
                                              size = 3, alpha = 1, width = 0.2) +
                                  labs(x = "Cluster") +
                                  # annotate("text", x = 1.5, y = max(data_patient_var$variance) + 0.01, size = 7) +
                                  scale_color_manual(values = coloursBarPlot[1:3], 
                                                     name = "Cluster" ) +
                                  theme_bw() + 
                                  theme(axis.title =  element_text(face = "bold"), 
                                        aspect.ratio = 1, 
                                        legend.position = "right", 
                                        text = element_text(size = 20), 
                                        panel.grid.major = element_blank(), 
                                        panel.grid.minor = element_blank(), 
                                        panel.background = element_rect(colour = "black", size=1.5), 
                                        axis.text.x = element_text(face = "bold"), 
                                        axis.text.y = element_text(face = "bold"))
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggplotPoint_MedianBetaValue_Cluster_", 
                             ImageName, ".", PNGorPDF))
      
      # Plotting median beta value for patients in each cluster, seperated by institute
      # Setting the number of comparisons for ggpubr
      if(length(unique(data_patient_var$institute)) == 2) {
        ComparisonOptionsInstitute <- list(names(table(data_patient_var$institute))[1:2])
      }else if(length(unique(data_patient_var$institute)) == 3) {
        ComparisonOptionsInstitute <- list(names(table(data_patient_var$institute))[1:2], 
                                           names(table(data_patient_var$institute))[2:3],
                                           names(table(data_patient_var$institute))[c(1, 3)])
      }else if(length(unique(data_patient_var$institute)) == 4) {
        ComparisonOptionsInstitute <- list(names(table(data_patient_var$institute))[1:2], 
                                           names(table(data_patient_var$institute))[c(1, 3)], 
                                           names(table(data_patient_var$institute))[c(1, 4)], 
                                           names(table(data_patient_var$institute))[c(2, 3)], 
                                           names(table(data_patient_var$institute))[c(2, 4)], 
                                           names(table(data_patient_var$institute))[c(3, 4)])
      }
      
      p0_boxplot_median_institute <- ggpubr::ggboxplot(data_patient_var, 
                                                       x = "institute", 
                                                       y = "medians", 
                                                       fill = "Cluster",
                                                       add = "boxplot", 
                                                       xlab = "institute", 
                                                       ylab = "Median Beta Values",
                                                       palette = coloursBarPlot[1:3]) + 
                                                       stat_compare_means(comparisons = ComparisonOptionsInstitute) + 
                                                       # Add pairwise comparisons p-value
                                                       stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggboxplot_MedianBetaValueAcrossClusters_Institute_", 
                                                      ImageName, ".", PNGorPDF))
      
      
      # Plot variance of beta value for patients in each cluster 
      p0_variance <- ggplot2::ggplot(data_patient_var, 
                                     aes(x = Cluster, y = variance, color = factor(Cluster))) +
                                     #geom_point(size=3)+
                                     scale_y_continuous(name = "Varaince of Beta values") +
                                     geom_jitter(aes(color = Cluster), 
                                                 size = 3, alpha = 1, width = 0.2) +
                                     labs(x = "Cluster") +
                                     #annotate("text", x = 1.5, y = max(data_patient_var$variance) + 0.01, size = 7) +
                                     theme_bw() + theme(text = element_text(size = 20), 
                                                       panel.grid.major = element_blank(),
                                                       panel.grid.minor = element_blank(), 
                                                       axis.text.x = element_text(face="bold"), 
                                                       axis.text.y = element_text(face="bold") ) +
                                     scale_color_manual(values = coloursBarPlot[1:3], name = "Cluster" ) +
                                     theme(aspect.ratio = 1,
                                          legend.position = "right", 
                                          panel.background = element_rect(colour = "black", size=1.5),
                                          axis.title =  element_text(face = "bold"))
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggplotPoint_Variance_Cluster_", 
                             ImageName, ".", PNGorPDF))
      

      p0_violin <- ggpubr::ggboxplot(data_patient_var, 
                                     x = "Cluster", 
                                     y = "variance", 
                                     fill = "Cluster",
                                     add = "boxplot", 
                                     xlab = "Cluster", 
                                     ylab = "Variance",
                                     palette = coloursBarPlot[1:3]) + 
                                     stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                                     # Add pairwise comparisons p-value
                                     stat_compare_means(paired = FALSE)
      ggsave(paste0(pathNow, "/img/21_ggboxplot_Variance_Violin_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      p0_scatter <- ggpubr::ggscatter(data_patient_var, 
                                      x = "Cluster", 
                                      y = "variance",
                                      color = "Cluster", 
                                      palette = coloursBarPlot[1:3])+ 
                                      stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                                      # Add pairwise comparisons p-value
                                      stat_compare_means(paired = FALSE)
      ggsave(paste0(pathNow, "/img/21_ggplotPoint_Variance_Scatter_Cluster_", 
                    ImageName, ".", PNGorPDF))
    }
    
    if (FigureGenerate == "Yes" && max(unique(ClusterLabels)) == 4) {
      # All probes and samples only corresponding to those specified in ClinicalCategoryToVisualize
      Cluster1Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
      Cluster2Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
      Cluster3Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      Cluster4Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 4)]))
      
      Cluster1Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
      Cluster2Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
      Cluster3Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      Cluster4Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 4)]))      
      
      Cluster1Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 1)]
      Cluster2Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 2)]
      Cluster3Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 3)]
      Cluster4Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 4)]
      
      Cluster1Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 1)]
      Cluster2Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 2)]
      Cluster3Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 3)]
      Cluster4Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 4)]
      
      # Create data frame
      data_patient_var <- data.frame(Cluster = c(rep("one", length(Cluster1Betas_Var)),
                                                 rep("two", length(Cluster2Betas_Var)),
                                                 rep("three", length(Cluster3Betas_Var)),
                                                 rep("four", length(Cluster4Betas_Var))),
                                     variance = c(Cluster1Betas_Var, 
                                                  Cluster2Betas_Var, 
                                                  Cluster3Betas_Var,
                                                  Cluster4Betas_Var),
                                     medians = c(Cluster1Betas_Median, 
                                                 Cluster2Betas_Median, 
                                                 Cluster3Betas_Median,
                                                 Cluster4Betas_Median),
                                     institute = c(Cluster1Institute, 
                                                   Cluster2Institute, 
                                                   Cluster3Institute,
                                                   Cluster4Institute),
                                     samples = c(Cluster1Samples, 
                                                 Cluster2Samples, 
                                                 Cluster3Samples,
                                                 Cluster4Samples))
      data_patient_var$Cluster <- factor(data_patient_var$Cluster, 
                                         levels = c("one", "two", "three", "four"), 
                                         labels = c("one", "two", "three", "four"))
      
      # Plot median beta value for patients in each cluster 
      p0_median <- ggplot2::ggplot(data_patient_var, 
                                   aes(x = Cluster, 
                                       y = medians, color=factor(Cluster))) +
                                   # geom_point(size=3)
                                   scale_y_continuous(name = "Median Beta values") +
                                   geom_jitter(aes(color = Cluster),
                                              size = 3, alpha = 1, width = 0.2) +
                                   labs(x = "Cluster") +
                                   # annotate("text", x = 1.5, y = max(data_patient_var$variance) + 0.01, size = 7) +
                                   scale_color_manual(values = coloursBarPlot[1:4], 
                                                     name = "Cluster" ) +
                                   theme_bw() + 
                                   theme(axis.title =  element_text(face = "bold"), 
                                         aspect.ratio = 1, 
                                         legend.position = "right", 
                                         text = element_text(size = 20), 
                                         panel.grid.major = element_blank(), 
                                         panel.grid.minor = element_blank(), 
                                         panel.background = element_rect(colour = "black", size=1.5), 
                                         axis.text.x = element_text(face = "bold"), 
                                         axis.text.y = element_text(face = "bold"))
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggplotPoint_MedianBetaValue_Cluster_", 
                             ImageName, ".", PNGorPDF))
      
      # Plotting median beta value for patients in each cluster, seperated by institute
      # Setting the number of comparisons for ggpubr
      if(length(unique(data_patient_var$institute)) == 2) {
        ComparisonOptionsInstitute <- list(names(table(data_patient_var$institute))[1:2])
      }else if(length(unique(data_patient_var$institute)) == 3) {
        ComparisonOptionsInstitute <- list(names(table(data_patient_var$institute))[1:2], 
                                           names(table(data_patient_var$institute))[2:3],
                                           names(table(data_patient_var$institute))[c(1, 3)])
      }else if(length(unique(data_patient_var$institute)) == 4) {
        ComparisonOptionsInstitute <- list( names(table(data_patient_var$institute))[1:2], 
                                            names(table(data_patient_var$institute))[c(1, 3)], 
                                            names(table(data_patient_var$institute))[c(1, 4)], 
                                            names(table(data_patient_var$institute))[c(2, 3)], 
                                            names(table(data_patient_var$institute))[c(2, 4)], 
                                            names(table(data_patient_var$institute))[c(3, 4)])
      }
      
      p0_boxplot_median_institute <- ggpubr::ggboxplot(data_patient_var, 
                                                       x = "institute", 
                                                       y = "medians", 
                                                       fill = "Cluster",
                                                       add = "boxplot", 
                                                       xlab = "institute", 
                                                       ylab = "Median Beta Values",
                                                       palette = coloursBarPlot[1:4]) + 
                                                       stat_compare_means(comparisons = ComparisonOptionsInstitute) + 
                                                       # Add pairwise comparisons p-value
                                                       stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggboxplot_MedianBetaValueAcrossClusters_Institute_", 
                    ImageName, ".", PNGorPDF))
      
      
      # Plot variance of beta value for patients in each cluster 
      p0_variance <- ggplot2::ggplot(data_patient_var, 
                                     aes(x = Cluster, y = variance, color=factor(Cluster))) +
                                     # geom_point(size=3) +
                                     scale_y_continuous(name = "Varaince of Beta values") +
                                     geom_jitter(aes(color = Cluster), size = 3, alpha = 1, width = 0.2) +
                                     labs(x = "Cluster") +
                                     #annotate("text", x = 1.5, y = max(data_patient_var$variance) + 0.01, size = 7) +
                                     theme_bw() + theme(text = element_text(size=20), 
                                                       panel.grid.major = element_blank(),
                                                       panel.grid.minor = element_blank(), 
                                                       axis.text.x = element_text(face="bold"), 
                                                       axis.text.y = element_text(face="bold") ) +
                                     scale_color_manual(values = coloursBarPlot[1:4], name = "Cluster" ) +
                                     theme(aspect.ratio = 1,
                                          legend.position = "right", 
                                          panel.background = element_rect(colour = "black", size=1.5),
                                          axis.title =  element_text(face = "bold"))
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggplotPoint_Variance_Cluster_", 
                             ImageName, ".", PNGorPDF))
      

      p0_violin <- ggpubr::ggboxplot(data_patient_var, 
                                     x = "Cluster", 
                                     y = "variance", 
                                     fill = "Cluster",
                                     add = "boxplot", 
                                     xlab = "Cluster", 
                                     ylab = "Variance",
                                     palette = coloursBarPlot[1:4]) + 
                                     stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                                     # Add pairwise comparisons p-value
                                     stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggboxplot_Variance_Violin_Cluster_", 
                             ImageName, ".", PNGorPDF))
      
      p0_scatter <- ggpubr::ggscatter(data_patient_var, 
                                      x = "Cluster", 
                                      y = "variance",
                                      color = "Cluster", 
                                      palette = coloursBarPlot[1:4])+ 
                                      stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                                      # Add pairwise comparisons p-value
                                      stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggplotPoint_Variance_Scatter_Cluster_", 
                             ImageName, ".", PNGorPDF))
    }
    
    if (FigureGenerate == "Yes" && max(unique(ClusterLabels)) == 5) {
      # All probes and samples only corresponding to those specified in ClinicalCategoryToVisualize
      Cluster1Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
      Cluster2Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
      Cluster3Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      Cluster4Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 4)]))
      Cluster5Betas_Var <- colVars(as.matrix(BetaMatrix[ , which(ClusterLabels == 5)]))
      
      Cluster1Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 1)]))
      Cluster2Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 2)]))
      Cluster3Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 3)]))
      Cluster4Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 4)]))     
      Cluster5Betas_Median <- colMedians(as.matrix(BetaMatrix[ , which(ClusterLabels == 5)]))      
      
      
      Cluster1Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 1)]
      Cluster2Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 2)]
      Cluster3Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 3)]
      Cluster4Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 4)]
      Cluster5Institute <- ClinicalFile$INSTITUTION[which(ClusterLabels == 5)]
      
      Cluster1Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 1)]
      Cluster2Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 2)]
      Cluster3Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 3)]
      Cluster4Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 4)]
      Cluster5Samples <- ClinicalFile$SAMPLE_ID[which(ClusterLabels == 5)]
      
      # Create data frame
      data_patient_var <- data.frame(Cluster = c(rep("one", length(Cluster1Betas_Var)),
                                                 rep("two", length(Cluster2Betas_Var)),
                                                 rep("three", length(Cluster3Betas_Var)),
                                                 rep("four", length(Cluster4Betas_Var)),
                                                 rep("five", length(Cluster5Betas_Var))),
                                     variance = c(Cluster1Betas_Var, 
                                                  Cluster2Betas_Var, 
                                                  Cluster3Betas_Var,
                                                  Cluster4Betas_Var,
                                                  Cluster5Betas_Var),
                                     medians = c(Cluster1Betas_Median, 
                                                 Cluster2Betas_Median, 
                                                 Cluster3Betas_Median,
                                                 Cluster4Betas_Median,
                                                 Cluster5Betas_Median),
                                     institute = c(Cluster1Institute, 
                                                   Cluster2Institute, 
                                                   Cluster3Institute,
                                                   Cluster4Institute,
                                                   Cluster5Institute),
                                     samples = c(Cluster1Samples, 
                                                 Cluster2Samples, 
                                                 Cluster3Samples,
                                                 Cluster4Samples,
                                                 Cluster5Samples))
      data_patient_var$Cluster <- factor(data_patient_var$Cluster, 
                                         levels = c("one", "two", "three", "four", "five"), 
                                         labels = c("one", "two", "three", "four", "five"))
      
      # Plot median beta value for patients in each cluster 
      p0_median <- ggplot2::ggplot(data_patient_var, 
                                   aes(x = Cluster, 
                                       y = medians, color=factor(Cluster))) +
                                   # geom_point(size=3)
                                   scale_y_continuous(name = "Median Beta values") +
                                   geom_jitter(aes(color = Cluster),
                                              size = 3, alpha = 1, width = 0.2) +
                                   labs(x = "Cluster") +
                                   # annotate("text", x = 1.5, y = max(data_patient_var$variance) + 0.01, size = 7) +
                                   scale_color_manual(values = coloursBarPlot[1:5], 
                                                     name = "Cluster" ) +
                                   theme_bw() + 
                                   theme(axis.title =  element_text(face = "bold"), 
                                        aspect.ratio = 1, 
                                        legend.position = "right", 
                                        text = element_text(size = 20), 
                                        panel.grid.major = element_blank(), 
                                        panel.grid.minor = element_blank(), 
                                        panel.background = element_rect(colour = "black", size=1.5), 
                                        axis.text.x = element_text(face = "bold"), 
                                        axis.text.y = element_text(face = "bold"))
      ggplot2::ggsave(paste0(pathNow,"/img/21_ggplotPoint_MedianBetaValue_Cluster_", 
                             ImageName, ".", PNGorPDF))
      
      # Plotting median beta value for patients in each cluster, seperated by institute
      # Setting the number of comparisons for ggpubr
      if(length(unique(data_patient_var$institute)) == 2) {
        ComparisonOptionsInstitute <- list(names(table(data_patient_var$institute))[1:2])
      } else if(length(unique(data_patient_var$institute)) == 3) {
        ComparisonOptionsInstitute <- list(names(table(data_patient_var$institute))[1:2], 
                                           names(table(data_patient_var$institute))[2:3],
                                           names(table(data_patient_var$institute))[c(1, 3)])
      } else if(length(unique(data_patient_var$institute)) == 4) {
        ComparisonOptionsInstitute <- list( names(table(data_patient_var$institute))[1:2], 
                                            names(table(data_patient_var$institute))[c(1, 3)], 
                                            names(table(data_patient_var$institute))[c(1, 4)], 
                                            names(table(data_patient_var$institute))[c(2, 3)], 
                                            names(table(data_patient_var$institute))[c(2, 4)], 
                                            names(table(data_patient_var$institute))[c(3, 4)])
      }
      
      p0_boxplot_median_institute <- ggpubr::ggboxplot(data_patient_var, 
                                                       x = "institute", 
                                                       y = "medians", 
                                                       fill = "Cluster",
                                                       add = "boxplot", 
                                                       xlab = "institute", 
                                                       ylab = "Median Beta Values",
                                                       palette = coloursBarPlot[1:5]) + 
                                                       stat_compare_means(comparisons = ComparisonOptionsInstitute) + 
                                                       # Add pairwise comparisons p-value
                                                       stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow,"/img/21_ggboxplot_MedianBetaValueAcrossClusters_Institute_", 
                    ImageName, ".", PNGorPDF))
      
      
      # Plot variance of beta value for patients in each cluster 
      p0_variance <- ggplot2::ggplot(data_patient_var, 
                                     aes(x = Cluster, y = variance, color=factor(Cluster))) +
                                     # geom_point(size=3) +
                                     scale_y_continuous(name = "Varaince of Beta values") +
                                     geom_jitter(aes(color = Cluster), size = 3, alpha = 1, width = 0.2) +
                                     labs(x = "Cluster") +
                                     # annotate("text", x = 1.5, y = max(data_patient_var$variance) + 0.01, size = 7) +
                                     theme_bw() + theme(text = element_text(size=20), 
                                                       panel.grid.major = element_blank(),
                                                       panel.grid.minor = element_blank(), 
                                                       axis.text.x = element_text(face="bold"), 
                                                       axis.text.y = element_text(face="bold") ) +
                                     scale_color_manual(values = coloursBarPlot[1:5], name = "Cluster" ) +
                                     theme(aspect.ratio = 1,
                                          legend.position = "right", 
                                          panel.background = element_rect(colour = "black", size=1.5),
                                          axis.title =  element_text(face = "bold"))
      ggplot2::ggsave(paste0(pathNow,"/img/21_ggplotPoint_Variance_Cluster_", 
                             ImageName, ".", PNGorPDF))
      
      
      p0_violin <- ggpubr::ggboxplot(data_patient_var, 
                                     x = "Cluster", 
                                     y = "variance", 
                                     fill = "Cluster",
                                     add = "boxplot", 
                                     xlab = "Cluster", 
                                     ylab = "Variance",
                                     palette = coloursBarPlot[1:5]) + 
                                     stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                                     # Add pairwise comparisons p-value
                                     stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow,"/img/21_ggboxplot_Variance_Violin_Cluster_", 
                             ImageName, ".", PNGorPDF))
      
      p0_scatter <- ggpubr::ggscatter(data_patient_var, 
                                      x = "Cluster", 
                                      y = "variance",
                                      color = "Cluster", 
                                      palette = coloursBarPlot[1:5])+ 
                                      stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                                      # Add pairwise comparisons p-value
                                      stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow,"/img/21_ggplotPoint_Variance_Scatter_Cluster_", 
                    ImageName, ".", PNGorPDF))
    }
    
    
    # Calculating values between 0.0 to 0.2 (low); >0.2 to <=0.6 (medium); and > 0.6 <= 1.0 (high)
    Cluster1Betas_0to0.2 <- which(data.table::between(Cluster1Betas, 
                                                      lower = 0, upper = 0.2) == TRUE)
    Cluster2Betas_0to0.2 <- which(data.table::between(Cluster2Betas, 
                                                      lower = 0, upper = 0.2) == TRUE)
    # Medium range is between 0.2 and 0.6 as defined by Landau et al., 2014
    Cluster1Betas_0.21to0.6 <- which(data.table::between(Cluster1Betas, 
                                                         lower = 0.20000000000001, upper = 0.6) == TRUE)
    Cluster2Betas_0.21to0.6 <- which(data.table::between(Cluster2Betas, 
                                                         lower = 0.20000000000001, upper = 0.6) == TRUE)
    # High range is between 0.6 and 1.0 as defined by Landau et al., 2014
    Cluster1Betas_0.61to1 <- which(data.table::between(Cluster1Betas, 
                                                       lower = 0.60000000000001, upper = 1.0) == TRUE)
    Cluster2Betas_0.61to1 <- which(data.table::between(Cluster2Betas, 
                                                       lower = 0.60000000000001, upper = 1.0) == TRUE)
    
    # Assign proportions
    Cluster1Betas_proportion <- vector(mode = "character", 
                                       length = length(Cluster1Betas))
    Cluster1Betas_proportion[Cluster1Betas_0to0.2] <- "low"
    Cluster1Betas_proportion[Cluster1Betas_0.21to0.6] <- "mid"
    Cluster1Betas_proportion[Cluster1Betas_0.61to1] <- "high"
    
    Cluster2Betas_proportion <- vector(mode = "character", 
                                       length = length(Cluster2Betas))
    Cluster2Betas_proportion[Cluster2Betas_0to0.2] <- "low"
    Cluster2Betas_proportion[Cluster2Betas_0.21to0.6] <- "mid"
    Cluster2Betas_proportion[Cluster2Betas_0.61to1] <- "high"
    
    # Drawing bar graph
    probes <- c(rownames(BetaMatrix), rownames(BetaMatrix))
    average_probes <- c(Cluster1Betas, Cluster2Betas)
    ranges_probes <- c(Cluster1Betas_proportion, Cluster2Betas_proportion)
    class <- c(as.character(rep("Cluster1", length(Cluster1Betas))), 
               as.character(rep("Cluster2", length(Cluster2Betas))))
    
    average <- as.data.frame(cbind(probes, average_probes, ranges_probes, class))
    average$average_probes <- as.numeric(as.character(average$average_probes))
    
    # If three clusters present
    if (max(unique(ClusterLabels)) == 3) {
      # Low range is between 0 and 0.2 as defined by Landau et al., 2014
      Cluster3Betas_0to0.2 <- which(data.table::between(Cluster3Betas, 
                                                        lower = 0, upper = 0.2) == TRUE)
      # Medium range is between 0.2 and 0.6 as defined by Landau et al., 2014
      Cluster3Betas_0.21to0.6 <- which(data.table::between(Cluster3Betas, 
                                                           lower = 0.20000000000001, upper = 0.6) == TRUE)
      # High range is between 0.6 and 1.0 as defined by Landau et al., 2014
      Cluster3Betas_0.61to1 <- which(data.table::between(Cluster3Betas, 
                                                         lower = 0.60000000000001, upper = 1.0) == TRUE)
      
      # Assign proportions
      Cluster3Betas_proportion <- vector(mode = "character", 
                                         length = length(Cluster3Betas))
      Cluster3Betas_proportion[Cluster3Betas_0to0.2] <- "low"
      Cluster3Betas_proportion[Cluster3Betas_0.21to0.6] <- "mid"
      Cluster3Betas_proportion[Cluster3Betas_0.61to1] <- "high"
      
      probes = c(rownames(BetaMatrix), rownames(BetaMatrix), rownames(BetaMatrix))
      average_probes = c(Cluster1Betas, Cluster2Betas, Cluster3Betas)
      ranges_probes = c(Cluster1Betas_proportion, 
                        Cluster2Betas_proportion, 
                        Cluster3Betas_proportion)
      class <- c(as.character(rep("Cluster1", length(Cluster1Betas))), 
                 as.character(rep("Cluster2", length(Cluster2Betas))), 
                 as.character(rep("Cluster3", length(Cluster3Betas))) )
      
      average <- as.data.frame(cbind(probes, average_probes, ranges_probes, class))
      average$average_probes <- as.numeric(as.character(average$average_probes))
      
    }
    
    # If four clusters present
    if (max(unique(ClusterLabels)) == 4) {
      
      Cluster3Betas_0to0.2 <- which(data.table::between(Cluster3Betas, 
                                                        lower = 0, upper = 0.2) == TRUE)
      # Medium range is between 0.2 and 0.6 as defined by Landau et al., 2014
      Cluster3Betas_0.21to0.6 <- which(data.table::between(Cluster3Betas, 
                                                           lower = 0.20000000000001, upper = 0.6) == TRUE)
      # High range is between 0.6 and 1.0 as defined by Landau et al., 2014
      Cluster3Betas_0.61to1 <- which(data.table::between(Cluster3Betas, 
                                                         lower = 0.60000000000001, upper = 1.0) == TRUE)
      
      # Assign proportions
      Cluster3Betas_proportion <- vector(mode = "character", 
                                         length = length(Cluster3Betas))
      Cluster3Betas_proportion[Cluster3Betas_0to0.2] <- "low"
      Cluster3Betas_proportion[Cluster3Betas_0.21to0.6] <- "mid"
      Cluster3Betas_proportion[Cluster3Betas_0.61to1] <- "high"
      
      
      
      # Low range is between 0 and 0.2 as defined by Landau et al., 2014
      Cluster4Betas_0to0.2 <- which(data.table::between(Cluster4Betas, 
                                                        lower = 0, upper = 0.2) == TRUE)
      # Medium range is between 0.2 and 0.6 as defined by Landau et al., 2014
      Cluster4Betas_0.21to0.6 <- which(data.table::between(Cluster4Betas, 
                                                           lower = 0.20000000000001, upper = 0.6) == TRUE)
      # High range is between 0.6 and 1.0 as defined by Landau et al., 2014
      Cluster4Betas_0.61to1 <- which(data.table::between(Cluster4Betas, 
                                                         lower = 0.60000000000001, upper = 1.0) == TRUE)
      
      # Assign proportions
      Cluster4Betas_proportion <- vector(mode = "character", 
                                         length = length(Cluster4Betas))
      Cluster4Betas_proportion[Cluster4Betas_0to0.2] <- "low"
      Cluster4Betas_proportion[Cluster4Betas_0.21to0.6] <- "mid"
      Cluster4Betas_proportion[Cluster4Betas_0.61to1] <- "high"
      
      probes <- c(rownames(BetaMatrix), rownames(BetaMatrix), 
                 rownames(BetaMatrix), rownames(BetaMatrix))
      average_probes = c(Cluster1Betas, Cluster2Betas, Cluster3Betas, Cluster4Betas)
      ranges_probes = c(Cluster1Betas_proportion, 
                        Cluster2Betas_proportion, 
                        Cluster3Betas_proportion,
                        Cluster4Betas_proportion)
      class <- c(as.character(rep("Cluster1", length(Cluster1Betas))), 
                 as.character(rep("Cluster2", length(Cluster2Betas))), 
                 as.character(rep("Cluster3", length(Cluster3Betas))),
                 as.character(rep("Cluster4", length(Cluster3Betas))))
      
      average <- as.data.frame(cbind(probes, average_probes, ranges_probes, class))
      average$average_probes <- as.numeric(as.character(average$average_probes))
      
    }
    
    # If five clusters present
    if (max(unique(ClusterLabels)) == 5) {
      
      Cluster3Betas_0to0.2 <- which(data.table::between(Cluster3Betas, 
                                                        lower = 0, upper = 0.2) == TRUE)
      # Medium range is between 0.2 and 0.6 as defined by Landau et al., 2014
      Cluster3Betas_0.21to0.6 <- which(data.table::between(Cluster3Betas, 
                                                           lower = 0.20000000000001, upper = 0.6) == TRUE)
      # High range is between 0.6 and 1.0 as defined by Landau et al., 2014
      Cluster3Betas_0.61to1 <- which(data.table::between(Cluster3Betas, 
                                                         lower = 0.60000000000001, upper = 1.0) == TRUE)
      
      # Assign proportions
      Cluster3Betas_proportion <- vector(mode = "character", 
                                         length = length(Cluster3Betas))
      Cluster3Betas_proportion[Cluster3Betas_0to0.2] <- "low"
      Cluster3Betas_proportion[Cluster3Betas_0.21to0.6] <- "mid"
      Cluster3Betas_proportion[Cluster3Betas_0.61to1] <- "high"
      
      
      
      # Low range is between 0 and 0.2 as defined by Landau et al., 2014
      Cluster4Betas_0to0.2 <- which(data.table::between(Cluster4Betas, 
                                                        lower = 0, upper = 0.2) == TRUE)
      # Medium range is between 0.2 and 0.6 as defined by Landau et al., 2014
      Cluster4Betas_0.21to0.6 <- which(data.table::between(Cluster4Betas, 
                                                           lower = 0.20000000000001, upper = 0.6) == TRUE)
      # High range is between 0.6 and 1.0 as defined by Landau et al., 2014
      Cluster4Betas_0.61to1 <- which(data.table::between(Cluster4Betas, 
                                                         lower = 0.60000000000001, upper = 1.0) == TRUE)
      # Assign proportions
      Cluster4Betas_proportion <- vector(mode = "character", 
                                         length = length(Cluster4Betas))
      Cluster4Betas_proportion[Cluster4Betas_0to0.2] <- "low"
      Cluster4Betas_proportion[Cluster4Betas_0.21to0.6] <- "mid"
      Cluster4Betas_proportion[Cluster4Betas_0.61to1] <- "high"
      
      
      # Low range is between 0 and 0.2 as defined by Landau et al., 2014
      Cluster5Betas_0to0.2 <- which(data.table::between(Cluster5Betas, 
                                                        lower = 0, upper = 0.2) == TRUE)
      # Medium range is between 0.2 and 0.6 as defined by Landau et al., 2014
      Cluster5Betas_0.21to0.6 <- which(data.table::between(Cluster5Betas, 
                                                           lower = 0.20000000000001, upper = 0.6) == TRUE)
      # High range is between 0.6 and 1.0 as defined by Landau et al., 2014
      Cluster5Betas_0.61to1 <- which(data.table::between(Cluster5Betas, 
                                                         lower = 0.60000000000001, upper = 1.0) == TRUE)
      
      
      # Assign proportions
      Cluster5Betas_proportion <- vector(mode = "character", 
                                         length = length(Cluster5Betas))
      Cluster5Betas_proportion[Cluster5Betas_0to0.2] <- "low"
      Cluster5Betas_proportion[Cluster5Betas_0.21to0.6] <- "mid"
      Cluster5Betas_proportion[Cluster5Betas_0.61to1] <- "high"
      
      probes <- c(rownames(BetaMatrix), rownames(BetaMatrix), 
                 rownames(BetaMatrix), rownames(BetaMatrix),
                 rownames(BetaMatrix))
      average_probes <- c(Cluster1Betas, Cluster2Betas, Cluster3Betas, Cluster4Betas, Cluster5Betas)
      ranges_probes <- c(Cluster1Betas_proportion, 
                        Cluster2Betas_proportion, 
                        Cluster3Betas_proportion,
                        Cluster4Betas_proportion,
                        Cluster5Betas_proportion)
      class <- c(as.character(rep("Cluster1", length(Cluster1Betas))), 
                 as.character(rep("Cluster2", length(Cluster2Betas))), 
                 as.character(rep("Cluster3", length(Cluster3Betas))),
                 as.character(rep("Cluster4", length(Cluster3Betas))),
                 as.character(rep("Cluster5", length(Cluster3Betas))))
      
      average <- as.data.frame(cbind(probes, average_probes, ranges_probes, class))
      average$average_probes <- as.numeric(as.character(average$average_probes))
      
    }
    
    if (FigureGenerate == "Yes") {
      
      names(average) <- c("Probe", "Value", "Range", "Cluster")
      average$Range <- factor(average$Range, levels = c("high", "mid", "low"))
    
      
      PLOT_average <- average%>% dplyr::group_by(Cluster, Range) %>%
                      dplyr::tally()  %>%
                      ggplot2::ggplot(aes(fill = Range, y = n, x = Cluster)) + 
                      ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
                       # or:
                      # geom_bar(position = position_fill(), stat = "identity")
                      scale_y_continuous("Average methylation using Beta values", 
                                         labels = scales::percent) +
                      # scale_x_discrete(as.character(CategoryToVisualize)) +
                      # adding colours
                      scale_fill_manual(values = c("red", "black", "darkgreen")) +
                      theme_bw() + theme( panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank()) +
                      theme(aspect.ratio = 1, text = element_text(size = 15))
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Average_CLUSTER_", ImageName, ".", PNGorPDF))
    }  
    
    
    
    # Getting proportion and variance between values 0.2 and 0.8 per patient
    Proportion0.2to0.8 <- Variance0.2to0.8 <- vector() # a vector to save proportion of each patient
    for (i in 1:ncol(BetaMatrix)) {
      Proportion0.2to0.8[i] <- length(BetaMatrix[, i][between(BetaMatrix[, i], 
                                      0.2000000000000, 0.8000000000000)]) / nrow(BetaMatrix)
      Variance0.2to0.8[i]  <- var(BetaMatrix[, i][between(BetaMatrix[, i], 
                                                   0.2000000000000, 0.8000000000000)])
    }
    
    ClusterLabels_words <- ClusterLabels
    ClusterLabels_words[which(ClusterLabels_words == 1)] <- "one"
    ClusterLabels_words[which(ClusterLabels_words == 2)] <- "two"
    if (max(ClusterLabels) == 3) {
      ClusterLabels_words[which(ClusterLabels_words == 3)] <- "three"
      ClusterLabels_words = factor(ClusterLabels_words, levels = c("one", "two", "three"))
    }
    if (max(ClusterLabels) == 4) {
      ClusterLabels_words[which(ClusterLabels_words == 3)] <- "three"
      ClusterLabels_words[which(ClusterLabels_words == 4)] <- "four"
      ClusterLabels_words = factor(ClusterLabels_words, levels = c("one", "two", "three", "four"))
    }
    if (max(ClusterLabels) == 5) {
      ClusterLabels_words[which(ClusterLabels_words == 3)] <- "three"
      ClusterLabels_words[which(ClusterLabels_words == 4)] <- "four"
      ClusterLabels_words[which(ClusterLabels_words == 5)] <- "five"
      ClusterLabels_words = factor(ClusterLabels_words, levels = c("one", "two", "three", "four", "five"))
    }
    DataPatients0.2to0.8 <- data.frame(cluster = ClusterLabels_words,
                                          proportion = Proportion0.2to0.8,
                                          variance = Variance0.2to0.8)
    
    # If only two groups, do t tests
    if (length(unique(DataPatients0.2to0.8$cluster)) == 2) {
      # compare the mean of two independent groups.
      # http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r
      
      # Assumption 1: Are the two samples independent?
      # Yes, since samples from patients are not related.
      
      # Assumption 2. Are the data from each of the 2 groups follow a normal distribution?
      # Use Shapiro-Wilk normality test as described at: Normality Test in R 
      # Null hypothesis: the data are normally distributed 
      # Alternative hypothesis: the data are not normally distributed
      
      # Shapiro-Wilk normality test for group 1
      group1_pvalue <- with(DataPatients0.2to0.8, 
                            shapiro.test(proportion[cluster == unique(DataPatients0.2to0.8$cluster)[1]]))$p.value 
      # p-value = 0.3022, cannot reject null
      
      # Shapiro-Wilk normality test for group 2
      group2_pvalue <- with(DataPatients0.2to0.8, 
                            shapiro.test(proportion[cluster == unique(DataPatients0.2to0.8$cluster)[2]]))$p.value 
      # p-value = 0.6608, cannot reject null
      
      # Only proceed if both p-values are not significant, which means that the distribution of the data
      # are not significantly different from the normal distribution. Thus, can assume normality.
      if( (group1_pvalue > 0.05) && (group2_pvalue > 0.05) ) {
        
        # Assumption 3. Do the two populations have the same variances?
        res.ftest <- var.test(proportion ~ cluster, data = DataPatients0.2to0.8)$p.value
        
        # Only proceed if p-value not significant, which means there is no significant difference between 
        # the variances of the sets of data. Therefore, can use the classic t-test witch assume 
        # equality of the variances among data sets.
        
        if (res.ftest > 0.05) {
          
          # # Compute t-test 
          test_results <- t.test(proportion ~ cluster, 
                                 data = DataPatients0.2to0.8, var.equal = TRUE)
          test_pvalue <- t.test(proportion ~ cluster, 
                                data = DataPatients0.2to0.8, var.equal = TRUE)$p.value 
          cat("\n T-test performed using proportions calculated based on Beta values.\n")
        }
      } else {
        test_pvalue = NA
      }
      
    } else if (length(unique(DataPatients0.2to0.8$cluster)) > 2) {
      # If more than two groups to compare, use One-Way ANOVA Test
      # http://www.sthda.com/english/wiki/one-way-anova-test-in-r
      # This is an extension of independent two-samples t-test for comparing 
      # means in a situation where there are more than two groups
      
      # compute One-Way ANOVA 
      test_results <- summary(aov(proportion ~ cluster, data = DataPatients0.2to0.8))
      test_pvalue <- summary(aov(proportion ~ cluster, data = DataPatients0.2to0.8))[[1]][5][[1]][1]
      cat("\n One-Way ANOVA performed using proportions calculated based on Beta values.\n")
      
      # In one-way ANOVA test, a significant p-value indicates that some of the group means are different
      # To determine if the mean difference between specific pairs of group are statistically significant, 
      # perform multiple pairwise-comparison between the means of groups.
      if (test_pvalue < 0.05) {
        test_results <- stats::TukeyHSD(aov(proportion ~ cluster, data = DataPatients0.2to0.8))
      }
    }
    
    
    # Plotting
    if (FigureGenerate == "Yes" && length(unique(DataPatients0.2to0.8$cluster)) == 2) {
      
      p3 <- ggplot2::ggplot(DataPatients0.2to0.8, 
                            aes(x = cluster, y = proportion, color=factor(cluster))) +
                            scale_y_continuous(name = "Proportion of beta values between 0.2 and 0.8") +
                            geom_jitter(aes(color = cluster), size = 3, alpha = 1, width = 0.2) +
                            labs(x = "Cluster") +
                            annotate("text", x = 1.5, 
                                     y = max(DataPatients0.2to0.8$proportion) + 0.1, 
                                     label = paste("P.Value = ", round(test_pvalue,5)), size = 7) +
                            scale_color_manual(values = coloursBarPlot[1:2], 
                                               name = "Cluster" ) +
                            theme_bw() + theme(axis.title =  element_text(face = "bold"), 
                                                      aspect.ratio = 1, 
                                                      legend.position = "right", 
                                                      text = element_text(size = 20), 
                                                      panel.grid.major = element_blank(), 
                                                      panel.grid.minor = element_blank(), 
                                                      panel.background = element_rect(colour = "black", size = 1.5), 
                                                      axis.text.x = element_text(face = "bold"), 
                                                      axis.text.y = element_text(face = "bold"))
      pathNow <- getwd()
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_ggplotPoint_Cluster_", 
                    ImageName, ".", PNGorPDF))
      

      p4 <- ggpubr::ggboxplot(DataPatients0.2to0.8,
                              x = "cluster", 
                              y = "proportion", 
                              fill = "cluster",
                              add = "boxplot", 
                              palette = coloursBarPlot[1:2],
                              ylab = "Proportion between 0.2 and 0.8", 
                              font.label = list(size = 20, color = "black")) +
                              ggtitle("Cluster vs. Proportion between 0.2 and 0.8") +
                              stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                              # Add pairwise comparisons p-value
                              stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggboxplot_Proportion_Btw0.2and0.8_Boxplot_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      p5 <- ggpubr::ggscatter(DataPatients0.2to0.8, 
                              x = "cluster", y = "proportion",
                              color = "cluster", 
                              palette = coloursBarPlot[1:2],
                              ylab = "Proportion between 0.2 and 0.8") + 
                              stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                              # Add pairwise comparisons p-value
                              stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_ScatterPlot_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      
      p7 <- ggpubr::ggdensity(DataPatients0.2to0.8, x = "proportion",
                              add = "mean", rug = TRUE,
                              palette = coloursBarPlot[1:2],
                              color = "cluster", fill = "cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_DensityPlot_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      p8 <- ggpubr::ggdotchart(DataPatients0.2to0.8, 
                               x = "cluster", y = "proportion",
                               group = "cluster", color = "cluster", 
                               ylab = "Proportion between 0.2 and 0.8",
                               xlab = "Cluster",
                               sorting = "descending",
                               #ggtheme = theme_bw(),
                               palette = coloursBarPlot[1:2]) +
                               stat_compare_means(comparisons = ComparisonOptionsCluster)
                        
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_DotChart_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      p9 <- ggpubr::ggboxplot(DataPatients0.2to0.8, 
                              x = "cluster", y = "variance", 
                              fill = "cluster",
                              add = "boxplot",  
                              palette = coloursBarPlot[1:2],
                              ylab = "Variance of values between 0.2 and 0.8",
                              font.label = list(size = 20, color = "black")) +
                              ggtitle("Cluster vs. Variance") +
                              stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                              # Add pairwise comparisons p-value
                              stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggboxplot_Variance_Btw0.2and0.8_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
    } else if(
      FigureGenerate == "Yes" && length(unique(DataPatients0.2to0.8$cluster)) > 2) {
      
      p3 <- ggplot2::ggplot(DataPatients0.2to0.8, 
                            aes(x = cluster, y = proportion, color = factor(cluster))) +
                            # geom_point(size=3)+
                            theme_bw() + theme( panel.grid.major = element_blank(), 
                                                panel.grid.minor = element_blank()) +
                            scale_y_continuous(name = "Proportion of beta values between 0.2 and 0.8") +
                            geom_jitter(aes(color = cluster), size = 3, alpha = 1, width = 0.2) +
                            labs(x = paste0(ClinicalCategory)) +
                            annotate("text", x = 2, y = max(DataPatients0.2to0.8$proportion) + 0.2, 
                                     label= paste("P.Value = ",round(test_pvalue, 5)), size=7) +
                            theme(text = element_text(size=20), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), 
                                  axis.text.y = element_text(face = "bold")) +
                            scale_color_manual(values = coloursBarPlot[1:5], 
                                               name=paste0(ClinicalCategory)) +
                            theme(aspect.ratio = 1, legend.position = "right", 
                                  panel.background = element_rect(colour = "black", size=1.5),  
                                  axis.title =  element_text(face = "bold"))
      pathNow <- getwd()
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_ggplotPoint_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      p4 <- ggpubr::ggboxplot(DataPatients0.2to0.8,
                              x = "cluster", 
                              y = "proportion", 
                              fill = "cluster",
                              add = "boxplot", 
                              palette = coloursBarPlot[1:5],
                              ylab = "Proportion between 0.2 and 0.8", 
                              font.label = list(size = 20, color = "black")) +
                              ggtitle("Cluster vs. Proportion between 0.2 and 0.8") +
                              stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                              # Add pairwise comparisons p-value
                              stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggboxplot_Proportion_Btw0.2and0.8_Boxplot_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      p5 <- ggpubr::ggscatter(DataPatients0.2to0.8, 
                              x = "cluster", y = "proportion",
                              color = "cluster", 
                              palette = coloursBarPlot[1:5],
                              ylab = "Proportion between 0.2 and 0.8") + 
                              stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                              # Add pairwise comparisons p-value
                              stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_ScatterPlot_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      p7 <- ggpubr::ggdensity(DataPatients0.2to0.8, 
                              x = "proportion",
                              add = "mean", rug = TRUE,
                              palette = coloursBarPlot[1:5],
                              color = "cluster", fill = "cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_DensityPlot_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      p8 <- ggpubr::ggdotchart(DataPatients0.2to0.8, x = "cluster", y = "proportion",
                               group = "cluster", color = "cluster", 
                               ylab = "Proportion between 0.2 and 0.8",
                               xlab = "Cluster",
                               sorting = "descending",
                               #ggtheme = theme_bw(),
                               palette = coloursBarPlot[1:5]) +
                               stat_compare_means(comparisons = ComparisonOptionsCluster)
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_DotChart_Cluster_", 
                    ImageName, ".", PNGorPDF))
      
      p9 <- ggpubr::ggboxplot(DataPatients0.2to0.8, 
                              x = "cluster", 
                              y = "variance", 
                              fill = "cluster",
                              add = "boxplot",  
                              palette = coloursBarPlot[1:5],
                              ylab = "Variance of values between 0.2 and 0.8",
                              font.label = list(size = 20, color = "black")) +
                              ggtitle("Cluster vs. Variance") +
                              stat_compare_means(comparisons = ComparisonOptionsCluster) + 
                              # Add pairwise comparisons p-value
                              stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggboxplot_Variance_Btw0.2and0.8_DotChart_Cluster_", 
                    ImageName, ".", PNGorPDF))
    }
    
    
  } else {
    ############################################ #
    # If cluster labels are not provided 
    
    # Matching rownames in methylation file with corresponding location in annotation file
    match_ids_methylation_beta <- match(rownames(BetaMatrix), AnnotationFile$V1)
    
    # Bringing annotation file in the same order as BetaMatrix
    AnnotationFile2 <- data.frame(AnnotationFile[match_ids_methylation_beta, ])
    
    
    # Identify column number matching the CategoryToVisualize
    columnNumber <- as.numeric(which(colnames(AnnotationFile2) == CategoryToVisualize))
    
    # Only take in the non empty entries of the selected column
    if(CategoryToVisualize == "UCSC_RefGene_Group") {
      geneRegion <- sub("\\;.*", "", AnnotationFile2$UCSC_RefGene_Group)
      non_empty_entries <- which(geneRegion != "")
    } else {
      non_empty_entries <- which(AnnotationFile2[ , columnNumber] != "")
    }
    
    # Defining functions
    # if PlotWithinCategories is "Yes"
    eachplot_Yes <- function(AnnotationCategory, 
                             BetaMatrix, 
                             non_empty_entries, 
                             ClinicalFile, 
                             ClinicalCategory) {
      # AnnotationCategory is the category within the annotation file that is being analyzed
      # E.g. table(AnnotationFile2[non_empty_entries,columnNumber]) gives Island  N_Shore OpenSea N_Shelf S_Shelf S_Shore
      # Then set AnnotationCategory="Island"
      
      # rows corresponding to the current AnnotationCategory, e.g. Island from the list Island  N_Shore OpenSea N_Shelf S_Shelf S_Shore
      if(CategoryToVisualize == "UCSC_RefGene_Group") {
        geneRegion <- sub("\\;.*", "", AnnotationFile2$UCSC_RefGene_Group)
        correspondingRows <- which(geneRegion == AnnotationCategory)
      } else {
        correspondingRows <- which(AnnotationFile2[non_empty_entries, columnNumber] == AnnotationCategory)
      }
      
      # Initialize average_cat3
      average_cat3 <- 0
      
      # Looking at averge for all probes belonging to ClinicalCategory
      # ClinicalCategory options "SEX","STAGE","SITE_BIOPSY","TYPE_BIOPSY","COO","TRANSLOC_14_18","INSTITUTION","EPIC_QC","TYPE"
      vectorNames <- unique(na.omit(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]))
      # Saving row means for first ClinicalCategory 
      average_cat1 <- rowMeans(as.matrix(BetaMatrix[correspondingRows,
                                                    which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)), 
                                                                                 which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[1])]))
      # Saving row means for second ClinicalCategory 
      average_cat2 <- rowMeans(as.matrix(BetaMatrix[correspondingRows, 
                                                    which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)),
                                                                                 which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[2])]))
      if (length(vectorNames) > 2) {
        # Saving row means for third ClinicalCategory 
        average_cat3 <- rowMeans(as.matrix(BetaMatrix[correspondingRows, 
                                                      which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)), 
                                                                                   which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[3])]))
      }
      if (length(vectorNames) > 3) {
        stop("\nCode only support upto 3 Clinical Categories Please alter code");
      }
      
      # Looking at probe ranges
      # Low range is between 0 and 0.2 as defined by Landau et al., 2014
      # Saving row means for first ClinicalCategory 
      average_cat1_0to0.2 <- which(data.table::between(average_cat1, lower = 0, upper = 0.2) == TRUE)
      # Saving row means for second ClinicalCategory 
      average_cat2_0to0.2 <- which(data.table::between(average_cat2, lower = 0, upper = 0.2) == TRUE)
      
      # Medium range is between 0.2 and 0.6 as defined by Landau et al., 2014
      average_cat1_0.21to0.6 <- which(data.table::between(average_cat1, 
                                                          lower = 0.20000000000001, 
                                                          upper = 0.6) == TRUE)
      average_cat2_0.21to0.6 <- which(data.table::between(average_cat2, 
                                                          lower = 0.20000000000001, 
                                                          upper = 0.6) == TRUE)
      
      # High range is between 0.6 and 1.0 as defined by Landau et al., 2014
      average_early_0.61to1 <- which(data.table::between(average_cat1, 
                                                         lower = 0.60000000000001, 
                                                         upper = 1.0) == TRUE)
      average_cat2_0.61to1 <- which(data.table::between(average_cat2, 
                                                        lower = 0.60000000000001, 
                                                        upper = 1.0) == TRUE)
      
      # Assign proportions
      average_cat1_proportion <- vector(mode = "character", 
                                        length = length(average_cat1))
      average_cat1_proportion[average_cat1_0to0.2] <- "low"
      average_cat1_proportion[average_cat1_0.21to0.6] <- "mid"
      average_cat1_proportion[average_early_0.61to1] <- "high"
      # Length(average_DMR_lateNo_0to0.2)+length(average_DMR_lateNo_0.21to0.6)+length(average_DMR_lateNo_0.61to1) # 14661; agrees!
      average_cat2_proportion <- vector(mode = "character", 
                                        length = length(average_cat2))
      average_cat2_proportion[average_cat2_0to0.2] <- "low"
      average_cat2_proportion[average_cat2_0.21to0.6] <- "mid"
      average_cat2_proportion[average_cat2_0.61to1] <- "high"
      
      # Drawing bar graph for islands
      probes <- c(names(average_cat1), names(average_cat2))
      average_probes <- c(average_cat1, average_cat2)
      ranges_probes <- c(average_cat1_proportion, average_cat2_proportion)
      class <- c(as.character(rep(vectorNames[1], length(average_cat1))), 
                 as.character(rep(vectorNames[2], length(average_cat2))))
      
      # Forming a data frame with all the info 
      average <- as.data.frame(cbind(probes, 
                                     as.numeric(average_probes), 
                                     ranges_probes, class))
      
      if (length(vectorNames) == 3) {
        # Low range is between 0 and 0.2 as defined by Landau et al., 2014
        average_cat3_0to0.2 <- which(data.table::between(average_cat3, lower = 0, 
                                                         upper = 0.2) == TRUE)
        # Medium range is between 0.2 and 0.6 as defined by Landau et al., 2014
        average_cat3_0.21to0.6 <- which(data.table::between(average_cat3, 
                                                            lower = 0.20000000000001, 
                                                            upper = 0.6) == TRUE)
        # High range is between 0.6 and 1.0 as defined by Landau et al., 2014
        average_cat3_0.61to1 <- which(data.table::between(average_cat3, 
                                                          lower = 0.60000000000001, 
                                                          upper = 1.0) == TRUE)
        
        # Assign proportions
        average_cat3_proportion <- vector(mode = "character", length = length(average_cat3))
        average_cat3_proportion[average_cat3_0to0.2] <- "low"
        average_cat3_proportion[average_cat3_0.21to0.6] <- "mid"
        average_cat3_proportion[average_cat3_0.61to1] <- "high"
        
        # Drawing bar graph for islands
        probes <- c(names(average_cat1), names(average_cat2), names(average_cat3))
        average_probes <- c(as.numeric(average_cat1), as.numeric(average_cat2), 
                            as.numeric(average_cat3))
        ranges_probes <- c(average_cat1_proportion, average_cat2_proportion, 
                           average_cat3_proportion)
        class <- c(as.character(rep(vectorNames[1], length(average_cat1))), 
                   as.character(rep(vectorNames[2], length(average_cat2))), 
                   as.character(rep(vectorNames[3], length(average_cat3))) )
        
        average <- as.data.frame(cbind(probes, as.numeric(average_probes), 
                                       ranges_probes, class))
       
      }
      
      names(average) <- c("Probe", "Value", "Range", "Class")
      average$Value <- as.numeric(as.character(average$Value))
      # Releveling the range
      average$Range <- factor(average$Range, levels = c("high", "mid", "low"))
      # cat("\n Printing",as.character(AnnotationCategory), "plot \n")
      if (FigureGenerate == "Yes") {
        
        # Violin plot 
       # calculating values across probes corresponding the the variable
       colMeanVariable1 <- colMeans(as.matrix(BetaMatrix[correspondingRows, 
                                                         which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)), 
                                                                                      which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[1])]))
       colMeanVariable2 <- colMeans(as.matrix(BetaMatrix[correspondingRows, 
                                                         which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)), 
                                                                                      which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[2])]))
       
       colMedianVariable1 <- colMedians(as.matrix(BetaMatrix[correspondingRows, 
                                                         which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)), 
                                                                                      which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[1])]))
       
       colMedianVariable2 <- colMedians(as.matrix(BetaMatrix[correspondingRows, 
                                                         which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)), 
                                                                                      which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[2])]))
       
       
       colMeanTable <- data.frame(cbind(Class = c(rep(vectorNames[1], length(colMeanVariable1)), 
                                                  rep(vectorNames[2], length(colMeanVariable2))) , 
                                        Mean = c(colMeanVariable1, colMeanVariable2),
                                        Median = c(colMedianVariable1, colMedianVariable2)))
       colMeanTable$Mean <- as.numeric(as.character(colMeanTable$Mean))
       colMeanTable$Median <- as.numeric(as.character(colMeanTable$Median))
       
       if (length(vectorNames) == 3) {
         colMeanVariable3 <- colMeans(as.matrix(BetaMatrix[correspondingRows, 
                                                           which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)), 
                                                                                        which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[3])]))
         
         colMedianVariable3 <- colMedians(as.matrix(BetaMatrix[correspondingRows, 
                                                           which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)), 
                                                                                        which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[3])]))
         
         colMeanTable <- data.frame(cbind (
                                    Class = c(rep(vectorNames[1], length(colMeanVariable1)), 
                                              rep(vectorNames[2], length(colMeanVariable2)), 
                                              rep(vectorNames[3], length(colMeanVariable3))), 
                                    Mean = c(colMeanVariable1, colMeanVariable2, colMeanVariable3),
                                    Median = c(colMedianVariable1, colMedianVariable2, colMedianVariable3)))
         colMeanTable$Mean <- as.numeric(as.character(colMeanTable$Mean))
         colMeanTable$Median <- as.numeric(as.character(colMeanTable$Median))
         
       }
       
       # Setting the number of comparisons for ggpubr
       if(length(unique(vectorNames)) == 2) {
         ComparisonOptions <- list(names(table(vectorNames))[1:2])
       } else if(length(unique(vectorNames)) == 3) {
         ComparisonOptions <- list(names(table(vectorNames))[1:2], 
                                   names(table(vectorNames))[2:3],
                                   names(table(vectorNames))[c(1, 3)])
       } else if(length(unique(vectorNames)) == 4) {
         ComparisonOptions <- list( names(table(vectorNames))[1:2], 
                                    names(table(vectorNames))[c(1, 3)], 
                                    names(table(vectorNames))[c(1, 4)], 
                                    names(table(vectorNames))[c(2, 3)], 
                                    names(table(vectorNames))[c(2, 4)], 
                                    names(table(vectorNames))[c(3, 4)])
       }
       
       
       if(ClinicalCategory == "STAGE") {
         colour = c("#762a83", "#c2a5cf")
       } else if (ClinicalCategory == "TYPE") {
         colour = c("#d6604d", "#66bd63", "#4575b4")
       } else if (ClinicalCategory == "SEX") { 
         colour = c("#b35806", "#fdb863")
       } else if (ClinicalCategory == "CLUSTER") {
         colour = c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                    '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                    '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                    '#000075', '#808080')
       } else if (ClinicalCategory == "TRANSLOC_14_18") { 
         colour = c("#e0e0e0", "#878787")
       } else if (ClinicalCategory == "TYPE_BIOPSY") { 
         colour = c("#a6dba0", "#878787")
       } else if (ClinicalCategory == "SITE_BIOPSY") { 
         colour = c("#f1b6da", "#c51b7d")
       }


       

       ViolinPlot1 <- ggpubr::ggboxplot(colMeanTable, x = "Class", y = "Mean", fill = "Class",
                                          add = "boxplot", 
                                          ylab = " Mean Beta Value", 
                                          title = paste0(AnnotationCategory),
                                          palette = colour) + 
                                          stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                                          stat_compare_means(paired = FALSE)
       ggsave(paste0(pathNow,"/img/21_ggboxplot_MeanBetaValueAcrossProbes_ANOVA_", 
                     AnnotationCategory, "_", ClinicalCategory, ".", PNGorPDF))

       
       ViolinPlot2 <- ggpubr::ggboxplot(colMeanTable, x = "Class", y = "Median", fill = "Class",
                                        add = "boxplot", 
                                        ylab = " Median Beta Value", 
                                        title = paste0(AnnotationCategory),
                                        palette = colour) + 
                                        stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                                        stat_compare_means(paired = FALSE)
       ggplot2::ggsave(paste0(pathNow,"/img/21_ggboxplot_MedianBetaValueAcrossProbes_ANOVA_", 
                     AnnotationCategory, "_", ClinicalCategory, ".", PNGorPDF))
      
        
        # Bar plots
        PLOT_average <- average %>% dplyr::group_by(Class, Range) %>%
                        dplyr::tally()  %>%
                        ggplot(aes(fill = Range, y = n, x = Class)) + 
                        ggplot2::geom_bar(position="fill", stat="identity", width = 0.7) +
                        # or:
                        # geom_bar(position = position_fill(), stat = "identity")
                        scale_y_continuous("Average methylation using Beta values", 
                                           labels = scales::percent) +
                        scale_x_discrete(as.character(AnnotationCategory)) +
                        # adding colours
                        scale_fill_manual(values = c("red", "black", "darkgreen")) +
                        theme_bw() + 
                        theme(axis.title =  element_text(face = "bold"), 
                              aspect.ratio = 1, 
                              legend.position = "right", 
                              text = element_text(size = 20), 
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(), 
                              panel.background = element_rect(colour = "black", size=1.5), 
                              axis.text.x = element_text(face = "bold"), 
                              axis.text.y = element_text(face = "bold"))
        ggsave(paste0(pathNow, "/img/21_Proportion_Average_", AnnotationCategory, 
                      "_", ClinicalCategory, ".", PNGorPDF))
      }
    }
    
    # if PlotWithinCategories is "No"
    eachplot_No <- function(CategoryToVisualize, 
                            BetaMatrix, 
                            non_empty_entries, 
                            ClinicalFile, 
                            ClinicalCategory) {
      
      # initialize average_cat3
      average_cat3 <- 0
      
      # looking at averge for all probes belonging to DMR  (average based on # of patietns)
      vectorNames <- unique(na.omit(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]))
      average_cat1 <- rowMeans(as.matrix(BetaMatrix[, which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)),
                                                                                   which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[1])]))
      average_cat2 <- rowMeans(as.matrix(BetaMatrix[, which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)),
                                                                                   which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[2])]))
      if (length(vectorNames) > 2) {
        average_cat3 <- rowMeans(as.matrix(BetaMatrix[, which(stri_trim(ClinicalFile[c(which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)),
                                                                                    which(colnames(ClinicalFile) == ClinicalCategory)]) == vectorNames[3])]))
      }
      
      # Looking at probes
      # Assign Proportions
      # Low range is between 0 and 0.2 as defined by Landau et al., 2014
      average_cat1_0to0.2 <- which(data.table::between(average_cat1, 
                                                       lower = 0, 
                                                       upper = 0.2) == TRUE)
      average_cat2_0to0.2 <- which(data.table::between(average_cat2, 
                                                       lower = 0, 
                                                       upper = 0.2) == TRUE)
      # Medium range is between 0.2 and 0.6 as defined by Landau et al., 2014
      average_cat1_0.21to0.6 <- which(data.table::between(average_cat1, 
                                                          lower = 0.20000000000001, 
                                                          upper = 0.6) == TRUE)
      average_cat2_0.21to0.6 <- which(data.table::between(average_cat2, 
                                                          lower = 0.20000000000001, 
                                                          upper = 0.6) == TRUE)
      # High range is between 0.6 and 1.0 as defined by Landau et al., 2014
      average_cat1_0.61to1 <- which(data.table::between(average_cat1, 
                                                        lower = 0.60000000000001, 
                                                        upper = 1.0) == TRUE)
      average_cat2_0.61to1 <- which(data.table::between(average_cat2, 
                                                        lower = 0.60000000000001, 
                                                        upper = 1.0) == TRUE)
      
      # Assign proportions
      average_cat1_proportion <- vector(mode = "character", 
                                        length = length(average_cat1))
      average_cat1_proportion[average_cat1_0to0.2] <- "low"
      average_cat1_proportion[average_cat1_0.21to0.6] <- "mid"
      average_cat1_proportion[average_cat1_0.61to1] <- "high"
      
      average_cat2_proportion <- vector(mode = "character", 
                                        length = length(average_cat2))
      average_cat2_proportion[average_cat2_0to0.2] <- "low"
      average_cat2_proportion[average_cat2_0.21to0.6] <- "mid"
      average_cat2_proportion[average_cat2_0.61to1] <- "high"
      
      # Drawing bar graph
      probes <- c(names(average_cat1),names(average_cat2))
      average_probes <- c(average_cat1,average_cat2)
      ranges_probes <- c(average_cat1_proportion, average_cat2_proportion)
      class <- c(as.character(rep(vectorNames[1], length(average_cat1))), 
                 as.character(rep(vectorNames[2], length(average_cat2))))
      
      average <- as.data.frame(cbind(probes, average_probes, ranges_probes, class))
      
      # If three variables present
      if (length(vectorNames) == 3) {
        # Low range is between 0 and 0.2 as defined by Landau et al., 2014
        average_cat3_0to0.2 <- which(data.table::between(average_cat3, 
                                                         lower = 0, 
                                                         upper = 0.2) == TRUE)
        # Medium range is between 0.2 and 0.6 as defined by Landau et al., 2014
        average_cat3_0.21to0.6 <- which(data.table::between(average_cat3, 
                                                            lower = 0.20000000000001, 
                                                            upper = 0.6) == TRUE)
        # High range is between 0.6 and 1.0 as defined by Landau et al., 2014
        average_cat3_0.61to1 <- which(data.table::between(average_cat3, 
                                                          lower = 0.60000000000001, 
                                                          upper = 1.0) == TRUE)
        
        # Assign proportions
        average_cat3_proportion <- vector(mode = "character", length = length(average_cat3))
        average_cat3_proportion[average_cat3_0to0.2] <- "low"
        average_cat3_proportion[average_cat3_0.21to0.6] <- "mid"
        average_cat3_proportion[average_cat3_0.61to1] <- "high"
        
        probes = c(names(average_cat1), names(average_cat2), names(average_cat3))
        average_probes = c(average_cat1, average_cat2, average_cat3)
        ranges_probes = c(average_cat1_proportion, average_cat2_proportion, 
                          average_cat3_proportion)
        class <- c(as.character(rep(vectorNames[1], 
                                    length(average_cat1))), 
                   as.character(rep(vectorNames[2], 
                                    length(average_cat2))), 
                   as.character(rep(vectorNames[3], 
                                    length(average_cat3))) )
        
        average <- as.data.frame(cbind(probes, average_probes, ranges_probes, class))
        
      }
      
      if (FigureGenerate == "Yes") {
        names(average) <- c("Probe", "Value", "Range", "Class")
        average$Range <- factor(average$Range, levels = c("high", "mid", "low"))
        
        PLOT_average <- average %>% dplyr::group_by(Class, Range) %>%
                        dplyr::tally()  %>%
                        ggplot(aes(fill = Range, y = n, x = Class)) + 
                        geom_bar(position="fill", stat="identity", width = 0.7) +
                        # or:
                        # geom_bar(position = position_fill(), stat = "identity")
                        scale_y_continuous("Average methylation using Beta values", 
                                           labels = function(x) paste0(x*100, "%")) +
                        scale_x_discrete(as.character(CategoryToVisualize)) +
                        # adding colours
                        scale_fill_manual(values = c("red","black", "darkgreen")) +
                        theme_bw() + 
                        theme(axis.title =  element_text(face = "bold"), 
                              aspect.ratio = 1, 
                              legend.position = "right", 
                              text = element_text(size = 20), 
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(), 
                              panel.background = element_rect(colour = "black", size=1.5), 
                              axis.text.x = element_text(face = "bold"), 
                              axis.text.y = element_text(face = "bold"))
        ggsave(paste0(pathNow, "/img/21_Proportion_Average_", ClinicalCategory, ".", PNGorPDF))
      }  
      
    }
    
    # Plotting
    if (PlotWithinCategories == "Yes" && FigureGenerate == "Yes") {
      if(CategoryToVisualize == "UCSC_RefGene_Group"){
        geneRegion <- sub("\\;.*", "", AnnotationFile2$UCSC_RefGene_Group)
        CategoriesToPlot <- unique(geneRegion)[which(unique(geneRegion) != "")]
      } else {
        CategoriesToPlot <- unique(AnnotationFile2[non_empty_entries, columnNumber])
      }
      
      plots <- lapply(CategoriesToPlot, eachplot_Yes, BetaMatrix = BetaMatrix,
                      non_empty_entries = non_empty_entries,
                      ClinicalFile = ClinicalFile,
                      ClinicalCategory = ClinicalCategory)
    } else if (PlotWithinCategories == "No" && FigureGenerate == "Yes") {
      plots <- eachplot_No(CategoryToVisualize = CategoryToVisualize,
                           BetaMatrix = BetaMatrix,
                           non_empty_entries = non_empty_entries,
                           ClinicalFile = ClinicalFile,
                           ClinicalCategory = ClinicalCategory)
    }
    
    ####
    
    # Not run for now
    ProportionBetween0.2AND0.8 <- function() {
      

    # Getting proportion between 0.2 and 0.8 per patient
    Proportion0.2to0.8 <- Variance0.2to0.8 <- vector() # a vector to save proportion of each patient
    for (i in 1:ncol(BetaMatrix)) {
      Proportion0.2to0.8[i] <- length(BetaMatrix[ , i][between(BetaMatrix[ , i], 
                                                               0.2000000000000, 0.8000000000000)]) / nrow(BetaMatrix)
      Variance0.2to0.8[i]  <- var(BetaMatrix[ , i][between(BetaMatrix[ , i], 
                                                          0.2000000000000, 0.8000000000000)])
      
    }
    
    if (length(which(is.na(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)) > 0) {
      DataPatients0.2to0.8 <- 
        data.frame(names = ClinicalFile[- which(is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE), 
                                                                 which(colnames(ClinicalFile) == ClinicalCategory)],
                   proportion = Proportion0.2to0.8[- which(is.na(ClinicalFile[, which(colnames(ClinicalFile) == ClinicalCategory)]) == TRUE)])
    } else {
      DataPatients0.2to0.8 <- 
        data.frame(names = ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategory)],
                   proportion = Proportion0.2to0.8)
    }
    
    ####
    # Running t-tests
    
    # If only two groups, do t tests
    if (length(unique(DataPatients0.2to0.8$names)) == 2) {
      # compare the mean of two independent groups.
      # http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r
      
      # Assumption 1: Are the two samples independent?
      # Yes, since samples from patients are not related.
      
      # Assumption 2. Are the data from each of the 2 groups follow a normal distribution?
      # Use Shapiro-Wilk normality test as described at: Normality Test in R 
      # Null hypothesis: the data are normally distributed 
      # Alternative hypothesis: the data are not normally distributed
      
      # Shapiro-Wilk normality test for group 1
      group1_pvalue <- with(DataPatients0.2to0.8, 
                            shapiro.test(proportion[names == 
                                                      unique(DataPatients0.2to0.8$names)[1]]))$p.value 
      # p-value = 0.3022, cannot reject null
      
      # Shapiro-Wilk normality test for group 2
      group2_pvalue <- with(DataPatients0.2to0.8, 
                            shapiro.test(proportion[names == 
                                                      unique(DataPatients0.2to0.8$names)[2]]))$p.value 
      # p-value = 0.6608, cannot reject null
      
      # Only proceed if both p-values are not significant, which means that the distribution of the data
      # are not significantly different from the normal distribution. Thus, can assume normality.
      if( (group1_pvalue > 0.05) && (group2_pvalue > 0.05) ) {
        
        # Assumption 3. Do the two populations have the same variances?
        res.ftest <- var.test(proportion ~ names, data = DataPatients0.2to0.8)$p.value
        
        # Only proceed if p-value not significant, which means there is no significant difference between 
        # the variances of the sets of data. Therefore, can use the classic t-test witch assume 
        # equality of the variances among data sets.
        
        if (res.ftest > 0.05) {
          
          # # Compute t-test 
          test_results <- t.test(proportion ~ names, 
                                 data = DataPatients0.2to0.8, var.equal = TRUE)
          test_pvalue <- t.test(proportion ~ names, 
                                data = DataPatients0.2to0.8, var.equal = TRUE)$p.value 
          cat("\n T-test performed using proportions calculated based on Beta values.\n")
        }
      }
    } else if (length(unique(DataPatients0.2to0.8$names)) > 2) {
      # If more than two groups to compare, use One-Way ANOVA Test
      # http://www.sthda.com/english/wiki/one-way-anova-test-in-r
      # This is an extension of independent two-samples t-test for comparing means in a situation where there are more than two groups
      
      # compute One-Way ANOVA 
      test_results <- summary(aov(proportion ~ names, data = DataPatients0.2to0.8))
      test_pvalue <- summary(aov(proportion ~ names, data = DataPatients0.2to0.8))[[1]][5][[1]][1]
      cat("\n One-Way ANOVA performed using proportions calculated based on Beta values.\n")
      
      # In one-way ANOVA test, a significant p-value indicates that some of the group means are different
      # To determine if the mean difference between specific pairs of group are statistically significant, 
      # perform multiple pairwise-comparison between the means of groups.
      if (test_pvalue < 0.05) {
        test_results <- stats::TukeyHSD(aov(proportion ~ names, data = DataPatients0.2to0.8))
      }
    }
    
    # Plotting
    if (FigureGenerate == "Yes" && length(unique(DataPatients0.2to0.8$names)) == 2) {
      
      p3 <- ggplot2::ggplot(DataPatients0.2to0.8, aes(x = names, y = proportion, color=factor(names))) +
                            # geom_point(size=3)+
                            theme_bw() + theme( panel.grid.major = element_blank(), 
                                                panel.grid.minor = element_blank()) +
                            scale_y_continuous(name = "Proportion of beta values between 0.2 and 0.8") +
                            geom_jitter(aes(color = names), size = 3, alpha = 1, width = 0.2) +
                            labs(x = paste0(ClinicalCategory)) +
                            annotate("text", x = 1.5, y = max(DataPatients0.2to0.8$proportion) + 0.1, 
                                     label = paste("P.Value = ", round(test_pvalue, 5)), size = 7) +
                            theme(text = element_text(size = 20), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), 
                                  axis.text.x = element_text(face = "bold"), 
                                  axis.text.y = element_text(face = "bold") ) +
                            scale_color_manual(values = c("#fddbc7", "#b2182b", "red"), 
                                               name = paste0(ClinicalCategory)) +
                            theme(aspect.ratio = 1, legend.position = "right", 
                                  panel.background = element_rect(colour = "black", size = 1.5),  
                                  axis.title =  element_text(face = "bold"))
                          
      pathNow <- getwd()
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_ggplotPoint_", 
                    ClinicalCategory, ".", PNGorPDF))
      
      # Setting the number of comparisons for ggpubr
      ComparisonOptions <- list(names(table(DataPatients0.2to0.8$names))[1:2])
      
      p4 <- ggpubr::ggboxplot(DataPatients0.2to0.8, 
                              x = "names", y = "proportion",
                              color = "names", palette = "jco")+ 
                              stat_compare_means(comparisons = ComparisonOptions) + 
                              # Add pairwise comparisons p-value
                              stat_compare_means(label.y = 0.5, paired = FALSE)     
                             # Add global p-value
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggboxplot_Proportion_Btw0.2and0.8_Boxplot_",
                    ClinicalCategory, ".", PNGorPDF))
      
      p5 <- ggscatter(DataPatients0.2to0.8, x = "names", y = "proportion",
                      color = "names", palette = "jco")+ 
                      stat_compare_means(comparisons = ComparisonOptions) + # Add pairwise comparisons p-value
                      stat_compare_means(label.y = 0.5, paired = FALSE)
                    
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_ScatterPlot_",
                    ClinicalCategory, ".", PNGorPDF))
      
      
      p7 <- ggdensity(DataPatients0.2to0.8, x = "proportion",
                      add = "mean", rug = TRUE,
                      color = "names", fill = "names")
      
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_DensityPlot_",
                    ClinicalCategory, ".", PNGorPDF))
      
      p8 <- ggdotchart(DataPatients0.2to0.8, x = "names", y = "proportion",
                       group = "names", color = "names",
                       sorting = "descending",
                       ggtheme = theme_bw(),
                       y.text.col = TRUE )
      
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_DotChart_",
                    ClinicalCategory, ".", PNGorPDF))
      
    } else if(FigureGenerate == "Yes" && length(unique(DataPatients0.2to0.8$names)) > 2) {
      
      p3 <- ggplot2::ggplot(DataPatients0.2to0.8, aes(x = names, y = proportion, color = factor(names))) +
                            # geom_point(size=3)+
                            theme_bw() + theme( panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank()) +
                            scale_y_continuous(name = "Proportion of beta values between 0.2 and 0.8") +
                            geom_jitter(aes(color = names), size = 3, alpha = 1, width = 0.2) +
                            labs(x = paste0(ClinicalCategory)) +
                            annotate("text", x = 2, y = max(DataPatients0.2to0.8$proportion) + 0.2, 
                                     label= paste("P.Value = ",round(test_pvalue, 5)), size = 7) +
                            theme(text = element_text(size = 20), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), 
                                  axis.text.x = element_text(face = "bold"), 
                                  axis.text.y = element_text(face = "bold") ) +
                            scale_color_manual(values=c("#fddbc7", "#b2182b", "red"), 
                                               name = paste0(ClinicalCategory)) +
                            theme(aspect.ratio = 1, legend.position = "right", 
                                  panel.background = element_rect(colour = "black", size = 1.5), 
                                  axis.title =  element_text(face = "bold"))
                          
      pathNow <- getwd()
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_ggplotPoint_", 
                    ClinicalCategory, ".", PNGorPDF))
      
      # Setting the number of comparisons for ggpubr
      # https://www.r-bloggers.com/add-p-values-and-significance-levels-to-ggplots/
      if(length(unique(DataPatients0.2to0.8$names)) == 3) {
        ComparisonOptions <- list(names(table(DataPatients0.2to0.8$names))[1:2],
                               names(table(DataPatients0.2to0.8$names))[2:3],
                               names(table(DataPatients0.2to0.8$names))[c(1, 3)])
      } else if(length(unique(DataPatients0.2to0.8$names)) == 4) {
        ComparisonOptions <- list( names(table(DataPatients0.2to0.8$names))[1:2], 
                                names(table(DataPatients0.2to0.8$names))[c(1, 3)], 
                                names(table(DataPatients0.2to0.8$names))[c(1, 4)], 
                                names(table(DataPatients0.2to0.8$names))[c(2, 3)], 
                                names(table(DataPatients0.2to0.8$names))[c(2, 4)], 
                                names(table(DataPatients0.2to0.8$names))[c(3, 4)])
      }
      
      p4 <- ggpubr::ggboxplot(DataPatients0.2to0.8, x = "names", y = "proportion",
                      color = "names", palette = "jco")+ 
                      stat_compare_means(comparisons = ComparisonOptions)+ # Add pairwise comparisons p-value
                      stat_compare_means(paired = FALSE)     # Add global p-value
      ggplot2::ggsave(paste0(pathNow, "/img/21_ggboxplot_Proportion_Btw0.2and0.8_Boxplot_", 
                    ClinicalCategory, ".", PNGorPDF))
      
      p5 <- ggpubr::ggscatter(DataPatients0.2to0.8, x = "names", y = "proportion",
                      color = "names", palette = "jco")+ 
                      stat_compare_means(comparisons = ComparisonOptions)+ # Add pairwise comparisons p-value
                      stat_compare_means(paired = FALSE)
                    
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_ScatterPlot_", 
                    ClinicalCategory, ".", PNGorPDF))
      
      p7 <- ggdensity(DataPatients0.2to0.8, x = "proportion",
                      add = "mean", rug = TRUE,
                      color = "names", fill = "names")
      
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_DensityPlot_", 
                    ClinicalCategory, ".", PNGorPDF))
      
      p8 <- ggdotchart(DataPatients0.2to0.8, x = "names", y = "proportion",
                       group = "names", color = "names",
                       sorting = "descending",
                       ggtheme = theme_bw(),
                       y.text.col = TRUE )
      ggplot2::ggsave(paste0(pathNow, "/img/21_Proportion_Btw0.2and0.8_DotChart_", 
                    ClinicalCategory, ".", PNGorPDF))
      }
    
    RESULTS <- list(Proportion0.35and0.6DataFrame = DataPatients0.2to0.8,
                    TestResults = test_results)
    
    class(RESULTS) <- "ProportionVisualization_ASilva"
    return(RESULTS)
    }
    
  }
  return(NULL)
}




