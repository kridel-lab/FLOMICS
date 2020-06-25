# Updated 26 Feb 2019
# Function: Generate tSNE plots
# Author: Anjali Silva

# Input:
# tSNEOutput: Output produced from running Rtsne() function from Rtsne package
# Perplexity: Value used for perplexity parameter
# tSNEDataFrame: The dataframe used to produce the ggplot of tSNE output
# ClusterLabels: Provide cluster labels from another method for comparison purposes
# Filter: If to filter by patient type, options = "FL".
# FigureGenerate: Produce images or not, options = "Yes" or "No"; default "Yes" 
# PNGorPDF: Output format of the image, options = "png" or "pdf"; default "png"
# ImageName: Unique name for the image

# Output
# Visuals saved to img folder
# tSNE results 13_tSNE_perplexity=X.p*

tSNEPlotGeneration <- function(BetaMatrix, 
                               PerplexityParameter, 
                               ClinicalFile,
                               ClusterLabels = NA, 
                               Filter = NA, 
                               FigureGenerate = "Yes", 
                               PNGorPDF = "png", 
                               ImageName) {
    
  # Code adapted based on Alberto (OICR) Sent 25 Feb 2019
  
  # Loading needed packages
  # LoadCheckPkg(RegularPckgs=c("ggplot2","Rtsne","tidyr","plotly"))
  library(ggplot2)
  library(Rtsne)
  library(tidyr)
  library(plotly)

  set.seed(1234)
  # T-distributed Stochastic Neighbor Embedding implementation
  tsne_out <- Rtsne::Rtsne(t(BetaMatrix),
                           perplexity = PerplexityParameter,
                           check_duplicates = F)
  
  facna_stage <- addNA(ClinicalFile$STAGE)
  levels(facna_stage) <- c(levels(facna_stage), "NoInformation")
  if(is.na(Filter) == FALSE) {
    facna_stage <- addNA(ClinicalFile$STAGE[which(substr(colnames(BetaMatrix), 4, 5) == Filter)])
    levels(facna_stage) <- c(levels(facna_stage), "NoInformation")
  }
  
  facna_sex <- addNA(ClinicalFile$SEX)
  levels(facna_sex) <- c(levels(facna_sex), "NoInformation")
  if(is.na(Filter) == FALSE) {
    facna_sex <- addNA(ClinicalFile$SEX[which(substr(colnames(BetaMatrix), 4, 5) == Filter)])
    levels(facna_sex) <- c(levels(facna_sex), "NoInformation")
  }
  
  facna_translocation <- addNA(ClinicalFile$TRANSLOC_14_18)
  levels(facna_translocation) <- c(levels(factor(facna_translocation)), "NoInformation")
  if(is.na(Filter) == FALSE) {
    facna_translocation <- addNA(ClinicalFile$TRANSLOC_14_18)[which(substr(colnames(BetaMatrix), 4, 5) == Filter)]
    levels(facna_translocation) <- c(levels(factor(facna_translocation)), "NoInformation")
  }
  
  df.tsne <- data.frame(colnames(BetaMatrix), 
                        as.data.frame(tsne_out$Y), 
                        factor(facna_stage), 
                        factor(facna_sex), 
                        factor(facna_translocation))
  colnames(df.tsne) <- c("samples", "D1", "D2", "Stage", "Sex", "Translocation")
  
  if(is.na(ClusterLabels)[1] != TRUE) {
    df.tsne <- data.frame(colnames(BetaMatrix), 
                          as.data.frame(tsne_out$Y), 
                          factor(facna_stage), 
                          factor(facna_sex), 
                          factor(facna_translocation), 
                          factor(ClusterLabels))
    colnames(df.tsne) <- c("samples", "D1", "D2", "Stage", "Sex", "Translocation", "ClusterLabels")
  }

  df.tsne$samples <- as.character(df.tsne$samples)
  df.tsne$samples[grepl("FL", df.tsne$sample)] <- "FL"
  df.tsne$samples[grepl("DLC", df.tsne$sample)] <- "DLC"
  df.tsne$samples[grepl("RLN", df.tsne$sample)] <- "RLN"
  
  if(FigureGenerate == "Yes") {
    # Getting the path of the file, which should contain a folder called "img"
    pathNow <- getwd()
    
    ggplot2::ggplot(as.data.frame(df.tsne), 
                    aes(x = D1, 
                        y = D2, 
                        color = samples, 
                        sample = samples)) + 
                    geom_point(size = 2) +
                    theme_bw() + 
                    theme(panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank()) + 
                    scale_color_manual(values = c("black", "black", "black")) + 
                    theme(aspect.ratio = 1, text = element_text(size = 15))
    
    # Getting the path of the file, which should contain a folder called "img"
    pathNow <- getwd()
    ggplot2::ggsave(paste0(pathNow, "/img/13_tSNE_none_perplexity=",
                           PerplexityParameter, "_", ImageName, ".", PNGorPDF))
    
    
    ggplot2::ggplot(as.data.frame(df.tsne), aes(x = D1, 
                                                y = D2, 
                                                color = samples,
                                                sample = samples)) + 
                  geom_point(size = 2) +
                  theme_bw() + 
                  theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank()) + 
                  scale_color_manual(values = c("#8c510a", "#a6dba0", "#9970ab")) + 
                  theme(aspect.ratio = 1, text = element_text(size = 15))
    
    ggplot2::ggsave(paste0(pathNow, "/img/13_tSNE_TYPE_perplexity=",
                           PerplexityParameter,"_", ImageName, ".", PNGorPDF))
    
    
    
    ggplot2::ggplot(as.data.frame(df.tsne), 
                    aes(x = D1, 
                        y = D2, 
                        color = Stage, 
                        sample = samples)) + 
                    geom_point(size = 2) +
                    theme_bw() + 
                    theme(panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank()) + 
                    scale_color_manual(values = c("#c51b7d", "#a6dba0", "#4d4d4d")) + 
                    theme(aspect.ratio = 1, text = element_text(size = 15))
    ggplot2::ggsave(paste0(pathNow,"/img/13_tSNE_STAGE_perplexity=",
                           PerplexityParameter,"_", ImageName, ".", PNGorPDF))
    
  
    ggplot2::ggplot(as.data.frame(df.tsne),
                    aes(x = D1, 
                        y = D2, 
                        color = Sex, 
                        sample = samples)) + 
                    geom_point(size = 2) +
                    theme_bw() + 
                    theme(panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank()) + 
                    scale_color_manual(values = c("#c51b7d", "#a6dba0", "#4d4d4d")) + 
                    theme(aspect.ratio = 1, text = element_text(size = 15))
    ggplot2::ggsave(paste0(pathNow,"/img/13_tSNE_SEX_perplexity=",
                           PerplexityParameter,"_",ImageName,".",PNGorPDF))
    
    
    ggplot2::ggplot(as.data.frame(df.tsne), 
                    aes(x = D1, 
                        y = D2, 
                        color = Translocation, 
                        sample = samples)) + 
                    geom_point(size = 2) +
                    theme_bw() + 
                    theme(panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank()) + 
                    scale_color_manual(values = c("#c51b7d", "#a6dba0", "#4d4d4d")) + 
                    theme(aspect.ratio = 1, text = element_text(size = 15))
    ggplot2::ggsave(paste0(pathNow,"/img/13_tSNE_TRANSLOCATION_perplexity=",
                  PerplexityParameter,"_",ImageName,".",PNGorPDF))
    
    
    # Define colours:
    coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                        '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080')
    
    if (length(levels(ClusterLabels)) == 2 || max(ClusterLabels) == 2) {
      # Cluster <- relevel(Cluster, ref = "2")
      # setting colour palette to be two colours
      colourPalette = coloursBarPlot[1:2]
      # provide reference group as cluster three
    } else if (length(levels(ClusterLabels)) == 3 || max(ClusterLabels) == 3) {
      # Cluster[which(Cluster == 3)] <- "three"
      # setting colour palette to be three colours
      colourPalette = coloursBarPlot[1:3]
      # Cluster <- factor(Cluster, levels = c("one", "three", "two"))
      # Cluster <- relevel(Cluster, ref = "three")
      #Cluster <- relevel(Cluster, ref = "3")
    } else if (length(levels(ClusterLabels)) == 4 || max(ClusterLabels) == 4) {
      # Cluster[which(Cluster == 4)] <- "four"
      # setting colour palette to be four colours
      colourPalette = coloursBarPlot[1:4]
      # Cluster <- factor(Cluster, levels = c("one", "two", "three", "four"))
      # Cluster <- relevel(Cluster, ref = "four")
      # Cluster <- relevel(Cluster, ref = "4")
    } else if (length(levels(ClusterLabels)) == 5 || max(ClusterLabels) == 5) {
      colourPalette = coloursBarPlot[1:5]
    } else {
      stop("\n Only upto 6 clusters are supported.")
    }
    
    if(is.na(ClusterLabels)[1] != TRUE) {
      ggplot2::ggplot(as.data.frame(df.tsne), 
                      aes(x = D1, 
                          y = D2, 
                          color = ClusterLabels, 
                          sample = ClusterLabels)) + 
                      geom_point(size = 2) +
                      theme_bw() + 
                      theme(panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank()) + 
                      scale_color_manual(values = colourPalette) + 
                      theme(aspect.ratio = 1, text = element_text(size = 15))
      ggplot2::ggsave(paste0(pathNow, "/img/13_tSNE_ExternalClusterLabs_perplexity=",
                             PerplexityParameter, "_", ImageName, ".", PNGorPDF))
    }
  }
  
  RESULTS <- list(tSNEOutput = tsne_out,
                  Perplexity = PerplexityParameter,
                  tSNEDataFrame = as.data.frame(df.tsne))
  
  class(RESULTS) <- "tSNEPlotGeneration_ASilva"
  return(RESULTS)
}




