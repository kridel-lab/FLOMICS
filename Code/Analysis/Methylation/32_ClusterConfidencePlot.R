# Updated 10 Feb 2020
# Function: Plot the probabilities of clusters for each sample. 
#           Only support upto 5 cluster labels.
# Author: Anjali Silva

# Input:
# ProbabilityMatrix: A matrix of probabilities of size patients x clusters, with patients 
#                    in rows and each cluster along each column. Each row should sum to 1.
# FigureGenerate: Character string specifying if image should be produced or not. 
#                 Options = "Yes" or "No"; default "Yes" 
# PNGorPDF: Output format of the image, options = "png" or "pdf"
# ImageName: Character string specifying a unique name for the image to be saved. 


# Output: 

# Visuals saved to img folder
# 31_MeanSDPlot_ImageName.p*

ClusterConfidencePlot <- function(ProbabilityMatrix, 
                                  FigureGenerate = "Yes", 
                                  PNGorPDF ="png", 
                                  ImageName = "Sample1") {
    
  # Loading needed packages
  # LoadCheckPkg(RegularPckgs = c("dplyr","ggplot2","magrittr"))
  library(dplyr)
  library(ggplot2)
  library(magrittr)
  library(mclust)
  library(reshape2)
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  # checking 
  if (! all.equal(sum(rowSums(ProbabilityMatrix)), nrow(ProbabilityMatrix))) {
    stop("\n The probabilities may not add upto 1 for each row of ProbabilityMatrix.")
  }
  if (ncol(ProbabilityMatrix) > 5) {
    stop("\n 32_ClusterConfidencePlot.R only supports upto 5 cluster labels. 
         Update the code.")
  }
  
   
  TableProbabilities <- data.frame(as.numeric(mclust::map(ProbabilityMatrix)), 
                                    round(as.matrix(ProbabilityMatrix), 5))
  
  TableProbabilities <- cbind(colnames(BetaMatrix_T1), TableProbabilities)
  
  names(TableProbabilities) <- c("Sample", "Cluster", 
                                 paste0("Probability of ", c(1:ncol(ProbabilityMatrix))))
  TableProbabilities$Cluster <- factor(TableProbabilities$Cluster)
  
  if (length(unique(mclust::map(ProbabilityMatrix))) == 5) {
    for (i in 3:7) {
      TableProbabilities[, i] <- factor(TableProbabilities[, i]) 
    }
    
    TableProbabilitiesMelt <- reshape::melt(TableProbabilities, 
                                            id.vars = c("Sample", "Cluster"), 
                                            measure.vars = c(3:7))
  
  } else if (length(unique(mclust::map(ProbabilityMatrix))) == 4) {
    for (i in 3:6) {
      TableProbabilities[, i] <- factor(TableProbabilities[, i]) 
    }
    
    TableProbabilitiesMelt <- reshape::melt(TableProbabilities, 
                                            id.vars = c("Sample", "Cluster"), 
                                            measure.vars = c(3:6))
    
  } else if (length(unique(mclust::map(ProbabilityMatrix))) == 3) {
    for (i in 3:5) {
      TableProbabilities[, i] <- factor(TableProbabilities[, i]) 
    }
    
    TableProbabilitiesMelt <- reshape::melt(TableProbabilities, 
                                            id.vars = c("Sample", "Cluster"), 
                                            measure.vars = c(3:5))
    
  } else if (length(unique(mclust::map(ProbabilityMatrix))) == 2) {
    for (i in 3:4) {
      TableProbabilities[, i] <- factor(TableProbabilities[, i]) 
    }
    
    TableProbabilitiesMelt <- reshape::melt(TableProbabilities, 
                                            id.vars = c("Sample", "Cluster"), 
                                            measure.vars = c(3:4))
  } else {
    stop("\n 32_ClusterConfidencePlot only supports up to 5 groups in ClusterLabels.
         Update the code.")
  }
  
  
  TableProbabilitiesMelt$value <- as.numeric(levels(TableProbabilitiesMelt$value))[TableProbabilitiesMelt$value]
  
  
  # Define colours
  coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                      '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                      '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                      '#000075', '#808080')

  # Stacked + percent
  ggplot2::ggplot(TableProbabilitiesMelt, aes(fill = variable, y = value, x = Sample,)) +
                  ggplot2::geom_bar(position = "fill", stat = "identity") +
                  ggplot2::theme(axis.text.x = element_text(angle = 90)) + 
                  ggplot2::coord_cartesian(ylim = c(0, 1)) +
                  labs(fill = "Cluster") +
                  scale_fill_manual(values = coloursBarPlot[sort(unique(TableProbabilitiesMelt$Cluster))]) +
                  # labs(y = "Posterior probability") +
                  # y axis tick mark lables
                  ggplot2::scale_y_continuous(name = "Probability", limits = c(0: 1)) 
  ggplot2::ggsave(paste0(pathNow, "/img/32_ClusterConfidencePlot_", ImageName, ".", PNGorPDF))
  
  
  
  # Plotting median beta value for patients in each cluster, seperated by stage
  # Setting the number of comparisons for ggpubr
  if(length(unique(Cluster)) == 2) {
    ComparisonOptionsClusters <- list(names(table(TableProbabilitiesMelt$Cluster))[1:2])
  } 
  
  p0_cluster_stage <- ggpubr::ggboxplot(TableProbabilitiesMelt, 
                                        x = "Cluster", 
                                        y = "value", 
                                        fill = "variable",
                                        add = "boxplot", 
                                        xlab = "Cluster", 
                                        ylab = "Probability",
                                        palette =  coloursBarPlot[sort(unique(TableProbabilitiesMelt$Cluster))]) + 
                                        stat_compare_means(comparisons = ComparisonOptionsClusters) + 
                                        # Add pairwise comparisons p-value
                                        stat_compare_means(paired = FALSE)
  ggplot2::ggsave(paste0(pathNow, "/img/32_ClusterConfidencePlot2_", ImageName, ".", PNGorPDF))
  
   
  return(NULL) # Only the plots are returned
}
  


