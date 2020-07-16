# Updated 10 Feb 2020
# Function: Plot the mean vs standard deviation values for probes seperated by cluster labels.
#           Only support upto 5 cluster labels.
# Author: Introduced by Robert Kridel and adapted to a function by Anjali Silva.

# Input:
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients 
#             as columns. The number of probes should (typically) equal to the number of probes
#             that was used for clustering. The number of patients MUST equal 
#               length(ClusterLabels).               
# ClinicalFile: File of patients and corresponding clinical category; matrix of size 
#               patients x clinical categories. The number of patients MUST equal 
#               length(ClusterLabels).
# ClusterLabels: A vector of integers indicating cluster membership. Only support upto
#                5 cluster labels.
# FigureGenerate: Character string specifying if image should be produced or not. 
#                 Options = "Yes" or "No"; default "Yes" 
# PNGorPDF: Output format of the image, options = "png" or "pdf"
# ImageName: Character string specifying a unique name for the image to be saved. 


# Output: 

# Visuals saved to img folder
# 31_MeanSDPlot_ImageName.p*

MeanSDPlot <- function(BetaMatrix, 
                       ClinicalFile, 
                       ClusterLabels, 
                       FigureGenerate = "Yes", 
                       PNGorPDF ="png", 
                       ImageName = "Sample1") {
    
  # Loading needed packages
  # LoadCheckPkg(RegularPckgs = c("dplyr","ggplot2","magrittr"))
  library(dplyr)
  library(ggplot2)
  library(magrittr)
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  if (max(ClusterLabels) > 5) {
    stop("\n 31_MeanSDPlot.R only supports upto 5 cluster labels. 
         Update the code.")
  }
  

  B.matrix.sd.mean1 <- t(BetaMatrix[, which(ClusterLabels == 1)]) %>%
                       data.frame() %>%
                       dplyr::mutate(SAMPLE_ID = row.names(.)) %>%
                       dplyr::left_join(ClinicalFile[,c("SAMPLE_ID", "TYPE")], by = "SAMPLE_ID") %>%
                       tidyr::pivot_longer(cols = -c(SAMPLE_ID, TYPE), names_to = "cg") %>%
                       dplyr::group_by(cg, TYPE) %>%
                       dplyr::summarize(mean = mean(value), sd = sd(value)) %>%
                       dplyr::mutate(matrix = "C1")
                    
  B.matrix.sd.mean2 <- t(BetaMatrix[, which(ClusterLabels == 2)]) %>%
                       data.frame() %>%
                       dplyr::mutate(SAMPLE_ID = row.names(.)) %>%
                       dplyr::left_join(ClinicalFile[,c("SAMPLE_ID", "TYPE")], by = "SAMPLE_ID") %>%
                       tidyr::pivot_longer(cols = -c(SAMPLE_ID, TYPE), names_to = "cg") %>%
                       dplyr::group_by(cg, TYPE) %>%
                       dplyr::summarize(mean = mean(value), sd = sd(value)) %>%
                       dplyr::mutate(matrix = "C2")
  
  if (length(unique(ClusterLabels )) == 2) {
    
    CombinedMatrix <- B.matrix.sd.mean1 %>%
                      rbind(B.matrix.sd.mean2) %>%
                      dplyr::mutate(matrix = factor(matrix, levels = c("C1", "C2")))
  } else if (length(unique(ClusterLabels)) == 3) {
    
    B.matrix.sd.mean3 <- t(BetaMatrix[, which(ClusterLabels == 3)]) %>%
                         data.frame() %>%
                         dplyr::mutate(SAMPLE_ID = row.names(.)) %>%
                         dplyr::left_join(ClinicalFile[,c("SAMPLE_ID", "TYPE")], by = "SAMPLE_ID") %>%
                         tidyr::pivot_longer(cols = -c(SAMPLE_ID, TYPE), names_to = "cg") %>%
                         dplyr::group_by(cg, TYPE) %>%
                         dplyr::summarize(mean = mean(value), sd = sd(value)) %>%
                         dplyr::mutate(matrix = "C3")
    
    CombinedMatrix <- B.matrix.sd.mean1 %>%
                      rbind(B.matrix.sd.mean2, B.matrix.sd.mean3) %>%
                      dplyr::mutate(matrix = factor(matrix, levels = c("C1", "C2", "C3")))
  } else if (length(unique(ClusterLabels)) == 4) {
    
    B.matrix.sd.mean3 <- t(BetaMatrix[, which(ClusterLabels == 3)]) %>%
      data.frame() %>%
      dplyr::mutate(SAMPLE_ID = row.names(.)) %>%
      dplyr::left_join(ClinicalFile[,c("SAMPLE_ID", "TYPE")], by = "SAMPLE_ID") %>%
      tidyr::pivot_longer(cols = -c(SAMPLE_ID, TYPE), names_to = "cg") %>%
      dplyr::group_by(cg, TYPE) %>%
      dplyr::summarize(mean = mean(value), sd = sd(value)) %>%
      dplyr::mutate(matrix = "C3")
    
    B.matrix.sd.mean4 <- t(BetaMatrix[, which(ClusterLabels == 4)]) %>%
                         data.frame() %>%
                         dplyr::mutate(SAMPLE_ID = row.names(.)) %>%
                         dplyr::left_join(ClinicalFile[,c("SAMPLE_ID", "TYPE")], by = "SAMPLE_ID") %>%
                         tidyr::pivot_longer(cols = -c(SAMPLE_ID, TYPE), names_to = "cg") %>%
                         dplyr::group_by(cg, TYPE) %>%
                         dplyr::summarize(mean = mean(value), sd = sd(value)) %>%
                         dplyr::mutate(matrix = "C4")
    
    CombinedMatrix <- B.matrix.sd.mean1 %>%
                      rbind(B.matrix.sd.mean2, B.matrix.sd.mean3, B.matrix.sd.mean4) %>%
                      dplyr::mutate(matrix = factor(matrix, levels = c("C1", "C2", "C3", "C4")))
  } else if (length(unique(ClusterLabels)) == 5) {
    
    B.matrix.sd.mean3 <- t(BetaMatrix[, which(ClusterLabels == 3)]) %>%
      data.frame() %>%
      dplyr::mutate(SAMPLE_ID = row.names(.)) %>%
      dplyr::left_join(ClinicalFile[,c("SAMPLE_ID", "TYPE")], by = "SAMPLE_ID") %>%
      tidyr::pivot_longer(cols = -c(SAMPLE_ID, TYPE), names_to = "cg") %>%
      dplyr::group_by(cg, TYPE) %>%
      dplyr::summarize(mean = mean(value), sd = sd(value)) %>%
      dplyr::mutate(matrix = "C3")
    
    B.matrix.sd.mean4 <- t(BetaMatrix[, which(ClusterLabels == 4)]) %>%
      data.frame() %>%
      dplyr::mutate(SAMPLE_ID = row.names(.)) %>%
      dplyr::left_join(ClinicalFile[,c("SAMPLE_ID", "TYPE")], by = "SAMPLE_ID") %>%
      tidyr::pivot_longer(cols = -c(SAMPLE_ID, TYPE), names_to = "cg") %>%
      dplyr::group_by(cg, TYPE) %>%
      dplyr::summarize(mean = mean(value), sd = sd(value)) %>%
      dplyr::mutate(matrix = "C4")
    
    B.matrix.sd.mean5 <- t(BetaMatrix[, which(ClusterLabels == 5)]) %>%
      data.frame() %>%
      dplyr::mutate(SAMPLE_ID = row.names(.)) %>%
      dplyr::left_join(ClinicalFile[,c("SAMPLE_ID", "TYPE")], by = "SAMPLE_ID") %>%
      tidyr::pivot_longer(cols = -c(SAMPLE_ID, TYPE), names_to = "cg") %>%
      dplyr::group_by(cg, TYPE) %>%
      dplyr::summarize(mean = mean(value), sd = sd(value)) %>%
      dplyr::mutate(matrix = "C5")
    
    CombinedMatrix <- B.matrix.sd.mean1 %>%
      rbind(B.matrix.sd.mean2, B.matrix.sd.mean3, B.matrix.sd.mean4, B.matrix.sd.mean5) %>%
      dplyr::mutate(matrix = factor(matrix, levels = c("C1", "C2", "C3", "C4", "C5")))
  } else {
    stop("\n 31_MeanSDPlot only supports up to 5 groups in ClusterLabels. Update the code.")
  }
  
  if (FigureGenerate == "Yes") {
    CombinedMatrix %>%
      ggplot2::ggplot(aes(x = mean, y = sd)) +
      ggplot2::geom_point(alpha = 0.1) +
      ggplot2::facet_grid(cols = vars(TYPE), rows = vars(matrix)) +
      ggplot2::guides(fill = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     text = element_text(size = 15)) 
    ggplot2::ggsave(paste0(pathNow, "/img/31_MeanSDPlot_", ImageName, ".", PNGorPDF))
  }
  
  return(NULL) # Only the plots are returned
}
  


