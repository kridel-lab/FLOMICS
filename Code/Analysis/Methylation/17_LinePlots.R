# Updated 3 Aug 2021
# Updated 26 April 2019
# Function: Generate Line Plot (still under construction) of provided values and mean (in yellow).
# Author: Anjali Silva

# Input:
# Dataset: Numericmatrix of values for probes/ ENSEMBL IDs x patients. Smaller the better. 
# ClusterMembershipVector: vector of length(patients) indicating the cluster membership
# FigureGenerate: Produce images or not, options = "Yes" or "No"; default "Yes" 
# PNGorPDF: Output format of the image, options = "png" or "pdf"; default "png"


LinePlot <- function(Dataset, 
                     ClusterMembershipVector, 
                     FigureGenerate = "Yes",
                     PNGorPDF = "png") {
  
  # Checking/ Loading needed packages
  # LoadCheckPkg(RegularPckgs=c("pheatmap","gplots","RColorBrewer","MASS"))
  library(pheatmap)
  library(gplots)
  library(RColorBrewer)
  library(MASS)
  
  # Obtaining path to save images
  pathNow <- getwd()
  
  # Saving cluster membership for each observation
  DataPlusLabs <- rbind(Dataset, ClusterMembershipVector)
  ordervector <- anothervector <- list()
  
  for (i in 1:max(ClusterMembershipVector)) {
    ordervector[[i]] <- which(DataPlusLabs[nrow(Dataset) + 1, ] == i)
    anothervector[[i]] <- rep(i, length(which(DataPlusLabs[nrow(Dataset) + 1, ] == i)))
  }
  
  vec <- unlist(ordervector)
  colorsvector <- unlist(anothervector)
  
  # Define colours
  coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                      '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                      '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                      '#000075', '#808080')

  
  # Line Plots
  if (PNGorPDF == "png") {
    grDevices::png(paste0(pathNow, "/img/17_LinePlots.", PNGorPDF))
  }
  if (PNGorPDF == "pdf") {
    grDevices::pdf(paste0(pathNow, "/img/17_LinePlots.", PNGorPDF))
  }
  
  # Layout
  if (length(unique(ClusterMembershipVector)) == 2){
    graphics::par(mfrow = c(2, 1)) # setting space for two plots
  } else{
    graphics::par(mfrow = c(2, length(unique(ClusterMembershipVector))))
  }

  
  for(cluster in unique(ClusterMembershipVector)) {
    # Save how many observations below to each cluster size, given by 'cluster'
    toplot1 <- t(as.matrix(DataPlusLabs[c(1:nrow(Dataset)),
                                         which(DataPlusLabs[nrow(Dataset) + 1, ] == cluster)], 
                            ncol = ncol(Dataset)))
    # Save column mean in last row
    toplot1 <- cbind(toplot1, rowMeans(toplot1 + 1))
    

    graphics::matplot(toplot1, type = "l", pch = 1,
                      col = c(rep(coloursBarPlot[cluster], nrow(toplot1)), 7),
                      xlab = "Samples", ylab = "Value", cex = 1,
                      lty = c(rep(2, nrow(toplot1)), 1),
                      lwd = c(rep(3, nrow(toplot1)), 4),
                      xaxt = "n", xlim = c(1, ncol(toplot1)),
                      main = paste("Cluster ", cluster))
    axis(1, at = c(1:(nrow(toplot1) - 1)), labels = rownames(toplot1)[- (nrow(toplot1))])
    grDevices::dev.off()
  }
  
  
  return(NULL)
}
# [END]
