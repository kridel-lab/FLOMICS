# Updated 12 Feb 2019
# Function: Create the heat plot using beta values
# Author: Anjali Silva

# Input:
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns
# ClinicalFile: File of patients and corresponding clinical category; matrix of size patients x clinical categories
# CategoryToVisualize: Category to visualize from the Clinical File;
#                      Options: "SAMPLE_ID", SAMPLE_ID_TGL", "INCLUDE", "TIME_POINT", "RES_ID","LY_FL_ID", "OTHER_ID", "SEX","SITE_BIOPSY"          
#                      "SITE_EXTRANODAL", "TYPE_BIOPSY", "TRANSLOC_14_18", "INSTITUTION", "TYPE", "STAGE", "COO", "EPIC_QC"
# PNGorPDF: Output format of the image, options = "png" or "pdf"

# Output: None 

# Visuals saved to img folder
# 6_HeatPlot_.p*

HeatPlot <- function(BetaMatrix, 
                     ClinicalFile, 
                     CategoryToVisualize, 
                     PNGorPDF) {
  
  # Loading needed packages
  # LoadCheckPkg(RegularPckgs=c("tidyverse","ggplot2"))
  library(tidyverse)
  library(ggplot2)
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  # Obtaining the unique categories within CategoryToVisualize
  NamesWithinCategoryToVisualize <- names(table(ClinicalFile[ , which(colnames(ClinicalFile) == CategoryToVisualize)]))
  
  # Obtaining rowMeans from BetaMatrix
  BetaMatrix_ProbeMeans <- as.data.frame(sapply(c(1:length(NamesWithinCategoryToVisualize)), 
                                                function(x) rowMeans(BetaMatrix[ , which(ClinicalFile[ , which(colnames(ClinicalFile) == CategoryToVisualize)] 
                                                                                         == NamesWithinCategoryToVisualize[x])])))
  colnames(BetaMatrix_ProbeMeans) <- NamesWithinCategoryToVisualize
  
  # Remove NA values
  if(length(which(is.na(BetaMatrix_ProbeMeans) == TRUE)) > 0) {
    BetaMatrix_ProbeMeans <- BetaMatrix_ProbeMeans[which(is.na(BetaMatrix_ProbeMeans) == TRUE), ]
  }

  # Convert aes names to symbol with sym (from rlang) and do the evaluation
  # https://stackoverflow.com/questions/50960339/create-ggplot2-function-and-specify-arguments-as-variables-in-data-as-per-ggplot
  x.var <- rlang::sym(NamesWithinCategoryToVisualize[1])
  y.var <- rlang::sym(NamesWithinCategoryToVisualize[2])
  
  # Plotting based on NamesWithinCategoryToVisualize
  HeatPlot <- ggplot(BetaMatrix_ProbeMeans, aes(x = !! x.var, y= !! y.var)) +
                     stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
                     scale_fill_distiller(palette= "Spectral", direction=1) +
                     scale_x_continuous(expand = c(0, 0)) +
                     scale_y_continuous(expand = c(0, 0)) +
                     theme(legend.position='none') +
                     geom_abline(slope=1,intercept=0) +
                     theme(aspect.ratio=1)
      ggsave(paste0(pathNow, "/img/11_HeatPlot_", CategoryToVisualize, ".",PNGorPDF))
      
  if(length(NamesWithinCategoryToVisualize) == 3) {
    
    # Convert aes names to symbol with sym (from rlang) and do the evaluation
    z.var <- rlang::sym(NamesWithinCategoryToVisualize[3])
    
    # Plot
    HeatPlot2 <- ggplot(BetaMatrix_ProbeMeans, aes(x = !! x.var, y = !! z.var)) +
                        stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
                        scale_fill_distiller(palette = "Spectral", direction = 1) +
                        scale_x_continuous(expand = c(0, 0)) +
                        scale_y_continuous(expand = c(0, 0)) +
                        theme(legend.position='none') +
                        geom_abline(slope = 1, intercept = 0) +
                        theme(aspect.ratio = 1)
      ggsave(paste0(pathNow, "/img/11_HeatPlot_", CategoryToVisualize,"_2.", PNGorPDF))
      
    HeatPlot3 <- ggplot(BetaMatrix_ProbeMeans, aes(x = !! y.var, y = !! z.var)) +
                        stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
                        scale_fill_distiller(palette = "Spectral", direction=1) +
                        scale_x_continuous(expand = c(0, 0)) +
                        scale_y_continuous(expand = c(0, 0)) +
                        theme(legend.position = 'none') +
                        geom_abline(slope = 1, intercept = 0) +
                        theme(aspect.ratio = 1)
      ggsave(paste0(pathNow, "/img/11_HeatPlot_", CategoryToVisualize, "_3.", PNGorPDF))
  }
  
  if(length(NamesWithinCategoryToVisualize)>3){
    cat("\n*** Warning: The Argument CategoryToVisualize =",CategoryToVisualize ,"has more than 3 categories.
        \nFunction can only intake in upto 3 different categories within CategoryToVisualize.\n")
  }
      # Developed by Anjali Silva
  return(NULL)
}



