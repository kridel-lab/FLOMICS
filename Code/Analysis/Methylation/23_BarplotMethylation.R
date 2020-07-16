# Date: 29 May 2019
# Function: Generate bar plots of sample number given clinical category. Cluster labels maybe provided
#           by adding CLUSTER column to ClinicalFile. 
# Author: Anjali Silva

# Input
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns.
# ClinicalFile: File of patients and corresponding clinical category; matrix of size patients x clinical categories.
# ClinicalCategory: Specify column of ClinicalFile to use. Options "STAGE", "SEX", "SITE_BIOPSY", "TYPE_BIOPSY",
#                   "INSTITUTION", "COO", "TYPE", "TRANSLOC_14_18".

# Output: None, plots are saved
# 23_ggplotBar_Mean_.*

BarplotMethylation <- function(BetaMatrix, 
                               ClinicalFile, 
                               ClinicalCategory, 
                               FigureGenerate = "Yes", 
                               PNGorPDF = "png") {

  # Loading needed packages
  # LoadCheckPkg(RegularPckgs = c("ggplot2","stringi","data.table","ggpubr"))
  # "stringi" used to remove extra white space, if present, e.g. "EN " to "EN"
  library(ggplot2)
  library(stringi)
  library(data.table)
  library(ggpubr)
  library(grDevices)
  
  # Remove empty entries in the clinical file
  ClinicalFile[ClinicalFile == ""] <- NA
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  dataframeToPlot <- data.frame(values = colMeans(BetaMatrix),
                                category = as.character(factor(unlist(ClinicalFile[ , which(colnames(ClinicalFile) == ClinicalCategory)]))))
  # remove NA values if present
  if(length(which(is.na(dataframeToPlot$category)) == TRUE) > 0) {
    dataframeToPlot <- dataframeToPlot[! is.na(dataframeToPlot$category), ]
  }
  if(FigureGenerate == "Yes") {
    if(PNGorPDF == "pdf") {   
      pdf(paste0(pathNow, "/img/23_Barplot_", ClinicalCategory,".pdf"), 
          width = 60, height = 60, pointsize = 50)
    } else {
      png(paste0(pathNow, "/img/23_Barplot_", ClinicalCategory,".png"))
    }
    
    # Define colours
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
    
    
    
    p <- dataframeToPlot %>% dplyr::group_by(category) %>%
          dplyr::tally() %>%
          ggplot(aes(x = category, y = n, fill = category)) +
          geom_bar(stat = "identity", width = 0.7, fill = colour) +
          labs(y = "Number of samples", x = paste(ClinicalCategory)) +
          theme_bw() + theme( panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank()) +
          theme(aspect.ratio = 1, text = element_text(size = 15))
    
    
    ggplot2::ggsave(paste0(pathNow, "/img/23_ggplotBar_Mean_", 
                           ClinicalCategory, ".", PNGorPDF))
    grDevices::dev.off()
  }
  
  
  return(NULL)
}
