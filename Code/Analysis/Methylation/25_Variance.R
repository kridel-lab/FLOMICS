# Updated 9 September 2019
# Function: Visualize varaince of data between Beta values 0.2 to 0.8 falling within a category. 
#           If PlotWithinCategories is set to "Yes", then categories within CategoryToVisualize 
#           will be ploted. Otherwise, only probes corresponding to CategoryToVisualize will be
#           plotted. 
# Author: Anjali Silva

# Input:
# CategoryToVisualize: Specifies which column needs to be visualized from AnnotationFile.
#                      Need exact spelling and caps/noCaps for CategoryToVisualize, e.g.: 
#                      CategoryToVisualize = "DMR". Options "Relation_to_Island", "chr", 
#                     "Regulatory_Feature_Group", "X450k_Enhancer", "Phantom4_Enhancers", 
#                     "Phantom5_Enhancers".
# ClusterLabels: If proportions are being visualized for clusters, provide a vector of integers 
#                indicating cluster membership.
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns
# AnnotationFile: Matrix of annotations for all the probes found in BetaMatrix. It is of size 
#                 probes x annotations.
# ClinicalFile:  File of patients and corresponding clinical category; matrix of size patients x clinical categories.
# ClinicalCategory: Specify column of ClinicalFile to use. Options "STAGE", "SEX", "SITE_BIOPSY",
#                   "TYPE_BIOPSY", "INSTITUTION", "COO", "TYPE", "TRANSLOC_14_18".
# PlotWithinCategories: Should indicate "Yes" or "No", as to whether categories within each category 
#                       to visualize (CategoryToVisualize) will to be investigated.
# PNGorPDF: Output format of the image, options = "png" or "pdf".

# Output: 
# Proportion0.35and0.6: A data frame of proportion of data values that fell between 0.35 and 0.6 for each person. 

# Visuals saved to img folder
# 25_Varaince_*",ClinicalCategory.p*

Varaince <- function(CategoryToVisualize = NA, 
                     ClusterLabels = NA, 
                     BetaMatrix, 
                     AnnotationFile = NA, 
                     ClinicalFile = NA, 
                     ClinicalCategory = NA, 
                     PlotWithinCategories = NA, 
                     FigureGenerate = "No", 
                     PNGorPDF = "png", 
                     ImageName, 
                     Lower = 0.2, 
                     Upper = 0.8) {
  
  # Loading needed packages
  library(ggplot2)
  library(stringi)
  library(data.table)
  library(ggpubr)
  
  ImageName <- paste0(ImageName, "_range_", Lower, "_to_", Upper)
  
  # Remove empty entries in the clinical file
  ClinicalFile[ClinicalFile == ""] <- NA
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  # If cluster labels are provided 
  if (sum(is.na(ClusterLabels)) == 0 ) {
    
    # all probes and samples only corresponding to those specified in ClinicalCategoryToVisualize
    Cluster1Betas <- as.matrix(BetaMatrix[ , ClusterLabels == 1])
    Cluster2Betas <- as.matrix(BetaMatrix[ , ClusterLabels == 2])
    if (max(unique(ClusterLabels)) > 2) {
      Cluster3Betas <- as.matrix(BetaMatrix[ , ClusterLabels == 3])
    }
    if (max(unique(ClusterLabels)) > 3) {
      stop("\nCode only support upto 3 clusters. Please alter code");
    }
    
    Cluster1_Var <- sapply(c(1:ncol(Cluster1Betas)),
                           function(a) var(Cluster1Betas[ data.table::between(Cluster1Betas[ , a], 
                                                                              lower = Lower, upper = Upper), a]))
    Cluster2_Var <- sapply(c(1:ncol(Cluster2Betas)), 
                           function(a) var(Cluster2Betas[ data.table::between(Cluster2Betas[ , a], 
                                                                              lower = Lower, upper = Upper), a]))
    if (max(unique(ClusterLabels)) > 2) {
      Cluster3_Var <- sapply(c(1:ncol(Cluster3Betas)),
                             function(a) var(Cluster3Betas[ data.table::between(Cluster3Betas[ , a], 
                                                                                lower = Lower, upper = Upper), a]))
    }
    
  # Plotting variability
  if (FigureGenerate == "Yes" && max(unique(ClusterLabels)) == 2) {
 
    data_patient_var <- data.frame(names = c(rep("one", length(Cluster1_Var)), 
                                            rep("two", length(Cluster2_Var)),
                                            rep("three", length(Cluster3_Var))),
                                   variance = c(Cluster1_Var, Cluster2_Var, Cluster3_Var))
    data_patient_var$names <- factor(data_patient_var$names, levels = c("one", "two", "three"))
    Variance_plot <- ggplot2::ggplot(data_patient_var, aes(x = names, y = as.numeric(variance), color = factor(names))) +
                            #geom_point(size=3)+
                            theme_bw() + 
                            scale_y_continuous(name = "Variance") +
                            geom_jitter(aes(color = names), size = 3, alpha = 1, width = 0.2) +
                            labs(x = "Cluster") +
                            #annotate("text", x = 1.5, y = max(data_patient_var$variance) + 0.01, size = 7) +
                            theme(text = element_text(size = 20), 
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), 
                                  axis.text.x = element_text(face = "bold"),
                                  axis.text.y = element_text(face = "bold")) +
                            scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
                            theme(aspect.ratio = 1, 
                                  legend.position = "right", 
                                  panel.background = element_rect(colour = "black", size = 1.5),  
                                  axis.title =  element_text(face = "bold")) 
    Variance_plot + ggtitle("Variance based on Beta between 0.3 & 0.7, per sample") 
    ggplot2::ggsave(paste0(pathNow, "/img/25_ggplotPoint_Varaince_Cluster_", ImageName, ".", PNGorPDF))
    
    # Setting the number of comparisons for ggpubr
    my_comparisons <- list(c("one","two"), c("one","three"), c("two","three"))
    p0_violin_median <- ggpubr::ggviolin(data_patient_var, x = "names", y = "variance", fill = "names",
                                 add = "boxplot") + 
                                 stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
                                 stat_compare_means( paired = FALSE)
    ggplot2::ggsave(paste0(pathNow, "/img/25_ggplotPoint_Variance_Violin_Cluster_", ImageName, ".", PNGorPDF))

  }
    
    
  } else {
    stop("\nCode only support cluster results for now. Please alter code");
  }
    
    RESULTS <- list(Variance0.2and0.8DataFrame = data_patient_var)
    
    class(RESULTS) <- "Variance_ASilva"
    return(RESULTS)
}
  





