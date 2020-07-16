# 24 Feb 2020
# Function: Function to implement similarity network fusion (SNF) clustering. 
# Author: Anjali Silva


# Input:
# RNAseqCountMatrix: A matrix of counts (type "integer") size features (e.g., genes) x patients, with 
#                    features in rows and patients as columns. The ncol(RNAseqCountMatrix)
#                    should match ncol(MvalueMatrix), ncol(BetaMatrix), nrow(ClinicalFile) and 
#                    nrow(QCMatrix).
# MvalueMatrix: A matrix of M values (type "double") size probes x patients, with 
#               probes in rows and patients as columns. The ncol(MvalueMatrix)
#               should match ncol(RNAseqCountMatrix), ncol(BetaMatrix), nrow(ClinicalFile) and 
#               nrow(QCMatrix).
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns.
#            The ncol(BetaMatrix) should match ncol(RNAseqCountMatrix), ncol(MvalueMatrix), nrow(ClinicalFile)
#            and nrow(QCMatrix).
# AnnotationFile: A matrix of annotations for all the probes found in MvalueMatrix/BetaMatrix. It is of 
#                 size probes x annotations.
# QCMatrix: A data frame of QC metrics, that has size patients x QC metric categories. The
#           nrow(QCMatrix) should match ncol(RNAseqCountMatrix).
# TumorPurity: A vector of tumor purity values of length ncol(BetaMatrix), ncol(RNAseqCountMatrix), 
#              ncol(MvalueMatrix), nrow(ClinicalFile) and nrow(QCMatrix).
# ClinicalFile: A data frame of clinical categories, that has size patients x clinical categories. The
#               nrow(ClinicalFile) should match ncol(RNAseqCountMatrix), ncol(MvalueMatrix), ncol(BetaMatrix),
#               and nrow(QCMatrix).
# SurvivalFile: A data frame of survival data,  with size patients x survival categories.
# NumberOfClusters: An integer, less than 6, specifying the number of clusters for which plots should be 
#                   generated.
# ImageName: A character vector specifying a unique name for the images to be produced. 
# PNGorPDF: Output format of the image, options = "png" or "pdf". Default: png.



# Output
# NumberOfClusters: NumberOfClusters defined by the user. 
# SNFClusterLabels: A data frame of size samples x 2, where the two columns define "Cluster" and 
#                   "SampleID". Cluster will contain a vector of integers specifying the cluster of
#                    of each sample defined in "SampleID".
# RNASeqCountMatrixMatchedFiltered: A matrix of of size transcripts x samples for RNAseq data, where the 
#                                   transcripts and samples have been filtered using methods specified. 
# MvalueMatrixMatchedFiltered: A matrix of size probes x samples, where the probes and samples have 
#                              been filtered using methods specified. 
# ClinicalFileMatchedFiltered: A data frame of size samples x clinical categories, where the samples 
#                              have been filtered using methods specified. 
# BetaMatrixMatchedFiltered: A matrix of size probes x samples, where the probes and samples have 
#                              been filtered using methods specified. 
# FinalFusedMatrix: A matrix of size samples x samples providing the final fused matrix generated via
#                   similarity network fusion (SNF).
# EstimatedNumberOfClusters: Number of clusters selected by model selection methods provided by 
#                            similarity network fusion (SNF) method: Eigen-gap and Rotation cost. 
#                            Generated from SNFtool::estimateNumberOfClustersGivenGraph().
# MethylationDensityOutput: Output obatined from running 3_DensityPlotMethylation function for
#                           clusters defined by SNF. 
# ProportionVisualizationOutput: Output obtained from running 21_ProportionVisualization function 
#                                for clusters defined by SNF. 
# DifferentialExpressionRNAseqOutput: Output obtained from running 36_DifferentialExpressionRNAseq
#                                     function for clusters defined by SNF. 
# DiffMethylatedRegionsOutput: Output obtained from running 9_DifferentiallyMethylatedRegions
#                              function for clusters defined by SNF. 

# Visuals saved to img folder
# 33_SNF_BetaMatrix_Image_.*
# 33_SNF_MValueMatrix_Image_.*
# 33_SNF_FusedMatrix_Image_.*
# 33_SNF_RNAseqMatrix_Image_.*
# 33_SNF_plotAlluvial_Image_.*

# Run and generate images from few other functions including:
# 3_DensityBeta_.*
# 9_DifferentiallyMethylatedRegions_.*
# 14_medianBeta_Translocation_.*
# 21_ggboxplot_MedianBetaValueAcrossClusters_Institute_.*
# 21_ggboxplot_Proportion_Btw0.2and0.8_Boxplot_.*
# 21_ggboxplot_Variance_Violin_.*
# 21_Proportion_Average_.*
# 21_Proportion_Btw0.2and0.8_DotChart_.*
# 31_MeanSDPlot_33_MeanSDPlot_BetaMatrix_Image_.*
# 33_PurityPlot_Image_.*
# 33_medianBeta_.*
# 36_EnhancedVolcano_Contrast_.*


SNFClustering <- function(RNAseqCountMatrix, 
                          MvalueMatrix, 
                          BetaMatrix, 
                          AnnotationFile,
                          QCMatrix,
                          TumorPurity,
                          ClinicalFile, 
                          SurvivalFile,
                          NumberOfClusters,
                          ImageName = paste0("Image_", date()),
                          PNGorPDF = "png",
                          ProduceImages = "Yes") {
  
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  library("SNFtool")
  library("edgeR")
  library("EnhancedVolcano")
  
  # Checking user input - sizes should be the same
  if(ncol(RNAseqCountMatrix) != ncol(MvalueMatrix)) {
    stop("\n ncol(RNAseqCountMatrix) and ncol(MvalueMatrix) should match.")
  } else if (ncol(RNAseqCountMatrix) != ncol(BetaMatrix)) {
    stop("\n ncol(RNAseqCountMatrix) and ncol(BetaMatrix) should match.")
  } else if (ncol(RNAseqCountMatrix) != nrow(ClinicalFile)) {
    stop("\n ncol(RNAseqCountMatrix) and nrow(ClinicalFile) should match.")
  } else if (ncol(RNAseqCountMatrix) != nrow(QCMatrix)) {
    stop("\n ncol(RNAseqCountMatrix) and nrow(QCMatrix) should match.")
  }
  
  # Checking user input - sample names should be the same
  if(! identical(colnames(RNAseqCountMatrix), colnames(MvalueMatrix))) {
    stop("\n Sample names in RNAseqCountMatrix and MvalueMatrix should match.")
  } else if (! identical(colnames(RNAseqCountMatrix), colnames(BetaMatrix))) {
    stop("\n Sample names in RNAseqCountMatrix and BetaMatrix should match.")
  } else if (! identical(colnames(RNAseqCountMatrix), ClinicalFile$SAMPLE_ID)) {
    stop("\n Sample names in RNAseqCountMatrix and ClinicalFile should match.")
  } else if (! identical(colnames(RNAseqCountMatrix), QCMatrix$SAMPLE_ID)) {
    stop("\n  Sample names in RNAseqCountMatrix and QCMatrix should match.")
  }
  
  
  if(typeof(RNAseqCountMatrix) != "integer") {
    stop("\n RNAseqCountMatrix should be a count matrix of non-normaized values.")
  }
  
  
  # Methylation feature selection
  set.seed(1234)
  
  # Methylation normalization
  Data1MvalueFilteredNorm <- SNFtool::standardNormalization(t(MvalueMatrix))
  # dim(Data1MvalueFilteredNorm) # 132 5000
  Data1BetaFilteredNorm <- SNFtool::standardNormalization(t(BetaMatrix))
  # dim(Data1BetaFilteredNorm) # 132 5000
  
  # RNAseq Normalization
  # As per Dr. Bo Wang's response from 27 Oct 2018
  # "I would recommend to first take transformation of  (X = log2(1 + counts)), 
  # and then apply standardNormalization()"
  Data2RNAseqFilteredNorm <- SNFtool::standardNormalization(log2(1 + t(RNAseqCountMatrix)))
  
  
  # Setting number of Parameters for SNF
  K <- round(nrow(Data1MvalueFilteredNorm) / 10); # number of neighbors; Number of patients / Number of clusters ;
  cat("\n Based on sample size, SNF recommended number of clusters/ neigbours is", K)
  #       If C is unknown, set K = N/10 usually (10~30);
  #       In our case N = 122/10 = 12
  alpha <- 0.5; # hyperparameter, usually (0.3~0.8)
  T <- 20; # Number of Iterations, usually (10~20)
  
  
  ## Calculate the Pairwise squared Euclidean distances
  # Computes the squared Euclidean distances between all pairs of data point given
  # If the data is continuous, we recommend to use the function "dist2" as follows
  Dist1Data1MvalueFilteredNorm <- (SNFtool::dist2(as.matrix(Data1MvalueFilteredNorm), 
                                                  as.matrix(Data1MvalueFilteredNorm)))^(1/2)
  # dim(Dist1Data1MvalueFilteredNorm) # sample x sample
  Dist2Data2RNAseqFilteredNorm <- (SNFtool::dist2(as.matrix(Data2RNAseqFilteredNorm), 
                                                  as.matrix(Data2RNAseqFilteredNorm)))^(1/2)
  # dim(Dist2Data2RNAseqFilteredNorm) # sample x sample
  Dist2Data2BetaFilteredNorm <- (SNFtool::dist2(as.matrix(Data1BetaFilteredNorm), 
                                                  as.matrix(Data1BetaFilteredNorm)))^(1/2)
  # dim(Dist2Data2BetaFilteredNorm) # sample x sample
  
  
  # Construct similarity graphs
  # Returns an affinity matrix that represents the neighborhood graph of the data points
  W1Data1MvalueFilteredNorm <- SNFtool::affinityMatrix(Dist1Data1MvalueFilteredNorm, K, alpha)
  # dim(W1Data1MvalueFilteredNorm) # sample x sample
  W2Data2RNAseqFilteredNorm <- SNFtool::affinityMatrix(Dist2Data2RNAseqFilteredNorm, K, alpha)
  # dim(W2Data2RNAseqFilteredNorm) # sample x sample
  W2Data2BetaFilteredNorm <- SNFtool::affinityMatrix(Dist2Data2BetaFilteredNorm, K, alpha)
  
  
  
  
  ## Fuse all the graphs
  ## then the overall matrix can be computed by similarity network fusion(SNF):
  # Similarity Network Fusion takes multiple views of a network and fuses them
  # together to construct an overall status matrix. The input to our algorithm
  # can be feature vectors, pairwise distances, or pairwise similarities. 
  # The learned status matrix can then be used for retrieval, clustering, and 
  # classification.
  set.seed(1234)
  WCombinedFiltered <- SNFtool::SNF(Wall = list(W1Data1MvalueFilteredNorm, 
                                                W2Data2RNAseqFilteredNorm), 
                                    K = K, t = T)
  
  ## Here we provide two ways to estimate the number of clusters. Note that,
  ## these two methods cannot guarantee the accuracy of esstimated number of
  ## clusters, but just to offer two insights about the datasets.
  estimationClusters <- SNFtool::estimateNumberOfClustersGivenGraph(W = WCombinedFiltered, 
                                                                    NUMC = 2:K)
  
  
  
  
  if(NumberOfClusters == 1) {
    RESULTS <- list(NumberOfClusters = NumberOfClusters,
                    RNASeqCountMatrixNormalized = Data2RNAseqFilteredNorm,
                    MvalueMatrixNormalized = Data1MvalueFilteredNorm,
                    BetaMatrixNormalized = Data1BetaFilteredNorm,
                    estimateNumberOfClustersGivenGraph = estimationClusters)
    #SurvivalAnalysisOutput = SurvivalAnalysisOutput)
  } else if (NumberOfClusters > 1) {
    # Running spectral clustering for specified number of clusters
    ClusterLabelsSNF <- list() # A list for saving results
    for (Clusters in NumberOfClusters:NumberOfClusters) {
      
      ClusterLabelsSNF[[Clusters]] <- SNFtool::spectralClustering(affinity = WCombinedFiltered, 
                                                                  K = Clusters) 
      cat("\n Tabulation based on", Clusters, "clusters \n")
      print(table(ClusterLabelsSNF[[Clusters]]))
      
      ClusterLabels =  ClusterLabelsSNF[[Clusters]]
      
      cat("\n Results based on", Clusters, "clusters for TYPE \n")
      print(table(ClusterLabelsSNF[[Clusters]], ClinicalFile$TYPE))
      print(chisq.test(table(ClusterLabelsSNF[[Clusters]], ClinicalFile$TYPE))) 
      # X-squared = 
      
      # printing labels against clinical classes
      # table(rpmmClass_numbers,ClinicalFile$TYPE)
      cat("\n Results based on", Clusters, "clusters for stage \n")
      print(table(ClusterLabelsSNF[[Clusters]], ClinicalFile$STAGE))
      # Chi-squared Test of Independence
      # Null: cluster membership is independent of stage at .05 significance level.
      print(chisq.test(table(ClusterLabelsSNF[[Clusters]], ClinicalFile$STAGE))) 
      
      cat("\n Results based on", Clusters, "clusters for sex \n")
      print(table(ClusterLabelsSNF[[Clusters]], ClinicalFile$SEX))
      print(chisq.test(table(ClusterLabelsSNF[[Clusters]], ClinicalFile$SEX))) 
      
      cat("\n Results based on", Clusters, "clusters for site biopsy \n")
      print(table(ClusterLabelsSNF[[Clusters]], trimws(ClinicalFile$SITE_BIOPSY)))
      print(chisq.test(table(ClusterLabelsSNF[[Clusters]], trimws(ClinicalFile$SITE_BIOPSY)))) 
      
      cat("\n Results based on blcTree for type biopsy \n")
      print(table(ClusterLabelsSNF[[Clusters]], trimws(ClinicalFile$TYPE_BIOPSY)))
      print(chisq.test(table(ClusterLabelsSNF[[Clusters]], trimws(ClinicalFile$TYPE_BIOPSY)))) 
      
      
      
      cat("\n Results based on", Clusters, "clusters for institution \n")
      print(table(ClusterLabelsSNF[[Clusters]], ClinicalFile$INSTITUTION))
      print(chisq.test(table(ClusterLabelsSNF[[Clusters]], ClinicalFile$INSTITUTION))) 
      
      cat("\n Results based on", Clusters, "clusters for 14_18 translocation \n")
      NA_entries <- which(is.na(ClinicalFile$TRANSLOC_14_18) == TRUE)
      print(table(ClusterLabelsSNF[[Clusters]][- NA_entries], 
                  ClinicalFile$TRANSLOC_14_18[- NA_entries]) )
      print(chisq.test((table(ClusterLabelsSNF[[Clusters]][-NA_entries], 
                              ClinicalFile$TRANSLOC_14_18[- NA_entries]) ))) 
      
      
      cat("\n Results based on blcTree for grade \n")
      matchedIDs <- match(substr(colnames(BetaMatrix[, which(substr(colnames(BetaMatrix), 4, 5) == "FL")]), 1, 9), 
                          SurvivalFile$LY_FL_ID)
      if (sum(is.na(matchedIDs) ) > 0) {
        # if NAs are present
        SurvivalFileOrderedbyBetaMatrixPatients <- SurvivalFile[matchedIDs[- which(is.na(matchedIDs) == TRUE)], ]
      } else {
        SurvivalFileOrderedbyBetaMatrixPatients <- SurvivalFile[matchedIDs, ]
      }
    }
    
    # recategorize grade
    # response on SLACK from Robert on 12 May 2020
    # It's best to regroup 1 and 2 into 1-2 category,
    # 3 and 3A into 3A category, and in situ into 1-2 category.
    # In situ means that a lymph node is not pathologically enlarged, nor otherwise abnormal, 
    # other than incidental findings of BCL2 positive B cells within germinal centre. However, 
    # if lymph node, or any lymph node was enlarged, than should be called low-grade FL, 
    # hence suggestion to choose 1-2.
    recatGRADE <- SurvivalFileOrderedbyBetaMatrixPatients$GRADE
    recatGRADE[which(SurvivalFileOrderedbyBetaMatrixPatients$GRADE == "1")] <- "1-2"
    recatGRADE[which(SurvivalFileOrderedbyBetaMatrixPatients$GRADE == "2")] <- "1-2"
    recatGRADE[which(SurvivalFileOrderedbyBetaMatrixPatients$GRADE == "IN_SITU")] <- "1-2"
    recatGRADE[which(SurvivalFileOrderedbyBetaMatrixPatients$GRADE == "3")] <- "3A"
    
    
    print(table(ClusterLabelsSNF[[Clusters]][which(substr(colnames(BetaMatrix), 4, 5) == "FL")], 
                trimws(recatGRADE)))
    print(chisq.test(table(ClusterLabelsSNF[[Clusters]][which(substr(colnames(BetaMatrix), 4, 5) == "FL")], 
                           trimws(recatGRADE)))) 
    
    
    
    
    pathNow <- getwd() # getting the path
    
    
    
    # Printing SNF fused matrix plots
    if(ProduceImages == "Yes") {
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/33_SNF_FusedMatrix_", ImageName, "_only.pdf"),
            width = 60, height = 60, pointsize = 50)
      } else {
        grDevices::png(paste0(pathNow, "/img/33_SNF_FusedMatrix_", ImageName, "_only.png"))
      }
      par(mfrow = c(1, 1))
      SNFFused <- SNFtool::displayClusters(WCombinedFiltered, ClusterLabelsSNF[[Clusters]])
      dev.off()
      
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/33_SNF_MValueMatrix_", ImageName, "_only.pdf"),
            width = 60, height = 60, pointsize = 50)
      } else {
        grDevices::png(paste0(pathNow, "/img/33_SNF_MValueMatrix_", ImageName, "_only.png"))
      }
      par(mfrow = c(1, 1))
      SNFtool::displayClusters(W1Data1MvalueFilteredNorm, ClusterLabelsSNF[[Clusters]])
      dev.off()
      
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/33_SNF_RNAseqMatrix_", ImageName, "_only.pdf"),
            width = 60, height = 60, pointsize = 50)
      } else {
        grDevices::png(paste0(pathNow, "/img/33_SNF_RNAseqMatrix_", ImageName, "_only.png"))
      }
      par(mfrow = c(1, 1))
      SNFtool::displayClusters(W2Data2RNAseqFilteredNorm, ClusterLabelsSNF[[Clusters]])
      dev.off()
      
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/33_SNF_BetaMatrix_", ImageName, "_only.pdf"),
            width = 60, height = 60, pointsize = 50)
      } else {
        grDevices::png(paste0(pathNow, "/img/33_SNF_BetaMatrix_", ImageName, "_only.png"))
      }
      par(mfrow = c(1, 1))
      SNFtool::displayClusters(W2Data2BetaFilteredNorm, ClusterLabelsSNF[[Clusters]])
      dev.off()
      
      
      
      
      # plotAlluvial
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/33_SNF_plotAlluvial_", ImageName, "_only.pdf"),
            width = 60, height = 60, pointsize = 50)
      } else {
        grDevices::png(paste0(pathNow, "/img/33_SNF_plotAlluvial_", ImageName, "only.png"))
      }
      graphics::par(mfrow = c(2, 1))
      # TYPE
      SNFtool::plotAlluvial(W = WCombinedFiltered, 
                            clust.range = c(1, Clusters),
                            color.vect = factor(ClinicalFile$TYPE, 
                                                levels = c("DLBCL", "FL", "RLN"), 
                                                labels = c("#d6604d", "#66bd63", "#4575b4")))
      
      # STAGE
      STAGE <- factor(ClinicalFile$STAGE[! is.na(ClinicalFile$STAGE)],
                      levels = c("LIMITED", "ADVANCED"), 
                      labels = c("#c2a5cf", "#762a83"))
      
      
      SNFtool::plotAlluvial(W = WCombinedFiltered[! is.na(ClinicalFile$STAGE), ! is.na(ClinicalFile$STAGE)], 
                            clust.range = c(1, Clusters),
                            color.vect = STAGE)
      grDevices::dev.off()
      
      
      # Plot clusters by purity
      par(mfrow = c(1, 1))
      PurityCluster <- data.frame(purity = TumorPurity, cluster = ClusterLabelsSNF[[Clusters]])
      
      # Define comparisons
      if(length(unique(ClusterLabelsSNF[[Clusters]])) == 2) {
        ClusterComparisonOptions <- list(c("1", "2"))
      } else if(length(unique(ClusterLabelsSNF[[Clusters]])) == 3) {
        ClusterComparisonOptions <- list(c("1", "2"),
                                         c("1", "3"), 
                                         c("2", "3"))
      } else if(length(unique(ClusterLabelsSNF[[Clusters]])) == 4) {
        ClusterComparisonOptions <- list(c("1", "2"), 
                                         c("1", "3"), 
                                         c("1", "4"), 
                                         c("2", "3"), 
                                         c("2", "4"), 
                                         c("3", "4"))
      } else if(length(unique(ClusterLabelsSNF[[Clusters]])) == 5) {
        ClusterComparisonOptions <- list(c("1", "2"), 
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
      
      
      PurityPlot <- ggpubr::ggboxplot(PurityCluster, x = "cluster", y = "purity", fill = "cluster",
                                      add = "boxplot", ylab = "Tumor Purity", 
                                      font.label = list(size = 20, color = "black"), 
                                      palette = coloursBarPlot[sort(unique(ClusterLabelsSNF[[Clusters]]))]) +
                                      ggtitle("SNF cluster vs. tumor purity") +
                                      stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                      # Add pairwise comparisons p-value
                                      stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_PurityPlot_", ImageName, ".png"))
    }
      
    
    cat("\n Mean SD plot for", Clusters, "clusters: \n")
    if(ProduceImages == "Yes") {
      MeanSDPlot(BetaMatrix = WCombinedFiltered, 
                 ClinicalFile = ClinicalFile, 
                 ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                 FigureGenerate = "Yes", 
                 PNGorPDF ="png", 
                 ImageName = paste0("33_MeanSDPlot_FusedMatrix_", ImageName))
      
      MeanSDPlot(BetaMatrix = BetaMatrix, 
                 ClinicalFile = ClinicalFile, 
                 ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                 FigureGenerate = "Yes", 
                 PNGorPDF ="png", 
                 ImageName = paste0("33_MeanSDPlot_BetaMatrix_", ImageName))
      
      MeanSDPlot(BetaMatrix = MvalueMatrix, 
                 ClinicalFile = ClinicalFile, 
                 ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                 FigureGenerate = "Yes", 
                 PNGorPDF ="png", 
                 ImageName = paste0("33_MeanSDPlot_MvalueMatrix_", ImageName))
      
      MeanSDPlot(BetaMatrix = RNAseqCountMatrix, 
                 ClinicalFile = ClinicalFile, 
                 ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                 FigureGenerate = "Yes", 
                 PNGorPDF ="png", 
                 ImageName = paste0("33_MeanSDPlot_RNAseq_", ImageName))
    } 
    
    
    cat("\n Methylation density plot for", Clusters, "clusters:\n")
    if(ProduceImages == "Yes") {
      MethylationDensityOutput <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                         ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                                                         PlotWithinAnnotationCategories = "No", 
                                                         ClinicalCategoryToVisualize = "CLUSTER", 
                                                         BetaMatrix = BetaMatrix, 
                                                         AnnotationFile = AnnotationFile, 
                                                         ClinicalFile = ClinicalFile, 
                                                         SampleSheet = NA, 
                                                         FigureGenerate = "Yes", 
                                                         ImageName = paste0("MethylationPlot_AllProbes", ImageName), 
                                                         PNGorPDF = "png")
    } else {
      MethylationDensityOutput <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                         ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                                                         PlotWithinAnnotationCategories = "No", 
                                                         ClinicalCategoryToVisualize = "CLUSTER", 
                                                         BetaMatrix = BetaMatrix, 
                                                         AnnotationFile = AnnotationFile, 
                                                         ClinicalFile = ClinicalFile, 
                                                         SampleSheet = NA, 
                                                         FigureGenerate = "No")
    }
    
    
    cat("\n Proportion visualization for", Clusters, "clusters:\n")
    if(ProduceImages == "Yes") {
      ProportionVisualizationOutput <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", 
                                                               ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                                                               BetaMatrix = BetaMatrix, 
                                                               AnnotationFile = AnnotationFile, 
                                                               ClinicalFile = ClinicalFile, 
                                                               ClinicalCategory = "TYPE", 
                                                               PlotWithinCategories = "Yes", 
                                                               FigureGenerate = "Yes", 
                                                               PNGorPDF = "png", 
                                                               ImageName = paste0("AllProbes_", ImageName))
    } else{
      ProportionVisualizationOutput <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", 
                                                               ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                                                               BetaMatrix = BetaMatrix, 
                                                               AnnotationFile = AnnotationFile, 
                                                               ClinicalFile = ClinicalFile, 
                                                               ClinicalCategory = "TYPE", 
                                                               PlotWithinCategories = "Yes", 
                                                               FigureGenerate = "No")
    }
    
    
    cat("\n Differentially expression analysis for", Clusters, "clusters:\n")
    if(ProduceImages == "Yes") {
      DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = RNAseqCountMatrix, 
                                                                         ContrastColumnName = "CLUSTER", 
                                                                         ClinicalFile = ClinicalFile, 
                                                                         ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                                                                         ProduceImages = "Yes", 
                                                                         PNGorPDF = "png") 
    } else {
      DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = RNAseqCountMatrix, 
                                                                         ContrastColumnName = "CLUSTER", 
                                                                         ClinicalFile = ClinicalFile, 
                                                                         ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                                                                         ProduceImages = "No") 
    }
    
    cat("\n Differentially methylated regions analysis for", Clusters, "clusters:\n")
    if(ProduceImages == "Yes") {
      DiffMethylatedRegionsOutput <- DiffMethylatedRegions(Method = "DMRcate", 
                                                           BetaMatrix = BetaMatrix, 
                                                           MvalueMatrix = MvalueMatrix, 
                                                           ContrastColumnName = "CLUSTER", 
                                                           ClinicalFile = ClinicalFile, 
                                                           ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                                                           AnnotationFile = AnnotationFile, 
                                                           ProduceImages = "Yes", 
                                                           DMR = 1, 
                                                           PNGorPDF = "png",
                                                           ExpressionFile = NA)
    } else {
      DiffMethylatedRegionsOutput <- DiffMethylatedRegions(Method = "DMRcate", 
                                                           BetaMatrix = BetaMatrix, 
                                                           MvalueMatrix = MvalueMatrix, 
                                                           ContrastColumnName = "CLUSTER", 
                                                           ClinicalFile = ClinicalFile, 
                                                           ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                                                           AnnotationFile = AnnotationFile, 
                                                           ProduceImages = "No", 
                                                           DMR = 1, 
                                                           ExpressionFile = NA)
    }
    
    
    cat("\n Survival analysis for", Clusters, "clusters:\n")
    if(ProduceImages == "NotRun") {
      SurvivalAnalysisOutput <- SurvivalAnalysis(BetaMatrix = BetaMatrix, 
                                                 ClinicalFile = ClinicalFile, 
                                                 SurvivalFile = SurvivalFile, 
                                                 ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                                                 FigureGenerate = "Yes", 
                                                 PNGorPDF = "png")
    }
    
    # PCA analysis
    cat("\n PCA analysis for", Clusters, "clusters:\n")
    if(ProduceImages == "Yes") {  
      set.seed(1234)
      res.pcaFiltered <- list()
      # using selected beta probes 
      res.pcaFiltered[[1]] <- prcomp(t(WCombinedFiltered), 
                                center = TRUE, 
                                scale = TRUE)
      
      res.pcaFiltered[[2]] <- prcomp(t(RNAseqCountMatrix), 
                                     center = TRUE, 
                                     scale = TRUE)
      
      res.pcaFiltered[[3]] <- prcomp(t(BetaMatrix), 
                                     center = TRUE, 
                                     scale = TRUE)
      
      res.pcaFiltered[[4]] <- prcomp(t(MvalueMatrix), 
                                     center = TRUE, 
                                     scale = TRUE)
      
      ListOfCategories <- c("WCombinedFiltered", "RNAseqCountMatrix", "BetaMatrix", "MvalueMatrix")
      
      for (category in 1:4) {
        
      
        if (PNGorPDF == "pdf") {
          grDevices::pdf(paste0(pathNow, "/img/33_PCA_Cluster_", ListOfCategories[[category]], ImageName, ".pdf"),
                         width = 60, height = 60, pointsize = 50)
        } else {
          grDevices::png(paste0(pathNow,"/img/33_PCA_Cluster_", ListOfCategories[[category]], ImageName, ".png"))
        }
        par(mfrow = c(1, 1))
        # cluster
        p0Filtered <- factoextra::fviz_pca_ind(res.pcaFiltered[[category]], 
                                               label = "none", 
                                               habillage = factor(ClusterLabelsSNF[[Clusters]]),
                                               addEllipses = TRUE, ellipse.level = 0.95,  
                                               palette = c("#4363d8", "#f58231")) +
                                               theme_minimal()
        grDevices::dev.off()
        
        
        if (PNGorPDF == "pdf") {
          grDevices::pdf(paste0(pathNow, "/img/33_PCA_TYPE_", ListOfCategories[[category]], ImageName, ".pdf"),
                         width = 60, height = 60, pointsize = 50)
        } else {
          grDevices::png(paste0(pathNow,"/img/33_PCA_TYPE_", ListOfCategories[[category]], ImageName, ".png"))
        }
        par(mfrow = c(1, 1))
        # type
        p1Filtered <- factoextra::fviz_pca_ind(res.pcaFiltered[[category]], 
                                               label = "none", 
                                               habillage = factor(ClinicalFile$TYPE),
                                               addEllipses = TRUE, ellipse.level = 0.95,  
                                               palette = c("#d6604d", "#66bd63", "#4575b4")) +
                                               theme_minimal()
        grDevices::dev.off()
        
        
        if (PNGorPDF == "pdf") {
          grDevices::pdf(paste0(pathNow, "/img/33_PCA_STAGE_", ListOfCategories[[category]], ImageName, ".pdf"),
                         width = 60, height = 60, pointsize = 50)
        } else {
          grDevices::png(paste0(pathNow,"/img/33_PCA_STAGE_", ListOfCategories[[category]], ImageName, ".png"), 
                         height = 15.41,
                         width = 30.24,
                         units = "cm",
                         res = 250)
        }
        par(mfrow = c(1, 1))
        # stage
        p2Filtered <- factoextra::fviz_pca_ind(res.pcaFiltered[[category]], 
                                               label = "none", 
                                               habillage = factor(ClinicalFile$STAGE),
                                               addEllipses = TRUE, ellipse.level = 0.95,  
                                               palette = c("#762a83", "#c2a5cf")) +
                                               theme_minimal()
        grDevices::dev.off()
        
        
        if (PNGorPDF == "pdf") {
          grDevices::pdf(paste0(pathNow, "/img/33_PCA_translocation_", ListOfCategories[[category]], ImageName, ".pdf"),
                         width = 60, height = 60, pointsize = 50)
        } else {
          grDevices::png(paste0(pathNow,"/img/33_PCA_translocation_", ListOfCategories[[category]], ImageName, ".png"), 
                         height = 15.41,
                         width = 30.24,
                         units = "cm",
                         res = 250)
        }
        par(mfrow = c(1, 1))
        # translocation status
        p3Filtered <- factoextra::fviz_pca_ind(res.pcaFiltered[[category]], 
                                               label = "none", 
                                               habillage = factor(ClinicalFile$TRANSLOC_14_18),
                                               addEllipses = TRUE, ellipse.level = 0.95,  
                                               palette = c("#e0e0e0", "#878787")) +
                                               theme_minimal()
        grDevices::dev.off()
      }
    }
    
    cat("\n Beta heatmap for", Clusters, "clusters:\n")
    cat("\n Heatmaps only support upto 3 clusters. Otherwise alter code.\n")
    BetaHeatmap <- function() {
      ########################## Plot heatmap
      ClusterLabels =  ClusterLabelsSNF[[Clusters]]
      
      matchRowNames <- match(rownames(BetaMatrix), AnnotationFile$V1)
      geneRegion <- sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Group)
      geneRegion[which(geneRegion == "")] <- "NA"
      geneRegionBreakdownIsland <- AnnotationFile$Relation_to_Island[matchRowNames]
      geneRegionBreakdownGeneregion <- geneRegion[matchRowNames]
      
      orderRows <- c(which(geneRegionBreakdownIsland == "Island"),
                     which(geneRegionBreakdownIsland == "N_Shelf"),
                     which(geneRegionBreakdownIsland == "N_Shore"),
                     which(geneRegionBreakdownIsland == "OpenSea"),
                     which(geneRegionBreakdownIsland == "S_Shelf"),
                     which(geneRegionBreakdownIsland == "S_Shore"))
      
      orderRowsGeneRegions <- c(which(geneRegionBreakdownGeneregion == "1stExon"),
                                which(geneRegionBreakdownGeneregion == "3'UTR"),
                                which(geneRegionBreakdownGeneregion == "5'UTR"),
                                which(geneRegionBreakdownGeneregion == "Body"),
                                which(geneRegionBreakdownGeneregion == "ExonBnd"),
                                which(geneRegionBreakdownGeneregion == "TSS1500"),
                                which(geneRegionBreakdownGeneregion == "TSS200"),
                                which(geneRegionBreakdownGeneregion == "NA"))
      
      annotation_row <- data.frame (
        GeneRegion = factor(geneRegionBreakdownGeneregion)[orderRowsGeneRegions],
        RelationToIsland = factor(geneRegionBreakdownIsland[orderRowsGeneRegions]))
      
      rownames(annotation_row) = rownames(BetaMatrix[orderRowsGeneRegions, ])
      
      
      
      # order patients by cluster 
      orderColumns <- c(which(ClusterLabels == 1),
                        which(ClusterLabels == 2),
                        which(ClusterLabels == 3))
      
      annotation_col <- data.frame(
        Cluster = factor(ClusterLabels[orderColumns]),
        Disease = factor(ClinicalFile$TYPE[orderColumns]),
        Stage = factor(ClinicalFile$STAGE[orderColumns]),
        Sex = factor(ClinicalFile$SEX[orderColumns]),
        Translocation = factor(ClinicalFile$TRANSLOC_14_18[orderColumns]),
        TypeBiopsy =  ClinicalFile$TYPE_BIOPSY[orderColumns],
        SiteBiopsy =  ClinicalFile$SITE_BIOPSY[orderColumns])
      rownames(annotation_col) = colnames(BetaMatrix[, orderColumns])
      
      # use http://colorbrewer2.org/#type=diverging&scheme=RdGy&n=8 to decide colors 
      
      
      ann_colors <- list(
        RelationToIsland = c("Island" = "#a6cee3", "N_Shelf" = "#1f78b4",
                             "N_Shore" = "#fb9a99", "OpenSea" = "#fdbf6f", 
                             "S_Shelf" = "#ff7f00", "S_Shore" = "#cab2d6"),
        GeneRegion = c("1stExon" = "#8dd3c7", "3'UTR" = "#ffffb3", 
                       "5'UTR" = "#bebada", "Body" = "#fb8072",
                       "ExonBnd" = "#80b1d3", "NA" = "#fdb462", 
                       "TSS1500" = "#fccde5", "TSS200" = "#d9d9d9"),
        Cluster = c("1" = "#4363d8", "2" = "#f58231", "3" = '#911eb4'), 
        Disease = c(DLBCL = "#d6604d", FL = "#66bd63", RLN = "#4575b4"), #option 5
        Stage = c(ADVANCED = "#762a83", LIMITED = "#c2a5cf"),
        Sex = c("F"="#b35806", "M"="#fdb863"),
        Translocation = c("0"="#e0e0e0","1"="#878787"),
        TypeBiopsy = c("TISSUE" = "#a6dba0", "CORE"="#878787"),
        SiteBiopsy = c("LN" = "#f1b6da" , "EN" = "#c51b7d"))
      
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/33_Beta_Heatmap_", ImageName, ".pdf"),
            width = 60, height = 60, pointsize = 50)
      } else {
        grDevices::png(paste0(pathNow,"/img/33_Beta_Heatmap_", ImageName, ".png"), 
                       height = 30.41,
                       width = 60.24,
                       units = "cm",
                       res = 250)
      }
      par(mfrow = c(1, 1))
      heatmap_all <- pheatmap::pheatmap(as.matrix(BetaMatrix[orderRowsGeneRegions, orderColumns]), 
                                        show_colnames = T, 
                                        show_rownames = F, 
                                        fontface = "italic", 
                                        legend = T, scale ="row", 
                                        annotation_colors = ann_colors, 
                                        border_color = "black", 
                                        cluster_row = FALSE, 
                                        cluster_cols = FALSE,
                                        annotation_col = annotation_col, 
                                        annotation_row = annotation_row,
                                        color =  rev(redgreen(1000)))
      grDevices::dev.off()
      
      # recategorize grade
      # response on SLACK from Robert on 12 May 2020
      # It's best to regroup 1 and 2 into 1-2 category,
      # 3 and 3A into 3A category, and in situ into 1-2 category.
      # In situ means that a lymph node is not pathologically enlarged, nor otherwise abnormal, 
      # other than incidental findings of BCL2 positive B cells within germinal centre. However, 
      # if lymph node, or any lymph node was enlarged, than should be called low-grade FL, 
      # hence suggestion to choose 1-2.
      matchedIDs <- match(substr(colnames(RNAseqCountMatrix), 1, 9), 
                          SurvivalFile$LY_FL_ID)
      recatGRADE <- SurvivalFile$GRADE[matchedIDs]
      recatGRADE[which(SurvivalFile$GRADE[matchedIDs] == "1")] <- "1-2"
      recatGRADE[which(SurvivalFile$GRADE[matchedIDs] == "2")] <- "1-2"
      recatGRADE[which(SurvivalFile$GRADE[matchedIDs] == "IN_SITU")] <- "1-2"
      recatGRADE[which(SurvivalFile$GRADE[matchedIDs] == "3")] <- "3A"
      
      par(mfrow = c(1, 1))
      PurityCluster <- data.frame(purity = TumorPurity, 
                                  cluster = ClusterLabels,
                                  Clusters = as.character(ClusterLabels),
                                  stage = ClinicalFile$STAGE,# [c(11:165), ]
                                  translocation = ClinicalFile$TRANSLOC_14_18, # [c(11:165), ]
                                  type = ClinicalFile$TYPE, # [c(11:165), ]
                                  sex = ClinicalFile$SEX,
                                  sitebiopsy = ClinicalFile$SITE_BIOPSY,
                                  institution = ClinicalFile$INSTITUTION,
                                  typebiopsy = ClinicalFile$TYPE_BIOPSY,
                                  grade = trimws(recatGRADE),
                                  meanBeta = rowMeans(t(BetaMatrix)),
                                  medianBeta = rowMedians(t(BetaMatrix)),
                                  sdBeta = apply(t(BetaMatrix), 1, sd, na.rm = TRUE))
      
      # Define comparisons
      if(length(unique(ClusterLabels)) == 2) {
        ClusterComparisonOptions <- list(c("1","2"))
      } else if(length(unique(ClusterLabels)) == 3) {
        ClusterComparisonOptions <- list(c("1", "2"),
                                         c("1", "3"), 
                                         c("2", "3"))
      } else if(length(unique(ClusterLabels)) == 4) {
        ClusterComparisonOptions <- list(c("1", "2"), 
                                         c("1", "3"), 
                                         c("1", "4"), 
                                         c("2", "3"), 
                                         c("2", "4"), 
                                         c("3", "4"))
      } else if(length(unique(ClusterLabels)) == 5) {
        ClusterComparisonOptions <- list(c("1", "2"), 
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
      
      
      grDevices::dev.off()
      meanBetaPlot <- ggpubr::ggboxplot(PurityCluster, 
                                        x = "cluster", 
                                        y = "meanBeta", 
                                        fill = "cluster",
                                        add = "boxplot",
                                        xlab = "Cluster",
                                        ylab = "Mean Beta Value", 
                                        font.label = list(size = 20, color = "black"), 
                                        palette = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                        ggtitle("Cluster vs. meanBeta values") +
                                        stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                        # Add pairwise comparisons p-value
                                        stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_meanBeta_", ImageName, ".png"))
      
      
      medianBetaPlot <- ggpubr::ggboxplot(PurityCluster, 
                                          x = "cluster", 
                                          y = "medianBeta", 
                                          fill = "cluster",
                                          add = "boxplot", 
                                          xlab = "Cluster", 
                                          ylab = "Median Beta Value", 
                                          font.label = list(size = 20, color = "black"), 
                                          palette = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                          ggtitle("Cluster vs. medianBeta values") +
                                          stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                          # Add pairwise comparisons p-value
                                          stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianBeta_", ImageName, ".png"))  
      
      
      # Plotting median beta value for patients in each cluster, seperated by type
      # Setting the number of comparisons for ggpubr
      if(length(unique(PurityCluster$type)) == 2) {
        ComparisonOptionsType <- list(names(table(PurityCluster$type))[1:2])
      } else if(length(unique(PurityCluster$type)) == 3) {
        ComparisonOptionsType <- list(names(table(PurityCluster$type))[1:2], 
                                      names(table(PurityCluster$type))[2:3],
                                      names(table(PurityCluster$type))[c(1, 3)])
      } else if(length(unique(PurityCluster$type)) == 4) {
        ComparisonOptionsType <- list( names(table(PurityCluster$type))[1:2], 
                                       names(table(PurityCluster$type))[c(1, 3)], 
                                       names(table(PurityCluster$type))[c(1, 4)], 
                                       names(table(PurityCluster$type))[c(2, 3)], 
                                       names(table(PurityCluster$type))[c(2, 4)], 
                                       names(table(PurityCluster$type))[c(3, 4)])
      }
      
      p0ClusterType <- ggpubr::ggboxplot(PurityCluster, 
                                           x = "type", 
                                           y = "medianBeta", 
                                           fill = "Clusters",
                                           add = "boxplot", 
                                           xlab = "Type", 
                                           ylab = "Median Beta Values",
                                           palette = coloursBarPlot[1:4]) + 
                                           stat_compare_means(comparisons = ComparisonOptionsType) + 
                                           # Add pairwise comparisons p-value
                                           stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianBeta_TypeCluster", ImageName, ".png")) 
      
      # Plotting median beta value for patients in each cluster, seperated by stage
      # Setting the number of comparisons for ggpubr
      if(length(unique(PurityCluster$stage[-which(is.na(PurityCluster$stage))])) == 2) {
        ComparisonOptionsStage <- list(names(table(PurityCluster$stage))[1:2])
      } 
      
      p0ClusterStage <- ggpubr::ggboxplot(PurityCluster[- which(is.na(PurityCluster$stage)),], 
                                            x = "stage", 
                                            y = "medianBeta", 
                                            fill = "Clusters",
                                            add = "boxplot", 
                                            xlab = "Stage", 
                                            ylab = "Median Beta Values",
                                            palette = coloursBarPlot[1:4]) + 
                                            stat_compare_means(comparisons = ComparisonOptionsStage) + 
                                            # Add pairwise comparisons p-value
                                            stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianBeta_StageCluster", ImageName, ".png")) 
      
      
      # Plotting median beta value for patients in each cluster, seperated by translocation
      # Setting the number of comparisons for ggpubr
      if(length(unique(PurityCluster$stage[- which(is.na(PurityCluster$translocation))])) == 2) {
        ComparisonOptionsTranslocation <- list(names(table(PurityCluster$translocation))[1:2])
      } 
      
      p0ClusterTranslocation <- ggpubr::ggboxplot(PurityCluster[- which(is.na(PurityCluster$translocation)),], 
                                                    x = "translocation", 
                                                    y = "medianBeta", 
                                                    fill = "Clusters",
                                                    add = "boxplot", 
                                                    xlab = "Translocation Status", 
                                                    ylab = "Median Beta Values",
                                                    palette = coloursBarPlot[1:4]) + 
                                                    stat_compare_means(comparisons = ComparisonOptionsTranslocation) + 
                                                    # Add pairwise comparisons p-value
                                                    stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/14_medianBeta_TranslocationCluster", ImageName, ".png")) 
      
      
      if(length(unique(PurityCluster$sex[- which(is.na(PurityCluster$sex))])) == 2) {
        ComparisonOptionsSex <- list(names(table(PurityCluster$sex))[1:2])
      } 
      p0ClusterSex <- ggpubr::ggboxplot(PurityCluster[- which(is.na(PurityCluster$sex)),], 
                                        x = "sex", 
                                        y = "medianBeta", 
                                        fill = "Clusters",
                                        add = "boxplot", 
                                        xlab = "Sex", 
                                        ylab = "Median Beta Values",
                                        palette = coloursBarPlot[1:4]) + 
                                        stat_compare_means(comparisons = ComparisonOptionsSex) + 
                                        # Add pairwise comparisons p-value
                                        stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianBeta_SexCluster", ImageName, ".png")) 
      
      
      if(length(unique(PurityCluster$sitebiopsy[- which(is.na(PurityCluster$sitebiopsy))])) == 2) {
        ComparisonOptionsSitebiopsy <- list(names(table(PurityCluster$sitebiopsy))[1:2])
      } 
      p0ClusterSitebiopsy <- ggpubr::ggboxplot(PurityCluster[- which(is.na(PurityCluster$sitebiopsy)),], 
                                               x = "sitebiopsy", 
                                               y = "medianBeta", 
                                               fill = "Clusters",
                                               add = "boxplot", 
                                               xlab = "Site biopsy", 
                                               ylab = "Median Beta Values",
                                               palette = coloursBarPlot[1:4]) + 
                                               stat_compare_means(comparisons = ComparisonOptionsSitebiopsy) + 
                                               # Add pairwise comparisons p-value
                                               stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianBeta_SitebiopsyCluster", ImageName, ".png")) 
      
      
      
      if(length(unique(PurityCluster$institution)) == 3) {
        ComparisonOptionsInstitution <- list(names(table(PurityCluster$institution))[1:3])
      } 
      p0ClusterInstitute <- ggpubr::ggboxplot(PurityCluster, 
                                              x = "institution", 
                                              y = "medianBeta", 
                                              fill = "Clusters",
                                              add = "boxplot", 
                                              xlab = "Institute", 
                                              ylab = "Median Beta Values",
                                              palette = coloursBarPlot[1:4]) + 
                                              stat_compare_means(comparisons = ComparisonOptionsInstitution) + 
                                              # Add pairwise comparisons p-value
                                              stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianBeta_InstituteCluster", ImageName, ".png")) 
      
      
      if(length(unique(PurityCluster$typebiopsy)) == 2) {
        ComparisonOptionsTypebiopsy <- list(names(table(PurityCluster$typebiopsy))[1:2])
      } 
      p0ClusterTypebiopsy <- ggpubr::ggboxplot(PurityCluster, 
                                              x = "typebiopsy", 
                                              y = "medianBeta", 
                                              fill = "Clusters",
                                              add = "boxplot", 
                                              xlab = "Type biopsy", 
                                              ylab = "Median Beta Values",
                                              palette = coloursBarPlot[1:4]) + 
                                              stat_compare_means(comparisons = ComparisonOptionsTypebiopsy) + 
                                              # Add pairwise comparisons p-value
                                              stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianBeta_TypebiopsyCluster", ImageName, ".png")) 
      
      
      if(length(unique(PurityCluster$grade[- which(is.na(PurityCluster$grade))])) == 2) {
        ComparisonOptionsGrade <- list(names(table(PurityCluster$grade))[1:2])
      } 
      p0ClusterSitebiopsy <- ggpubr::ggboxplot(PurityCluster[- which(is.na(PurityCluster$grade)),], 
                                               x = "grade", 
                                               y = "medianBeta", 
                                               fill = "Clusters",
                                               add = "boxplot", 
                                               xlab = "Grade", 
                                               ylab = "Median Beta Values",
                                               palette = coloursBarPlot[1:4]) + 
                                               stat_compare_means(comparisons = ComparisonOptionsGrade) + 
                                               # Add pairwise comparisons p-value
                                               stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianBeta_GradeCluster", ImageName, ".png")) 
      
      
      
    }
    if(ProduceImages == "Yes") {
      BetaHeatmap()
    }
    
    cat("\n RNAseq heatmap for", Clusters, "clusters:\n")
    RNAseqHeatmap <- function() {
      # RNAseq data counts
      ClusterLabels =  ClusterLabelsSNF[[Clusters]]
      
      orderColumnsRNAseq <- c(which(ClusterLabels == 1),
                              which(ClusterLabels == 2),
                              which(ClusterLabels == 3))
      
      Cluster = factor(ClusterLabels)
      Disease = factor(ClinicalFile$TYPE)
      Stage = factor(ClinicalFile$STAGE)
      Sex = factor(ClinicalFile$SEX)
      Translocation = factor(ClinicalFile$TRANSLOC_14_18)
      TypeBiopsy =  factor(ClinicalFile$TYPE_BIOPSY)
      SiteBiopsy =  factor(ClinicalFile$SITE_BIOPSY)
      
      annotation_colRNAseq <- data.frame(
        Cluster = factor(Cluster[orderColumnsRNAseq]),
        Disease = factor(Disease[orderColumnsRNAseq]),
        Stage = factor(Stage[orderColumnsRNAseq]),
        Sex = factor(Sex[orderColumnsRNAseq]),
        Translocation = factor(Translocation[orderColumnsRNAseq]),
        TypeBiopsy = factor(TypeBiopsy[orderColumnsRNAseq]),
        SiteBiopsy = factor(SiteBiopsy[orderColumnsRNAseq]))
      rownames(annotation_colRNAseq) <- colnames(RNAseqCountMatrix[, orderColumnsRNAseq])
      
      ann_colorsRNAseq <- list(
        Cluster = c("1" = "#4363d8", "2" = "#f58231", "3" = '#911eb4'), 
        Disease = c(DLBCL = "#d6604d", FL = "#66bd63", RLN = "#4575b4"), #option 5
        Stage = c(ADVANCED = "#762a83", LIMITED = "#c2a5cf"),
        Sex = c("F" = "#b35806", "M" = "#fdb863"),
        Translocation = c("0" = "#e0e0e0","1" = "#878787"),
        TypeBiopsy = c("TISSUE" = "#a6dba0", "CORE" = "#878787"),
        SiteBiopsy = c("LN" = "#f1b6da" , "EN" = "#c51b7d"))
      
      
      
      ## Get some nicer colours
      mypalette <- RColorBrewer::brewer.pal(11, "RdYlBu")
      morecols <- grDevices::colorRampPalette(mypalette)
      
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/33_LogRNAseq_Heatmap_", ImageName, ".pdf"),
                       width = 60, height = 60, pointsize = 50)
      } else {
        grDevices::png(paste0(pathNow,"/img/33_LogRNAseq_Heatmap_", ImageName, ".png"), 
                       height = 30.41,
                       width = 60.24,
                       units = "cm",
                       res = 250)
      }
      par(mfrow = c(1, 1))
      # Plot the heatmap
      pheatmap::pheatmap((log(RNAseqCountMatrix + 0.01))[, orderColumnsRNAseq],
                         show_colnames = T, 
                         show_rownames = F,
                         fontface = "italic", 
                         legend = T,
                         annotation_colors = ann_colorsRNAseq, 
                         border_color = "black", 
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,
                         annotation_col = annotation_colRNAseq, 
                         color =  rev(morecols(50)))
      grDevices::dev.off()
      
      
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/33_SNFStandardNormalizeRNAseq_Heatmap_", ImageName, ".pdf"),
                       width = 60, height = 60, pointsize = 50)
      } else {
        grDevices::png(paste0(pathNow,"/img/33_SNFStandardNormalizeRNAseq_Heatmap_", ImageName, ".png"), 
                       height = 30.41,
                       width = 60.24,
                       units = "cm",
                       res = 250)
      }
      par(mfrow = c(1, 1))
      # Plot the heatmap of normalized counts 
      pheatmap::pheatmap(t(Data2RNAseqFilteredNorm)[, orderColumnsRNAseq],
                         show_colnames = T, 
                         show_rownames = F,
                         fontface = "italic", 
                         legend = T,
                         annotation_colors = ann_colorsRNAseq, 
                         border_color = "black", 
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,
                         annotation_col = annotation_colRNAseq, 
                         color =  rev(morecols(50)))
      grDevices::dev.off()
      
      
      # generate edgeR normalized counts
      ClinicalFileT1Cluster <- data.frame(ClinicalFile, "CLUSTER" = ClusterLabels)
      # dim(ClinicalFileT1Cluster) # 132  26
      exprsRNAseq <- edgeR::DGEList(counts = RNAseqCountMatrix, 
                                    group = factor(ClinicalFileT1Cluster$CLUSTER))
      exprsRNAseq$counts <- exprsRNAseq$counts
      
      # perform the TMM normalization and display the normalization factors
      Data2RNAseqFilterededgeRNormFactors <- edgeR::calcNormFactors(exprsRNAseq)
      # dim(Data2RNAseqFilterededgeRNormFactors) # 45266   132
      # plotMDS(Data2RNAseqFilterededgeRNorm)
      # Generate matrix
      Data2RNAseqFilterededgeRNorm <- edgeR::cpm(Data2RNAseqFilterededgeRNormFactors, 
                                                 normalized.lib.sizes = TRUE, log = TRUE)
      # dim(Data2RNAseqFilterededgeRNorm) # 45266   132
      
      
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/33_edgeRNormalizeRNAseq_Heatmap_", ImageName, ".pdf"),
                       width = 60, height = 60, pointsize = 50)
      } else {
        grDevices::png(paste0(pathNow,"/img/33_edgeRNormalizeRNAseq_Heatmap_", ImageName, ".png"), 
                       height = 30.41,
                       width = 60.24,
                       units = "cm",
                       res = 250)
      }
      par(mfrow = c(1, 1))
      # Plot the heatmap of normalized counts 
      pheatmap::pheatmap(Data2RNAseqFilterededgeRNorm[, orderColumnsRNAseq],
                         show_colnames = T, 
                         show_rownames = F,
                         fontface = "italic", 
                         legend = T,
                         annotation_colors = ann_colorsRNAseq, 
                         border_color = "black", 
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,
                         annotation_col = annotation_colRNAseq, 
                         color =  rev(morecols(50)))
      grDevices::dev.off()
      
      
      # recategorize grade
      # response on SLACK from Robert on 12 May 2020
      # It's best to regroup 1 and 2 into 1-2 category,
      # 3 and 3A into 3A category, and in situ into 1-2 category.
      # In situ means that a lymph node is not pathologically enlarged, nor otherwise abnormal, 
      # other than incidental findings of BCL2 positive B cells within germinal centre. However, 
      # if lymph node, or any lymph node was enlarged, than should be called low-grade FL, 
      # hence suggestion to choose 1-2.
      matchedIDs <- match(substr(colnames(RNAseqCountMatrix), 1, 9), 
                          SurvivalFile$LY_FL_ID)
      recatGRADE <- SurvivalFile$GRADE[matchedIDs]
      recatGRADE[which(SurvivalFile$GRADE[matchedIDs] == "1")] <- "1-2"
      recatGRADE[which(SurvivalFile$GRADE[matchedIDs] == "2")] <- "1-2"
      recatGRADE[which(SurvivalFile$GRADE[matchedIDs] == "IN_SITU")] <- "1-2"
      recatGRADE[which(SurvivalFile$GRADE[matchedIDs] == "3")] <- "3A"
      
      par(mfrow = c(1, 1))
      RNAseqCluster <- data.frame(purity = TumorPurity, 
                                  cluster = ClusterLabels,
                                  Clusters = as.character(ClusterLabels),
                                  stage = ClinicalFile$STAGE,# [c(11:165), ]
                                  translocation = ClinicalFile$TRANSLOC_14_18, # [c(11:165), ]
                                  type = ClinicalFile$TYPE, # [c(11:165), ]
                                  sex = ClinicalFile$SEX,
                                  sitebiopsy = ClinicalFile$SITE_BIOPSY,
                                  institution = ClinicalFile$INSTITUTION,
                                  typebiopsy = ClinicalFile$TYPE_BIOPSY,
                                  grade = trimws(recatGRADE),
                                  meanCount = rowMeans(t(as.matrix(RNAseqCountMatrix))),
                                  medianCount = rowMedians(t(as.matrix(RNAseqCountMatrix))),
                                  sdCount = apply(t(as.matrix(RNAseqCountMatrix)),1, sd, na.rm = TRUE))
      
      # Define comparisons
      if(length(unique(ClusterLabels)) == 2) {
        ClusterComparisonOptions <- list(c("1", "2"))
      } else if(length(unique(ClusterLabels)) == 3) {
        ClusterComparisonOptions <- list(c("1", "2"),
                                         c("1", "3"), 
                                         c("2", "3"))
      } else if(length(unique(ClusterLabels)) == 4) {
        ClusterComparisonOptions <- list(c("1", "2"), 
                                         c("1", "3"), 
                                         c("1", "4"), 
                                         c("2", "3"), 
                                         c("2", "4"), 
                                         c("3", "4"))
      } else if(length(unique(ClusterLabels)) == 5) {
        ClusterComparisonOptions <- list(c("1", "2"), 
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
      
      
      meanCountPlot <- ggpubr::ggboxplot(RNAseqCluster, x = "cluster", y = "meanCount", fill = "cluster",
                                         add = "boxplot", ylab = "Mean Count Value", 
                                         font.label = list(size = 20, color = "black"), 
                                         palette = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                         ggtitle("Cluster vs. Mean Count Value") +
                                         stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                         # Add pairwise comparisons p-value
                                         stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_meanCount_", ImageName, ".png"))      
      
      
      medianCountPlot <- ggpubr::ggboxplot(RNAseqCluster, x = "cluster", y = "medianCount", fill = "cluster",
                                           add = "boxplot", ylab = "Median Count Value", 
                                           font.label = list(size = 20, color = "black"), 
                                           palette = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                           ggtitle("Cluster vs. Median Count Value") +
                                           stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                           # Add pairwise comparisons p-value
                                           stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianCount_", ImageName, ".png"))    
      
      
      # Plotting median beta value for patients in each cluster, seperated by type
      # Setting the number of comparisons for ggpubr
      if(length(unique(RNAseqCluster$type)) == 2) {
        ComparisonOptionsType <- list(names(table(RNAseqCluster$type))[1:2])
      } else if(length(unique(RNAseqCluster$type)) == 3) {
        ComparisonOptionsType <- list(names(table(RNAseqCluster$type))[1:2], 
                                      names(table(RNAseqCluster$type))[2:3],
                                      names(table(RNAseqCluster$type))[c(1, 3)])
      } else if(length(unique(RNAseqCluster$type)) == 4) {
        ComparisonOptionsType <- list( names(table(RNAseqCluster$type))[1:2], 
                                       names(table(RNAseqCluster$type))[c(1, 3)], 
                                       names(table(RNAseqCluster$type))[c(1, 4)], 
                                       names(table(RNAseqCluster$type))[c(2, 3)], 
                                       names(table(RNAseqCluster$type))[c(2, 4)], 
                                       names(table(RNAseqCluster$type))[c(3, 4)])
      }
      
      p0ClusterType <- ggpubr::ggboxplot(RNAseqCluster, 
                                           x = "type", 
                                           y = "medianCount", 
                                           fill = "Clusters",
                                           add = "boxplot", 
                                           xlab = "Type", 
                                           ylab = "Median Count Values",
                                           palette = coloursBarPlot[1:4]) + 
                                           stat_compare_means(comparisons = ComparisonOptionsType) + 
                                           # Add pairwise comparisons p-value
                                           stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianCount_TypeCluster", ImageName, ".png")) 
      
      # Plotting median beta value for patients in each cluster, seperated by stage
      # Setting the number of comparisons for ggpubr
      if(length(unique(RNAseqCluster$stage[-which(is.na(RNAseqCluster$stage))])) == 2) {
        ComparisonOptionsStage <- list(names(table(RNAseqCluster$stage))[1:2])
      } 
      
      p0ClusterStage <- ggpubr::ggboxplot(RNAseqCluster[- which(is.na(RNAseqCluster$stage)),], 
                                          x = "stage", 
                                          y = "medianCount", 
                                          fill = "Clusters",
                                          add = "boxplot", 
                                          xlab = "Stage", 
                                          ylab = "Median Count Values",
                                          palette = coloursBarPlot[1:4]) + 
                                          stat_compare_means(comparisons = ComparisonOptionsStage) + 
                                          # Add pairwise comparisons p-value
                                          stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianCount_StageCluster", ImageName, ".png")) 
      
      
      # Plotting median beta value for patients in each cluster, seperated by translocation
      # Setting the number of comparisons for ggpubr
      if(length(unique(RNAseqCluster$stage[- which(is.na(RNAseqCluster$translocation))])) == 2) {
        ComparisonOptionsTranslocation <- list(names(table(RNAseqCluster$translocation))[1:2])
      } 
      
      p0ClusterTranslocation <- ggpubr::ggboxplot(RNAseqCluster[- which(is.na(RNAseqCluster$translocation)),], 
                                                    x = "translocation", 
                                                    y = "medianCount", 
                                                    fill = "Clusters",
                                                    add = "boxplot", 
                                                    xlab = "Translocation Status", 
                                                    ylab = "Median Count Values",
                                                    palette = coloursBarPlot[1:4]) + 
                                                    stat_compare_means(comparisons = ComparisonOptionsTranslocation) + 
                                                    # Add pairwise comparisons p-value
                                                    stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianCount_TranslocationCluster", ImageName, ".png")) 
      
      
      if(length(unique(RNAseqCluster$sex[- which(is.na(RNAseqCluster$sex))])) == 2) {
        ComparisonOptionsSex <- list(names(table(RNAseqCluster$sex))[1:2])
      } 
      p0ClusterSex <- ggpubr::ggboxplot(RNAseqCluster[- which(is.na(RNAseqCluster$sex)),], 
                                                    x = "sex", 
                                                    y = "medianCount", 
                                                    fill = "Clusters",
                                                    add = "boxplot", 
                                                    xlab = "Sex", 
                                                    ylab = "Median Count Values",
                                                    palette = coloursBarPlot[1:4]) + 
                                                    stat_compare_means(comparisons = ComparisonOptionsSex) + 
                                                    # Add pairwise comparisons p-value
                                                    stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianCount_SexCluster", ImageName, ".png")) 
      
      
      if(length(unique(RNAseqCluster$sitebiopsy[- which(is.na(RNAseqCluster$sitebiopsy))])) == 2) {
        ComparisonOptionsSitebiopsy <- list(names(table(RNAseqCluster$sitebiopsy))[1:2])
      } 
      p0ClusterSitebiopsy <- ggpubr::ggboxplot(RNAseqCluster[- which(is.na(RNAseqCluster$sitebiopsy)),], 
                                        x = "sitebiopsy", 
                                        y = "medianCount", 
                                        fill = "Clusters",
                                        add = "boxplot", 
                                        xlab = "Site biopsy", 
                                        ylab = "Median Count Values",
                                        palette = coloursBarPlot[1:4]) + 
                                        stat_compare_means(comparisons = ComparisonOptionsSitebiopsy) + 
                                        # Add pairwise comparisons p-value
                                        stat_compare_means(paired = FALSE)
      ggplot2::ggsave(paste0(pathNow, "/img/33_medianCount_SitebiopsyCluster", ImageName, ".png")) 
      
      
      
      if(length(unique(RNAseqCluster$institution)) == 3) {
        ComparisonOptionsInstitution <- list(names(table(RNAseqCluster$institution))[1:3])
      } 
      p0ClusterInstitute <- ggpubr::ggboxplot(RNAseqCluster, 
                                               x = "institution", 
                                               y = "medianCount", 
                                               fill = "Clusters",
                                               add = "boxplot", 
                                               xlab = "Institute", 
                                               ylab = "Median Count Values",
                                               palette = coloursBarPlot[1:4]) + 
                                               stat_compare_means(comparisons = ComparisonOptionsInstitution) + 
                                               # Add pairwise comparisons p-value
                                               stat_compare_means(paired = FALSE)
       ggplot2::ggsave(paste0(pathNow, "/img/33_medianCount_InstituteCluster", ImageName, ".png")) 
      
       
       if(length(unique(RNAseqCluster$typebiopsy)) == 2) {
         ComparisonOptionsTypebiopsy <- list(names(table(RNAseqCluster$typebiopsy))[1:2])
       } 
       p0ClusterTypebiopsy<- ggpubr::ggboxplot(RNAseqCluster, 
                                                x = "typebiopsy", 
                                                y = "medianCount", 
                                                fill = "Clusters",
                                                add = "boxplot", 
                                                xlab = "Type biopsy", 
                                                ylab = "Median Count Values",
                                                palette = coloursBarPlot[1:4]) + 
                                                stat_compare_means(comparisons = ComparisonOptionsTypebiopsy) + 
                                                # Add pairwise comparisons p-value
                                                stat_compare_means(paired = FALSE)
       ggplot2::ggsave(paste0(pathNow, "/img/33_medianCount_TypebiopsyCluster", ImageName, ".png")) 
       
       
       if(length(unique(RNAseqCluster$grade[- which(is.na(RNAseqCluster$grade))])) == 2) {
         ComparisonOptionsGrade <- list(names(table(RNAseqCluster$grade))[1:2])
       } 
       p0ClusterSitebiopsy <- ggpubr::ggboxplot(RNAseqCluster[- which(is.na(RNAseqCluster$grade)),], 
                                                x = "grade", 
                                                y = "medianCount", 
                                                fill = "Clusters",
                                                add = "boxplot", 
                                                xlab = "Grade", 
                                                ylab = "Median Count Values",
                                                palette = coloursBarPlot[1:4]) + 
                                                stat_compare_means(comparisons = ComparisonOptionsGrade) + 
                                                # Add pairwise comparisons p-value
                                                stat_compare_means(paired = FALSE)
       ggplot2::ggsave(paste0(pathNow, "/img/33_medianCount_GradeCluster", ImageName, ".png")) 
       
       
       ContamCluster <- data.frame(rrnacontam = QCMatrix$rRNAcontam, 
                                   cluster = ClusterLabels)
       ContamPlot <- ggpubr::ggboxplot(ContamCluster, x = "cluster", y = "rrnacontam", fill = "cluster",
                                       add = "boxplot", ylab=" rRNA contamination level",
                                       font.label = list(size = 20, color = "black"), 
                                       palette = coloursBarPlot[sort(unique(ClusterLabelsSNF[[Clusters]]))]) +
                                       ggtitle("SNF cluster vs. rRNA contamination") +
                                       stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                       # Add pairwise comparisons p-value
                                       stat_compare_means(paired = FALSE)
       # kruskal.test(ContamCluster$rrnacontam ~ ContamCluster$cluster, data = ContamCluster)
       ggsave(paste0(pathNow, "/img/33_rRNAContaminationPlot_", ImageName, ".png"))
       
       
       
       uniqReadCluster <- data.frame(uniqReads = QCMatrix$Uniquely.mapped, 
                                     cluster = ClusterLabelsSNF[[Clusters]])
       uniqReadPlot <- ggpubr::ggboxplot(uniqReadCluster, x = "cluster", y = "uniqReads", fill = "cluster",
                                         add = "boxplot", ylab = "uniquely mapped read counts",
                                         font.label = list(size = 20, color = "black"), 
                                         palette = coloursBarPlot[sort(unique(ClusterLabelsSNF[[Clusters]]))]) +
                                         ggtitle("SNF cluster vs. uniquely mapped reads") +
                                         stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                         # Add pairwise comparisons p-value
                                         stat_compare_means(paired = FALSE)
       ggsave(paste0(pathNow, "/img/33_uniquelyMappedReads_", ImageName, ".png"))

    }
    if(ProduceImages == "Yes") {
      RNAseqHeatmap()
    }
    # tSNE
    cat("\n tSNE for", Clusters, "clusters:\n")
    if(ProduceImages == "Yes") {  
      tSNEOutputFilteredBeta <- tSNEPlotGeneration(BetaMatrix = BetaMatrix, 
                                                   PerplexityParameter = 6, 
                                                   ClinicalFile = ClinicalFile, 
                                                   ClusterLabels = ClusterLabels, 
                                                   Filter = NA, 
                                                   FigureGenerate = "Yes", 
                                                   PNGorPDF = "png", 
                                                   ImageName = "SNFClust_Beta_AllSamples")
      
      tSNEOutputFilteredMvalue <- tSNEPlotGeneration(BetaMatrix = MvalueMatrix, 
                                                     PerplexityParameter = 6, 
                                                     ClinicalFile = ClinicalFile, 
                                                     ClusterLabels = ClusterLabels, 
                                                     Filter = NA, 
                                                     FigureGenerate = "Yes", 
                                                     PNGorPDF = "png", 
                                                     ImageName = "SNFClust_MvalueMatrix_AllSamples")
      
      tSNEOutputFilteredRNAseq <- tSNEPlotGeneration(BetaMatrix = RNAseqCountMatrix, 
                                                     PerplexityParameter = 6, 
                                                     ClinicalFile = ClinicalFile, 
                                                     ClusterLabels = ClusterLabels, 
                                                     Filter = NA, 
                                                     FigureGenerate = "Yes", 
                                                     PNGorPDF = "png", 
                                                     ImageName = "SNFClust_RNAseq_AllSamples")
      
      
      tSNEOutputFiltered <- tSNEPlotGeneration(BetaMatrix = WCombinedFiltered, 
                                               PerplexityParameter = 6, 
                                               ClinicalFile = ClinicalFile, 
                                               ClusterLabels = ClusterLabels, 
                                               Filter = NA, 
                                               FigureGenerate = "Yes", 
                                               PNGorPDF = "png", 
                                               ImageName = "SNFClust_Fused_AllSamples")
    
    }
    
    cat("\n Saving results")
    
    RESULTS <- list(NumberOfClusters = NumberOfClusters,
                    SNFClusterLabels = data.frame(Cluster = ClusterLabelsSNF[[Clusters]],
                                                  SampleID = colnames(WCombinedFiltered)),
                    RNASeqCountMatrixNormalized = Data2RNAseqFilteredNorm,
                    MvalueMatrixNormalized = Data1MvalueFilteredNorm,
                    BetaMatrixNormalized = Data1BetaFilteredNorm,
                    estimateNumberOfClustersGivenGraph = estimationClusters,
                    FinalFusedMatrix = WCombinedFiltered,
                    EstimatedNumberOfClusters = estimationClusters,
                    MethylationDensityOutput = MethylationDensityOutput,
                    ProportionVisualizationOutput = ProportionVisualizationOutput,
                    DifferentialExpressionRNAseqOutput = DifferentialExpressionRNAseqOutput,
                    DiffMethylatedRegionsOutput = DiffMethylatedRegionsOutput)
    
    #SurvivalAnalysisOutput = SurvivalAnalysisOutput)
  } else {
    stop("\n NumberOfClusters  is not specified.")
  }
  
  
  class(RESULTS) <- "SNFClusteringFunction_ASilva"
  return(RESULTS)
  
}
