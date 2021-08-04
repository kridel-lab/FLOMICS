# Updated 3 Aug 2021
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
# TargetedSeqMatrix: A matrix of 0 and 1 values for gene x patients, with genes in rows and patients as columns.
#            The ncol(TargetedSeqMatrix) should match ncol(RNAseqCountMatrix), ncol(MvalueMatrix), nrow(ClinicalFile),
#            nrow(QCMatrix) and ncol(BetaMatrix). Default to NA, if this matrix is not provided. 
# MethylationAnnotationFile: A matrix of annotations for all the probes found in MvalueMatrix/BetaMatrix. It is of 
#                 size probes x annotations.
# RNAseqAnnotationFile: A matrix of annotations for all the ENSEMBL IDs found in RNAseq matrix. It is of 
#                 size ENSEMBL IDs x annotations.
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


SNFClustering33 <- function(RNAseqCountMatrix, 
                            NormalizeRNAseqCountMatrix = "Yes",
                            MvalueMatrix = NA, 
                            BetaMatrix = NA, 
                            TargetedSeqMatrix = NA,
                            MethylationAnnotationFile,
                            RNAseqAnnotationFile = NA,
                            QCMatrix,
                            TumorPurity,
                            ClinicalFile, 
                            SurvivalFile,
                            NumberOfClusters = NA,
                            ImageName = paste0("Image_", date()),
                            PNGorPDF = "PNG",
                            ProduceImages = "No") {
  
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  library("SNFtool")
  library("edgeR")
  library("EnhancedVolcano")
  library("reshape2")
  library("fpc")
  
  
  # Checking user input - are NA values present within matrices?
  if(all(is.na(MvalueMatrix)) != TRUE) {
    if (length(which(is.na(rowSums(MvalueMatrix)) == TRUE)) > 0) {
      stop("\n MvalueMatrix has NA values. Fix and rerun.")
    }
  }
  
  if(all(is.na(RNAseqCountMatrix)) != TRUE) {
    if (length(which(is.na(rowSums(RNAseqCountMatrix)) == TRUE)) > 0) {
      stop("\n RNAseqCountMatrix has NA values. Fix and rerun.")
    }
  } 
  
  if(all(is.na(TargetedSeqMatrix)) != TRUE) {
    if (length(which(is.na(rowSums(TargetedSeqMatrix)) == TRUE)) > 0) {
      stop("\n TargetedSeqMatrix has NA values. Fix and rerun.")
    }
    if(length(which(rowSums(TargetedSeqMatrix) == 0)) > 0) {
      stop("\n TargetedSeqMatrix has genes with 0 mutations. Fix and rerun.")
    } 
  } 
  
  
  
  # Checking user input - sizes should be the same
  if(all(is.na(MvalueMatrix)) != TRUE & all(is.na(RNAseqCountMatrix)) != TRUE) {
    if(ncol(RNAseqCountMatrix) != ncol(MvalueMatrix)) {
      stop("\n ncol(RNAseqCountMatrix) and ncol(MvalueMatrix) should match.")
    } 
    
    # Checking user input - sample names should be the same
    if(! identical(colnames(RNAseqCountMatrix), colnames(MvalueMatrix))) {
      stop("\n Sample names in RNAseqCountMatrix and MvalueMatrix should match.")
    } 
  } 
  
  if(all(is.na(BetaMatrix)) != TRUE & all(is.na(RNAseqCountMatrix)) != TRUE) {

    if (ncol(RNAseqCountMatrix) != ncol(BetaMatrix)) {
      stop("\n ncol(RNAseqCountMatrix) and ncol(BetaMatrix) should match.")
    } 
    
    if (! identical(colnames(RNAseqCountMatrix), colnames(BetaMatrix))) {
      stop("\n Sample names in RNAseqCountMatrix and BetaMatrix should match.")
    }
  }
  
  if(all(is.na(RNAseqCountMatrix)) != TRUE) {
    if (ncol(RNAseqCountMatrix) != nrow(ClinicalFile)) {
      stop("\n ncol(RNAseqCountMatrix) and nrow(ClinicalFile) should match.")
    } 
    
    if (ncol(RNAseqCountMatrix) != nrow(QCMatrix)) {
      stop("\n ncol(RNAseqCountMatrix) and nrow(QCMatrix) should match.")
    } 
  }
  
  
  if(all(is.na(TargetedSeqMatrix)) != TRUE & all(is.na(RNAseqCountMatrix)) != TRUE) {
   if (ncol(TargetedSeqMatrix) != ncol(RNAseqCountMatrix)) {
    stop("\n ncol(TargetedSeqMatrix) and ncol(RNAseqCountMatrix) should match.")
   }
  }
  
  if(all(is.na(RNAseqCountMatrix)) != TRUE) {
    # Checking user input - sample names should be the same
    if (! identical(colnames(RNAseqCountMatrix), ClinicalFile$SAMPLE_ID)) {
      stop("\n Sample names in RNAseqCountMatrix and ClinicalFile should match.")
    } 
    
    if (! identical(colnames(RNAseqCountMatrix), QCMatrix$SAMPLE_ID)) {
      stop("\n  Sample names in RNAseqCountMatrix and QCMatrix should match.")
    } 
  }
  
  if(all(is.na(TargetedSeqMatrix)) != TRUE & all(is.na(RNAseqCountMatrix)) != TRUE) {
    if (! identical(colnames(TargetedSeqMatrix), QCMatrix$SAMPLE_ID)) {
      stop("\n Sample names in TargetedSeqMatrix and QCMatrix should match.")
    }
  }
  
  if(all(is.na(RNAseqCountMatrix)) != TRUE) {
    if(NormalizeRNAseqCountMatrix == "Yes") {
      if(all(is.na(RNAseqCountMatrix)) != TRUE) {
        if(is.integer(RNAseqCountMatrix) != TRUE) {
          stop("\n RNAseqCountMatrix should be a count matrix of non-normaized values.")
        }
      }
    }
  }
  
  # Normalization of data
  set.seed(1234)
  
  if(all(is.na(MvalueMatrix)) != TRUE) {
    # Methylation normalization
    Data1MvalueFilteredNorm <- SNFtool::standardNormalization(t(MvalueMatrix))
    # Data1MvalueFilteredNorm<- t(MvalueMatrix)
    # dim(Data1MvalueFilteredNorm) # 132 5000
  }
  
  if(all(is.na(BetaMatrix)) != TRUE) {
    Data1BetaFilteredNorm <- SNFtool::standardNormalization(t(BetaMatrix))
    # dim(Data1BetaFilteredNorm) # 132 5000
  }
  
  # RNAseq Normalization
  # As per Dr. Bo Wang's response from 27 Oct 2018
  # "I would recommend to first take transformation of  (X = log2(1 + counts)), 
  # and then apply standardNormalization()"
  if(all(is.na(RNAseqCountMatrix)) != TRUE) {
    if(NormalizeRNAseqCountMatrix == "Yes") {
      if(all(is.na(RNAseqCountMatrix)) != TRUE) {
        Data2RNAseqFilteredNorm <- SNFtool::standardNormalization(log2(1 + t(RNAseqCountMatrix)))
      }
    } else {
      Data2RNAseqFilteredNorm <- t(RNAseqCountMatrix)
      # Data2RNAseqFilteredNorm <- SNFtool::standardNormalization(log2(1 + t(RNAseqCountMatrix)))
    }
  }
  
  # Targeted DNAseq data were not normalized as NA values result - doesn't work as NAs result
  # So proceed without normalization
  # Data3MutationFilteredNorm <- SNFtool::standardNormalization(t(TargetedSeqMatrix))
  # Data3MutationFilteredNorm <- (t(TargetedSeqMatrix))
  # range(Data3MutationFilteredNorm) # -1.272405 10.296080
  # range(Data3MutationFilteredNorm) # -1.272405 10.296080
  Data3MutationFilteredNorm <- SNFtool::standardNormalization(log2(1 + t(TargetedSeqMatrix)))
  # Data3MutationFilteredNorm <- SNFtool::standardNormalization(as.matrix(t(TargetedSeqMatrix)))
  # Data3MutationFilteredNorm <- as.matrix(t(TargetedSeqMatrix))
  
  
  
  # Setting number of Parameters for SNF
  if(all(is.na(RNAseqCountMatrix)) != TRUE) {
    K <- round(ncol(RNAseqCountMatrix) / 10); # number of neighbors; Number of patients / Number of clusters ;
  } else if(all(is.na(BetaMatrix)) != TRUE) {
    K <- round(ncol(BetaMatrix) / 10); # number of neighbors; Number of patients / Number of clusters ;
  } else if(all(is.na(TargetedSeqMatrix)) != TRUE) {
    K <- round(ncol(TargetedSeqMatrix) / 10); # number of neighbors; Number of patients / Number of clusters ;
  }
  
  
  cat("\n Based on sample size, SNF recommended number of clusters/neigbours is", K)
  #       If C is unknown, set K = N/10 usually (10~30);
  #       In our case N = 122/10 = 12
  alpha <- 0.5; # hyperparameter, usually (0.3~0.8)
  T <- 20; # Number of Iterations, usually (10~20)
  
  
  ## Calculate the Pairwise squared Euclidean distances
  # Computes the squared Euclidean distances between all pairs of data point given
  # If the data is continuous, we recommend to use the function "dist2" as follows
  
  if(all(is.na(MvalueMatrix)) != TRUE) {
    Dist1Data1MvalueFilteredNorm <- (SNFtool::dist2(as.matrix(Data1MvalueFilteredNorm), 
                                                    as.matrix(Data1MvalueFilteredNorm))) ^ (1 / 2)
    # dim(Dist2Data2RNAseqFilteredNorm) # sample x sample
  }
  
  if(all(is.na(BetaMatrix)) != TRUE) {
    Dist2Data2BetaFilteredNorm <- (SNFtool::dist2(as.matrix(Data1BetaFilteredNorm), 
                                                    as.matrix(Data1BetaFilteredNorm))) ^ (1 / 2)
    # dim(Dist2Data2BetaFilteredNorm) # sample x sample
  }
  # dim(Dist1Data1MvalueFilteredNorm) # sample x sample
  
  if(all(is.na(RNAseqCountMatrix)) != TRUE) {
    Dist2Data2RNAseqFilteredNorm <- (SNFtool::dist2(as.matrix(Data2RNAseqFilteredNorm), 
                                                    as.matrix(Data2RNAseqFilteredNorm))) ^ (1 / 2)
  }
  
  if(all(is.na(TargetedSeqMatrix)) != TRUE) {
    # Error - NAs result
    Dist3Data3DNAseqFilteredNorm <- (SNFtool::dist2(as.matrix(Data3MutationFilteredNorm), 
                                                    as.matrix(Data3MutationFilteredNorm))) ^ (1 / 2)
    
    # dim(Dist3Data3BetaFilteredNorm) # sample x sample
    # range(Dist3Data3DNAseqFilteredNorm, na.rm = TRUE) # NA ; 0.00000 21.89594
    
    # Dist3Data3DNAseqFilteredNorm <- (SNFtool::dist2(as.matrix(t(TargetedSeqMatrix)), 
     #                                               as.matrix(t(TargetedSeqMatrix))))^(1/2)
    # dim(Dist3Data3DNAseqFilteredNorm) # sample x sample
    # range(Dist3Data3DNAseqFilteredNorm, na.rm = TRUE) # NA ;  0.000000 4.582576
  }
  

  # Construct similarity graphs
  # Returns an affinity matrix that represents the neighborhood graph of the data points
  if(all(is.na(MvalueMatrix)) != TRUE) {
    W1Data1MvalueFilteredNorm <- SNFtool::affinityMatrix(Dist1Data1MvalueFilteredNorm, K, alpha)
    # dim(W1Data1MvalueFilteredNorm) # sample x sample
  }
  
  if(all(is.na(BetaMatrix)) != TRUE) {
    W2Data2BetaFilteredNorm <- SNFtool::affinityMatrix(Dist2Data2BetaFilteredNorm, K, alpha)
    # dim(W2Data2BetaFilteredNorm) # sample x sample
  }
  
  if(all(is.na(RNAseqCountMatrix)) != TRUE) {
    W2Data2RNAseqFilteredNorm <- SNFtool::affinityMatrix(Dist2Data2RNAseqFilteredNorm, K, alpha)
  }
  # dim(W2Data2RNAseqFilteredNorm) # sample x sample  # 1.436719e-09 2.738448e-05
  
  if(all(is.na(TargetedSeqMatrix)) != TRUE) {
    W3Data3DNAseqFilteredNorm <- SNFtool::affinityMatrix(Dist3Data3DNAseqFilteredNorm, K, alpha)
  }
  # dim(W3Data3DNAseqFilteredNorm) # must have a positive length 

  
  ## Fuse all the graphs
  ## then the overall matrix can be computed by similarity network fusion(SNF):
  # Similarity Network Fusion takes multiple views of a network and fuses them
  # together to construct an overall status matrix. The input to our algorithm
  # can be feature vectors, pairwise distances, or pairwise similarities. 
  # The learned status matrix can then be used for retrieval, clustering, and 
  # classification.
  set.seed(1234)
  if(all(is.na(MvalueMatrix)) == TRUE) { # If M matrix is not present 
    cat("\n Mvalue matrix is not present ")
    WCombinedFiltered <- SNFtool::SNF(Wall = list(W2Data2RNAseqFilteredNorm, 
                                                  W3Data3DNAseqFilteredNorm), 
                                      K = K, t = T)
    
  } else if(all(is.na(TargetedSeqMatrix)) == TRUE) {  # if targeted seq not present
    cat("\n targSeq matrix is not present ")
    WCombinedFiltered <- SNFtool::SNF(Wall = list(W1Data1MvalueFilteredNorm, 
                                                  W2Data2RNAseqFilteredNorm), 
                                      K = K, t = T)
  } else if(all(is.na(RNAseqCountMatrix)) == TRUE) { # if RNAaseq not present
    cat("\n RNAseq matrix is not present ")
    WCombinedFiltered <- SNFtool::SNF(Wall = list(W1Data1MvalueFilteredNorm, 
                                                  W3Data3DNAseqFilteredNorm), 
                                      K = K, t = T)
  } else {
    cat("\n All matrices present ")
    WCombinedFiltered <- SNFtool::SNF(Wall = list(W1Data1MvalueFilteredNorm, 
                                                  W2Data2RNAseqFilteredNorm, 
                                                  W3Data3DNAseqFilteredNorm), 
                                      K = K, t = T)
  }
  # write.csv(WCombinedFiltered, file = "WCombinedFiltered_10Feb2021.csv")

  
  ## Here we provide two ways to estimate the number of clusters. Note that,
  ## these two methods cannot guarantee the accuracy of esstimated number of
  ## clusters, but just to offer two insights about the datasets.
  estimationClusters <- SNFtool::estimateNumberOfClustersGivenGraph(W = WCombinedFiltered, 
                                                                    NUMC = 2:K)
  
  print(estimationClusters)
  
  
  EigenGapValue <- function() {
  # calculate eigen gap
  # install.packages("Spectrum")
  # library(Spectrum)
    EigenGapValue <- Spectrum::Spectrum(data = WCombinedFiltered, 
                                fontsize = 12, dotsize = 5, 
                                method = 2,
                                maxk = K)
                                #showpca = TRUE)
    
    if(ProduceImages == "Yes") {  
      pathNow <- getwd()
      ggplot2::ggplot(data.frame(EigenGapValue$eigenvector_analysis), 
                      aes(x = K, y = Z)) +
                      geom_point() +
                      labs(x = "Cluster (k)", y = "Eigen Gap Value") +
                      theme_bw() + theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         text = element_text(size=15)) +
                      scale_x_discrete(limits = c(1:(K+1)))
      ggplot2::ggsave(paste0(pathNow, "/img/33_EigenGapPlot", ImageName, ".", PNGorPDF))
      
    }
  }
  NumberOfClusters <- estimationClusters$`Eigen-gap best`
  
  # Perform other analyses
  if(NumberOfClusters == 1) {
    RESULTS <- list(NumberOfClusters = NumberOfClusters,
                    RNASeqCountMatrixNormalized = Data2RNAseqFilteredNorm,
                    MvalueMatrixNormalized = Data1MvalueFilteredNorm,
                    BetaMatrixNormalized = Data1BetaFilteredNorm,
                    estimateNumberOfClustersGivenGraph = estimationClusters)
    # SurvivalAnalysisOutput = SurvivalAnalysisOutput)
  } else if (NumberOfClusters > 1) {
    # Running spectral clustering for specified number of clusters
    ClusterLabelsSNF <- list() # A list for saving results
    
    for (Clusters in NumberOfClusters:NumberOfClusters) {
      
      ClusterLabelsSNF[[Clusters]] <- SNFtool::spectralClustering(affinity = WCombinedFiltered, 
                                                                  K = Clusters) 
      
      
      cat("\n Tabulation based on", Clusters, "clusters \n")
      print(table(ClusterLabelsSNF[[Clusters]]))
      
      ClusterLabels <- ClusterLabelsSNF[[Clusters]]
      
      cat("\n Results based on", Clusters, "clusters for TYPE \n")
      print(table(ClusterLabels, ClinicalFile$TYPE))
      print(chisq.test(table(ClusterLabels, ClinicalFile$TYPE))) 
      # X-squared = 
      
      # printing labels against clinical classes
      # table(rpmmClass_numbers,ClinicalFile$TYPE)
      cat("\n Results based on", Clusters, "clusters for stage \n")
      print(table(ClusterLabels, ClinicalFile$STAGE))
      # Chi-squared Test of Independence
      # Null: cluster membership is independent of stage at .05 significance level.
      print(chisq.test(table(ClusterLabels, ClinicalFile$STAGE))) 
      print(fisher.test(table(ClusterLabels, ClinicalFile$STAGE))) 
      
      
      cat("\n Results based on", Clusters, "clusters for sex \n")
      print(table(ClusterLabels, ClinicalFile$SEX))
      print(chisq.test(table(ClusterLabels, ClinicalFile$SEX))) 
      
      
      cat("\n Results based on", Clusters, "clusters for site biopsy \n")
      print(table(ClusterLabels, trimws(ClinicalFile$SITE_BIOPSY)))
      print(chisq.test(table(ClusterLabels, trimws(ClinicalFile$SITE_BIOPSY)))) 
      
      
      cat("\n Results based on blcTree for type biopsy \n")
      print(table(ClusterLabels, trimws(ClinicalFile$TYPE_BIOPSY)))
      print(chisq.test(table(ClusterLabels, trimws(ClinicalFile$TYPE_BIOPSY)))) 
      
      
      cat("\n Results based on", Clusters, "clusters for institution \n")
      print(table(ClusterLabels, ClinicalFile$INSTITUTION))
      print(chisq.test(table(ClusterLabels, ClinicalFile$INSTITUTION))) 
    
      
      cat("\n Results based on", Clusters, "clusters for 14_18 translocation \n")
      if(length(which(is.na(ClinicalFile$TRANSLOC_14_18) == TRUE)) > 0) {
        NAentries <- which(is.na(ClinicalFile$TRANSLOC_14_18) == TRUE)
        print(table(ClusterLabels[- NAentries], 
                    ClinicalFile$TRANSLOC_14_18[- NAentries]))
        print(chisq.test((table(ClusterLabels[-NAentries], 
                                ClinicalFile$TRANSLOC_14_18[- NAentries])))) 
      } else {
        print(table(ClusterLabels, ClinicalFile$TRANSLOC_14_18) )
        print(chisq.test((table(ClusterLabels, ClinicalFile$TRANSLOC_14_18)))) 
      }
      
      # if(all(is.na(TargetedSeqMatrix)) != TRUE) {
      #   print(table(ClusterLabels, 
      #     TargetedSeqMatrix[which(rownames(TargetedSeqMatrix) == "BCL2_BA_consensus"), ]))
      #   print(chisq.test(table(ClusterLabels, 
      #                          TargetedSeqMatrix[which(rownames(TargetedSeqMatrix) == "BCL2_BA_consensus"), ])))
      # }
      
      cat("\n Results based on blcTree for grade \n")
      if(all(is.na(RNAseqCountMatrix)) != TRUE) {
        matchedIDs <- match(substr(colnames(RNAseqCountMatrix[, which(substr(colnames(RNAseqCountMatrix), 4, 5) == "FL")]), 1, 9), 
                            SurvivalFile$LY_FL_ID)
      } else if(all(is.na(BetaMatrix)) != TRUE) {
        matchedIDs <- match(substr(colnames(BetaMatrix[, which(substr(colnames(BetaMatrix), 4, 5) == "FL")]), 1, 9), 
                            SurvivalFile$LY_FL_ID)        
      }
      if (sum(is.na(matchedIDs)) > 0) {
        # if NAs are present
        SurvivalFileOrderedbyRNAseqPatients <- SurvivalFile[matchedIDs[- which(is.na(matchedIDs) == TRUE)], ]
      } else {
        SurvivalFileOrderedbyRNAseqPatients <- SurvivalFile[matchedIDs, ]
      }
      
      # recategorize grade
      # response on SLACK from Robert on 12 May 2020
      # It's best to regroup 1 and 2 into 1-2 category,
      # 3 and 3A into 3A category, and in situ into 1-2 category.
      # In situ means that a lymph node is not pathologically enlarged, nor otherwise abnormal, 
      # other than incidental findings of BCL2 positive B cells within germinal centre. However, 
      # if lymph node, or any lymph node was enlarged, than should be called low-grade FL, 
      # hence suggestion to choose 1-2.
      recatGRADE <- SurvivalFileOrderedbyRNAseqPatients$GRADE
      recatGRADE[which(SurvivalFileOrderedbyRNAseqPatients$GRADE == "1")] <- "1-2"
      recatGRADE[which(SurvivalFileOrderedbyRNAseqPatients$GRADE == "2")] <- "1-2"
      recatGRADE[which(SurvivalFileOrderedbyRNAseqPatients$GRADE == "IN_SITU")] <- "1-2"
      recatGRADE[which(SurvivalFileOrderedbyRNAseqPatients$GRADE == "3")] <- "3A"
      
      Clusters <- NumberOfClusters
      if(all(is.na(RNAseqCountMatrix)) != TRUE) {
        print(table(ClusterLabels[which(substr(colnames(RNAseqCountMatrix), 4, 5) == "FL")], 
                    trimws(recatGRADE)))
        print(chisq.test(table(ClusterLabels[which(substr(colnames(RNAseqCountMatrix), 4, 5) == "FL")], 
                               trimws(recatGRADE)))) 
      } else if (all(is.na(BetaMatrix)) != TRUE) {
        print(table(ClusterLabels[which(substr(colnames(BetaMatrix), 4, 5) == "FL")], 
                    trimws(recatGRADE)))
        print(chisq.test(table(ClusterLabels[which(substr(colnames(BetaMatrix), 4, 5) == "FL")], 
                               trimws(recatGRADE)))) 
      }
      
    }
  }
    
  SNFNMI33(MvalueMatrix = MvalueMatrix,
           TargetedSeqMatrix = TargetedSeqMatrix,
           RNAseqCountMatrix = RNAseqCountMatrix,
           W1Data1MvalueFilteredNorm = W1Data1MvalueFilteredNorm, 
           W3Data3DNAseqFilteredNorm = W3Data3DNAseqFilteredNorm,
           W2Data2RNAseqFilteredNorm = W2Data2RNAseqFilteredNorm,
           WCombinedFiltered = WCombinedFiltered,
           NumberOfClusters = NumberOfClusters,
           ProduceImages = ProduceImages,
           ImageName = "NMI",
           PNGorPDF = PNGorPDF)
  
  SNFFusedMatrixPlot33(ClusterLabels = ClusterLabels,
                       MvalueMatrix = MvalueMatrix,
                       BetaMatrix = BetaMatrix,
                       RNAseqCountMatrix = RNAseqCountMatrix,
                       TargetedSeqMatrix = TargetedSeqMatrix,
                       W1Data1MvalueFilteredNorm = W1Data1MvalueFilteredNorm,
                       W3Data3DNAseqFilteredNorm = W3Data3DNAseqFilteredNorm,
                       W2Data2RNAseqFilteredNorm = W2Data2RNAseqFilteredNorm,
                       W2Data2BetaFilteredNorm = W2Data2BetaFilteredNorm,
                       WCombinedFiltered = WCombinedFiltered,
                       ProduceImages = ProduceImages,
                       PNGorPDF = PNGorPDF,
                       ImageName = "FusedMatrix")
    
    
  cat("\n Saving results")
  
  if(all(is.na(RNAseqCountMatrix)) == TRUE) {
    RNASeqCountMatrixNormalized = NA
  } else {
      RNASeqCountMatrixNormalized = Data2RNAseqFilteredNorm
 }
    
  RESULTS <- list(NumberOfClusters = NumberOfClusters,
                  SNFClusterLabels = data.frame(Cluster = ClusterLabelsSNF[[Clusters]],
                                                  SampleID = colnames(WCombinedFiltered)),
                  RNASeqCountMatrixNormalized = RNASeqCountMatrixNormalized,
                  MvalueMatrixNormalized = Data1MvalueFilteredNorm,
                  BetaMatrixNormalized = Data1BetaFilteredNorm,
                  estimateNumberOfClustersGivenGraph = estimationClusters,
                  FinalFusedMatrix = WCombinedFiltered,
                  EstimatedNumberOfClusters = estimationClusters)
  class(RESULTS) <- "SNFClustering33"
  return(RESULTS)
}   
    
    
    # Calculate cluster validation statistics using fused matrix
    SNFClusterValidation33 <- function(WCombinedFiltered,
                                       ClusterLabels) {
    
    
    
    set.seed(1234)
    WCombinedFilteredDistance <- dist(WCombinedFiltered)
    ClusterStatsSummary <- fpc::cluster.stats(d = WCombinedFilteredDistance, 
                                              clustering = ClusterLabels)
    cat("\n Average silhouette width: ", round(ClusterStatsSummary$avg.silwidth, 3))
    cat("\n Entropy of the distribution of cluster memberships: ", round(ClusterStatsSummary$entropy, 3)) 
    cat("\n Dunn index: ", round(ClusterStatsSummary$dunn, 3)) 
    cat("\n Separation index: ", round(ClusterStatsSummary$sindex, 3)) 
    cat("\n Calinski and Harabasz index: ", round(ClusterStatsSummary$ch, 3))
    
    RESULTS <- list(AverageSilhouetteWidth = round(ClusterStatsSummary$avg.silwidth, 3),
                    Entropy = round(ClusterStatsSummary$entropy, 3),
                    DunnIndex = round(ClusterStatsSummary$dunn, 3),
                    SeparationIndex = round(ClusterStatsSummary$sindex, 3),
                    CalinskiHarabaszIndex = round(ClusterStatsSummary$ch, 3))
    
    class(RESULTS) <- "SNFClusterValidation33"
    return(RESULTS)
    }
    
    # Concordance Network NMI calculation
    SNFNMI33 <- function(MvalueMatrix = NA,
                         TargetedSeqMatrix = NA,
                         RNAseqCountMatrix = NA,
                         W1Data1MvalueFilteredNorm, 
                         W3Data3DNAseqFilteredNorm,
                         W2Data2RNAseqFilteredNorm,
                         WCombinedFiltered,
                         NumberOfClusters,
                         ProduceImages = "No",
                         PNGorPDF,
                         ImageName) {
      
      # Concordance Network NMI calculation matrix
      if(all(is.na(MvalueMatrix)) == TRUE) { # If M matrix is not present 
        SNFdata <- list(RNAseq = W2Data2RNAseqFilteredNorm,
                        DNAseq = W3Data3DNAseqFilteredNorm,
                        Fused = WCombinedFiltered)
        
      } else if(all(is.na(TargetedSeqMatrix)) == TRUE) {  # if targeted seq not present
        SNFdata <- list("DNA methylation" = W1Data1MvalueFilteredNorm, 
                        RNAseq = W2Data2RNAseqFilteredNorm,
                        Fused = WCombinedFiltered)
      } else if(all(is.na(RNAseqCountMatrix)) == TRUE) { # if RNAaseq not present
        SNFdata <- list("DNA methylation" = W1Data1MvalueFilteredNorm, 
                        DNAseq = W3Data3DNAseqFilteredNorm,
                        Fused = WCombinedFiltered)
      } else {
        SNFdata <- list("DNA methylation" = W1Data1MvalueFilteredNorm, 
                        RNAseq = W2Data2RNAseqFilteredNorm,
                        DNAseq = W3Data3DNAseqFilteredNorm,
                        Fused = WCombinedFiltered)
      }
      
      
      ConcordanceMatrix <- SNFtool::concordanceNetworkNMI(Wall = SNFdata, C = NumberOfClusters)
      colnames(ConcordanceMatrix) <- rownames(ConcordanceMatrix) <- names(SNFdata)
      
      if(ProduceImages == "Yes") {
        
        # http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
        # Get lower triangle of the correlation matrix
        getLowerTri <- function(cormat) {
          cormat[upper.tri(cormat)] <- NA
          return(cormat)
        }
        # Get upper triangle of the correlation matrix
        getUpperTri <- function(cormat) {
          cormat[lower.tri(cormat)] <- NA
          return(cormat)
        }
        
        upperTri <- round(getUpperTri(cormat = ConcordanceMatrix), 4)
        
        meltedCormat <- reshape2::melt(upperTri, na.rm = TRUE)
        # dim(melted_cormat) # 10  3
        
        # Heatmap of NMI
        pathNow <- getwd()
        NMIHeatMap <- ggplot2::ggplot(data = meltedCormat, 
                                      aes(Var2, Var1, fill = value)) +
                                      geom_tile(color = "white") +
                                      geom_text(aes(label = format(round(value, 2), nsmall = 2)), size = 6) + 
                                      labs(x = "", y = "") +
                                      # scale_fill_gradient2(low = "white", high = "brown", mid = "red", 
                                      scale_fill_gradient2(low = "white", high = "#4d4d4d", mid = "#878787", 
                                                           midpoint = 0.5, limit = c(0.000, 1.000), space = "Lab",
                                                           name = "NMI") +
                                      theme_bw(base_rect_size = 1) + 
                                      theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            legend.text = element_text(size = 20),
                                            #legend.key.size = unit(1, "in"),
                                            axis.text.x = element_text(angle = 0, 
                                                                       #vjust = 1, 
                                                                       #hjust = 1, 
                                                                       color = "black",
                                                                       size = 20),
                                            axis.text.y = element_text(color = "black", 
                                                                       size = 20)) +
                                      coord_fixed()
        ggplot2::ggsave(paste0(pathNow, "/img/33_ConcordanceNetworkNMI", ImageName, ".", PNGorPDF),
                        width = 10.5, height = 10.5)
        NMIHeatMap
        RESULTS <- list(NMI = upperTri)
        class(RESULTS) <- "SNFNMI33"
        return(RESULTS)
      }
    }
    
    
    # Printing SNF fused matrix plots
    SNFFusedMatrixPlot33 <- function(ClusterLabels,
                                     MvalueMatrix,
                                     BetaMatrix,
                                     RNAseqCountMatrix,
                                     TargetedSeqMatrix,
                                     W1Data1MvalueFilteredNorm,
                                     W3Data3DNAseqFilteredNorm,
                                     W2Data2RNAseqFilteredNorm,
                                     W2Data2BetaFilteredNorm,
                                     WCombinedFiltered,
                                     ProduceImages = "No",
                                     PNGorPDF = "PNG",
                                     ImageName) {
      
    pathNow <- getwd() # getting the path
    
    if(ProduceImages == "Yes") {
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/33_SNF_FusedMatrix_", ImageName, "_only.pdf"),
            width = 80, height = 80, pointsize = 200)
      } else {
        grDevices::png(paste0(pathNow, "/img/33_SNF_FusedMatrix_", ImageName, "_only.png"))
      }
      par(mfrow = c(1, 1))
      SNFFused <- SNFtool::displayClusters(WCombinedFiltered, group = ClusterLabels)
      grDevices::dev.off()
      
      # write.csv(WCombinedFiltered, file = "WCombinedFiltered_6Jan2021.csv",  row.names = TRUE)
      
      if(all(is.na(MvalueMatrix)) != TRUE) {
        if (PNGorPDF == "pdf") {
          grDevices::pdf(paste0(pathNow, "/img/33_SNF_MValueMatrix_", ImageName, "_only.pdf"),
                         width = 80, height = 80, pointsize = 200)
        } else {
          grDevices::png(paste0(pathNow, "/img/33_SNF_MValueMatrix_", ImageName, "_only.png"))
        }
        par(mfrow = c(1, 1))
        SNFtool::displayClusters(W1Data1MvalueFilteredNorm, group = ClusterLabels)
        grDevices::dev.off()
      }
      
      if(all(is.na(RNAseqCountMatrix)) != TRUE) {
        if (PNGorPDF == "pdf") {
          grDevices::pdf(paste0(pathNow, "/img/33_SNF_RNAseqMatrix_", ImageName, "_only.pdf"),
                         width = 80, height = 80, pointsize = 200)
        } else {
          grDevices::png(paste0(pathNow, "/img/33_SNF_RNAseqMatrix_", ImageName, "_only.png"))
        }
        par(mfrow = c(1, 1))
        SNFtool::displayClusters(W2Data2RNAseqFilteredNorm, ClusterLabels)
        grDevices::dev.off()
      }
      
      
      if(all(is.na(BetaMatrix)) != TRUE) {
        if (PNGorPDF == "pdf") {
          grDevices::pdf(paste0(pathNow, "/img/33_SNF_BetaMatrix_", ImageName, "_only.pdf"),
                         width = 80, height = 80, pointsize = 200)
        } else {
          grDevices::png(paste0(pathNow, "/img/33_SNF_BetaMatrix_", ImageName, "_only.png"))
        }
        par(mfrow = c(1, 1))
        SNFtool::displayClusters(W2Data2BetaFilteredNorm, ClusterLabels)
        grDevices::dev.off()
      }
      
      if(all(is.na(TargetedSeqMatrix)) != TRUE) {
        if (PNGorPDF == "pdf") {
          grDevices::pdf(paste0(pathNow, "/img/33_SNF_DNAseq_", ImageName, "_only.pdf"),
                         width = 80, height = 80, pointsize = 200)
        } else {
          grDevices::png(paste0(pathNow, "/img/33_SNF_DNAseq_", ImageName, "_only.png"))
        }
        par(mfrow = c(1, 1))
        SNFtool::displayClusters(W3Data3DNAseqFilteredNorm, ClusterLabels)
        grDevices::dev.off()
      }
      
      return(invisible(NULL))
      }
    }
      
    
    # Printing SNF Alluvial plots
    SNFAlluvialPlot33 <- function(ClusterLabels,
                                  ClinicalFile,
                                  WCombinedFiltered,
                                  ProduceImages = "No",
                                  PNGorPDF = "PNG",
                                  ImageName) {
      
      if(ProduceImages == "Yes") {
        cat("\n plotAlluvial printing. \n")
        
        Clusters <- max(ClusterLabels)
        
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
      }
      return(invisible(NULL))
    }
      
    
    # Plot purity 
    SNFPurityPlot33 <- function(ClusterLabels,
                                TumorPurity,
                                ProduceImages = "No",
                                PNGorPDF = "PNG",
                                ImageName) {
                                
    if(ProduceImages == "Yes") {
      cat("\n Plot clusters by purity. \n")
      # Plot clusters by purity
      par(mfrow = c(1, 1))
      PurityCluster <- data.frame(purity = TumorPurity, cluster = ClusterLabels)
      
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
      
      
      PurityPlot <- ggpubr::ggboxplot(PurityCluster, x = "cluster", y = "purity", fill = "cluster",
                                      add = "boxplot", ylab = "Tumor Purity", 
                                      font.label = list(size = 20, color = "black"), 
                                      palette = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                      ggtitle("SNF cluster vs. tumor purity") +
                                      stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                    # Add pairwise comparisons p-value
                                     stat_compare_means(paired = FALSE)
      pathNow <- getwd()
      ggplot2::ggsave(paste0(pathNow, "/img/33_PurityPlot_", ImageName, ".", PNGorPDF))
    }
      
    return(invisible(null))
      
    }
      
    # MeanSDPlotS
    SNFMeanSDPlot33 <- function(ClusterLabels,
                                TumorPurity,
                                ClinicalFile,
                                WCombinedFiltered,
                                BetaMatrix,
                                MvalueMatrix,
                                RNAseqCountMatrix,
                                ProduceImages = "No",
                                PNGorPDF = "PNG",
                                ImageName) {
      cat("\n Plot Mean SD plot. \n")
      cat("\n Mean SD plot for", Clusters, "clusters: \n")
      if(ProduceImages == "Yes") {
        MeanSDPlot(BetaMatrix = WCombinedFiltered, 
                   ClinicalFile = ClinicalFile, 
                   ClusterLabels = ClusterLabels, 
                   FigureGenerate = ProduceImages, 
                   PNGorPDF ="png", 
                   ImageName = paste0("33_MeanSDPlot_FusedMatrix_", ImageName))
        
        if(all(is.na(BetaMatrix)) != TRUE) {
          if(nrow(BetaMatrix) < 10000) {
            MeanSDPlot(BetaMatrix = BetaMatrix, 
                       ClinicalFile = ClinicalFile, 
                       ClusterLabels = ClusterLabels, 
                       FigureGenerate = ProduceImages, 
                       PNGorPDF ="png", 
                       ImageName = paste0("33_MeanSDPlot_BetaMatrix_", ImageName)) 
          }
          
          if(nrow(BetaMatrix) < 10000) {
            MeanSDPlot(BetaMatrix = MvalueMatrix, 
                       ClinicalFile = ClinicalFile, 
                       ClusterLabels = ClusterLabels, 
                       FigureGenerate = ProduceImages, 
                       PNGorPDF ="png", 
                       ImageName = paste0("33_MeanSDPlot_MvalueMatrix_", ImageName))
          }
        }
        MeanSDPlot(BetaMatrix = RNAseqCountMatrix, 
                   ClinicalFile = ClinicalFile, 
                   ClusterLabels = ClusterLabels, 
                   FigureGenerate = ProduceImages, 
                   PNGorPDF ="png", 
                   ImageName = paste0("33_MeanSDPlot_RNAseq_", ImageName))
      } 
      return(invisible(NULL))
    }
    
    
    # MethylationDensityPlot
    SNFMethylationDensityPlot33 <- function(ClusterLabels,
                                            TumorPurity,
                                            ClinicalFile,
                                            WCombinedFiltered,
                                            BetaMatrix,
                                            MvalueMatrix,
                                            RNAseqCountMatrix,
                                            MethylationAnnotationFile,
                                            ProduceImages = "No",
                                            PNGorPDF = "PNG") {
    MethylationDensityOutput <- NA
    if(all(is.na(BetaMatrix)) != TRUE) {
      cat("\n Plot methylation density plot. \n")
      cat("\n Methylation density plot for", Clusters, "clusters:\n")
      if(ProduceImages == "Yes") {
        MethylationDensityOutput <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                           ClusterLabels = ClusterLabels, 
                                                           PlotWithinAnnotationCategories = "No", 
                                                           ClinicalCategoryToVisualize = "CLUSTER", 
                                                           BetaMatrix = BetaMatrix, 
                                                           AnnotationFile = MethylationAnnotationFile, 
                                                           ClinicalFile = ClinicalFile, 
                                                           SampleSheet = NA, 
                                                           FigureGenerate = ProduceImages, 
                                                           ImageName = paste0("MethylationPlot_AllProbes", ImageName), 
                                                           PNGorPDF = PNGorPDF)
      } else {
        MethylationDensityOutput <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                           ClusterLabels = ClusterLabels, 
                                                           PlotWithinAnnotationCategories = "No", 
                                                           ClinicalCategoryToVisualize = "CLUSTER", 
                                                           BetaMatrix = BetaMatrix, 
                                                           AnnotationFile = MethylationAnnotationFile, 
                                                           ClinicalFile = ClinicalFile, 
                                                           SampleSheet = NA, 
                                                           FigureGenerate = ProduceImages,
                                                           ImageName = paste0("MethylationPlot_AllProbes", ImageName),
                                                           PNGorPDF = PNGorPDF)
      }
    }
    return(invisible(NULL))
  }
    
    
    # ProportionVisualization
    SNFProportionVisualization33 <- function(ClusterLabels,
                                             TumorPurity,
                                             ClinicalFile,
                                             WCombinedFiltered,
                                             BetaMatrix,
                                             MvalueMatrix,
                                             RNAseqCountMatrix,
                                             MethylationAnnotationFile,
                                             ProduceImages = "No",
                                             PNGorPDF = "PNG", 
                                             ImageName) {
    ProportionVisualizationOutput <- NA
    if(all(is.na(BetaMatrix)) != TRUE) {
      cat("\n Proportion visualization for", Clusters, "clusters:\n")
      if(ProduceImages == "Yes") {
        ProportionVisualizationOutput <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", 
                                                                 ClusterLabels = ClusterLabels, 
                                                                 BetaMatrix = BetaMatrix, 
                                                                 AnnotationFile = MethylationAnnotationFile, 
                                                                 ClinicalFile = ClinicalFile, 
                                                                 ClinicalCategory = "TYPE", 
                                                                 PlotWithinCategories = "Yes", 
                                                                 FigureGenerate = "Yes", 
                                                                 PNGorPDF = "png", 
                                                                 ImageName = paste0("AllProbes_", ImageName))
      } else {
        ProportionVisualizationOutput <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", 
                                                                 ClusterLabels = ClusterLabels, 
                                                                 BetaMatrix = BetaMatrix, 
                                                                 AnnotationFile = MethylationAnnotationFile, 
                                                                 ClinicalFile = ClinicalFile, 
                                                                 ClinicalCategory = "TYPE", 
                                                                 PlotWithinCategories = "Yes", 
                                                                 FigureGenerate = "No")
      }
    }
    return(invisible(NULL))
  }
    
    
    
    # DifferentialExpressionRNAseq
    SNFDifferentialExpressionRNAseq33 <- function(ClusterLabels,
                                             RNAseqCountMatrix,
                                             RNAseqAnnotationFile,
                                             ClinicalFile,
                                             ProduceImages = "No",
                                             PNGorPDF = "PNG") {
    cat("\n Differentially expression analysis for", Clusters, "clusters:\n")
    DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = RNAseqCountMatrix, 
                                                                         ContrastColumnName = "CLUSTER", 
                                                                         AnnotationFileEnsemblToGenes = RNAseqAnnotationFile,
                                                                         ClinicalFile = ClinicalFile, 
                                                                         ClusterLabels = ClusterLabels, 
                                                                         ProduceImages = ProduceImages, 
                                                                         PNGorPDF = PNGorPDF) 
    return(DifferentialExpressionRNAseqOutput)
  }
    
    
    # DifferentialExpressionRNAseq
    SNFDiffMethylatedRegions33 <- function(Clusters,
                                                 ClusterLabels,
                                                 BetaMatrix, 
                                                 MvalueMatrix,
                                                 ClinicalFile,
                                                 MethylationAnnotationFile,
                                                 ProduceImages = "No",
                                                 PNGorPDF = "PNG") {
    DiffMethylatedRegionsOutput <- NA
    if(all(is.na(BetaMatrix)) != TRUE) {
      cat("\n Differentially methylated regions analysis for", Clusters, "clusters:\n")
      if(ProduceImages == "Yes") {
        DiffMethylatedRegionsOutput <- DiffMethylatedRegions(Method = "DMRcate", 
                                                             BetaMatrix = BetaMatrix, 
                                                             MvalueMatrix = MvalueMatrix, 
                                                             ContrastColumnName = "CLUSTER", 
                                                             ClinicalFile = ClinicalFile, 
                                                             ClusterLabels = ClusterLabels, 
                                                             AnnotationFile = MethylationAnnotationFile, 
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
                                                             ClusterLabels = ClusterLabels, 
                                                             AnnotationFile = MethylationAnnotationFile, 
                                                             ProduceImages = "No", 
                                                             DMR = 1, 
                                                             ExpressionFile = NA)
      }
    }
    return(DiffMethylatedRegionsOutput)
  }
  
    
    # SurvivalAnalysis
    SNFSurvivalAnalysis33 <- function(Clusters,
                                      ClusterLabels,
                                      BetaMatrix, 
                                      ClinicalFile,
                                      ProduceImages = "No",
                                      PNGorPDF = "PNG") {
    SurvivalAnalysisOutput <- NA
    if(all(is.na(BetaMatrix)) != TRUE) {
      cat("\n Survival analysis for", Clusters, "clusters:\n")
      if(ProduceImages == ProduceImages) {
        SurvivalAnalysisOutput <- SurvivalAnalysis15(BetaMatrix = BetaMatrix, 
                                                   ClinicalFile = ClinicalFile, 
                                                   SurvivalFile = SurvivalFile, 
                                                   ClusterLabels = ClusteringSNFNoTranslocation52Mutations$SNFClusterLabels$Cluster, 
                                                   FigureGenerate = "Yes", 
                                                   PNGorPDF = PNGorPDF)
      }
    }
    return(invisible(NULL))
  }
    
    
    # PCA Analysis
    SNFPCAAnalysis33 <- function(Clusters,
                                 ClusterLabels,
                                 WCombinedFiltered, 
                                 W2Data2RNAseqFilteredNorm,
                                 W1Data1MvalueFilteredNorm,
                                 W3Data3DNAseqFilteredNorm,
                                 BetaMatrix,
                                 ProduceImages = "No",
                                 PNGorPDF = "PNG",
                                 ImageName) {
    # PCA analysis
    cat("\n PCA analysis for", Clusters, "clusters:\n")
    if(ProduceImages == "Yes") {  
      set.seed(1234)
      res.pcaFiltered <- list()
      rangePCA <- c(1:2)
      # using selected beta probes 
      res.pcaFiltered[[1]] <- tryCatch(prcomp(t(WCombinedFiltered), 
                                center = TRUE, 
                                scale = TRUE), 
                                error = function(err) NA)
      
      res.pcaFiltered[[2]] <- tryCatch(prcomp(t(W2Data2RNAseqFilteredNorm), # RNAseqCountMatrix
                                       center = TRUE, 
                                       scale = TRUE), 
                               error = function(err) NA)
      
      if(all(is.na(BetaMatrix)) != TRUE) {
        res.pcaFiltered[[3]] <- tryCatch(prcomp(t(BetaMatrix),  
                                         center = TRUE, 
                                         scale = TRUE), 
                                  error = function(err) NA)
      
        res.pcaFiltered[[4]] <- tryCatch(prcomp(t(W1Data1MvalueFilteredNorm), # MvalueMatrix
                                       center = TRUE, 
                                       scale = TRUE), 
                                  error = function(err) NA)
        
        res.pcaFiltered[[5]] <- tryCatch(prcomp(t(W3Data3DNAseqFilteredNorm), # TargetedSeqMatrix
                                                center = TRUE, 
                                                scale = FALSE), 
                                         error = function(err) NA)
        rangePCA <- c(1:5)
      }
      
      ListOfCategories <- c("WCombinedFiltered", 
                            "RNAseqCountMatrix", 
                            "BetaMatrix", 
                            "MvalueMatrix")
      
      
      for (category in rangePCA) {
        
        if(all(! is.na(res.pcaFiltered[[category]]))) {
      
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
                                                 habillage = factor(ClusterLabels),
                                                 addEllipses = TRUE, ellipse.level = 0.95,  
                                                 palette = c("#4363d8", "#f58231")) +
            theme_bw(base_rect_size = 1) + 
            labs(title = " ") +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(size = rel(1.5)),
                  axis.title = element_text(size = rel(1.5)),
                  legend.text = element_text(size = rel(1.5)),
                  legend.title = element_text(size = rel(1.5)),
                  aspect.ratio = 1)
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
    }
      return(invisible(NULL))
    }
    
    # Beta plots
    SNFBetaPlots33 <- function(TumorPurity,
                               ClusterLabels,
                               recatGRADE,
                               ClinicalFile,
                               BetaMatrix,
                               ProduceImages = "No",
                               PNGorPDF = "PNG",
                               ImageName) {
    
      
      par(mfrow = c(1, 1))
      PurityCluster <- data.frame(purity = TumorPurity, 
                                  cluster = ClusterLabels,
                                  Clusters = as.character(ClusterLabels),
                                  stage = ClinicalFile$STAGE,# [c(11:165), ]
                                  #translocation = as.vector(t(TargetedSeqMatrix["BCL2_BA_consensus", ])), # [c(11:165), ]
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
      
      pathNow <- getwd()
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
      if(length(unique(PurityCluster$type)) == 1) {
        ComparisonOptionsType <- list(names(table(PurityCluster$type))[1])
      } else if(length(unique(PurityCluster$type)) == 2) {
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
      if(length(unique(PurityCluster$stage[- which(is.na(PurityCluster$stage))])) == 2) {
        ComparisonOptionsStage <- list(names(table(PurityCluster$stage))[1:2])
        
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
      } else {
        ComparisonOptionsStage <- list(names(table(PurityCluster$stage))[1:2])
        
        p0ClusterStage <- ggpubr::ggboxplot(PurityCluster, 
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
      }
      

      
      
      # Plotting median beta value for patients in each cluster, seperated by translocation
      # Setting the number of comparisons for ggpubr
      if(length(unique(PurityCluster$stage[- which(is.na(PurityCluster$translocation))])) == 2) {
        ComparisonOptionsTranslocation <- list(names(table(PurityCluster$translocation))[1:2])
        
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
        
      } else {
        ComparisonOptionsTranslocation <- list(names(table(PurityCluster$translocation))[1:2])
        
        p0ClusterTranslocation <- ggpubr::ggboxplot(PurityCluster, 
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
      }
      
      
      
      
      if(length(unique(PurityCluster$sex[- which(is.na(PurityCluster$sex))])) == 2) {
        ComparisonOptionsSex <- list(names(table(PurityCluster$sex))[1:2])
        
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
      } 
      
      
      if(length(unique(PurityCluster$sitebiopsy[- which(is.na(PurityCluster$sitebiopsy))])) == 2) {
        ComparisonOptionsSitebiopsy <- list(names(table(PurityCluster$sitebiopsy))[1:2])
        
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
      } 
      
      
      
      if(length(unique(PurityCluster$institution)) == 3) {
        ComparisonOptionsInstitution <- list(names(table(PurityCluster$institution))[1:3])
        
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
        
      } 
      
      
      if(length(unique(PurityCluster$typebiopsy)) == 2) {
        ComparisonOptionsTypebiopsy <- list(names(table(PurityCluster$typebiopsy))[1:2])
        
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
        
      } 
            
      
      if(length(unique(PurityCluster$grade[- which(is.na(PurityCluster$grade))])) == 2) {
        ComparisonOptionsGrade <- list(names(table(PurityCluster$grade))[1:2])
        
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
    }
      
      
      
    # RNAseq plots
    SNFRNAseqPlots33 <- function(TumorPurity,
                               ClusterLabels,
                               recatGRADE,
                               ClinicalFile,
                               BetaMatrix,
                               Data2RNAseqFilteredNorm,
                               ProduceImages = "No",
                               PNGorPDF = "PNG",
                               ImageName) {

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
                                  meanCount = rowMeans((as.matrix(Data2RNAseqFilteredNorm))), 
                                  medianCount = rowMedians((as.matrix(Data2RNAseqFilteredNorm))),
                                  sdCount = apply((as.matrix(Data2RNAseqFilteredNorm)),1, sd, na.rm = TRUE))
      
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
      
      par(mfrow = c(1, 1))
      pathNow <- getwd()
      meanCountPlot <- ggpubr::ggboxplot(RNAseqCluster, x = "cluster", y = "meanCount", fill = "cluster",
                                         add = "boxplot", ylab = "Mean SNF Normalized Count Value", 
                                         font.label = list(size = 20, color = "black"), 
                                         palette = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                         ggtitle("Cluster vs. Mean SNF Normalized Count Value") +
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
      if(length(unique(RNAseqCluster$stage[- which(is.na(RNAseqCluster$stage))])) == 2) {
        ComparisonOptionsStage <- list(names(table(RNAseqCluster$stage))[1:2])
        
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
      } else {
        ComparisonOptionsStage <- list(names(table(RNAseqCluster$stage))[1:2])
        
        p0ClusterStage <- ggpubr::ggboxplot(RNAseqCluster, 
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
      }
      

      
      
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
       p0ClusterTypebiopsy <- ggpubr::ggboxplot(RNAseqCluster, 
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


    
    # tSNE plot
    SNFtSNEPlots33 <- function(BetaMatrix,
                               ClinicalFile,
                               ClusterLabels, 
                               W1Data1MvalueFilteredNorm,
                               W2Data2RNAseqFilteredNorm,
                               W3Data3DNAseqFilteredNorm,
                               WCombinedFiltered,
                               ProduceImages = "No",
                               PNGorPDF = "PNG",
                               ImageName) {
    
    cat("\n tSNE for", Clusters, "clusters:\n")
    if(ProduceImages == ProduceImages) {  
      
      if(all(is.na(BetaMatrix)) != TRUE) {
        tSNEOutputFilteredBeta <- tSNEPlotGeneration(InputMatrix = BetaMatrix, 
                                                     PerplexityParameter = NA, 
                                                     ClinicalFile = ClinicalFile, 
                                                     ClusterLabels = ClusterLabels, 
                                                     Filter = NA, 
                                                     FigureGenerate = "Yes", 
                                                     PNGorPDF = "png", 
                                                     ImageName = "SNFClust_Beta_AllSamples")
        
        tSNEOutputFilteredMvalue <- tSNEPlotGeneration(InputMatrix = W1Data1MvalueFilteredNorm, # MvalueMatrix, 
                                                       PerplexityParameter = NA, 
                                                       ClinicalFile = ClinicalFile, 
                                                       ClusterLabels = ClusterLabels, 
                                                       Filter = NA, 
                                                       FigureGenerate = "Yes", 
                                                       PNGorPDF = "png", 
                                                       ImageName = "SNFClust_MvalueMatrix_AllSamples")
      }
      tSNEOutputFilteredRNAseq <- tSNEPlotGeneration(InputMatrix = W2Data2RNAseqFilteredNorm, # RNAseqCountMatrix, 
                                                     PerplexityParameter = NA, 
                                                     ClinicalFile = ClinicalFile, 
                                                     ClusterLabels = ClusterLabels, 
                                                     Filter = NA, 
                                                     FigureGenerate = "Yes", 
                                                     PNGorPDF = "png", 
                                                     ImageName = "SNFClust_RNAseq_AllSamples")
      
      tSNEOutputFilteredMutation <- tSNEPlotGeneration(InputMatrix = W3Data3DNAseqFilteredNorm, # TargetedSeqMatrix, 
                                                     PerplexityParameter = NA, 
                                                     ClinicalFile = ClinicalFile, 
                                                     ClusterLabels = ClusterLabels, 
                                                     Filter = NA, 
                                                     FigureGenerate = "Yes", 
                                                     PNGorPDF = "png", 
                                                     ImageName = "SNFClust_Mutation_AllSamples")
      
      
      tSNEOutputFiltered <- tSNEPlotGeneration(InputMatrix = WCombinedFiltered, 
                                               PerplexityParameter = NA, 
                                               ClinicalFile = ClinicalFile, 
                                               ClusterLabels = ClusterLabels, 
                                               Filter = NA, 
                                               FigureGenerate = "Yes", 
                                               PNGorPDF = "png", 
                                               ImageName = "SNFClust_Fused_AllSamplesNA")
      
      
     combMatrix <- rbind(MvalueMatrix, 
                         RNAseqCountMatrix,
                         TargetedSeqMatrix)
     
     
     tSNEOutputFiltered <- tSNEPlotGeneration(BetaMatrix = combMatrix, 
                                              PerplexityParameter = NA, 
                                              ClinicalFile = ClinicalFile, 
                                              ClusterLabels = ClusterLabels, 
                                              Filter = NA, 
                                              FigureGenerate = "Yes", 
                                              PNGorPDF = "png", 
                                              ImageName = "SNFClust_combDataMatrix_AllSamples")
      
    }
   return(invisible(NULL))
  }
      

# For model selection, see estimateNumberOfClustersGivenGraph.R
# For heatmaps, use 33_SNFClusteringScript.R
# save.image("33_SNFClustering_10Feb2021.RData")
# [END]
