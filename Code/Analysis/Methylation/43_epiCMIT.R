# Updated 3 Aug 2021
# 11 November 2020
# Function: Generate epiCMIT scores provided BetaMatrix with probes in rows and samples
#           along columns. The function uses program originally developed by Duran
#           FerrerM at https://github.com/Duran-FerrerM/Pan-B-cell-methylome
#           Original R code: 
#           https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/Estimate.epiCMIT.R 
# Author: Anjali Silva


# Input:
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns.
# ClusterLabels: A dataframe of size patients x 2, where rows correspond to patient ID and two columns 
#                containing ID and cluster label. Column names should read 'ID' and 'Cluster'. 
# TumorPurity: A numeric vector of tumor purity values of length ncol(BetaMatrix), nrow(ClinicalFile),
#              and nlength(ClusterLabels).
# ClinicalFile: A data frame of clinical categories, that has size patients x clinical categories. The
#               nrow(ClinicalFile) should match ncol(BetaMatrix), nrow(TumorPurity) and nlength(ClusterLabels).
# SurvivalFile: A data frame of survival data, with size patients x survival categories.
# MuataionFile: A dataframe with patients x mutations, where rows correspond to patient IDs and columns contain
#               gene name. The data frame should be a 0 and 1 matrix, where 0 means no mutation and 1 means
#               presence of mutation. Currently columns 'BCL2_BA_consensus' and 'EZH2' and 'BCL2' will be used.
# ImageName: A character vector giving a unique name for the images being produced. Otherwise, Image_date() will 
#            be used. 
# PNGorPDF: Output format of the image, options = "png" or "pdf". Default: png.
# ProduceImages: Produce images or not, options = "Yes" or "No". Default: Yes

# Output:
# FLResultsDataFrame: A data frame containing patients x categories, categories containing "epiCMIT", "Type",
#                    "Stage", "Sex", "BCL2Translocation", "Clusters", "EZH2Mut", "BCL2Mut", "KMT2DMut", "CREBBPMut",
#                    "EP300Mut", "epiCMIThyper", "epiCMIThypo", "indicator"



epiCMITScoreCalculatorFL43 <- function(BetaMatrix,
                                       ClusterLabels,
                                       ClinicalFile,
                                       MuataionFile = NA,
                                       ImageName = paste0("Image_", date()),
                                       PNGorPDF = "png",
                                       ProduceImages = "Yes") {
  
  library(EnvStats)
  library(plotly)
  library(ggpubr)
  
  
  # Checking inputs
  # Check if RNAseqCountMatrix and QCMatrix has the same samples, if not stop
  if (length(which((colnames(BetaMatrix) %in% ClinicalFile$SAMPLE_ID) == FALSE)) != 0) {
    stop("\n BetaMatrix and ClinicalFile should have the same samples.") }
  
  if (length(which((colnames(BetaMatrix) %in% ClusterLabels$ID) == FALSE)) != 0) {
    stop("\n BetaMatrix and ClusterLabels should have the same samples.") }
  
  if(all(is.na(MuataionFile)) != TRUE) {
    if (length(which((colnames(BetaMatrix) %in% rownames(MuataionFile)) == FALSE)) != 0) {
      stop("\n BetaMatrix and MuataionFile should have the same samples.") }
  }
  
  # Adjuste Survival file size
  if(class(SurvivalFile)[1] != "logical") {
    SurvivalFileMatchedSampleFiltered <- SurvivalFile[match(substr(colnames(BetaMatrix), 1, 9), 
                                                            SurvivalFile$LY_FL_ID), ]
  }

  
  # Duran FerrerM is the author of Estimate.epiCMIT
  Estimate.epiCMIT <- function(betas = NULL, epiCMIT.annot = NULL, export = FALSE) {
    
    ##checkings
    stopifnot(!is.null(betas))
    stopifnot(!is.null(epiCMIT.annot))
    
    #separate CpGs
    epiCMIT.hyper.CpGs <- epiCMIT.annot$Name[grep("hyper",epiCMIT.annot$epiCMIT.class)]
    epiCMIT.hypo.CpGs <- epiCMIT.annot$Name[grep("hypo",epiCMIT.annot$epiCMIT.class)]
    
    ##check how many of epiCMIT are in betas
    epiCMIT.hyper.CpGs.inbetas <- epiCMIT.hyper.CpGs[which(epiCMIT.hyper.CpGs%in%rownames(betas))]
    epiCMIT.hypo.CpGs.inbetas <- epiCMIT.hypo.CpGs[which(epiCMIT.hypo.CpGs%in%rownames(betas))]
    
    ## how many epiCMIT-CpGs might not be present in beta matrix
    if(length(epiCMIT.hyper.CpGs)!=length(epiCMIT.hyper.CpGs.inbetas)){
      w <- epiCMIT.hyper.CpGs[!epiCMIT.hyper.CpGs%in%epiCMIT.hyper.CpGs.inbetas]
      message(paste0(length(w)," epiCMIT-hyper CpGs are not in your matrix!"))
    }
    if(length(epiCMIT.hypo.CpGs)!=length(epiCMIT.hypo.CpGs.inbetas)){
      w <- epiCMIT.hypo.CpGs[!epiCMIT.hypo.CpGs%in%epiCMIT.hypo.CpGs.inbetas]
      message(paste0(length(w)," epiCMIT-hypo CpGs are not in your matrix!"))
    }
    
    if(any(is.na(betas[c(epiCMIT.hypo.CpGs.inbetas,epiCMIT.hyper.CpGs.inbetas),]))){
      message("Your matrix contain NA's. CpGs with NA's are not considered to calcualte the epiCMIT score.")
    }
    epiCMIT.hypo.res <- 1-colMeans(as.matrix(betas[epiCMIT.hypo.CpGs.inbetas,]),na.rm = T)
    epiCMIT.hyper.res <- colMeans(as.matrix(betas[epiCMIT.hyper.CpGs.inbetas,]),na.rm = T)
    
    epiCMIT.res <- epiCMIT.hypo.res
    Wlower <- names(epiCMIT.hypo.res)[which(epiCMIT.hypo.res<epiCMIT.hyper.res)]
    epiCMIT.res[Wlower] <- epiCMIT.hyper.res[Wlower]
    res <- data.frame(epiCMIT=epiCMIT.res,
                      epiCMIT.hyper=epiCMIT.hyper.res,
                      epiCMIT.hypo=epiCMIT.hypo.res)
    if(export){
      if (!requireNamespace("xlsx")){
        install.packages("xlsx")
      }
      library(xlsx)
      write.xlsx2(x = res,file = "epiCMIT.results.xlsx")
    }
    return(res)
  }
  load("Estimate.epiCMIT.RData") # to get annotation file from epiCMIT
  
  
  # Applying to our data 
  set.seed(1234)
  FLResults <- Estimate.epiCMIT(betas = BetaMatrix, 
                                epiCMIT.annot = epiCMIT.annot,
                                export = FALSE)
  # 10 epiCMIT-hyper CpGs are not in your matrix!
  # 74 epiCMIT-hypo CpGs are not in your matrix!

  # check missing CpGs
  checkMissingCpGs <- function() {
    # check which CpGs are missing in BetaMatrix_T1
    matchProbes <- match(epiCMIT.annot$Name, rownames(BetaMatrix_T1))
    head(rownames(BetaMatrix_T1[matchProbes, ]))
    # missing CpGs
    dim(epiCMIT.annot[which(is.na(matchProbes) == TRUE), ]) # 425   7
    table(epiCMIT.annot[which(is.na(matchProbes) == TRUE), ]$Relation_to_Island)
    epiCMITEntriesMissing <- data.frame(table(epiCMIT.annot[which(is.na(matchProbes) == TRUE), ]$Relation_to_Island))
    epiCMITEntriesMissing <- data.frame(table(epiCMIT.annot[which(is.na(matchProbes) == TRUE), ]$epiCMIT.class))
  
    plotly::plot_ly(epiCMITEntriesMissing, 
                    labels = ~Var1, 
                    values = ~Freq, 
                    type = 'pie',
                    textposition = 'outside', 
                    textinfo = 'label+percent') %>%
      layout(title = '', showlegend = FALSE, 
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  }


  # Create own dataframe
  # If mutation data is not present
  if(all(is.na(MuataionFile)) != TRUE) {
    FLResultsDataFrame <- data.frame(SampleID = rownames(FLResults),
                                     epiCMIT = as.numeric(FLResults$epiCMIT),
                                     Type = factor(ClinicalFile$TYPE),
                                     Stage = factor(ClinicalFile$STAGE),
                                     Sex = factor(ClinicalFile$SEX),
                                     BCL2Translocation = factor(MuataionFile$BCL2_BA_consensus),
                                     Clusters = factor(ClusterLabels$Cluster), 
                                     EZH2Mut = factor(MuataionFile$EZH2),
                                     BCL2Mut = factor(MuataionFile$BCL2),
                                     KMT2DMut = factor(MuataionFile$KMT2D),
                                     CREBBPMut = factor(MuataionFile$CREBBP),
                                     EP300Mut = factor(MuataionFile$EP300),
                                     epiCMIThyper = as.numeric(FLResults$epiCMIT.hyper),
                                     epiCMIThypo = as.numeric(FLResults$epiCMIT.hypo),
                                     indicator = NA)
  } else {
    FLResultsDataFrame <- data.frame(SampleID = rownames(FLResults),
                                     epiCMIT = as.numeric(FLResults$epiCMIT),
                                     Type = factor(ClinicalFile$TYPE),
                                     Stage = factor(ClinicalFile$STAGE),
                                     Sex = factor(ClinicalFile$SEX),
                                     Clusters = factor(ClusterLabels$Cluster), 
                                     epiCMIThyper = as.numeric(FLResults$epiCMIT.hyper),
                                     epiCMIThypo = as.numeric(FLResults$epiCMIT.hypo),
                                     indicator = NA)
  }

  # write a function to see hyper or hypomethylated is max
  for (num in seq_along(1:length(FLResults$epiCMIT))) {
    result <- which(FLResultsDataFrame[num, c("epiCMIThyper", "epiCMIThypo")] == 
              max(FLResultsDataFrame[num, c("epiCMIThyper", "epiCMIThypo")]))
    FLResultsDataFrame[num, "indicator"] <- c("epiCMIThyper", "epiCMIThypo")[result]
  }
  
  FLResultsDataFrame$epiCMIThyper <- factor(FLResultsDataFrame$epiCMIThyper)
  FLResultsDataFrame$epiCMIThypo <- factor(FLResultsDataFrame$epiCMIThypo)
  

  # Plotting 
  if(ProduceImages == "Yes") {
    
    pathNow <- getwd()
    
    # A function for comparion options
    comparisonOptions <- function(variableofInterest) {
      if(length(unique(variableofInterest)) == 1) {
        ComparisonOptions <- list(names(table(variableofInterest))[1])
      } else if(length(unique(variableofInterest)) == 2) {
        ComparisonOptions <- list(names(table(variableofInterest))[1:2])
      } else if(length(unique(variableofInterest)) == 3) {
        ComparisonOptions <- list(names(table(variableofInterest))[1:2], 
                                      names(table(variableofInterest))[2:3],
                                      names(table(variableofInterest))[c(1, 3)])
      } else if(length(unique(variableofInterest)) == 4) {
        ComparisonOptions <- list( names(table(variableofInterest))[1:2], 
                                       names(table(variableofInterest))[c(1, 3)], 
                                       names(table(variableofInterest))[c(1, 4)], 
                                       names(table(variableofInterest))[c(2, 3)], 
                                       names(table(variableofInterest))[c(2, 4)], 
                                       names(table(variableofInterest))[c(3, 4)])
      } else if(length(unique(variableofInterest)) == 5) {
        ComparisonOptions <- list(names(table(variableofInterest))[1:2], 
                                   names(table(variableofInterest))[c(1, 3)], 
                                   names(table(variableofInterest))[c(1, 4)], 
                                   names(table(variableofInterest))[c(1, 5)], 
                                   names(table(variableofInterest))[c(2, 3)], 
                                   names(table(variableofInterest))[c(2, 4)], 
                                   names(table(variableofInterest))[c(2, 5)], 
                                   names(table(variableofInterest))[c(3, 4)],
                                   names(table(variableofInterest))[c(3, 5)],
                                   names(table(variableofInterest))[c(4, 5)])
      }
      return(ComparisonOptions)
    }
    
    # TYPE
    variableofInterest <- FLResultsDataFrame$Type
    ComparisonOptions <- comparisonOptions(variableofInterest = variableofInterest)
    
    # if multiple comparisons
    if(length(ComparisonOptions) > 1) {
      coloursBarPlot <- c("#d6604d", "#66bd63", "#4575b4")
      

      par(mfrow = c(1, 1))
      TypePlot1 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                     x = "Type", 
                                     y = "epiCMIT", 
                                     width = 0.8,
                                     size = 1,
                                     color = "Type",
                                     ylab = "epiCMIT score", 
                                     font.label = list(size = 20, color = "black"), 
                                     add = "jitter",
                                     legend = "none",
                                     palette = coloursBarPlot[sort(unique(unclass(factor(variableofInterest))))]) +
                                     ggtitle("Type vs. epiCMIT") +
                                     ggpubr::stat_compare_means(comparisons = ComparisonOptions) + 
                                     # Add pairwise comparisons p-value
                                     # ggpubr::stat_compare_means(paired = FALSE) +
                                     EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_TYPE_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      
      
      par(mfrow = c(1, 1))
      coloursBarPlot <- c("#d6604d", "#2166ac")
      TypePlot2 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                     x = "Type", 
                                     y = "epiCMIT",
                                     width = 0.8,
                                     size = 1,
                                     color = "indicator",
                                     ylab = "epiCMIT score", 
                                     font.label = list(size = 20, color = "black"), 
                                     add = "jitter",
                                     legend = "none",
                                     palette = coloursBarPlot[sort(unique(unclass(factor(variableofInterest))))]) +
                                     ggtitle("Type vs. epiCMIT") +
                                     stat_compare_means(comparisons = ComparisonOptions) + 
                                     # Add pairwise comparisons p-value
                                     # stat_compare_means(paired = FALSE)  +
                                     EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_TYPE_ByScore_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
    }
    
    # if only one comparison
    if(length(ComparisonOptions) == 1) {
      coloursBarPlot <- c("#d6604d", "#66bd63", "#4575b4")
      
      par(mfrow = c(1, 1))
      TypePlot1 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                     x = "Type", 
                                     y = "epiCMIT", 
                                     width = 0.8,
                                     size = 1,
                                     color = "Type",
                                     ylab = "epiCMIT score", 
                                     font.label = list(size = 20, color = "black"), 
                                     add = "jitter",
                                     legend = "none",
                                     palette = coloursBarPlot[sort(unique(unclass(factor(variableofInterest))))]) +
                                     ggtitle("Type vs. epiCMIT")  +
                                     stat_compare_means(comparisons = ComparisonOptions) + 
                                     EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_TYPE_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      
    
      par(mfrow = c(1, 1))
      coloursBarPlot <- c("#d6604d", "#2166ac")
      TypePlot2 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                     x = "Type", 
                                     y = "epiCMIT", 
                                     width = 0.8,
                                     size = 1,
                                     color = "indicator",
                                     ylab = "epiCMIT score", 
                                     font.label = list(size = 20, color = "black"),
                                     add = "jitter",
                                     legend = "none",
                                     palette = coloursBarPlot) +
                                     ggtitle("Type vs. epiCMIT") +
                                     stat_compare_means(comparisons = ComparisonOptions) + 
                                     EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_TYPE_ByScore_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
    
    
    # STAGE
    variableofInterest <- FLResultsDataFrame$Stage
    ComparisonOptions <- comparisonOptions(variableofInterest = variableofInterest)
    coloursBarPlot <- c("#762a83", "#c2a5cf")
    
    par(mfrow = c(1, 1))
    StagePlot1 <- ggpubr::ggboxplot(FLResultsDataFrame[! is.na(FLResultsDataFrame$Stage), ], 
                                   x = "Stage", 
                                   y = "epiCMIT",
                                   width = 0.8,
                                   size = 1,
                                   color = "Stage",
                                   ylab = "epiCMIT score", 
                                   font.label = list(size = 20, color = "black"), 
                                   add = "jitter",
                                   legend = "none",
                                   palette = coloursBarPlot[sort(unique(unclass(factor(variableofInterest))))]) +
                                   ggtitle("Stage vs. epiCMIT") +
                                   stat_compare_means(comparisons = ComparisonOptions) + 
                                   # Add pairwise comparisons p-value
                                   # stat_compare_means(paired = FALSE) +
                                   EnvStats::stat_n_text()
    ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_STAGE_", ImageName, ".", PNGorPDF),
                    width = 3.5, height = 3.5)
    
    
    par(mfrow = c(1, 1))
    coloursBarPlot <- c("#d6604d", "#2166ac")
    StagePlot2 <- ggpubr::ggboxplot(FLResultsDataFrame[! is.na(FLResultsDataFrame$Stage), ], 
                                   x = "Stage", 
                                   y = "epiCMIT", 
                                   width = 0.8,
                                   size = 1,
                                   color = "indicator",
                                   ylab = "epiCMIT score", 
                                   font.label = list(size = 20, color = "black"), 
                                   add = "jitter",
                                   legend = "none",
                                   palette = coloursBarPlot) +
                                   ggtitle("Stage vs. epiCMIT") +
                                   stat_compare_means(comparisons = ComparisonOptions) + 
                                   # Add pairwise comparisons p-value
                                   # stat_compare_means(paired = FALSE) +
                                   EnvStats::stat_n_text()
    ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_STAGE_ByScore_", ImageName, ".", PNGorPDF),
                    width = 3.5, height = 3.5)
    
    # CLUSTER
    variableofInterest <- FLResultsDataFrame$Clusters
    ComparisonOptions <- comparisonOptions(variableofInterest = variableofInterest)
    coloursBarPlot <- c("#4363d8", "#f58231")
    

    par(mfrow = c(1, 1))
    ClustPlot1 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                    x = "Clusters", 
                                    y = "epiCMIT", 
                                    width = 0.8,
                                    size = 1,
                                    color = "Clusters",
                                    ylab = "epiCMIT score", 
                                    font.label = list(size = 20, color = "black"), 
                                    add = "jitter",
                                    legend = "none",
                                    palette = coloursBarPlot[sort(unique(unclass(factor(variableofInterest))))]) +
                                    ggtitle("Clusters vs. epiCMIT") +
                                    stat_compare_means(comparisons = ComparisonOptions) + 
                                    # Add pairwise comparisons p-value
                                    # stat_compare_means(paired = FALSE) +
                                    EnvStats::stat_n_text()
    ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_CLUSTER_", ImageName, ".", PNGorPDF),
                    width = 3.5, height = 3.5)
    

    par(mfrow = c(1, 1))
    coloursBarPlot <- c("#d6604d", "#2166ac")
    ClustPlot2 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                    x = "Clusters", 
                                    y = "epiCMIT", 
                                    width = 0.8,
                                    size = 1,
                                    color = "indicator",
                                    ylab = "epiCMIT score", 
                                    font.label = list(size = 20, color = "black"), 
                                    add = "jitter",
                                    legend = "none",
                                    palette = coloursBarPlot) +
                                    ggtitle("Clusters vs. epiCMIT") +
                                    stat_compare_means(comparisons = ComparisonOptions) + 
                                    # Add pairwise comparisons p-value
                                    # stat_compare_means(paired = FALSE) +
                                    EnvStats::stat_n_text()
    ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_CLUSTER_ByScore_", ImageName, ".", PNGorPDF),
                    width = 3.5, height = 3.5)
    
    
    # For publication
    coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                        '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080')
    
    ClustPlot1 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                      x = "Clusters", 
                      y = "epiCMIT", 
                      width = 0.8,
                      size = 1,
                      ylab = "epiCMIT score", 
                      #font.label = list(size = 30, color = "black"),
                      color = "Clusters",
                      palette =  coloursBarPlot,
                      add = "jitter",
                      legend = "none") +
                      ggtitle("epiCMIT_CLUSTER") +
                      ggpubr::stat_compare_means() +
                      EnvStats::stat_n_text()
    ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_CLUSTER_Publication_", ImageName, ".", PNGorPDF),
                    width = 3.5, height = 3.5)
    
    
    # If mutation file is available
    if(all(is.na(MuataionFile)) != TRUE) {
      
      # TRANSLOCATAION
      variableofInterest <- FLResultsDataFrame[! is.na(FLResultsDataFrame$BCL2Translocation), ]$BCL2Translocation
      ComparisonOptions <- comparisonOptions(variableofInterest = variableofInterest)
      coloursBarPlot <- c("#e0e0e0", "#878787")
      

      par(mfrow = c(1, 1))
      BCL2TranslocationPlot1 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                                  x = "BCL2Translocation", 
                                                  y = "epiCMIT", 
                                                  width = 0.8,
                                                  size = 1,
                                                  ylab = "epiCMIT score", 
                                                  #font.label = list(size = 30, color = "black"),
                                                  color = "Clusters",
                                                  palette =  coloursBarPlot,
                                                  add = "jitter",
                                                  legend = "none") +
                                                  ggtitle("BCL2Translocation vs. epiCMIT") +
                                                  ggpubr::stat_compare_means()  +
                                                  EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_BCL2Translocation_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      

      par(mfrow = c(1, 1))
      coloursBarPlot <- c("#d6604d", "#2166ac")
      BCL2TranslocationPlot2 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                                  x = "BCL2Translocation", 
                                                  y = "epiCMIT", 
                                                  width = 0.8,
                                                  size = 1,
                                                  color = "indicator",
                                                  ylab = "epiCMIT score", 
                                                  font.label = list(size = 20, color = "black"), 
                                                  palette = coloursBarPlot,
                                                  add = "jitter",
                                                  legend = "none") +
                                                  ggtitle("BCL2Translocation vs. epiCMIT") +
                                                  stat_compare_means(comparisons = ComparisonOptions) + 
                                                  # Add pairwise comparisons p-value
                                                  # stat_compare_means(paired = FALSE) +
                                                  EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_BCL2Translocation_ByScore_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      
      
      
      # EZH2Mut
      variableofInterest <- FLResultsDataFrame[! is.na(FLResultsDataFrame$EZH2Mut), ]$EZH2Mut
      ComparisonOptions <- comparisonOptions(variableofInterest = variableofInterest)
      coloursBarPlot <- c("#e0e0e0", "#878787")
      

      
      par(mfrow = c(1, 1))
      EZH2MutPlot1 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                        x = "EZH2Mut", 
                                        y = "epiCMIT", 
                                        width = 0.8,
                                        size = 1,
                                        color = "EZH2Mut",
                                        ylab = "epiCMIT score", 
                                        add = "jitter",
                                        legend = "none",
                                        font.label = list(size = 20, color = "black"), 
                                        palette = coloursBarPlot[sort(unique(unclass(factor(variableofInterest))))]) +
                                        ggtitle("EZH2Mut vs. epiCMIT") +
                                        stat_compare_means(comparisons = ComparisonOptions) + 
                                        # Add pairwise comparisons p-value
                                        # stat_compare_means(paired = FALSE) +
                                        EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_EZH2Mut_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      
      

      par(mfrow = c(1, 1))
      coloursBarPlot <- c("#d6604d", "#2166ac")
      EZH2MutPlot2 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                        x = "EZH2Mut", 
                                        y = "epiCMIT",
                                        width = 0.8,
                                        size = 1,
                                        color = "indicator",
                                        ylab = "epiCMIT score", 
                                        add = "jitter",
                                        legend = "none",
                                        font.label = list(size = 20, color = "black"), 
                                        palette = coloursBarPlot) +
                                        ggtitle("EZH2Mut vs. epiCMIT") +
                                        stat_compare_means(comparisons = ComparisonOptions) + 
                                        # Add pairwise comparisons p-value
                                        # stat_compare_means(paired = FALSE) +
                                        EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_EZH2Mut_ByScore_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      
      
      
      # BCL2Mut
      variableofInterest <- FLResultsDataFrame[! is.na(FLResultsDataFrame$BCL2Mut), ]$BCL2Mut
      ComparisonOptions <- comparisonOptions(variableofInterest = variableofInterest)
      coloursBarPlot <- c("#e0e0e0", "#878787")

      par(mfrow = c(1, 1))
      BCL2MutPlot1 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                        x = "BCL2Mut", 
                                        y = "epiCMIT", 
                                        width = 0.8,
                                        size = 1,
                                        color = "BCL2Mut",
                                        ylab = "epiCMIT score", 
                                        add = "jitter",
                                        legend = "none",
                                        font.label = list(size = 20, color = "black"), 
                                        palette = coloursBarPlot[sort(unique(unclass(factor(variableofInterest))))]) +
                                        ggtitle("BCL2Mut vs. epiCMIT") +
                                        stat_compare_means(comparisons = ComparisonOptions) + 
                                        # Add pairwise comparisons p-value
                                        # stat_compare_means(paired = FALSE) # +
                                        EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_BCL2Mut_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      
      par(mfrow = c(1, 1))
      coloursBarPlot <- c("#d6604d", "#2166ac")
      BCL2MutPlot2 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                        x = "BCL2Mut", 
                                        y = "epiCMIT", 
                                        width = 0.8,
                                        size = 1,
                                        color = "indicator",
                                        add = "jitter",
                                        legend = "none",
                                        ylab = "epiCMIT score", 
                                        font.label = list(size = 20, color = "black"), 
                                        palette = coloursBarPlot) +
                                        ggtitle("BCL2Mut vs. epiCMIT") +
                                        stat_compare_means(comparisons = ComparisonOptions) + 
                                        # Add pairwise comparisons p-value
                                        # stat_compare_means(paired = FALSE) # +
                                        EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_BCL2Mut_ByScore_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      
      
      
      
      
      # KMT2DMut
      variableofInterest <- FLResultsDataFrame[! is.na(FLResultsDataFrame$KMT2DMut), ]$KMT2DMut
      ComparisonOptions <- comparisonOptions(variableofInterest = variableofInterest)
      coloursBarPlot <- c("#e0e0e0", "#878787")
      
      par(mfrow = c(1, 1))
      KMT2DMutPlot1 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                         x = "KMT2DMut", 
                                         y = "epiCMIT", 
                                         width = 0.8,
                                         size = 1,
                                         color = "KMT2DMut",
                                         ylab = "epiCMIT score", 
                                         add = "jitter",
                                         legend = "none",
                                         font.label = list(size = 20, color = "black"), 
                                         palette = coloursBarPlot[sort(unique(unclass(factor(variableofInterest))))]) +
                                         ggtitle("KMT2DMut vs. epiCMIT") +
                                         stat_compare_means(comparisons = ComparisonOptions) + 
                                         # Add pairwise comparisons p-value
                                         # stat_compare_means(paired = FALSE) # +
                                         #EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_KMT2DMut_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      

      par(mfrow = c(1, 1))
      coloursBarPlot <- c("#d6604d", "#2166ac")
      KMT2DMutPlot2 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                         x = "KMT2DMut", 
                                         y = "epiCMIT", 
                                         width = 0.8,
                                         size = 1,
                                         color = "indicator",
                                         ylab = "epiCMIT score", 
                                         font.label = list(size = 20, color = "black"), 
                                         add = "jitter",
                                         legend = "none",
                                         palette = coloursBarPlot) +
                                         ggtitle("KMT2DMut vs. epiCMIT") +
                                         stat_compare_means(comparisons = ComparisonOptions) + 
                                         # Add pairwise comparisons p-value
                                         # stat_compare_means(paired = FALSE) # +
                                         #EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_KMT2DMut_ByScore_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      
      
      
      # CREBBPMut
      variableofInterest <- FLResultsDataFrame[! is.na(FLResultsDataFrame$CREBBPMut), ]$CREBBPMut
      ComparisonOptions <- comparisonOptions(variableofInterest = variableofInterest)
      coloursBarPlot <- c("#e0e0e0", "#878787")
      
      par(mfrow = c(1, 1))
      CREBBPMutPlot1 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                          x = "CREBBPMut", 
                                          y = "epiCMIT", 
                                          width = 0.8,
                                          size = 1,
                                          color = "CREBBPMut",
                                          ylab = "epiCMIT score", 
                                          font.label = list(size = 20, color = "black"), 
                                          palette = coloursBarPlot[sort(unique(unclass(factor(variableofInterest))))]) +
                                          ggtitle("CREBBPMut vs. epiCMIT") +
                                          stat_compare_means(comparisons = ComparisonOptions) + 
                                          # Add pairwise comparisons p-value
                                          # stat_compare_means(paired = FALSE) +
                                          EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_CREBBPMut_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      

      par(mfrow = c(1, 1))
      coloursBarPlot <- c("#d6604d", "#2166ac")
      CREBBPMutPlot2 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                          x = "CREBBPMut", 
                                          y = "epiCMIT", 
                                          width = 0.8,
                                          size = 1,
                                          color = "indicator",
                                          ylab = "epiCMIT score", 
                                          font.label = list(size = 20, color = "black"), 
                                          add = "jitter",
                                          legend = "none",
                                          palette = coloursBarPlot) +
                                          ggtitle("CREBBPMut vs. epiCMIT") +
                                          stat_compare_means(comparisons = ComparisonOptions) + 
                                          # Add pairwise comparisons p-value
                                          # stat_compare_means(paired = FALSE) +
                                          EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_CREBBPMut_ByScore_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      
      
      
      
      # EP300Mut
      variableofInterest <- FLResultsDataFrame[! is.na(FLResultsDataFrame$EP300Mut), ]$EP300Mut
      ComparisonOptions <- comparisonOptions(variableofInterest = variableofInterest)
      coloursBarPlot <- c("#e0e0e0", "#878787")
      
      par(mfrow = c(1, 1))
      EP300MutPlot1 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                         x = "EP300Mut", 
                                         y = "epiCMIT", 
                                         width = 0.8,
                                         size = 1,
                                         color = "EP300Mut",
                                         ylab = "epiCMIT score", 
                                         add = "jitter",
                                         legend = "none",
                                         font.label = list(size = 20, color = "black"), 
                                         palette = coloursBarPlot[sort(unique(unclass(factor(variableofInterest))))]) +
                                         ggtitle("EP300Mut vs. epiCMIT") +
                                         stat_compare_means(comparisons = ComparisonOptions) + 
                                         # Add pairwise comparisons p-value
                                         # stat_compare_means(paired = FALSE) +
                                         EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_EP300Mut_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
      

      
      par(mfrow = c(1, 1))
      coloursBarPlot <- c("#d6604d", "#2166ac")
      EP300MutPlot2 <- ggpubr::ggboxplot(FLResultsDataFrame, 
                                         x = "EP300Mut", 
                                         y = "epiCMIT", 
                                         width = 0.8,
                                         size = 1,
                                         fill = "indicator",
                                         ylab = "epiCMIT score", 
                                         add = "jitter",
                                         legend = "none",
                                         font.label = list(size = 20, color = "black"), 
                                         palette = coloursBarPlot) +
                                         ggtitle("EP300Mut vs. epiCMIT") +
                                         stat_compare_means(comparisons = ComparisonOptions) + 
                                         # Add pairwise comparisons p-value
                                         # stat_compare_means(paired = FALSE) +
                                         EnvStats::stat_n_text()
      ggplot2::ggsave(paste0(pathNow, "/img/43_epiCMIT_EP300Mut_ByScore_", ImageName, ".", PNGorPDF),
                      width = 3.5, height = 3.5)
    } # End mutation data plots
  } 
  }# End plots
  
  RESULTS <- list(ResultsDataFrame = FLResultsDataFrame)
  class(RESULTS) <- "epiCMITScoreCalculatorFL"
  return(RESULTS)
} # End function


epiCMITExample <- function() {
  # Link https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/Estimate.epiCMIT.R
  # epiCMIT calculation example
  
  # Version v.1.0, last update 07/08/2020
  
  # download a file necessary data load it in R, and delete the file
  # download.file("https://github.com/Duran-FerrerM/Pan-B-cell-methylome/raw/master/Estimate.epiCMIT.RData", 
  #               destfile = "Estimate.epiCMIT.RData", method="libcurl")
  load("Estimate.epiCMIT.RData")
  
  
  # code take form https://stackoverflow.com/questions/43427958/how-to-check-name-of-file-rdata-loaded-in-current-session-r
  # help see objects in .RData
  list_of_vars <- sapply(list.files(pattern = "Estimate.epiCMIT.RData"),
                         function(f) {
                           e <- new.env(parent = emptyenv())
                           load(f, envir = e)
                           ls(envir = e)
                         }, simplify = FALSE)
  list_of_vars
  # "betas.example"    "epiCMIT.annot"    "Estimate.epiCMIT
  
  # file.remove("Estimate.epiCMIT.RData")
  
  # Description of the function Estimate.epiCMIT()
  # This function calculates the epiCMIT, epiCMIT-hyper and epiCMIT-hypo in B-cell tumors as published in Duran-Ferrer, M, 2020.
  
  #INPUT
  # Function arguments
  #betas: DNA methylation matrix of samples to calculate the epiCMIT. CpGs on rows, samples on columns. 
  #       Rownames of matrix should be names of CpGs. See betas.example.
  #epiCMIT.annot: Annotation necessary for calculation of epiCMIT. epiCMIT have been built with 450K aray data.
  #export: Whether to export or not results in the current directory. Default to FALSE
  
  
  # OUTPUT
  # A data.frame with epiCMIT,epiCMIT-hyper and epiCMIT-hypo in all the samples. 
  # If export=T, export an .xlsx file in yout current directory with epiCMIT results.
  
  # Example execution:
  dim(betas.example) # 1348  117
  colnames(betas.example)
  exampleResults <- Estimate.epiCMIT(betas = betas.example, 
                                     epiCMIT.annot = epiCMIT.annot,
                                     export = T)
  names(exampleResults) #  "epiCMIT"       "epiCMIT.hyper" "epiCMIT.hypo" 
  
  exampleResults$epiCMIT
  exampleResults$epiCMIT.hyper
  
  
  chromosomeEntries <- data.frame(table(epiCMIT.annot$chr))
  ReIslandEntries <- data.frame(table(epiCMIT.annot$Relation_to_Island))
  epiCMITEntries <- data.frame(table(epiCMIT.annot$epiCMIT.class))
  plotly::plot_ly(epiCMITEntries, 
                  labels = ~Var1, 
                  values = ~Freq, 
                  type = 'pie',
                  textposition = 'outside', 
                  textinfo = 'label+percent') %>%
    layout(title = '', showlegend = FALSE, 
           xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
}

# [END]
