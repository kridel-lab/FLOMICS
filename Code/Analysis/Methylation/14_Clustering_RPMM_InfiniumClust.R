# Updated 23 April 2020
# Function: Clustering of Methylation Data
# Author: Anjali Silva

# Input:
# TypeofClustering: A character vector indicating what type of clustering to be performed. Options include
#                   "all", "RPMM", "InfiniumClust", "Kmeans", "Medoids", "Hierarchical". Default "all". 
# TumorPurity: A numeric vector of length equalling number of samples, indicating the purity level of each sample. Only need to be provided 
#        if TypeofClustering = "InfiniumClust", "RPMM" or "all". Default is NA. 
# BetaMatrix: A matrix of beta values for probes (or genes) x patients, with probes (or genes) in rows and patients as columns
# MvalueMatrix: A matrix of M values for probes (or genes) x patients, with probes (or genes) in rows and patients as columns;
#               or NA (default) if the MvalueMatrix is not provided
# ListofProbes: List of probes (or genes) to serve as clustering input, that are also present in rownames of BetaMatrix and MvalueMatrix
# ClinicalFile: File with patient sample names, and categories. 
# FigureGenerate: Produce images or not, options = "Yes" or "No"; default "Yes" 
# PNGorPDF: Output format of the image, options = "png" or "pdf"; default "png"

# *** If using probes, all input files should be for probes. If using genes, then all input files should be for genes. 

# Output:
# BetaMatrixProbesUsedClustering: Input probe list with corresponding Beta values used for clustering,
# Kmeans_PValue_3: Kmeans clustering results from cutting tree at k=3 groups
# Kmeans_PValue_2: Kmeans clustering results from cutting tree at k=2 groups
# Medoids_PValue_3: Medoids clustering results from cutting tree at k=3 groups
# Medoids_PValue_2: Medoids clustering results from cutting tree at k=2 groups
# Hierarchical_PValue_3: Hierarchical clustering results from cutting tree at k=3 groups
# Hierarchical_PValue_2: Hierarchical clustering results from cutting tree at k=2 groups


# Visuals saved to img folder; must have this folder in working directory in advance of running this function

Clustering <- function(TypeofClustering = "all", 
                       BetaMatrix, 
                       MvalueMatrix = NA, 
                       AnnotationFile,
                       ListofProbes, 
                       ClinicalFile, 
                       TumorPurity,
                       FigureGenerate = "Yes", 
                       PNGorPDF = "png", 
                       ImageName) {
  
  # Loading needed packages
  # LoadCheckPkg(RegularPckgs=c("gplots","ggplot2","RPMM","RColorBrewer","stats","cluster","mclust"))
  library(gplots)
  library(ggplot2)
  library(RPMM)
  library(RColorBrewer)
  library(stats)
  library(cluster)
  
  library(mclust)
  
  set.seed(1) # altering set.seed() did not alter the results
  cat("\n Clustering...")
  # Using differentially expressed results provided as input, select a subset of probes/genes for analysis 
  
  # From DifferentialExpressionResults select top probes (arranged by p-value) with p values less than 0.01
  matchedIDs <- match(ListofProbes, rownames(BetaMatrix))
  BetaMatrix_OrderListofProbes <- as.matrix(BetaMatrix[matchedIDs, ])
  
  # Obtaining path to save images
  pathNow <- getwd()
  
  qual_col_pals <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual', ]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # Visualizations
  #  if(FigureGenerate == "Yes") {
    
  #   TypeLymphomVector <- c(which(substr(colnames(BetaMatrix), 4, 5 ) == "FL"), 
  #                          which(substr(colnames(BetaMatrix), 4, 5 ) == "DL"),
  #                          which(substr(colnames(BetaMatrix), 4, 5 ) == "RL"))
  #   ColVector <- c(rep(1, length(which(substr(colnames(BetaMatrix), 4, 5 ) == "FL"))), 
  #                  rep(2, length(which(substr(colnames(BetaMatrix), 4, 5 ) == "DL"))), 
  #                  rep(3, length(which(substr(colnames(BetaMatrix), 4, 5 ) == "RL"))))
    
  #   qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
  #   col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
  #   if (PNGorPDF == "png") {
  #     png(paste0(pathNow, "/img/14_No_Clustering_", ImageName, ".", PNGorPDF))
  #   }
  #   if (PNGorPDF == "pdf") {
  #     grDevices::pdf(paste0(pathNow, "/img/14_No_Clustering_", ImageName, ".", PNGorPDF))
  #   }
  #   cat("\n Dimension of BetaMatrix_OrderListofProbes[, TypeLymphomVector]: ", 
  #       dim(BetaMatrix_OrderListofProbes[ ,TypeLymphomVector]))
    
    # cat("\n Side colour vector length is: ", length(col_vector[c(ClinicalFile$TYPE[TypeLymphomVector])+3]))
  #   gplots::heatmap.2(BetaMatrix_OrderListofProbes[ ,TypeLymphomVector], 
  #             dendrogram = 'none', Rowv = FALSE, Colv = FALSE, trace = 'none', 
  #             col = rev(redgreen(75)), key = T, density.info="density", 
  #             ColSideColors = col_vector[as.numeric(factor(ClinicalFile$TYPE[TypeLymphomVector])) + 3], 
  #             main = "No Clustering")
    
  #   grDevices::dev.off()
    
  #   if (PNGorPDF == "png") {
  #     png(paste0(pathNow, "/img/14_Clustering_Probes_And_Samples_", ImageName, ".", PNGorPDF))
  #   }
  #   if (PNGorPDF == "pdf") { 
  #     grDevices::pdf(paste0(pathNow, "/img/14_Clustering_Probes_And_Samples_", ImageName, ".", PNGorPDF))
  #   }
    # heatmap(BetaMatrix_OrderListofProbes[1:50,], scale="n", col=colorRampPalette(c("yellow","black","blue"),space="Lab")(128))
  #   gplots::heatmap.2(BetaMatrix_OrderListofProbes, dendrogram = 'both', 
  #             Rowv = T, Colv = T,  trace = 'none', col = rev(redgreen(75)), 
  #             key = T, density.info = "density", main = "Clustering Probes And Samples", 
  #             ColSideColors = col_vector[as.numeric(factor(ClinicalFile$TYPE)) + 3])
  #   grDevices::dev.off()
    
  #   if (PNGorPDF == "png") {    
  #     png(paste0(pathNow, "/img/14_Clustering_Probes_", ImageName, ".", PNGorPDF))
  #   }
  #   if (PNGorPDF == "pdf") { 
  #     grDevices::pdf(paste0(pathNow, "/img/14_Clustering_Probes_", ImageName, ".", PNGorPDF))
  #   }
  #   gplots::heatmap.2(BetaMatrix_OrderListofProbes[ ,TypeLymphomVector], 
  #             dendrogram ="row", trace = "none", scale = "none", Rowv = T,  
  #             Colv = F, col = rev(redgreen(75)), key = T, density.info = "density", 
  #             ColSideColors = col_vector[as.numeric(factor(ClinicalFile$TYPE[TypeLymphomVector])) + 3],
  #             main = "Clustering Probes Only")
  #   grDevices::dev.off()    
    
  #   if (PNGorPDF == "png") {  
  #     png(paste0(pathNow,"/img/14_Clustering_Samples_",ImageName,".",PNGorPDF))
  #   }
  #   if (PNGorPDF == "pdf") { 
  #     grDevices::pdf(paste0(pathNow,"/img/14_Clustering_Samples_",ImageName,".",PNGorPDF))
  #    }
  #   gplots::heatmap.2(BetaMatrix_OrderListofProbes, dendrogram = "column", trace = "none", 
  #             scale = "none", Rowv = F,  col = rev(redgreen(75)), key = T, 
  #             density.info = "density", 
  #             ColSideColors = col_vector[as.numeric(factor(ClinicalFile$TYPE)) + 3],
  #             main = "Clustering Samples Only")
  #   grDevices::dev.off() 
  # }
  
  # Initialization of object names
  RPMM <- Hierarchical_PValue_4 <- Hierarchical_PValue_3 <- Hierarchical_PValue_2 <- NA
  Medoids_PValue_4 <- Medoids_PValue_3 <- Medoids_PValue_2 <- NA
  Kmeans_PValue_4 <- Kmeans_PValue_3 <- Kmeans_PValue_2 <- NA
  InfiniumClust_4 <- InfiniumClust_3 <- InfiniumClust_2 <- NA
  
  if (TypeofClustering == "all" || TypeofClustering == "RPMM") {
    
    cat("\n Start RPMM Clustering...")
    # Performs beta latent class modeling using recursively-partitioned mixture model
    rpmm <- blcTree(t(BetaMatrix_OrderListofProbes), verbose = 0)
    # Get weight matrix and show first few rows
    rpmmWeightMatrix <- blcTreeLeafMatrix(rpmm)
    
    # Visualizations
    if(FigureGenerate == "Yes") {
      if (PNGorPDF == "png") {
        png(paste0(pathNow, "/img/14_Clustering_Samples_RPMM_", ImageName, ".", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/14_Clustering_Samples_RPMM_", ImageName ,".", PNGorPDF))
      }
      par(mfrow = c(2, 1)) # setting space for two plots
      plot(rpmm) ; title("RPMM Profile")
      #plotTree.blcTree(rpmm) ; title("Dendrogram with Labels")
      plotTree.blcTree(rpmm, labelFunction = function(u, digits) table(as.character(ClinicalFile$TYPE[u$index])))
      # plotImage.blcTree(rpmm); same as plot(rpmm)
      title("Dendrogram with Tissue Counts")
      grDevices::dev.off()
    }
    
    # Get class assignments and compare with tissue
    rpmmClass <- blcTreeLeafClasses(rpmm) # classes are assigned based on terminal nodes
    rpmmClass_numbers <- c(rpmmClass) 
    table(rpmmClass_numbers)
    # write.csv(cbind(colnames(BetaMatrix_OrderListofProbes), rpmmClass_numbers), file = "265_point26SD_Probes_1ClusterModel_RPMM_Allsamples_Labels.csv")
    
    cat("\n Results based on blcTree for type \n")
    print(table(rpmmClass_numbers, ClinicalFile$TYPE))
    print(chisq.test(table(rpmmClass_numbers, ClinicalFile$TYPE))) 
    # X-squared = 
  
    # printing labels against clinical classes
    # table(rpmmClass_numbers,ClinicalFile$TYPE)
    cat("\n Results based on blcTree for stage \n")
    print(table(rpmmClass_numbers, ClinicalFile$STAGE))
    # Chi-squared Test of Independence
    # Null: cluster membership is independent of stage at .05 significance level.
    print(chisq.test(table(rpmmClass_numbers, ClinicalFile$STAGE))) 
    
    cat("\n Results based on blcTree for sex \n")
    print(table(rpmmClass_numbers, ClinicalFile$SEX))
    print(chisq.test(table(rpmmClass_numbers, ClinicalFile$SEX))) 
    
    cat("\n Results based on blcTree for site biopsy \n")
    print(table(rpmmClass_numbers, trimws(ClinicalFile$SITE_BIOPSY)))
    print(chisq.test(table(rpmmClass_numbers, trimws(ClinicalFile$SITE_BIOPSY)))) 
    
    cat("\n Results based on blcTree for institution \n")
    print(table(rpmmClass_numbers, ClinicalFile$INSTITUTION))
    print(chisq.test(table(rpmmClass_numbers, ClinicalFile$INSTITUTION))) 
    
    cat("\n Results based on blcTree for 14_18 translocation \n")
    NA_entries <- which(is.na(ClinicalFile$TRANSLOC_14_18) == TRUE)
    print(table(rpmmClass_numbers[- NA_entries], ClinicalFile$TRANSLOC_14_18[- NA_entries]) )
    print(chisq.test((table(rpmmClass_numbers[-NA_entries], ClinicalFile$TRANSLOC_14_18[- NA_entries]) ))) 
    cat("\n Done RPMM Clustering...")
  
    
    # Plot rpmm results
    if(FigureGenerate == "Yes") {
      cat("\n Start RPMM Visuals")
      if (PNGorPDF == "png") {
        png(paste0(pathNow, "/img/14_Clustering_Samples_RPMM_", ImageName, ".", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/14_Clustering_Samples_RPMM_", ImageName ,".", PNGorPDF))
      }
      par(mfrow = c(2, 1)) # setting space for two plots
      plot(rpmm) ; title("RPMM Profile")
      #plotTree.blcTree(rpmm) ; title("Dendrogram with Labels")
      plotTree.blcTree(rpmm, labelFunction = function(u, digits) table(as.character(ClinicalFile$TYPE[u$index])))
      # plotImage.blcTree(rpmm); same as plot(rpmm)
      title("Dendrogram with Tissue Counts")
      grDevices::dev.off()
      
      
      # if (PNGorPDF == "png") {
      #   png(paste0(pathNow,"/img/14_Clustering_Samples_RPMM_heatmap_", ImageName, ".",PNGorPDF))
      # }
      # if (PNGorPDF == "pdf") {
      #   grDevices::pdf(paste0(pathNow,"/img/14_Clustering_Samples_RPMM_heatmap_", ImageName, ".",PNGorPDF))
      # }
      # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
      # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      
      # if(length(unique(rpmmClass_numbers)) == 2) {
      #   SampleArrange_RPMM <- c(which(rpmmClass_numbers == 1), which(rpmmClass_numbers == 2))
      # } else if(length(unique(rpmmClass_numbers)) == 3) {
      #   SampleArrange_RPMM <- c(which(rpmmClass_numbers == 1), which(rpmmClass_numbers == 2), 
      #                           which(rpmmClass_numbers == 3))
      # } else if(length(unique(rpmmClass_numbers)) == 4) {
      #   SampleArrange_RPMM <- c(which(rpmmClass_numbers == 1), which(rpmmClass_numbers == 2), 
      #                           which(rpmmClass_numbers == 3), which(rpmmClass_numbers == 4))
      # } else if(length(unique(rpmmClass_numbers)) == 5) {
      #   SampleArrange_RPMM <- c(which(rpmmClass_numbers == 1), which(rpmmClass_numbers == 2), 
      #                           which(rpmmClass_numbers == 3), which(rpmmClass_numbers == 4), 
      #                           which(rpmmClass_numbers == 5))
      # }
      
      # gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange_RPMM],dendrogram='none', Rowv=F, Colv=F, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[rpmmClass_numbers[SampleArrange_RPMM]+3])
      # edited to cluster both row and columns
      # gplots::heatmap.2(BetaMatrix_OrderListofProbes[, SampleArrange_RPMM], dendrogram = 'row', 
      #           Rowv = TRUE, Colv = FALSE, trace = 'none', col = rev(redgreen(75)), key = T, 
      #           density.info = "density", ColSideColors = col_vector[rpmmClass_numbers[SampleArrange_RPMM] + 3])
      # graphics::legend("topright",      
      #        legend = paste0("C ",unique(rpmmClass_numbers[SampleArrange_RPMM])),
      #        col = unique(col_vector[rpmmClass_numbers[SampleArrange_RPMM] + 3]), 
      #        lty = 1,             
      #        lwd = 5,           
      #        cex =.7)
      # grDevices::dev.off()
      # cat("\n Done RPMM Visuals")
    }
    
    RPMM <- list(RPMMoutput = rpmmClass,RPMMoutputLabels = rpmmClass_numbers)

    # Alternate initialization
    # rpmm2 <- blcTree(t(BetaMatrix_OrderListofProbes[1:100,]), verbose=0, initFunctions=list(blcInitializeSplitEigen(), blcInitializeSplitFanny(nu=2.5)))
    # plotTree.blcTree(rpmm2, labelFunction=function(u,digits) table(as.character(ClinicalFile$TYPE[u$index])))
    # title("Dendrogram with Tissue Counts")
    
    # Alternate split criterion
      # rpmm3 <- blcTree(t(BetaMatrix_OrderListofProbes), verbose=0, maxlev=3, splitCriterion=blcSplitCriterionLevelWtdBIC)
    # Get class assignments and compare with tissue
      # rpmmClass3 <- blcTreeLeafClasses(rpmm3) # classes are assigned based on terminal nodes
      # rpmmClass_numbers3 <- c(rpmmClass3) 
    
    
    # Plotting rpmm3 results only if different clustering results from rpmmClass_numbers
      # if(FigureGenerate=="Yes" && sum(diag(table(rpmmClass_numbers, rpmmClass_numbers3))) != ncol(BetaMatrix_OrderListofProbes)){
      # if (PNGorPDF=="png"){
      #   png(paste0(pathNow,"/img/14_Clustering_Samples_RPMM_AlternateSplitCriterion_",ImageName,".",PNGorPDF))
      # }
      # if (PNGorPDF=="pdf"){
      #   pdf(paste0(pathNow,"/img/14_Clustering_Samples_RPMM_AlternateSplitCriterion_",ImageName,".",PNGorPDF))
      # }
      # par(mfrow=c(2,1))
      # plot(rpmm3) ; title("RPMM Profile with alternate split criterion")
      # plotTree.blcTree(rpmm3, labelFunction=function(u,digits) table(as.character(ClinicalFile$TYPE[u$index])))
      # title("Dendrogram with Tissue Counts")
      # grDevices::dev.off()    
      
      # if (PNGorPDF=="png"){
      #   png(paste0(pathNow,"/img/14_Clustering_Samples_RPMM_AlternateSplitCriterion_heatmap_",ImageName,".",PNGorPDF))
      # }
      # if (PNGorPDF=="pdf"){
      #   pdf(paste0(pathNow,"/img/14_Clustering_Samples_RPMM_AlternateSplitCriterion_heatmap_",ImageName,".",PNGorPDF))
      # }
      
      # if(length(unique(rpmmClass_numbers3))==2){
      #   SampleArrange_RPMM3 <- c(which(rpmmClass_numbers3==1), which(rpmmClass_numbers3==2))
      # }else if(length(unique(rpmmClass_numbers3))==3){
      #   SampleArrange_RPMM3 <- c(which(rpmmClass_numbers3==1), which(rpmmClass_numbers3==2), which(rpmmClass_numbers3==3))
      # }else if(length(unique(rpmmClass_numbers3))==4){
      #   SampleArrange_RPMM3 <- c(which(rpmmClass_numbers3==1), which(rpmmClass_numbers3==2), which(rpmmClass_numbers3==3), which(rpmmClass_numbers3==4))
      # }else if(length(unique(rpmmClass_numbers3))==5){
      #   SampleArrange_RPMM3 <- c(which(rpmmClass_numbers3==1), which(rpmmClass_numbers3==2), which(rpmmClass_numbers3==3), which(rpmmClass_numbers3==4), which(rpmmClass_numbers3==5))
      # }
      
      # gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange_RPMM3],dendrogram='row', Rowv=TRUE, Colv=FALSE, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[rpmmClass_numbers3[SampleArrange_RPMM3]+3])
      # legend("topright",      
             #        legend = paste0("C ",unique(rpmmClass_numbers3[SampleArrange_RPMM3])),
             #        col = unique(col_vector[rpmmClass_numbers3[SampleArrange_RPMM3]+3]), 
             #        lty= 1,             
             #        lwd = 5,           
             #        cex=.7
             #  )
             # grDevices::dev.off() 
      # }
    # }
    
    # rpmm4 <- blcTree(t(BetaMatrix_OrderListofProbes[1:100,]), verbose=0, maxlev=3, splitCriterion=blcSplitCriterionJustRecordEverything)
    # plotTree.blcTree(rpmm4, labelFunction=function(u,digits) table(as.character(ClinicalFile$TYPE[u$index])))
    # title("Dendrogram with Tissue Counts")
    
    # Performs Gaussian latent class modeling using recursively-partitioned mixture model
    # rpmm2<- glcTree(t(BetaMatrix_OrderListofProbes), verbose=0)
    # plotTree.glcTree(rpmm2,labelFunction=function(u,digits) table(as.character(ClinicalFile$TYPE[u$index])))
    # rpmm2Class <- glcTreeLeafClasses(rpmm2) # classes are assigned based on terminal nodes
    # rpmm2Class_numbers <- c(rpmm2Class) 
    
    # Plotting rpmm2 results only if different clustering results from rpmmClass_numbers
    # if(FigureGenerate=="Yes" && sum(diag(table(rpmmClass_numbers, rpmm2Class_numbers))) != ncol(BetaMatrix_OrderListofProbes)){
    
    #  cat("\n Results based on glcTree for stage \n")
    #  print(table(rpmm2Class_numbers,ClinicalFile$STAGE))
    
    #   if (PNGorPDF=="png"){
    #     png(paste0(pathNow,"/img/14_Clustering_Samples_GaussianRPMM_",ImageName,".",PNGorPDF))
    #   }
    #   if (PNGorPDF=="pdf"){
    #     pdf(paste0(pathNow,"/img/14_Clustering_Samples_GaussianRPMM_",ImageName,".",PNGorPDF))
    #   }
    #   par(mfrow=c(2,1))
    #   plot(rpmm2) ; title("Gaussian RPMM Profile")
    #   plotTree.glcTree(rpmm2, labelFunction=function(u,digits) table(as.character(ClinicalFile$TYPE[u$index])))
    #   title("Dendrogram with Tissue Counts")
    #   grDevices::dev.off()    
    
    #   if(length(unique(rpmm2Class_numbers))==2){
    #     SampleArrange_RPMM2 <- c(which(rpmm2Class_numbers==1), which(rpmm2Class_numbers==2))
    #   }else if(length(unique(rpmm2Class_numbers))==3){
    #     SampleArrange_RPMM2 <- c(which(rpmm2Class_numbers==1), which(rpmm2Class_numbers==2), which(rpmm2Class_numbers==3))
    #   }else if(length(unique(rpmm2Class_numbers))==4){
    #     SampleArrange_RPMM2 <- c(which(rpmm2Class_numbers==1), which(rpmm2Class_numbers==2), which(rpmm2Class_numbers==3), which(rpmm2Class_numbers==4))
    #   }else if(length(unique(rpmm2Class_numbers))==5){
    #     SampleArrange_RPMM2 <- c(which(rpmm2Class_numbers==1), which(rpmm2Class_numbers==2), which(rpmm2Class_numbers==3), which(rpmm2Class_numbers==4), which(rpmm2Class_numbers==5))
    #   }
    
    #   gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange_RPMM3],dendrogram='row', Rowv=TRUE, Colv=FALSE, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[rpmm2Class_numbers[SampleArrange_RPMM2]+3])
    #   legend("topright",      
    #          legend = paste0("C ",unique(rpmm2Class_numbers[SampleArrange_RPMM2])),
    #          col = unique(col_vector[rpmm2Class_numbers[SampleArrange_RPMM2]+3]), 
    #          lty= 1,             
    #          lwd = 5,           
    #          cex=.7
    #   )
    #   grDevices::dev.off() 
    # }
    
    
    # Plotting only set for RPMM 
    if(FigureGenerate == "Yes") {
      # Plot clusters by purity
      par(mfrow = c(1, 1))
      PurityCluster <- data.frame(purity = TumorPurity, cluster = rpmmClass_numbers)
      
      # Define comparisons
      if(length(unique(rpmmClass_numbers)) == 2) {
        ClusterComparisonOptions <- list(c("1","2"))
      } else if(length(unique(rpmmClass_numbers)) == 3) {
        ClusterComparisonOptions <- list(c("1", "2"),
                                         c("1", "3"), 
                                         c("2", "3"))
      } else if(length(unique(rpmmClass_numbers)) == 4) {
        ClusterComparisonOptions <- list(c("1", "2"), 
                                         c("1", "3"), 
                                         c("1", "4"), 
                                         c("2", "3"), 
                                         c("2", "4"), 
                                         c("3", "4"))
      } else if(length(unique(rpmmClass_numbers)) == 5) {
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
                                      palette = coloursBarPlot[sort(unique(rpmmClass_numbers))]) +
        ggtitle("Cluster vs. Tumor purity") +
        stat_compare_means(comparisons = ClusterComparisonOptions) + 
        # Add pairwise comparisons p-value
        stat_compare_means(paired = FALSE)
      ggsave(paste0(pathNow, "/img/33_PurityPlot_", ImageName, ".png"))
    }
  }
  
  if (TypeofClustering == "all" || TypeofClustering == "Kmeans") {
    
    cat("\n Start Kmeans Clustering")
    # Distance-based clustering using only the 100 first probes/genes
    # Kmeans
    kmeans_clustering_PValue_4 <- stats::kmeans(x = t(BetaMatrix_OrderListofProbes), 
                                         centers = 4, nstart = 10 )
    kmeans_clustering_PValue_3 <- stats::kmeans(x = t(BetaMatrix_OrderListofProbes), 
                                         centers = 3, nstart = 10 )
    kmeans_clustering_PValue_2 <- stats::kmeans(x = t(BetaMatrix_OrderListofProbes), 
                                         centers = 2, nstart = 10 )
    cat("\n Done Kmeans Clustering")
    # Comparing obtained labels of kmeans clustering vs type and stage
    # ARI_kmeans <- c(adjustedRandIndex(kmeans_clustering_PValue_3$cluster, c(ClinicalFile$TYPE)),
    #                 adjustedRandIndex(kmeans_clustering_PValue_2$cluster, c(ClinicalFile$STAGE)))
    
    if(FigureGenerate == "Yes") {
      
      cat("\n Start Kmeans Visuals")
      # Kmeans_PValue_2Clusters
      if (PNGorPDF == "png") {
        png(paste0(pathNow,"/img/14_Clustering_Samples_Kmeans_2Clusters_", ImageName, ".", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow,"/img/14_Clustering_Samples_Kmeans_2Clusters_", ImageName, ".", PNGorPDF))
      }
      SampleArrange <- c(which(as.vector(kmeans_clustering_PValue_2$cluster) == 1), 
                         which(as.vector(kmeans_clustering_PValue_2$cluster) == 2))
      
      # previous 
      # gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange],dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[as.vector(kmeans_clustering_PValue_2$cluster[SampleArrange])+3])
      # editing to add dendrograms for both row and columns on 15 May 2019
      gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange], dendrogram='row', 
                Rowv=T, Colv=F, trace='none', col = rev(redgreen(75)), key = T,
                density.info="density", 
                ColSideColors = col_vector[as.vector(kmeans_clustering_PValue_2$cluster[SampleArrange]) + 3])
      graphics::legend("topright",      
             legend = paste0("C ",unique(kmeans_clustering_PValue_2$cluster[SampleArrange])),
             col = unique(col_vector[as.vector(kmeans_clustering_PValue_2$cluster[SampleArrange]) + 3]), 
             lty = 1,             
             lwd = 5,           
             cex =.7)
      grDevices::dev.off()
      
      # Kmeans_PValue_3Clusters
      if (PNGorPDF == "png") {
        png(paste0(pathNow, "/img/14_Clustering_Samples_Kmeans_3Clusters_", ImageName,".", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/14_Clustering_Samples_Kmeans_3Clusters_", ImageName,".", PNGorPDF))
      }
      SampleArrange <- c(which(as.vector(kmeans_clustering_PValue_3$cluster) == 1), 
                         which(as.vector(kmeans_clustering_PValue_3$cluster) == 2), 
                         which(as.vector(kmeans_clustering_PValue_3$cluster) == 3))
      # gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange],dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[as.vector(kmeans_clustering_PValue_3$cluster[SampleArrange])+3])
      # editing to add dendrograms for both row and columns on 15 May 2019
      gplots::heatmap.2(BetaMatrix_OrderListofProbes[, SampleArrange], dendrogram = 'row', Rowv = T,
                Colv = F, trace='none', col = rev(redgreen(75)), key = T, density.info = "density", 
                ColSideColors = col_vector[as.vector(kmeans_clustering_PValue_3$cluster[SampleArrange]) + 3])
      graphics::legend("topright",      
             legend = paste0("C ",unique(kmeans_clustering_PValue_3$cluster[SampleArrange])),
             col = unique(col_vector[as.vector(kmeans_clustering_PValue_3$cluster[SampleArrange]) + 3]), 
             lty = 1,             
             lwd = 5,           
             cex =.7)
      grDevices::dev.off()
      
      # Kmeans_PValue_4Clusters
      if (PNGorPDF == "png") {
        png(paste0(pathNow, "/img/14_Clustering_Samples_Kmeans_4Clusters_", ImageName, ".", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/14_Clustering_Samples_Kmeans_4Clusters_", ImageName, ".", PNGorPDF))
      }
      SampleArrange <- c(which(as.vector(kmeans_clustering_PValue_4$cluster) == 1), 
                         which(as.vector(kmeans_clustering_PValue_4$cluster) == 2), 
                         which(as.vector(kmeans_clustering_PValue_4$cluster) == 3),
                         which(as.vector(kmeans_clustering_PValue_4$cluster) == 4))
      # gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange],dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[as.vector(kmeans_clustering_PValue_3$cluster[SampleArrange])+3])
      # editing to add dendrograms for both row and columns on 15 May 2019
      gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange], dendrogram = 'row', 
                Rowv = T, Colv = F, trace = 'none', col = rev(redgreen(75)), key = T, 
                density.info = "density", 
                ColSideColors = col_vector[as.vector(kmeans_clustering_PValue_4$cluster[SampleArrange]) + 3])
      graphics::legend("topright",      
             legend = paste0("C ",unique(kmeans_clustering_PValue_4$cluster[SampleArrange])),
             col = unique(col_vector[as.vector(kmeans_clustering_PValue_4$cluster[SampleArrange]) + 3]), 
             lty = 1,             
             lwd = 5,           
             cex = .7)
      grDevices::dev.off()
      cat("\n Done Kmeans Visuals")
    }
    Kmeans_PValue_4 = as.vector(kmeans_clustering_PValue_4$cluster)                  
    Kmeans_PValue_3 = as.vector(kmeans_clustering_PValue_3$cluster)
    Kmeans_PValue_2 = as.vector(kmeans_clustering_PValue_2$cluster)
  }

  if (TypeofClustering == "all" || TypeofClustering == "Medoids") { 
    
    cat("\n Start Medoids Clustering")
    # Medoids 
    medoids_clustering_PValue_4 <- cluster::pam(x = t(BetaMatrix_OrderListofProbes), k = 4)
    medoids_clustering_PValue_3 <- cluster::pam(x = t(BetaMatrix_OrderListofProbes), k = 3)
    medoids_clustering_PValue_2 <- cluster::pam(x = t(BetaMatrix_OrderListofProbes), k = 2)
    
    # Comparing obtained labels of Medoids clustering vs type and stage
    # ARI_medoids  <- c(adjustedRandIndex(medoids_clustering_PValue_3$cluster, c(ClinicalFile$TYPE)),
    #                 adjustedRandIndex(medoids_clustering_PValue_2$cluster, c(ClinicalFile$STAGE)))
    
    if(FigureGenerate == "Yes") {
      
      cat("\n Start Medoids Visuals")
      # Medoids_PValue_2Clusters
      if (PNGorPDF == "png") {
        png(paste0(pathNow,"/img/14_Clustering_Samples_Medoids_2Clusterss_",ImageName,".",PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow,"/img/14_Clustering_Samples_Medoids_2Clusters_",ImageName,".",PNGorPDF))
      }
      SampleArrange <- c(which(as.vector(medoids_clustering_PValue_2$cluster) == 1), 
                         which(as.vector(medoids_clustering_PValue_2$cluster) == 2))
      #gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange],dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[as.vector(medoids_clustering_PValue_2$cluster[SampleArrange])+3])
      # editing to add dendrograms for both row and columns on 15 May 2019
      gplots::heatmap.2(BetaMatrix_OrderListofProbes[, SampleArrange], dendrogram='row', Rowv=T, 
                Colv=F, trace='none', col = rev(redgreen(75)), key = T, density.info="density", 
                ColSideColors = col_vector[as.vector(medoids_clustering_PValue_2$cluster[SampleArrange]) + 3])
      graphics::legend("topright",      
             legend = paste0("C ",unique(medoids_clustering_PValue_2$cluster[SampleArrange])),
             col = unique(col_vector[as.vector(medoids_clustering_PValue_2$cluster[SampleArrange]) + 3]), 
             lty = 1,             
             lwd = 5,           
             cex =.7)
      grDevices::dev.off()
      
      # Medoids_PValue_3Clusters
      if (PNGorPDF == "png") {
        png(paste0(pathNow,"/img/14_Clustering_Samples_Medoids_3Clusters_",ImageName,".",PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow,"/img/14_Clustering_Samples_Medoids_3Clusters_",ImageName,".",PNGorPDF))
      }
      SampleArrange <- c(which(as.vector(medoids_clustering_PValue_3$cluster) == 1), 
                         which(as.vector(medoids_clustering_PValue_3$cluster) == 2), 
                         which(as.vector(medoids_clustering_PValue_3$cluster) == 3))
      #gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange],dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[as.vector(medoids_clustering_PValue_3$cluster[SampleArrange])+3])
      # editing to add dendrograms for both row and columns on 15 May 2019
      gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange], dendrogram='row', Rowv=T, 
                Colv=F, trace='none', col = rev(redgreen(75)), key = T, density.info="density", 
                ColSideColors = col_vector[as.vector(medoids_clustering_PValue_3$cluster[SampleArrange]) + 3])
      graphics::legend("topright",      
             legend = paste0("C ",unique(medoids_clustering_PValue_3$cluster[SampleArrange])),
             col = unique(col_vector[as.vector(medoids_clustering_PValue_3$cluster[SampleArrange]) + 3]), 
             lty = 1,             
             lwd = 5,           
             cex =.7)
      grDevices::dev.off()
      
      # Medoids_PValue_3Clusters
      if (PNGorPDF == "png") {
        png(paste0(pathNow,"/img/14_Clustering_Samples_Medoids_4Clusters_",ImageName,".",PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow,"/img/14_Clustering_Samples_Medoids_4Clusters_",ImageName,".",PNGorPDF))
      }
      SampleArrange <- c(which(as.vector(medoids_clustering_PValue_4$cluster) == 1), 
                         which(as.vector(medoids_clustering_PValue_4$cluster) == 2), 
                         which(as.vector(medoids_clustering_PValue_4$cluster) == 3), 
                         which(as.vector(medoids_clustering_PValue_4$cluster) == 4))
      #gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange],dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[as.vector(medoids_clustering_PValue_3$cluster[SampleArrange])+3])
      # editing to add dendrograms for both row and columns on 15 May 2019
      gplots::heatmap.2(BetaMatrix_OrderListofProbes[, SampleArrange], dendrogram='row', 
                Rowv=T, Colv=F, trace='none', col = rev(redgreen(75)), key = T, 
                density.info="density", 
                ColSideColors = col_vector[as.vector(medoids_clustering_PValue_4$cluster[SampleArrange]) + 3])
      graphics::legend("topright",      
             legend = paste0("C ",unique(medoids_clustering_PValue_4$cluster[SampleArrange])),
             col = unique(col_vector[as.vector(medoids_clustering_PValue_4$cluster[SampleArrange])+3]), 
             lty = 1,             
             lwd = 5,           
             cex =.7)
      grDevices::dev.off()
      cat("\n Done Medoids Visuals")
    }
    Medoids_PValue_4 = as.vector(medoids_clustering_PValue_4$cluster)
    Medoids_PValue_3 = as.vector(medoids_clustering_PValue_3$cluster)
    Medoids_PValue_2 = as.vector(medoids_clustering_PValue_2$cluster)
  }
  
  if (TypeofClustering == "all" || TypeofClustering == "Hierarchical") {
  
    cat("\n Start Hierarchical Clustering")
    # Hierarchical 
    hierarchical_clustering_PValue <- hclust(dist(t(BetaMatrix_OrderListofProbes)))
    hierarchical_cuttree_PValue_2 <- cutree(hierarchical_clustering_PValue, 2)
    hierarchical_cuttree_PValue_3 <- cutree(hierarchical_clustering_PValue, 3)
    hierarchical_cuttree_PValue_4 <- cutree(hierarchical_clustering_PValue, 4)
    
    cat("\n Done Hierarchical Clustering")
    
    if(FigureGenerate=="Yes") {
      
      cat("\n Start Hierarchical Visuals")
      # Hierarchical_PValue_2Clusters
      if (PNGorPDF == "png") {
        png(paste0(pathNow,"/img/14_Clustering_Samples_Hierarchical_2Clusters_",ImageName,".",PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow,"/img/14_Clustering_Samples_Hierarchical_2Clusters_",ImageName,".",PNGorPDF))
      }
      SampleArrange <- c(which(as.vector(hierarchical_cuttree_PValue_2) == 1), 
                         which(as.vector(hierarchical_cuttree_PValue_2) == 2))
      # gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange],dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[as.vector(hierarchical_cuttree_PValue_2[SampleArrange])+3])
      # editing to add dendrograms for both row and columns on 15 May 2019
      gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange],dendrogram='both', 
                Rowv=T, Colv=T, trace='none', col = rev(redgreen(75)), key = T, 
                density.info="density", 
                ColSideColors = col_vector[as.vector(hierarchical_cuttree_PValue_2[SampleArrange])+3])
      graphics::legend("topright",      
             legend = paste0("C ",unique(hierarchical_cuttree_PValue_2[SampleArrange])),
             col = unique(col_vector[as.vector(hierarchical_cuttree_PValue_2[SampleArrange])+3]), 
             lty = 1,             
             lwd = 5,           
             cex =.7)
      grDevices::dev.off()
      
      # Hierarchical_PValue_3Clusters
      if (PNGorPDF == "png") {
        png(paste0(pathNow,"/img/14_Clustering_Samples_Hierarchical_3Clusters_",ImageName,".",PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow,"/img/14_Clustering_Samples_Hierarchical_3Clusters_",ImageName,".",PNGorPDF))
      }
      SampleArrange <- c(which(as.vector(hierarchical_cuttree_PValue_3) == 1),
                         which(as.vector(hierarchical_cuttree_PValue_3) == 2), 
                         which(as.vector(hierarchical_cuttree_PValue_3) == 3))
      # gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange],dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[as.vector(hierarchical_cuttree_PValue_3[SampleArrange])+3])
      gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange], dendrogram = 'both', Rowv = T, 
                Colv = T, trace = 'none', col = rev(redgreen(75)), key = T, density.info = "density", 
                ColSideColors = col_vector[as.vector(hierarchical_cuttree_PValue_3[SampleArrange]) + 3])
      graphics::legend("topright",      
             legend = paste0("C ", unique(hierarchical_cuttree_PValue_3[SampleArrange])),
             col = unique(col_vector[as.vector(hierarchical_cuttree_PValue_3[SampleArrange]) + 3]), 
             lty = 1,             
             lwd = 5,           
             cex =.7)
      grDevices::dev.off()
      
      # Hierarchical_PValue_4Clusters
      if (PNGorPDF == "png") {
        png(paste0(pathNow,"/img/14_Clustering_Samples_Hierarchical_4Clusters_",ImageName,".",PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow,"/img/14_Clustering_Samples_Hierarchical_4Clusters_",ImageName,".",PNGorPDF))
      }
      SampleArrange <- c(which(as.vector(hierarchical_cuttree_PValue_4) == 1), 
                         which(as.vector(hierarchical_cuttree_PValue_4) == 2), 
                         which(as.vector(hierarchical_cuttree_PValue_4) == 3), 
                         which(as.vector(hierarchical_cuttree_PValue_4) == 4))
      # gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange],dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[as.vector(hierarchical_cuttree_PValue_3[SampleArrange])+3])
      gplots::heatmap.2(BetaMatrix_OrderListofProbes[,SampleArrange], dendrogram = 'both', Rowv=T, Colv=T, trace='none', col = rev(redgreen(75)), key = T, density.info="density", ColSideColors = col_vector[as.vector(hierarchical_cuttree_PValue_4[SampleArrange]) + 3])
      graphics::legend("topright",      
             legend = paste0("C ",unique(hierarchical_cuttree_PValue_4[SampleArrange])),
             col = unique(col_vector[as.vector(hierarchical_cuttree_PValue_4[SampleArrange]) + 3]), 
             lty = 1,             
             lwd = 5,           
             cex =.7)
      grDevices::dev.off()
      cat("\n Done Hierarchical Visuals")
    }
    Hierarchical_PValue_4 = as.vector(hierarchical_cuttree_PValue_4)
    Hierarchical_PValue_3 = as.vector(hierarchical_cuttree_PValue_3)
    Hierarchical_PValue_2 = as.vector(hierarchical_cuttree_PValue_2)
  }
  
  if (TypeofClustering == "all" || TypeofClustering == "InfiniumClust") {
    cat("\n Start InfiniumClust Clustering")    
    library(InfiniumPurify)
    data(iDMC) 
    probes <- iDMC[["DLBC"]]
    probes.true <- names(probes[probes == T])
    beta.sel <- BetaMatrix_T1[row.names(BetaMatrix) %in% probes.true,]
    purity_OICR_withNormal <- InfiniumPurify::getPurity(tumor.data = beta.sel, tumor.type = "DLBC")

    
    InfiniumClust_Clustering <- list() # create list to save data 
    
    
    InfiniumClust_Clustering[[4]] <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_OrderListofProbes, 
                                                                   purity = purity_OICR_withNormal, 
                                                                   K = 4, maxiter = 100, tol = 0.0001)
    InfiniumClust_Clustering[[4]] <- list(InfiniumClust_Clustering[[4]],
                                          mclust::map(InfiniumClust_Clustering[[4]]$Z))
    
    InfiniumClust_Clustering[[3]] <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_OrderListofProbes,
                                                                   purity = purity_OICR_withNormal, 
                                                                   K = 3, maxiter = 100, tol = 0.0001)
    InfiniumClust_Clustering[[3]] <- list(InfiniumClust_Clustering[[3]], 
                                          mclust::map(InfiniumClust_Clustering[[3]]$Z))
    
    set.seed(1234)
    InfiniumClust_Clustering[[2]] <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_OrderListofProbes, 
                                                                   purity = purity_OICR_withNormal, 
                                                                   K = 2, maxiter = 100, tol = 0.0001)
    
    InfiniumClust_Clustering[[2]] <-  list(InfiniumClust_Clustering[[2]], 
                                           mclust::map(InfiniumClust_Clustering[[2]]$Z))
    # saveRDS(InfiniumClust_Clustering, file = "InfiniumClustering2to4.rds")
    # Testing by varing purity levels
    # InfiniumClust_Clustering_PurityAltered <- list()
    # TumorPurity2 <- TumorPurity
    # # Scenario1
    # TumorPurity2[c(1:10)] <- 0.9
    # TumorPurity2[c(11:165)] <- 0.8
    # TumorPurity2[c(166:170)] <- 0.1
    # InfiniumClust_Clustering_PurityAltered[[2]] <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_OrderListofProbes, 
    #                                                                purity = TumorPurity2, 
    #                                                                K = 2, maxiter = 100, tol = 0.0001)
    # table(mclust::map(InfiniumClust_Clustering_PurityAltered[[2]]$Z))
    # table(mclust::map(InfiniumClust_Clustering_PurityAltered[[2]]$Z), ClinicalFile$TYPE)
    
    # # scenario2
    # TumorPurity3 <- TumorPurity
    # TumorPurity3[c(1:10)] <- 0.85
    # TumorPurity3[c(11:165)] <- 0.85
    # TumorPurity3[c(166:170)] <- 0.1
    # InfiniumClust_Clustering_PurityAltered[[3]] <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_OrderListofProbes, 
    #                                                                             purity = TumorPurity3, 
    #                                                                              K = 2, maxiter = 100, tol = 0.0001)
    # table(mclust::map(InfiniumClust_Clustering_PurityAltered[[3]]$Z))
    # table(mclust::map(InfiniumClust_Clustering_PurityAltered[[3]]$Z), ClinicalFile$TYPE)
    
    
    # # scenario3
    # TumorPurity4 <- TumorPurity
    # TumorPurity4[c(1:10)] <- 0.46
    # TumorPurity4[c(11:165)] <- 0.46
    # TumorPurity4[c(166:170)] <- 0.46
    # InfiniumClust_Clustering_PurityAltered[[4]] <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_OrderListofProbes, 
    #                                                                             purity = TumorPurity4, 
    #                                                                              K = 2, maxiter = 100, tol = 0.0001)
    # table(mclust::map(InfiniumClust_Clustering_PurityAltered[[4]]$Z))
    # table(mclust::map(InfiniumClust_Clustering_PurityAltered[[4]]$Z), ClinicalFile$TYPE)
    
    # # scenario 4
    # InfiniumClust_Clustering_PurityAltered[[5]] <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_OrderListofProbes[, c(1:165)], 
    #                                                                              purity = TumorPurity[c(1:165)], 
    #                                                                              K = 2, maxiter = 100, tol = 0.0001)
    # table(mclust::map(InfiniumClust_Clustering_PurityAltered[[5]]$Z))
    # table(mclust::map(InfiniumClust_Clustering_PurityAltered[[5]]$Z), ClinicalFile$TYPE[c(1:165)])
    
    
    
    
    # plotting sd vs mean/confidence
    if (FigureGenerate == "Yes") {
      i = 2
      ClusterLabels = InfiniumClust_Clustering[[i]][2][[1]]
      
      # results compared with clinical categories 
      cat("\n Results tabulation \n")
      print(table(ClusterLabels))
      
      
      cat("\n Results based on blcTree for TYPE \n")
      print(table(ClusterLabels, ClinicalFile$TYPE))
      print(chisq.test(table(ClusterLabels, ClinicalFile$TYPE))) 
      # X-squared = 
      
      
      
      # printing labels against clinical classes
      # table(rpmmClass_numbers,ClinicalFile$TYPE)
      cat("\n Results based on blcTree for stage \n")
      print(table(ClusterLabels, ClinicalFile$STAGE))
      # Chi-squared Test of Independence
      # Null: cluster membership is independent of stage at .05 significance level.
      print(chisq.test(table(ClusterLabels, ClinicalFile$STAGE))) 
      
      cat("\n Results based on blcTree for sex \n")
      print(table(ClusterLabels, ClinicalFile$SEX))
      print(chisq.test(table(ClusterLabels, ClinicalFile$SEX))) 
      
      cat("\n Results based on blcTree for site biopsy \n")
      print(table(ClusterLabels, trimws(ClinicalFile$SITE_BIOPSY)))
      print(chisq.test(table(ClusterLabels, trimws(ClinicalFile$SITE_BIOPSY)))) 
      
      cat("\n Results based on blcTree for type biopsy \n")
      print(table(ClusterLabels, trimws(ClinicalFile$TYPE_BIOPSY)))
      print(chisq.test(table(ClusterLabels, trimws(ClinicalFile$TYPE_BIOPSY)))) 
      
      
      cat("\n Results based on blcTree for grade \n")
      matchedIDs <- match(substr(colnames(BetaMatrix_OrderListofProbes[, which(substr(colnames(BetaMatrix_OrderListofProbes), 4, 5) == "FL")]), 1, 9), 
                          SurvivalFile$LY_FL_ID)
      if (sum(is.na(matchedIDs) ) > 0) {
        # if NAs are present
          SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile[matchedIDs[- which(is.na(matchedIDs) == TRUE)], ]
        } else {
          SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile[matchedIDs, ]
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
      recatGRADE <- SurvivalFile_OrderedbyBetaMatrixPatients$GRADE
      recatGRADE[which(SurvivalFile_OrderedbyBetaMatrixPatients$GRADE == "1")] <- "1-2"
      recatGRADE[which(SurvivalFile_OrderedbyBetaMatrixPatients$GRADE == "2")] <- "1-2"
      recatGRADE[which(SurvivalFile_OrderedbyBetaMatrixPatients$GRADE == "IN_SITU")] <- "1-2"
      recatGRADE[which(SurvivalFile_OrderedbyBetaMatrixPatients$GRADE == "3")] <- "3A"

      
      print(table(ClusterLabels[which(substr(colnames(BetaMatrix_OrderListofProbes), 4, 5) == "FL")], 
                  trimws(recatGRADE)))
      print(chisq.test(table(ClusterLabels[which(substr(colnames(BetaMatrix_OrderListofProbes), 4, 5) == "FL")], trimws(recatGRADE)))) 
      
      
      cat("\n Results based on blcTree for institution \n")
      print(table(ClusterLabels, ClinicalFile$INSTITUTION))
      print(chisq.test(table(ClusterLabels, ClinicalFile$INSTITUTION))) 
      
      cat("\n Results based on blcTree for 14_18 translocation \n")
      NA_entries <- which(is.na(ClinicalFile$TRANSLOC_14_18) == TRUE)
      print(table(ClusterLabels[- NA_entries], ClinicalFile$TRANSLOC_14_18[- NA_entries]) )
      print(chisq.test((table(ClusterLabels[-NA_entries], ClinicalFile$TRANSLOC_14_18[- NA_entries]) ))) 
      cat("\n Done InfiniumClust Clustering...")
      
      
      
      par(mfrow = c(1, 1))
      PurityCluster <- data.frame(purity = TumorPurity, 
                                  cluster = ClusterLabels,
                                  Clusters = as.character(ClusterLabels),
                                  stage = ClinicalFile$STAGE,# [c(11:165), ]
                                  translocation = ClinicalFile$TRANSLOC_14_18, # [c(11:165), ]
                                  type = ClinicalFile$TYPE, # [c(11:165), ]
                                  meanBeta = rowMeans(t(BetaMatrix_OrderListofProbes)),
                                  medianBeta = rowMedians(t(BetaMatrix_OrderListofProbes)),
                                  sdBeta = apply(t(BetaMatrix_OrderListofProbes),1, sd, na.rm = TRUE))
      
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
      
      
      PurityPlot <- ggpubr::ggboxplot(PurityCluster, x = "cluster", y = "purity", fill = "cluster",
                                      add = "boxplot", ylab = "Tumor Purity", 
                                      font.label = list(size = 20, color = "black"), 
                                      palette = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                      ggtitle("Cluster vs. Tumor purity") +
                                      stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                      # Add pairwise comparisons p-value
                                      stat_compare_means(paired = FALSE)
      ggsave(paste0(pathNow, "/img/14_PurityPlot_", ImageName, ".png"))
      

      
      meanBetaPlot <- ggpubr::ggboxplot(PurityCluster, x = "cluster", y = "meanBeta", fill = "cluster",
                                      add = "boxplot", ylab = "Mean Beta Value", 
                                      font.label = list(size = 20, color = "black"), 
                                      palette = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                      ggtitle("Cluster vs. meanBeta values") +
                                      stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                      # Add pairwise comparisons p-value
                                      stat_compare_means(paired = FALSE)
      ggsave(paste0(pathNow, "/img/14_meanBeta_", ImageName, ".png"))
        
      
      medianBetaPlot <- ggpubr::ggboxplot(PurityCluster, x = "cluster", y = "medianBeta", fill = "cluster",
                                        add = "boxplot", ylab = "Median Beta Value", 
                                        font.label = list(size = 20, color = "black"), 
                                        palette = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                        ggtitle("Cluster vs. medianBeta values") +
                                        stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                        # Add pairwise comparisons p-value
                                        stat_compare_means(paired = FALSE)
      ggsave(paste0(pathNow, "/img/14_medianBeta_", ImageName, ".png"))      

      
      
      
      # plotting sd vs mean plot - 10 Feb 2020
      # based on function 31
      # based on all probes in beta matrix
      MeanSDPlots <-  MeanSDPlot(BetaMatrix = BetaMatrix,
                                 ClinicalFile = ClinicalFile, 
                                 ClusterLabels = ClusterLabels, 
                                 FigureGenerate = "Yes", PNGorPDF ="png", 
                                 ImageName = "2ClusterModel_AllProbes")
      
      # based on selected probes in beta matrix
      MeanSDPlotsFiltered <-  MeanSDPlot(BetaMatrix = BetaMatrix_OrderListofProbes,
                                 ClinicalFile = ClinicalFile, 
                                 ClusterLabels = ClusterLabels, 
                                 FigureGenerate = "Yes", PNGorPDF ="png", 
                                 ImageName = "2ClusterModel_FilteredProbes")
  
      # plotting confidence in clustering 
      ClusterConfidencePlots  <- ClusterConfidencePlot(ProbabilityMatrix = InfiniumClust_Clustering[[i]][[1]]$Z, 
                                                       FigureGenerate = "Yes", PNGorPDF ="png", 
                                                       ImageName = "2ClusterModel")

    
      
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
      
      p0_cluster_type <- ggpubr::ggboxplot(PurityCluster, 
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
      ggsave(paste0(pathNow, "/img/14_medianBeta_TypeCluster", ImageName, ".png")) 
      
      # Plotting median beta value for patients in each cluster, seperated by stage
      # Setting the number of comparisons for ggpubr
      if(length(unique(PurityCluster$stage[-which(is.na(PurityCluster$stage))])) == 2) {
        ComparisonOptionsStage <- list(names(table(PurityCluster$stage))[1:2])
      } 
      
      p0_cluster_stage <- ggpubr::ggboxplot(PurityCluster[-which(is.na(PurityCluster$stage)),], 
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
      ggsave(paste0(pathNow, "/img/14_medianBeta_StageCluster", ImageName, ".png")) 
      
      
      # Plotting median beta value for patients in each cluster, seperated by translocation
      # Setting the number of comparisons for ggpubr
      if(length(unique(PurityCluster$stage[-which(is.na(PurityCluster$translocation))])) == 2) {
        ComparisonOptionsTranslocation <- list(names(table(PurityCluster$translocation))[1:2])
      } 
      
      p0_cluster_translocation <- ggpubr::ggboxplot(PurityCluster[-which(is.na(PurityCluster$translocation)),], 
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
      ggsave(paste0(pathNow, "/img/14_medianBeta_TranslocationCluster", ImageName, ".png")) 
      
      
      ########################## Plot heatmap
      ClusterLabels = InfiniumClust_Clustering[[i]][[2]]
      
      matchRowNames <- match(rownames(BetaMatrix_OrderListofProbes), AnnotationFile$V1)
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
        
      rownames(annotation_row) = rownames(BetaMatrix_OrderListofProbes[orderRowsGeneRegions, ])
      
      
      
      # order patients by cluster 
      orderColumns <- c(which(ClusterLabels == 1),
                 which(ClusterLabels == 2),
                 which(ClusterLabels == 3))
      
      annotation_col = data.frame(
                          Cluster = factor(ClusterLabels[orderColumns]),
                          Disease = factor(ClinicalFile$TYPE[orderColumns]),
                          Stage = factor(ClinicalFile$STAGE[orderColumns]),
                          Sex = factor(ClinicalFile$SEX[orderColumns]),
                          Translocation = factor(ClinicalFile$TRANSLOC_14_18[orderColumns]),
                          TypeBiopsy =  ClinicalFile$TYPE_BIOPSY[orderColumns],
                          SiteBiopsy =  ClinicalFile$SITE_BIOPSY[orderColumns])
      rownames(annotation_col) = colnames(BetaMatrix_OrderListofProbes[, orderColumns])
      
      # use http://colorbrewer2.org/#type=diverging&scheme=RdGy&n=8 to decide colors 
      
      ann_colors = list(
        RelationToIsland = c("Island" = "#a6cee3", "N_Shelf" = "#1f78b4",
                             "N_Shore" = "#fb9a99", "OpenSea" = "#fdbf6f", 
                             "S_Shelf" = "#ff7f00", "S_Shore" = "#cab2d6"),
        GeneRegion = c("1stExon" = "#8dd3c7", "3'UTR" = "#ffffb3", 
                       "5'UTR" = "#bebada", "Body" = "#fb8072",
                       "ExonBnd" = "#80b1d3", "NA" = "#fdb462", 
                       "TSS1500" = "#fccde5", "TSS200" = "#d9d9d9"),
        Cluster = c("1" = "#4363d8", "2" = "#f58231"), 
        Disease = c(DLBCL = "#d6604d", FL = "#66bd63", RLN = "#4575b4"), #option 5
        Stage = c(ADVANCED = "#762a83", LIMITED = "#c2a5cf"),
        Sex = c("F"="#b35806", "M"="#fdb863"),
        Translocation = c("0"="#e0e0e0","1"="#878787"),
        TypeBiopsy = c("TISSUE" = "#a6dba0", "CORE"="#878787"),
        SiteBiopsy = c("LN" = "#f1b6da" , "EN" = "#c51b7d"))
      
      
      heatmap_all <- pheatmap::pheatmap(as.matrix(BetaMatrix_OrderListofProbes[orderRowsGeneRegions, orderColumns]), 
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
  RNAseqHeatmap <- function() {
      # RNAseq data based on methylation probes
    
      # normalized RNAseq counts: t(Data2RNAseqFilterededgeRNorm)
      # raq counts: 
      #       rowNamesNormalized <- match(rownames(t(Data2RNAseqFilterededgeRNorm)), 
      # rownames(RNAseqCountsAllSamplesFilteredProbes))
      # RNAseqCountsAllSamplesFilteredProbes[rowNamesNormalized, ]
      
      # methylation gene names
      matchRowNames <- match(rownames(BetaMatrix_OrderListofProbes), AnnotationFile$V1)
      MethylationGenes <- geneRegion <- sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Name[matchRowNames])
      length(MethylationGenes) # 5000
      length(unique(MethylationGenes)) # 25097
      
      # RNAseq gene names
      ENSMBLid <- read.csv(file = paste0(RNAseqDirPath, "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv"))
      matchRNAseq <- match(substr(rownames(RNASeqCountMatrixMatched), 1, 15), ENSMBLid$gene)
      
      # match methylation genes with RNAseq genes
      matchRNAseqMethylation <- match(unique(MethylationGenes), ENSMBLid$name[matchRNAseq])
      RNAseqGenesFromMethylation <- ENSMBLid$name[matchRNAseq][matchRNAseqMethylation[! is.na(matchRNAseqMethylation)]]
      length(RNAseqGenesFromMethylation) # 1838
      
      # get RNAseq counts for methylation genes 
      match5000 <- match(RNAseqGenesFromMethylation, ENSMBLid$name)
      matchRNAseqMethylationGenes <- match(ENSMBLid$gene[match5000], 
                                           substr(rownames(RNASeqCountMatrixMatched), 1, 15))
      # matchRNAseqMethylationGenes[!is.na(matchRNAseqMethylationGenes)] # remove NA values
      dim(RNASeqCountMatrixMatched[matchRNAseqMethylationGenes[! is.na(matchRNAseqMethylationGenes)], ]) #  19620   132
      RNASeqCountMatrixMatchedMethylation <- RNASeqCountMatrixMatched[matchRNAseqMethylationGenes[! is.na(matchRNAseqMethylationGenes)], ]
      dim(RNASeqCountMatrixMatchedMethylation) # 1838   132
      # RNASeqCountMatrixOriginalFeatureFiltered2 <- RNASeqCountMatrixOriginalFeatureFiltered[matchRNAseqMethylationGenes[!is.na(matchRNAseqMethylationGenes)], ]
      # saveRDS(RNASeqCountMatrixMatchedMethylation, file = "RNASeqCountMatrixMatchedMethylationProbes.rds")
      # RNASeqCountMatrixMatchedMethylationNorm <- SNFtool::standardNormalization(log2(1 + RNASeqCountMatrixMatchedMethylation))
      # [,orderColumnsRNAseq]
      
      # order patients by cluster 
      
      #MatchRNAseqwith170 <- match(colnames(RNASeqCountMatrixMatchedMethylation), colnames(BetaMatrix_T1))
      # looking at Tier 2
      MatchRNAseqwith170 <- match(RNAseqQC18June2020T2Samples104$QCMatrixMatchedSampleFiltered$SAMPLE_ID, colnames(RNASeqCountMatrixMatchedMethylation))
      
      orderColumnsRNAseq <- c(which(ClusterLabels[MatchRNAseqwith170] == 1),
                              which(ClusterLabels[MatchRNAseqwith170] == 2))
      Cluster = factor(ClusterLabels[MatchRNAseqwith170])
      Disease = factor(ClinicalFile$TYPE[MatchRNAseqwith170])
      Stage = factor(ClinicalFile$STAGE[MatchRNAseqwith170])
      Sex = factor(ClinicalFile$SEX[MatchRNAseqwith170])
      Translocation = factor(ClinicalFile$TRANSLOC_14_18[MatchRNAseqwith170])
      TypeBiopsy =  factor(ClinicalFile$TYPE_BIOPSY[MatchRNAseqwith170])
      SiteBiopsy =  factor(ClinicalFile$SITE_BIOPSY[MatchRNAseqwith170])
      
      annotation_colRNAseq = data.frame(
        Cluster = factor(Cluster[orderColumnsRNAseq]),
        Disease = factor(Disease[orderColumnsRNAseq]),
        Stage = factor(Stage[orderColumnsRNAseq]),
        Sex = factor(Sex[orderColumnsRNAseq]),
        Translocation = factor(Translocation[orderColumnsRNAseq]),
        TypeBiopsy = factor(TypeBiopsy[orderColumnsRNAseq]),
        SiteBiopsy = factor(SiteBiopsy[orderColumnsRNAseq]))
      rownames(annotation_colRNAseq) = colnames(RNASeqCountMatrixMatchedMethylation[, orderColumnsRNAseq])
      
      ann_colorsRNAseq = list(
        Cluster = c("1" = "#4363d8", "2" = "#f58231"), 
        Disease = c(DLBCL = "#d6604d", FL = "#66bd63", RLN = "#4575b4"), #option 5
        Stage = c(ADVANCED = "#762a83", LIMITED = "#c2a5cf"),
        Sex = c("F" = "#b35806", "M" = "#fdb863"),
        Translocation = c("0" = "#e0e0e0","1" = "#878787"),
        TypeBiopsy = c("TISSUE" = "#a6dba0", "CORE" = "#878787"),
        SiteBiopsy = c("LN" = "#f1b6da" , "EN" = "#c51b7d"))
      

      
      ## Get some nicer colours
      mypalette <- brewer.pal(11,"RdYlBu")
      morecols <- colorRampPalette(mypalette)
      
      # Plot the heatmap
      pheatmap::pheatmap((log(RNASeqCountMatrixMatchedMethylation + 0.01))[, orderColumnsRNAseq],
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
      

      # generate normalized counts
      ClinicalFileT1Cluster <- data.frame(ClinicalFile, "CLUSTER" = ClusterLabels)
      ClinicalFileT1ClusterRNAseq <- ClinicalFileT1Cluster[match(colnames(RNASeqCountMatrixMatchedMethylation),
                                                                 ClinicalFileT1Cluster$SAMPLE_ID), ]
      # dim(ClinicalFileT1ClusterRNAseq) # 132  26
      exprsRNAseq <- edgeR::DGEList(counts = RNASeqCountMatrixMatchedMethylation, 
                                    group = factor(ClinicalFileT1ClusterRNAseq$CLUSTER))
      exprsRNAseq$counts <- exprsRNAseq$counts

      # perform the TMM normalization and display the normalization factors
      Data2RNAseqFilterededgeRNormFactors <- edgeR::calcNormFactors(exprsRNAseq)
      # dim(Data2RNAseqFilterededgeRNormFactors) # 19620   132
      # plotMDS(Data2RNAseqFilterededgeRNorm)
      # Generate matrix
      Data2RNAseqFilterededgeRNorm <- edgeR::cpm(Data2RNAseqFilterededgeRNormFactors, 
                                                 normalized.lib.sizes = TRUE, log = TRUE)
      # dim(Data2RNAseqFilterededgeRNorm) # 19620   132
      
      
      
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
      
      
      
      
      
      par(mfrow = c(1, 1))
      RNAseqCluster <- data.frame(purity = TumorPurity[MatchRNAseqwith170 [-28], ], 
                                  cluster = ClusterLabels[MatchRNAseqwith170][-28],
                                  Clusters = as.character(ClusterLabels[MatchRNAseqwith170][-28]),
                                  stage = ClinicalFile$STAGE[MatchRNAseqwith170][-28],# [c(11:165), ]
                                  translocation = ClinicalFile$TRANSLOC_14_18[MatchRNAseqwith170][-28], # [c(11:165), ]
                                  type = ClinicalFile$TYPE[MatchRNAseqwith170][-28], # [c(11:165), ]
                                  meanCount = rowMeans(t(as.matrix(RNASeqCountMatrixMatchedMethylation[, -28]))),
                                  medianCount = rowMedians(t(as.matrix(RNASeqCountMatrixMatchedMethylation[, -28]))),
                                  sdCount = apply(t(as.matrix(RNASeqCountMatrixMatchedMethylation[, -28])),1, sd, na.rm = TRUE))
      
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
      
      
   
    
      
      
      meanCountPlot <- ggpubr::ggboxplot(RNAseqCluster, x = "cluster", y = "meanCount", fill = "cluster",
                                        add = "boxplot", ylab = "Mean Count Value (normalized)", 
                                        font.label = list(size = 20, color = "black"), 
                                        palette = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                        ggtitle("Cluster vs. Mean Count Value") +
                                        stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                        # Add pairwise comparisons p-value
                                        stat_compare_means(paired = FALSE)
      ggsave(paste0(pathNow, "/img/14_meancount_normalized_", ImageName, ".png"))      
      
  
      medianCountPlot <- ggpubr::ggboxplot(RNAseqCluster, x = "cluster", y = "medianCount", fill = "cluster",
                                         add = "boxplot", ylab = "Median Count Value (normalized)", 
                                         font.label = list(size = 20, color = "black"), 
                                         palette = coloursBarPlot[sort(unique(ClusterLabels))]) +
                                         ggtitle("Cluster vs. Median Count Value") +
                                         stat_compare_means(comparisons = ClusterComparisonOptions) + 
                                         # Add pairwise comparisons p-value
                                         stat_compare_means(paired = FALSE)
      ggsave(paste0(pathNow, "/img/14_mediancount_", ImageName, ".png"))    
      
      
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
      
      p0_cluster_type <- ggpubr::ggboxplot(RNAseqCluster, 
                                           x = "type", 
                                           y = "medianCount", 
                                           fill = "Clusters",
                                           add = "boxplot", 
                                           xlab = "Type", 
                                           ylab = "Median Count Values (normalized)",
                                           palette = coloursBarPlot[1:4]) + 
                                           stat_compare_means(comparisons = ComparisonOptionsType) + 
                                           # Add pairwise comparisons p-value
                                           stat_compare_means(paired = FALSE)
      ggsave(paste0(pathNow, "/img/14_medianCount_TypeCluster", ImageName, ".png")) 
      
      # Plotting median beta value for patients in each cluster, seperated by stage
      # Setting the number of comparisons for ggpubr
      if(length(unique(RNAseqCluster$stage[-which(is.na(RNAseqCluster$stage))])) == 2) {
        ComparisonOptionsStage <- list(names(table(RNAseqCluster$stage))[1:2])
      } 
      
      p0_cluster_stage <- ggpubr::ggboxplot(RNAseqCluster[-which(is.na(RNAseqCluster$stage)),], 
                                            x = "stage", 
                                            y = "medianCount", 
                                            fill = "Clusters",
                                            add = "boxplot", 
                                            xlab = "Stage", 
                                            ylab = "Median Count Values (normalized)",
                                            palette = coloursBarPlot[1:4]) + 
                                            stat_compare_means(comparisons = ComparisonOptionsStage) + 
                                            # Add pairwise comparisons p-value
                                            stat_compare_means(paired = FALSE)
      ggsave(paste0(pathNow, "/img/14_medianCount_StageCluster", ImageName, ".png")) 
      
      
      # Plotting median beta value for patients in each cluster, seperated by translocation
      # Setting the number of comparisons for ggpubr
      if(length(unique(RNAseqCluster$stage[-which(is.na(RNAseqCluster$translocation))])) == 2) {
        ComparisonOptionsTranslocation <- list(names(table(RNAseqCluster$translocation))[1:2])
      } 
      
      p0_cluster_translocation <- ggpubr::ggboxplot(RNAseqCluster[-which(is.na(RNAseqCluster$translocation)),], 
                                                    x = "translocation", 
                                                    y = "meanCount", 
                                                    fill = "Clusters",
                                                    add = "boxplot", 
                                                    xlab = "Translocation Status", 
                                                    ylab = "Mean Count Values (normalized)",
                                                    palette = coloursBarPlot[1:4]) + 
                                                    stat_compare_means(comparisons = ComparisonOptionsTranslocation) + 
                                                    # Add pairwise comparisons p-value
                                                    stat_compare_means(paired = FALSE)
      ggsave(paste0(pathNow, "/img/14_medianBeta_TranslocationCluster", ImageName, ".png")) 
      
      
      
    } 

    }
  
  
  

    
  ImageName = "InfiniumClust"
  ClusterLabels = InfiniumClust_Clustering[[i]][[2]]
  if(FigureGenerate == "Yes") {
    cat("\n Methylation density plot:\n")
    # based on selected probes in beta matrix
    MethylationDensityOutputFiltered <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                         ClusterLabels = ClusterLabels, 
                                                         PlotWithinAnnotationCategories = "No", 
                                                         ClinicalCategoryToVisualize = "CLUSTER", 
                                                         BetaMatrix = BetaMatrix_OrderListofProbes, 
                                                         AnnotationFile = AnnotationFile, 
                                                         ClinicalFile = ClinicalFile, 
                                                         SampleSheet = NA, 
                                                         FigureGenerate = "Yes", 
                                                         ImageName = paste0("MethylationPlot_SelectedProbes", ImageName), 
                                                         PNGorPDF = "png")
    
    # based on all probes in beta matrix
    MethylationDensityOutput <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                       ClusterLabels = ClusterLabels, 
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
                                                         ClusterLabels = ClusterLabels, 
                                                         PlotWithinAnnotationCategories = "No", 
                                                         ClinicalCategoryToVisualize = "CLUSTER", 
                                                         BetaMatrix = BetaMatrix, 
                                                         AnnotationFile = AnnotationFile, 
                                                         ClinicalFile = ClinicalFile, 
                                                         SampleSheet = NA, 
                                                         FigureGenerate = "No")
    }
    
    
    
  if(FigureGenerate == "Yes") {
    cat("\n Proportion visualization:\n")
    # based on selected probes in beta matrix
    ProportionVisualizationOutputFiltered <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", 
                                                               ClusterLabels = ClusterLabels, 
                                                               BetaMatrix = BetaMatrix_OrderListofProbes, 
                                                               AnnotationFile = AnnotationFile, 
                                                               ClinicalFile = ClinicalFile, 
                                                               ClinicalCategory = "TYPE", 
                                                               PlotWithinCategories = "Yes", 
                                                               FigureGenerate = "Yes", 
                                                               PNGorPDF = "png", 
                                                               ImageName = paste0("AllProbes_", ImageName))
    
    # based on all probes in beta matrix
    ProportionVisualizationOutput <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", 
                                                             ClusterLabels = ClusterLabels, 
                                                             BetaMatrix = BetaMatrix, 
                                                             AnnotationFile = AnnotationFile, 
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
                                                               AnnotationFile = AnnotationFile, 
                                                               ClinicalFile = ClinicalFile, 
                                                               ClinicalCategory = "TYPE", 
                                                               PlotWithinCategories = "Yes", 
                                                               FigureGenerate = "No")
    }
    
    

    if(FigureGenerate == "Yes") {
      cat("\n Differentially methylated regions analysis:\n")
      # based on selected probes in beta matrix
      MvalueMatrix_OrderListofProbes <- MvalueMatrix[match(rownames(BetaMatrix_OrderListofProbes), rownames(MvalueMatrix)), ]
      DiffMethylatedRegionsOutputFiltered <- DiffMethylatedRegions(Method = "DMRcate", 
                                                           BetaMatrix = BetaMatrix_OrderListofProbes, 
                                                           MvalueMatrix = MvalueMatrix_OrderListofProbes, 
                                                           ContrastColumnName = "CLUSTER", 
                                                           ClinicalFile = ClinicalFile, 
                                                           ClusterLabels = ClusterLabels, 
                                                           AnnotationFile = AnnotationFile, 
                                                           ProduceImages = "Yes", 
                                                           DMR = 1, 
                                                           PNGorPDF = "png",
                                                           ExpressionFile = NA)
      
      # based on all probes in beta matrix
      DiffMethylatedRegionsOutput <- DiffMethylatedRegions(Method = "DMRcate", 
                                                           BetaMatrix = BetaMatrix, 
                                                           MvalueMatrix = MvalueMatrix, 
                                                           ContrastColumnName = "CLUSTER", 
                                                           ClinicalFile = ClinicalFile, 
                                                           ClusterLabels = ClusterLabels, 
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
                                                           ClusterLabels = ClusterLabels, 
                                                           AnnotationFile = AnnotationFile, 
                                                           ProduceImages = "No", 
                                                           DMR = 1, 
                                                           ExpressionFile = NA)
    }
    
    if(FigureGenerate == "Yes") {
      cat("\n Mean SD plot: \n")
      MeanSDPlot(BetaMatrix = BetaMatrix, 
                 ClinicalFile = ClinicalFile,
                 ClusterLabels = ClusterLabels, 
                 FigureGenerate = "Yes", 
                 PNGorPDF ="png", 
                 ImageName = paste0("MeanSDPlot_", ImageName))
    }
     
    if(FigureGenerate == "NotRun") {
      cat("\n Survival analysis: \n")    
      SurvivalAnalysisOutput <- SurvivalAnalysis(BetaMatrix = BetaMatrix, 
                                                 ClinicalFile = ClinicalFile, 
                                                 SurvivalFile = SurvivalFile, 
                                                 ClusterLabels = ClusterLabels, 
                                                 FigureGenerate = "Yes", 
                                                 PNGorPDF = "png")
    }
  
    if(FigureGenerate == "Yes") {  
      tSNEOutputFiltered <- tSNEPlotGeneration(BetaMatrix = BetaMatrix_OrderListofProbes, 
                                               PerplexityParameter = 6, 
                                               ClinicalFile = ClinicalFile, 
                                               ClusterLabels = ClusterLabels, 
                                               Filter = NA, 
                                               FigureGenerate = "Yes", 
                                               PNGorPDF = "png", 
                                               ImageName = "5000Probes_AllSamples")
      
      tSNEOutput <- tSNEPlotGeneration(BetaMatrix = BetaMatrix, 
                                       PerplexityParameter = 6, 
                                       ClinicalFile = ClinicalFile, 
                                       ClusterLabels = ClusterLabels, 
                                       Filter = NA, 
                                       FigureGenerate = "Yes", 
                                       PNGorPDF = "png", 
                                       ImageName = "AllProbes_AllSamples")
    }
  
  # Pie chart of probe distribution 
  if(FigureGenerate == "Yes") {  
    matchRowNames <- match(rownames(BetaMatrix_OrderListofProbes), AnnotationFile$V1)
    # (table(AnnotationFile$Relation_to_Island[matchRowNames])/sum(table(AnnotationFile$Relation_to_Island[matchRowNames])))*100
    geneRegion <- sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Group)
    geneRegion[which(geneRegion == "")] <- "NA"
    # table(geneRegion[matchRowNames])
    
    
    geneRegion_breakdown_up <- data.frame(table(AnnotationFile$Relation_to_Island[matchRowNames]))
    # Using plot_ly to create pie-chart for gene region
    geneRegionType <- plotly::plot_ly(geneRegion_breakdown_up, labels = ~Var1, values = ~Freq, 
                                      marker = list(colors = c("#a6cee3","#1f78b4","#fb9a99","#fdbf6f", "#ff7f00", "#cab2d6")),
                                      type = 'pie',
                                      textposition = 'outside',textinfo = 'label+percent') %>%
                                      layout(title = '', showlegend = FALSE, 
                                             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    
    
    
    geneRegion_breakdown_up <- data.frame(table(geneRegion[matchRowNames]))
    # Using plot_ly to create pie-chart for gene region
    geneRegionType <- plotly::plot_ly(geneRegion_breakdown_up, labels = ~Var1, values = ~Freq, 
                                      marker = list(colors = c("#8dd3c7", "#ffffb3", "#bebada","#fb8072", 
                                                               "#80b1d3", "#fdb462", "#fccde5", "#d9d9d9")),
                                      type = 'pie',
                                      textposition = 'outside',textinfo = 'label+percent') %>%
                                      layout(title = '', showlegend = FALSE, 
                                             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    
                                    
    
    # Create heatmap
    # Match significant probes with Beta Matrix
    # matchBetaSigProb_up <- match(rownames(top_P_Coef11[which((top_P_Coef11$logFC > 0) == TRUE), ]), rownames(BetaMatrix_T1))
    # write.csv(x = BetaMatrix_T1[matchBetaSigProb_up, ], 
    #           file = "matchBetaSigProb_up.csv")
    # heatmap(x = BetaMatrix_T1[matchBetaSigProb_up, ])
    
  }
  
  # PCA analysis
  if(FigureGenerate == "Yes") {  
    set.seed(1)
    
    # using selected beta probes 
    res.pcaFiltered <- prcomp(t(BetaMatrix_OrderListofProbes), 
                              center = TRUE, 
                              scale = TRUE)
    
    # cluster
    p0Filtered <- factoextra::fviz_pca_ind(res.pcaFiltered, 
                                           label = "none", 
                                           habillage = factor(ClusterLabels),
                                           addEllipses = TRUE, ellipse.level = 0.95,  
                                           palette = c("#4363d8", "#f58231")) +
                                           theme_minimal()

    
    # type
    p1Filtered <- factoextra::fviz_pca_ind(res.pcaFiltered, 
                                           label = "none", 
                                           habillage = factor(ClinicalFile$TYPE),
                                           addEllipses = TRUE, ellipse.level = 0.95,  
                                           palette = c("#d6604d", "#66bd63", "#4575b4")) +
                                           theme_minimal()
                
    # stage
    p2Filtered <- factoextra::fviz_pca_ind(res.pcaFiltered, 
                                           label = "none", 
                                           habillage = factor(ClinicalFile$STAGE),
                                           addEllipses = TRUE, ellipse.level = 0.95,  
                                           palette = c("#762a83", "#c2a5cf")) +
                                           theme_minimal()
    
    # translocation status
    p3Filtered <- factoextra::fviz_pca_ind(res.pcaFiltered, 
                                           label = "none", 
                                           habillage = factor(ClinicalFile$TRANSLOC_14_18),
                                           addEllipses = TRUE, ellipse.level = 0.95,  
                                           palette = c("#e0e0e0", "#878787")) +
                                           theme_minimal()
    
    
    
    # using all beta probes 
    res.pca <- prcomp(t(BetaMatrix), 
                      center = TRUE, 
                      scale = TRUE)
    
    # cluster
    p0 <- factoextra::fviz_pca_ind(res.pca, 
                                   label = "none", 
                                   habillage = factor(ClusterLabels),
                                   addEllipses = TRUE, ellipse.level = 0.95,  
                                   palette = c("#4363d8", "#f58231")) +
                                   theme_minimal()
    
    
    # type
    p1 <- factoextra::fviz_pca_ind(res.pca, 
                                   label = "none", 
                                   habillage = factor(ClinicalFile$TYPE),
                                   addEllipses = TRUE, ellipse.level = 0.95,  
                                   palette = c("#d6604d", "#66bd63", "#4575b4")) +
                                   theme_minimal()
    
    # stage
    p2 <- factoextra::fviz_pca_ind(res.pca, 
                                   label="none", 
                                   habillage = factor(ClinicalFile$STAGE),
                                   addEllipses = TRUE, ellipse.level = 0.95,  
                                   palette = c("#762a83", "#c2a5cf")) +
                                   theme_minimal()
    
    # translocation status
    p3 <- factoextra::fviz_pca_ind(res.pca, 
                                   label = "none", 
                                   habillage = factor(ClinicalFile$TRANSLOC_14_18),
                                   addEllipses = TRUE, ellipse.level = 0.95,  
                                   palette = c("#e0e0e0", "#878787")) +
                                   theme_minimal()
    
  }
  
  
  # cat("\n Saving results")
  RESULTS <- list(BetaMatrixProbesUsedClustering = BetaMatrix_OrderListofProbes,
                  Kmeans_PValue_4 = Kmeans_PValue_4,                  
                  Kmeans_PValue_3 = Kmeans_PValue_3,
                  Kmeans_PValue_2 = Kmeans_PValue_2,
                  Medoids_PValue_4 = Medoids_PValue_4,
                  Medoids_PValue_3 = Medoids_PValue_3,
                  Medoids_PValue_2 = Medoids_PValue_2,
                  Hierarchical_PValue_4 = Hierarchical_PValue_4,
                  Hierarchical_PValue_3 = Hierarchical_PValue_3,
                  Hierarchical_PValue_2 = Hierarchical_PValue_2,
                  RPMM = RPMM, 
                  InfiniumClust_Clustering = InfiniumClust_Clustering)
  
  class(RESULTS) <- "Clustering_ASilva"
  return(RESULTS)
}