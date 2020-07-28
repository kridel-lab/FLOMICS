# 24 July 2020
# Function: In take 1) RNAseq ENSEMBL IDs along with fold change (FC) and 2) methylation probes 
#           with FC OR genes falling within differentially methylated regions and their Stouffer
#           value. Use these to create a plot of expression FC against methylation FC. A probe/gene
#           is considered significant if P value is < 0.05. Gprofiler analysis is run for gene list falling 
#           in each of direction: NegMeth,NegExp; NegMeth,PosExp; PosMeth,NegExp; PosMeth,PosExp.
#           
# Author: Anjali Silva

# Input:
# RNAseqProbesFC: A dataframe such as the output from edgeR::topTags() or similar containing 
#                 column names logFC, PValue, FDR, and ensemblGeneId. 
# RNAseqAnnotationFile: A dataframe with column names gene (e.g., ENSG00000000003), name (e.g., TSPAN6),
#                       chr (e.g., chrX) and type (e.g., protein_coding).
# MethylationProbesFC: A dataframe such as the output from limma::topTable() or similar containing 
#                      column names logFC, P.Value and CpG Probe names as row names.
# MethylationAnnotationFile: Matrix of annotations for all the probes found in BetaMatrix. It is of size
#                 probes x annotations.
# MethylationRegionsFC: A dataframe such as the output from DMRcate::extractRanges() or similar containing
#                       "no.cpgs" (per region), "meandiff", "Stouffer", and "overlapping.genes". Deafult NA. 
# ProduceImages: Produce images or not, options = "Yes" or "No"
# PNGorPDF: Output format of the image, options = "png" or "pdf"

# Output: 
# AllResultsWithDirection: A dataframe containing logFC, P.Value and CpG Probe names as well as gene name, 
#         relation to island, UCSC ref group name along with corresponding RNAseq ENSEMBL ID, significance 
#         (if P value is < 0.05 in both RNAseq and methylation) and direction (e.g., PosMeth,PosExp).
# DirectionWRTUCSCRefGene: A table with columns containing direction (e.g., NegMeth,NegExp; NegMeth,PosExp;
#                          PosMeth,NegExp; PosMeth,PosExp) and rows containing number of probes or genes 
#                          in each UCSC ref gene group (e.g., 1stExon, Body, TSS1500, etc.).
# DirectionWRTRelationToIsland: A table with columns containing direction (e.g., NegMeth,NegExp; NegMeth,PosExp;
#                              PosMeth,NegExp; PosMeth,PosExp) and rows containing number of probes or genes 
#                              in each relation to island (e.g., 1stExon, Body, TSS1500, etc.).
  
# Visuals saved to img folder
# 39_RNAseqVsMethFC_SigNonsig_*, 
# 39_RNAseqVsMethFC_Direction_*
# 39_BarUCSCRefGeneGroup_*
# 39_BarRelationToIsland_*
# 26_Gprofiler() analysis plots

RNAseqVsMethylationFC <- function(RNAseqProbesFC, 
                                  RNAseqAnnotationFile,
                                  MethylationProbesFC, 
                                  MethylationRegionsFC = NA, 
                                  MethylationAnnotationFile,
                                  ProduceImages = "Yes", 
                                  PNGorPDF = "png") {
  
  library(gprofiler2)
  library(grDevices)
  library(ggplot2)
  
  
  # Comparing RNAseq expression FC against methylation probe FC
  if (all(is.na(MethylationRegionsFC)) == TRUE) {
    
    # Create methyaltion data table with MethylationProbesFC, probe name, gene name, relation 
    # to island UCSCRefGene
    matchMethProbesGenes <- match(rownames(MethylationProbesFC), MethylationAnnotationFile$V1)
    methylationResults <- data.table(MethylationProbesFC, 
                                     Probe = rownames(MethylationProbesFC),
                                     GeneName = sub("\\;.*", "", MethylationAnnotationFile$UCSC_RefGene_Name[matchMethProbesGenes]),
                                     RelationToIsland = sub("\\;.*", "", MethylationAnnotationFile$Relation_to_Island[matchMethProbesGenes]),
                                     UCSCRefGene = sub("\\;.*", "", MethylationAnnotationFile$UCSC_RefGene_Group[matchMethProbesGenes]))
    
    
    
    
    # Create RNAseq data table with RNAseqProbesFC
    RNAseqtoGene <- match(substr(RNAseqProbesFC$ensemblGeneId, 1, 15), RNAseqAnnotationFile$gene)
    RNAseqresults <- data.table(RNAseqProbesFC, 
                                GeneName = RNAseqAnnotationFile$name[RNAseqtoGene])
    # Remove entries with no GeneName, i.e., NA
    RNAseqresults <- RNAseqresults[- which(is.na(RNAseqresults$GeneName) == TRUE), ]
    
    
    
    
    # Match RNAseq and methyaltion probes via comparing gene name
    compareRNAseqMethylationGenes <- match(RNAseqresults$GeneName, 
                                           methylationResults$GeneName)
    # Get only the matching entries from methylation data
    matchedMethylationResults <- methylationResults[compareRNAseqMethylationGenes[! is.na(compareRNAseqMethylationGenes)], ]
    # Get only the matching entries from RNAseq data
    matchedRNAseqresults <- RNAseqresults[match(matchedMethylationResults$GeneName, RNAseqresults$GeneName), ]
    
    # Printing the total number of probes/ genes involved for the user
    if(nrow(matchedMethylationResults) == nrow(matchedRNAseqresults)) {
      cat("\n The total number of probes/genes involved: ", nrow(matchedRNAseqresults))
    }
    
    if(length(matchedMethylationResults$GeneName) != length(unique(matchedMethylationResults$GeneName))) {
      cat("\n Different methylation probes map to the same gene name. All probes are used.")
    }
    
    if(length(matchedRNAseqresults$GeneName) != length(unique(matchedRNAseqresults$GeneName))) {
      cat("\n Different RNAseq ENSEMBLIDs map to the same gene name. All ENSEMBLIDs are used.")
    }
    
    
    
    # Vector to determine if both pvalues are significant or not among RNAseq and methylation
    colVector <- rep("NoSig", nrow(matchedMethylationResults))
    for(i in 1:nrow(matchedMethylationResults)) {
      if((matchedMethylationResults$adj.P.Val[i] < 0.05) && (matchedRNAseqresults$PValue[i] < 0.05)) {
        colVector[i] <- "Sig"
      }
    }
    colVector <- relevel(factor(colVector), ref = "Sig")
    cat("\n The number of significant and non significant entries, respectively: \n", table(colVector))
    
    # Obtaining path to save images
    pathNow <- getwd()
    
    # Combine all results with significance 
    combinedResults <- data.table(matchedMethylationResults, matchedRNAseqresults, colVector)
    colnames(combinedResults)[1:10] <- paste0("Methylation", colnames(combinedResults)[1:10])
    colnames(combinedResults)[11:18] <- paste0("RNAseq", colnames(combinedResults)[11:18])
    
    
    
    if(ProduceImages == "Yes") {
      # plot of significant and non-siginficant entries 
      p2 <- ggplot2::ggplot(combinedResults, 
                            aes(x = MethylationlogFC, 
                                y = RNAseqlogFC, 
                                color = factor(colVector))) +
        geom_point(size = 2) +
        scale_y_continuous(name = "RNAseq logFC") +
        labs(x = "Methylation logFC")+
        theme_bw() + 
        theme(text = element_text(size=20), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(face = "bold"),
              axis.text.y = element_text(face = "bold") )+
        scale_color_manual(values = c("#b2182b", "#f0f0f0"),
                           name = "Significance") +
        theme(aspect.ratio = 1, 
              legend.position = "right", 
              panel.background = element_rect(colour = "black", size=1.5),  
              axis.title =  element_text(face = "bold"))
      
      ggplot2::ggsave(paste0(pathNow, "/img/39_RNAseqVsMethFC_SigNonsig.", PNGorPDF))

    }
    
    
    
    # Determine genes with inverse direction 
    inverseCorVector <- as.vector(colVector)
    for(i in 1:nrow(matchedMethylationResults)) {
      if(colVector[i] == "Sig") {
        if((matchedMethylationResults$logFC[i] > 0) && (matchedRNAseqresults$logFC[i] < 0)) {
          inverseCorVector[i] <- "PosMeth,NegExp"
        }
        if((matchedMethylationResults$logFC[i] < 0) && (matchedRNAseqresults$logFC[i] > 0)) {
          inverseCorVector[i] <- "NegMeth,PosExp"
        }
        if((matchedMethylationResults$logFC[i] < 0) && (matchedRNAseqresults$logFC[i] < 0)) {
          inverseCorVector[i] <- "NegMeth,NegExp"
        }
        if((matchedMethylationResults$logFC[i] > 0) && (matchedRNAseqresults$logFC[i] > 0)) {
          inverseCorVector[i] <- "PosMeth,PosExp"
        }
      }
    }
    
    # Combine all results with direction 
    combinedResultsCorrelation <- data.table(combinedResults, "direction" = inverseCorVector)
    
    if(ProduceImages == "Yes") {
      # Plot only significant entries with direction
      p3 <- ggplot2::ggplot(combinedResultsCorrelation[- which(inverseCorVector == "NoSig"), ], 
                            aes(x = MethylationlogFC, y = RNAseqlogFC, color = factor(direction))) +
        geom_point(size = 2) +
        scale_y_continuous(name = "RNAseq logFC") +
        labs(x = "Methylation logFC") +
        theme_bw() + 
        theme(text = element_text(size = 20), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(face = "bold"),
              axis.text.y = element_text(face = "bold") )+
        scale_color_manual(values = c("#4292c6", "#de2d26", "#41ab5d", "#9970ab"),
                           name = "Significance") +
        theme(aspect.ratio = 1, 
              legend.position = "right", 
              panel.background = element_rect(colour = "black", size = 1.5),  
              axis.title =  element_text(face = "bold"))
      ggplot2::ggsave(paste0(pathNow, "/img/39_RNAseqVsMethFC_Direction.", PNGorPDF))
      
      
      
      # Plot stacked bar chart of distribution by UCSC Ref Genome or Relation to Island
      
      TableToPlot <- data.frame(t(table(combinedResultsCorrelation$direction, 
                                        combinedResultsCorrelation$MethylationUCSCRefGene))[, -3])
      coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                          '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                          '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                          '#000075', '#808080')
      
      BarUCSCRefGeneGroup <- TableToPlot %>%
                             ggplot2::ggplot(aes(fill = Var1, y = Freq, x = Var2)) +
                             ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
                             labs(x = "UCSC RefGene Group",
                                  y = "Percentage WRT UCSC RefGene Group",
                                  fill = "UCSC RefGene Group") +
                             scale_y_continuous(labels = scales::percent) +
                             scale_fill_manual(values = coloursBarPlot[1:7]) +
                             theme_bw() + theme( panel.grid.major = element_blank(), 
                                                panel.grid.minor = element_blank(),
                                                axis.text.x = element_text(angle = 90)) +
                             theme(aspect.ratio = 1, text = element_text(size = 15))
      ggplot2::ggsave(paste0(pathNow, "/img/39_BarUCSCRefGeneGroup.", PNGorPDF))
      
      
      TableToPlot <- data.frame(t(table(combinedResultsCorrelation$direction, 
                                        combinedResultsCorrelation$MethylationRelationToIsland))[, -3])
      
      BarRelationToIsland <- TableToPlot %>%
                              ggplot2::ggplot(aes(fill = Var1, y = Freq, x = Var2)) +
                              ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
                              labs(x = "Relation To Island",
                                   y = "Percentage WRT Relation To Island",
                                   fill = "Relation To Island") +
                              scale_y_continuous(labels = scales::percent) +
                              scale_fill_manual(values = coloursBarPlot[1:7]) +
                              theme_bw() + theme( panel.grid.major = element_blank(), 
                                                  panel.grid.minor = element_blank(),
                                                  axis.text.x = element_text(angle = 90)) +
                              theme(aspect.ratio = 1, text = element_text(size = 15))
      ggplot2::ggsave(paste0(pathNow, "/img/39_BarRelationToIsland.", PNGorPDF))

    }
    
    
    # Running gprofiler analysis
    A <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,PosExp"), ]$RNAseqGeneName) 
    B <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,NegExp"), ]$RNAseqGeneName)
    C <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,NegExp"), ]$RNAseqGeneName)
    D <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName)
    E <- unfactor(combinedResults[which(inverseCorVector == "NoSig"), ]$RNAseqGeneName)
    
    gprofilerPosMethPosExp <- GProfilerAnalysis(GeneIDs = list(A),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = "positive", # with respect to expression
                                                ConditionName = "PosMeth,PosExp",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
    gprofilerPosMethNegExp <- GProfilerAnalysis(GeneIDs = list(B),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = "negative", # with respect to expression
                                                ConditionName = "PosMeth,NegExp",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
    gprofilerNegMethNegExp <- GProfilerAnalysis(GeneIDs = list(C),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = "negative", # with respect to expression
                                                ConditionName = "NegMeth,NegExp",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
    gprofilerNegMethPosExp <- GProfilerAnalysis(GeneIDs = list(D),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = "positive", # with respect to expression
                                                ConditionName = "NegMeth,PosExp",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
    gprofilerNoSig <- GProfilerAnalysis(GeneIDs = list(E),
                                        Organism = "hsapiens",
                                        OrderedQuery = TRUE,
                                        PvalAlphaLevel = 0.01,
                                        PositiveorNegFC = NA, # with respect to expression
                                        ConditionName = "NoSig",
                                        ProduceImages = "Yes", 
                                        PNGorPDF = "png")
    
    
    RESULTS <- list(AllResultsWithDirection = combinedResultsCorrelation,
                    DirectionWRTUCSCRefGene = t(table(combinedResultsCorrelation$direction, 
                                                      combinedResultsCorrelation$MethylationUCSCRefGene))[, -3], 
                    DirectionWRTRelationToIsland = t(table(combinedResultsCorrelation$direction, 
                                                           combinedResultsCorrelation$MethylationRelationToIsland))[, -3],
                    gprofilerPosMethPosExp = gprofilerPosMethPosExp,
                    gprofilerPosMethNegExp = gprofilerPosMethNegExp, 
                    gprofilerNegMethNegExp = gprofilerNegMethNegExp,
                    gprofilerNegMethPosExp = gprofilerNegMethPosExp, 
                    gprofilerNoSig = gprofilerNoSig)
    
    class(RESULTS) <- "RNAseqVsMethylationFC_ASilva"
    
    
  } else if (all(is.na(MethylationRegionsFC)) == FALSE) {
    
    # Comparing RNAseq expression FC against genes falling within DMRcate regions
    
    # First take overlapping regions from DMRcate analysis and save them with respective meandiff
    saveDetails <- matrix(NA, nrow = 1, ncol = 2)
    for(i in 1:nrow(MethylationRegionsFC@elementMetadata)) {
      if (MethylationRegionsFC@elementMetadata[i, ]$Stouffer < 0.05) {
        geneNames <- unlist(strsplit(MethylationRegionsFC@elementMetadata[i, ]$overlapping.genes, ","))
        Info <- cbind(geneNames, rep(MethylationRegionsFC@elementMetadata[i, ]$meandiff, length(geneNames)))
        saveDetails <- rbind(saveDetails, Info)
      }
    }
    saveDetails <- as.data.frame(saveDetails[- 1, ]) # Remove first row as it is empty 
    colnames(saveDetails) <- c("GeneName", "MeanDiff")
    saveDetails <- saveDetails[- which(is.na(saveDetails$GeneName) == TRUE), ] # Remove NA values
    
    
    
    # Create RNAseq data table with RNAseqProbesFC
    RNAseqtoGene <- match(substr(RNAseqProbesFC$ensemblGeneId, 1, 15), RNAseqAnnotationFile$gene)
    RNAseqresults <- data.table(RNAseqProbesFC, 
                                GeneName = RNAseqAnnotationFile$name[RNAseqtoGene])
    # Remove entries with no GeneName, i.e., NA
    RNAseqresults <- RNAseqresults[- which(is.na(RNAseqresults$GeneName) == TRUE), ]
    
    
    
    
    # Match RNAseq genes with genes from DMRcate via comparing gene names
    compareRNAseqDMRcateGenes <- match(unfactor(RNAseqresults$GeneName), saveDetails$GeneName)
    
    # Final matched DMRcate analysis
    matchedDMRcateResults <- saveDetails[compareRNAseqDMRcateGenes[! is.na(compareRNAseqDMRcateGenes)], ]

    # Final matched RNAseq analysis
    matchedRNAseqResults <- RNAseqresults[match(matchedDMRcateResults$Gene, RNAseqresults$GeneName), ]

    
    
    
    # Vector to determine if both pvalues are significant or not 
    # Remember for DMRcate only  < 0.05 are already selected, so here only looking at RNAseq
    colVectorDMRcate <- rep("NoSig", nrow(matchedDMRcateResults))
    for(i in 1:nrow(matchedDMRcateResults)) {
      if(matchedRNAseqResults$PValue[i] < 0.05) {
        colVectorDMRcate[i] <- "Sig"
      }
    }
    colVectorDMRcate <- relevel(factor(colVectorDMRcate), ref = "Sig")
    cat("\n The number of significant and non significant entries, respectively: \n", table(colVectorDMRcate))
    
    # Obtaining path to save images
    pathNow <- getwd()
    
    # Combine all results with significance 
    combinedResults <- data.table(matchedDMRcateResults, matchedRNAseqResults, colVectorDMRcate)
    colnames(combinedResults)[1:2] <- paste0("DMRcate", colnames(combinedResults)[1:2])
    colnames(combinedResults)[3:9] <- paste0("RNAseq", colnames(combinedResults)[3:9])
    # Save mean difference as numeric values
    combinedResults$DMRcateMeanDiff <- round(as.numeric(unfactor(combinedResults$DMRcateMeanDiff)), 4)
    
    
    
    
    if(ProduceImages == "Yes") {
      # Plot of significant and non-siginficant entries, all values for DMRcate mean diff > 0.01
      p2 <- ggplot2::ggplot(combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ], 
                            aes(x = DMRcateMeanDiff, y = RNAseqlogFC, color = factor(colVectorDMRcate))) +
                            geom_point(size = 2) +
                            xlim(-0.12, 0.24) + 
                            scale_y_continuous(name = "RNAseq logFC") +
                            labs(x = "DMRcate Mean Difference")+
                            theme_bw() + 
                            theme(text = element_text(size=20), 
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.text.x = element_text(face = "bold"),
                                  axis.text.y = element_text(face = "bold") ) +
                            scale_color_manual(values = c("#b2182b", "#f0f0f0"),
                                               name = "Significance")  +
                            theme(aspect.ratio = 1, 
                                  legend.position = "right", 
                                  panel.background = element_rect(colour = "black", size = 1.5),  
                                  axis.title =  element_text(face = "bold"))
                          
      
      ggplot2::ggsave(paste0(pathNow, "/img/39_RNAseqVsMethDMRcateFC_SigNonsig.", PNGorPDF))
    }
    
    
    # Determine genes with inverse correlation 
    # All values for DMRcate mean diff > 0.01
    inverseCorVectorDMRcate <- rep("NoSig", nrow(combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ]))
    for(i in 1:nrow(combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ])) {
      if(colVectorDMRcate[abs(combinedResults$DMRcateMeanDiff) > 0.01][i] == "Sig") {
        if((combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ]$DMRcateMeanDiff[i] > 0) && 
           (combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ]$RNAseqlogFC[i] < 0)) {
          inverseCorVectorDMRcate[i] <- "PosMeth,NegExp"
        }
        
        if((combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ]$DMRcateMeanDiff[i] < 0) && 
           (combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ]$RNAseqlogFC[i] > 0)) {
          inverseCorVectorDMRcate[i] <- "NegMeth,PosExp"
        }
        
        if((combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ]$DMRcateMeanDiff[i] < 0) && 
           (combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ]$RNAseqlogFC[i] < 0)) {
          inverseCorVectorDMRcate[i] <- "NegMeth,NegExp"
        }
        if((combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ]$DMRcateMeanDiff[i] > 0) &&
           (combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ]$RNAseqlogFC[i] > 0)) {
          inverseCorVectorDMRcate[i] <- "PosMeth,PosExp"
        }
      }
    }

    
    # Combine all results with direction 
    combinedResultsCorrelation <- data.table(combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ],
                                             "direction" = inverseCorVectorDMRcate)
    
    
    if(ProduceImages == "Yes") {
      # Plot only significant entries with direction
      
      p3 <- ggplot2::ggplot(combinedResultsCorrelation[- which(combinedResultsCorrelation$direction == "NoSig"), ], 
                           aes(x = DMRcateMeanDiff, y = RNAseqlogFC, color = factor(direction))) +
                           geom_point(size = 2) +
                           xlim(- 0.12, 0.24) + 
                           scale_y_continuous(name = "RNAseq logFC") +
                           labs(x = "DMRcate Mean Difference") +
                           theme_bw() + 
                           theme(text = element_text(size = 20), 
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.text.x = element_text(face = "bold"),
                                axis.text.y = element_text(face = "bold")) +
                           scale_color_manual(values = c("#4292c6", "#de2d26", "#41ab5d", "#9970ab"),
                                              name = "Significance") +
                           theme(aspect.ratio = 1, 
                                legend.position = "right", 
                                panel.background = element_rect(colour = "black", size = 1.5),  
                                axis.title =  element_text(face = "bold"))
      
      ggplot2::ggsave(paste0(pathNow, "/img/39_RNAseqVsMethDMRcateFC_Direction.", PNGorPDF))
    }
    
    # Running gprofiler analysis
    A <- unfactor(combinedResultsCorrelation[which(combinedResultsCorrelation$direction == "PosMeth,PosExp"), ]$RNAseqGeneName) 
    B <- unfactor(combinedResultsCorrelation[which(combinedResultsCorrelation$direction == "PosMeth,NegExp"), ]$RNAseqGeneName)
    C <- unfactor(combinedResultsCorrelation[which(combinedResultsCorrelation$direction == "NegMeth,NegExp"), ]$RNAseqGeneName)
    D <- unfactor(combinedResultsCorrelation[which(combinedResultsCorrelation$direction == "NegMeth,PosExp"), ]$RNAseqGeneName)
    E <-  unfactor(combinedResultsCorrelation[which(combinedResultsCorrelation$direction == "NoSig"), ]$RNAseqGeneName)
    
    gprofilerPosMethPosExp <- GProfilerAnalysis(GeneIDs = list(A),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = "positive", # with respect to expression
                                                ConditionName = "PosMeth,PosExp",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
    gprofilerPosMethNegExp <- GProfilerAnalysis(GeneIDs = list(B),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = "negative", # with respect to expression
                                                ConditionName = "PosMeth,NegExp",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
    gprofilerNegMethNegExp <- GProfilerAnalysis(GeneIDs = list(C),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = "negative", # with respect to expression
                                                ConditionName = "NegMeth,NegExp",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
    gprofilerNegMethPosExp <- GProfilerAnalysis(GeneIDs = list(D),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = "positive", # with respect to expression
                                                ConditionName = "NegMeth,PosExp",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
    gprofilerNoSig <- GProfilerAnalysis(GeneIDs = list(E),
                                        Organism = "hsapiens",
                                        OrderedQuery = TRUE,
                                        PvalAlphaLevel = 0.01,
                                        PositiveorNegFC = NA, # with respect to expression
                                        ConditionName = "NoSig",
                                        ProduceImages = "Yes", 
                                        PNGorPDF = "png")
    
    
    RESULTS <- list(AllResultsWithDirection = combinedResultsCorrelation,
                    gprofilerPosMethPosExp = gprofilerPosMethPosExp,
                    gprofilerPosMethNegExp = gprofilerPosMethNegExp, 
                    gprofilerNegMethNegExp = gprofilerNegMethNegExp,
                    gprofilerNegMethPosExp = gprofilerNegMethPosExp, 
                    gprofilerNoSig = gprofilerNoSig)
    
    class(RESULTS) <- "RNAseqVsMethylationDMRcateFC_ASilva"
    
  } else {
    stop("MethylationRegionsFC argument should either contain values or NA.")
  }
  
  return(RESULTS)
  
}
