# 24 July 2020
# Function: In take differentially expressed probes along with fold change (FC) and differentially   
#           methylated probes with FC OR genes falling within differentially methylated regions 
#           and their Stouffer value. Use these to create a plot of expression FC against methylation 
#           FC. A probe/gene is considered significant if P value is < 0.05. Gprofiler analysis is
#           run for gene list falling in each of direction: NegMeth,NegExp; NegMeth,PosExp;
#           PosMeth,NegExp; PosMeth,PosExp.
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
    E <-  unfactor(combinedResults[which(inverseCorVector == "NoSig"), ]$RNAseqGeneName)
    
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
    
    
    
    
    
  } else if (all(is.na(MethylationRegionsFC)) == FALSE) {
    # Comparing RNAseq expression FC against genes falling within DMRcate regions
    
  }
  
  
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
  return(RESULTS)
  
}










