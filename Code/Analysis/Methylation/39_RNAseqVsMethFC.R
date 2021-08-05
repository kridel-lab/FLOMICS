# Updated 3 Aug 2021
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
                                  ProduceImages = "No", 
                                  PNGorPDF = "png") {
  
  library(gprofiler2)
  library(grDevices)
  library(ggplot2)
  
  
  # Comparing RNAseq expression FC against methylation probe FC
  if (all(is.na(MethylationRegionsFC)) == TRUE) {
    
    # 1st Method - match all genes corresponding to probes of methylation with RNAseq genes
    
    # Create methyaltion data table with MethylationProbesFC, probe name, gene name, relation 
    # to island UCSCRefGene
    matchMethProbesGenes <- match(rownames(MethylationProbesFC), MethylationAnnotationFile$V1)
    methylationResults <- data.table(MethylationProbesFC, 
                                     Probe = rownames(MethylationProbesFC),
                                     GeneName = sub("\\;.*", "", MethylationAnnotationFile$UCSC_RefGene_Name[matchMethProbesGenes]),
                                     RelationToIsland = sub("\\;.*", "", MethylationAnnotationFile$Relation_to_Island[matchMethProbesGenes]),
                                     UCSCRefGene = sub("\\;.*", "", MethylationAnnotationFile$UCSC_RefGene_Group[matchMethProbesGenes]),
                                     Xsome = sub("\\;.*", "", MethylationAnnotationFile$chr[matchMethProbesGenes]))
    # if NA values present for gene name, remove those
    if(length(which(is.na(methylationResults$GeneName) == TRUE)) > 0) {
      methylationResults <- methylationResults[- which(is.na(methylationResults$GeneName) == TRUE), ]
      cat("\n Removing", length(which(is.na(methylationResults$GeneName) == TRUE)), 
          "entries from methylation due to no gene name.")
    }
     # length(unique(methylationResults$GeneName)) # 22560
    
    
    
    # Create RNAseq data table with RNAseqProbesFC
    RNAseqtoGene <- match(substr(RNAseqProbesFC$ensemblGeneId, 1, 15), RNAseqAnnotationFile$gene)
    RNAseqresults <- data.table(RNAseqProbesFC, 
                                GeneName = RNAseqAnnotationFile$name[RNAseqtoGene],
                                Xsome = RNAseqAnnotationFile$chr[RNAseqtoGene],
                                Type = RNAseqAnnotationFile$type[RNAseqtoGene])
    
    # Keep only protein coding genes 
    RNAseqresults <- RNAseqresults[which(RNAseqresults$Type == "protein_coding"), ]
    
    # if NA values present for GeneName, remove those
    if(length(which(is.na(RNAseqresults$GeneName) == TRUE)) > 0) {
      cat("\n Removing", length(which(is.na(RNAseqresults$GeneName) == TRUE)), 
          "entries from RNAseq due to no gene name.")
      RNAseqresults <- RNAseqresults[- which(is.na(RNAseqresults$GeneName) == TRUE), ]
    }
    # dim(RNAseqresults) #  14914     9
    
    
  
    cat("\n Method 1: match genes corresponding to all probes of methylation with genes of RNAseq.")
    cat("\n Total methylation probes/genes available for matching with RNAseq data is", nrow(methylationResults), ".
        Number of unqiue genes here is", length(unique(methylationResults$GeneName)), ".")
    cat("\n Total RNAseq protein-coding only genes available for matching with methylation data is", nrow(RNAseqresults), ".
        Number of unqiue protein-coding genes here is", length(unique(RNAseqresults$GeneName)), ".")
    
    # Match gene names that correspond to methylation probe data with RNAseq gene names
    compareRNAseqMethylationGenes <- match(methylationResults$GeneName, 
                                           RNAseqresults$GeneName)
    
    
    # Remove NA values from compareRNAseqMethylationGenes
    if (length(which(is.na(compareRNAseqMethylationGenes))) > 0) {
      compareRNAseqMethylationGenes <- compareRNAseqMethylationGenes[- which(is.na(compareRNAseqMethylationGenes))]
    }
    matchedRNAseqresults <- RNAseqresults[compareRNAseqMethylationGenes, ]
    # dim(matchedRNAseqresults) # 86230     9
    # Get only the matching entries from RNAseq data
    matchedMethylationResults <- methylationResults[match(matchedRNAseqresults$GeneName, 
                                                          methylationResults$GeneName), ]
    # dim(matchedMethylationResults) # 86230    11
    #identical(unfactor(matchedRNAseqresults$GeneName), matchedMethylationResults$GeneName) # TRUE
    
    # Printing the total number of probes/ genes involved for the user
    cat("\n The total number of matched genes between RNAseq and methylation: ", nrow(matchedRNAseqresults))

    
    
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
    cat("\n Number of unique genes matched across methylation and RNAseq", length(unique(combinedResults$MethylationGeneName)), ".")
    
    
    if(ProduceImages == "Yes") {
      # plot of significant and non-siginficant entries 
      p2 <- ggplot2::ggplot(combinedResults, 
                            aes(x = MethylationlogFC, 
                                y = RNAseqlogFC, 
                                color = factor(colVector))) +
        geom_point(size = 2) +
        scale_y_continuous(name = "RNAseq logFC") +
        labs(x = "Methylation logFC") +
        theme_bw(base_rect_size = 0.6) + 
        theme(text = element_text(size = 20), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(color = "black"),
              axis.text.y = element_text(color = "black"),
              axis.title.x = element_text(color = "black"),
              axis.title.y = element_text(color = "black"))+
        scale_color_manual(values = c("#b2182b", "#f0f0f0"),
                           name = "Significance") +
        theme(aspect.ratio = 1, 
              legend.position = "right", 
              panel.background = element_rect(colour = "black", size = 1.5),  
              axis.title =  element_text(face = "bold"))
      
      ggplot2::ggsave(paste0(pathNow, "/img/39_RNAseqVsMethFC_SigNonsig_Method1.", PNGorPDF),
                      width = 8.5, height = 8.5)

    }
    
    
    
    # Determine genes with inverse direction 
    inverseCorVector <- as.vector(combinedResults$colVector)
    for(i in 1:nrow(combinedResults)) {
      if(colVector[i] == "Sig") {
        if((combinedResults$MethylationlogFC[i] > 0) && (combinedResults$RNAseqlogFC[i] < 0)) {
          inverseCorVector[i] <- "PosMeth,NegExp"
        }
        if((combinedResults$MethylationlogFC[i] < 0) && (combinedResults$RNAseqlogFC[i] > 0)) {
          inverseCorVector[i] <- "NegMeth,PosExp"
        }
        if((combinedResults$MethylationlogFC[i] < 0) && (combinedResults$RNAseqlogFC[i] < 0)) {
          inverseCorVector[i] <- "NegMeth,NegExp"
        }
        if((combinedResults$MethylationlogFC[i] > 0) && (combinedResults$RNAseqlogFC[i] > 0)) {
          inverseCorVector[i] <- "PosMeth,PosExp"
        }
      }
    }
    
    # Combine all results with direction 
    combinedResultsCorrelation <- data.table(combinedResults, "direction" = inverseCorVector)
    cat("\n Summary of direction:")
    print(table(combinedResultsCorrelation$direction))
    cat("\n")
    
    if(ProduceImages == "Yes") {
      
      if(length(which(inverseCorVector == "NoSig")) > 0) {
        # Plot only significant entries with direction
        p3 <- ggplot2::ggplot(combinedResultsCorrelation[- which(inverseCorVector == "NoSig"), ], 
                              aes(x = MethylationlogFC, y = RNAseqlogFC, color = factor(direction))) +
          geom_point(size = 2) +
          # xlim(-2, 2) +
          # scale_y_continuous(name = "RNAseq logFC", breaks = seq (-9, 6, by = 2), limits = c(-9, 6)) +
          scale_y_continuous(name = "RNAseq logFC") +
          labs(x = "Methylation logFC") +
          theme_bw(base_rect_size = 0.6) + 
          theme(text = element_text(size = 20), 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black"),
                axis.title.x = element_text(color = "black"),
                axis.title.y = element_text(color = "black"),
                aspect.ratio = 1, 
                legend.position = "right", 
                legend.text = element_text(size = 20), 
                panel.background = element_rect(colour = "black", size = 0.6)) +
          scale_color_manual(values = c("#4292c6", "#de2d26", "#41ab5d", "#9970ab"),
                             name = "Direction") 
                #axis.title =  element_text(face = "bold"))
        
        ggplot2::ggsave(paste0(pathNow, "/img/39_RNAseqVsMethFC_Direction_Method1.", PNGorPDF),
                        width = 8.5, height = 8.5)
        
        
      } else if(length(which(inverseCorVector == "NoSig")) == 0) {
        # Plot only significant entries with direction
        p3 <- ggplot2::ggplot(combinedResultsCorrelation, 
                              aes(x = MethylationlogFC, y = RNAseqlogFC, color = factor(direction))) +
          geom_point(size = 2) +
          scale_y_continuous(name = "RNAseq logFC") +
          labs(x = "Methylation logFC") +
          theme_bw(base_rect_size = 0.6) + 
          theme(text = element_text(size = 20), 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black"),
                axis.title.x = element_text(color = "black"),
                axis.title.y = element_text(color = "black")) +
          scale_color_manual(values = c("#4292c6", "#de2d26", "#41ab5d", "#9970ab"),
                             name = "Significance") +
          theme(aspect.ratio = 1, 
                legend.position = "right", 
                panel.background = element_rect(colour = "black", size = 1.5))  
                # axis.title =  element_text(face = "bold"))
        ggplot2::ggsave(paste0(pathNow, "/img/39_RNAseqVsMethFC_Direction_Method1.", PNGorPDF),
                        width = 8.5, height = 8.5)
      }


      
      
      # Plot stacked bar chart of distribution by UCSC Ref Genome or Relation to Island
      coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                          '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                          '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                          '#000075', '#808080')
      
      RelationToIsland <- c("Island" = "#a6cee3", "N_Shelf" = "#1f78b4",
                           "N_Shore" = "#fb9a99", "OpenSea" = "#fdbf6f", 
                           "S_Shelf" = "#ff7f00", "S_Shore" = "#cab2d6")
      
      GeneRegion <- c("1stExon" = "#8dd3c7", "3'UTR" = "#ffffb3", 
                     "5'UTR" = "#bebada", "Body" = "#fb8072",
                     "ExonBnd" = "#80b1d3", "NA" = "#fdb462", 
                     "TSS1500" = "#fccde5", "TSS200" = "#d9d9d9")
    }
    
      # if only one unique entry for MethylationUCSCRefGene
      if (length(unique(combinedResultsCorrelation$MethylationUCSCRefGene)) == 1) {
        TableToPlotBarUCSCRefGeneGroupM1 <- data.frame(table(combinedResultsCorrelation$direction, 
                                        combinedResultsCorrelation$MethylationUCSCRefGene))[ - 3, ]
        
      if(ProduceImages == "Yes") {
        BarUCSCRefGeneGroup <- TableToPlotBarUCSCRefGeneGroupM1 %>%
                               ggplot2::ggplot(aes(fill = Var1, y = Freq, x = Var2)) +
                               ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
                               labs(x = "UCSC RefGene Group",
                                   y = "Percentage WRT Direction",
                                   fill = "Direction") +
                               scale_y_continuous(labels = scales::percent) +
                               scale_fill_manual(values = GeneRegion) +
                               theme_bw(base_rect_size = 1) + 
                               theme(axis.text.y = element_text(colour = "black"),
                                     axis.text.x = element_text(colour = "black", angle = 90),
                                     panel.grid.major = element_blank(), 
                                     panel.grid.minor = element_blank(),
                                     axis.title.x = element_text(color = "black"),
                                     axis.title.y = element_text(color = "black")) +
                               theme(aspect.ratio = 1, text = element_text(size = 15))
        ggplot2::ggsave(paste0(pathNow, "/img/39_BarUCSCRefGeneGroup_Method1.", PNGorPDF),
                        width = 8.5, height = 8.5)
      }
        
      } else { # if multiple entries for MethylationUCSCRefGene
        TableToPlotBarUCSCRefGeneGroupM1 <- data.frame(t(table(combinedResultsCorrelation$direction, 
                                          combinedResultsCorrelation$MethylationUCSCRefGene))[ , - 3])
        if(ProduceImages == "Yes") {
          BarUCSCRefGeneGroup <- TableToPlotBarUCSCRefGeneGroupM1 %>%
                                 ggplot2::ggplot(aes(fill = Var1, y = Freq, x = Var2)) +
                                 ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
                                 labs(x = "UCSC RefGene Group",
                                     y = "Percentage WRT UCSC RefGene Group",
                                     fill = "UCSC RefGene Group") +
                                 scale_y_continuous(labels = scales::percent) +
                                 scale_fill_manual(values = GeneRegion) +
                                 theme_bw(base_rect_size = 1) + theme(axis.text.y = element_text(colour = "black"),
                                                    axis.text.x = element_text(colour = "black", angle = 90),
                                                    panel.grid.major = element_blank(), 
                                                    panel.grid.minor = element_blank(),
                                                    axis.title.x = element_text(color = "black"),
                                                    axis.title.y = element_text(color = "black")) +
                                 theme(aspect.ratio = 1, text = element_text(size = 15))
          ggplot2::ggsave(paste0(pathNow, "/img/39_BarUCSCRefGeneGroup_Method1.", PNGorPDF),
                          width = 8.5, height = 8.5)
        }
      }
     

      
      # if only one unique entry for MethylationRelationToIsland
      if (length(unique(combinedResultsCorrelation$MethylationRelationToIsland)) == 1) {
        TableToPlotBarRelationToIslandM1 <- data.frame(table(combinedResultsCorrelation$direction, 
                                        combinedResultsCorrelation$MethylationRelationToIsland))[ - 3, ]
        if(ProduceImages == "Yes") {
        BarRelationToIsland <- TableToPlotBarRelationToIslandM1 %>%
                               ggplot2::ggplot(aes(fill = Var1, y = Freq, x = Var2)) +
                               ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
                               labs(x = "Relation To Island",
                                     y = "Percentage WRT Direction",
                                     fill = "Direction") +
                               scale_y_continuous(labels = scales::percent) +
                               scale_fill_manual(values = RelationToIsland) +
                               theme_bw(base_rect_size = 1) + theme(axis.text.y = element_text(colour = "black"),
                                                  axis.text.x = element_text(colour = "black", angle = 90),
                                                  panel.grid.major = element_blank(), 
                                                  panel.grid.minor = element_blank(),
                                                  axis.title.x = element_text(color = "black"),
                                                  axis.title.y = element_text(color = "black")) +
                               theme(aspect.ratio = 1, text = element_text(size = 15))
        ggplot2::ggsave(paste0(pathNow, "/img/39_BarRelationToIsland_Method1.", PNGorPDF),
                        width = 8.5, height = 8.5)
        }

      } else { # if multiple entries for MethylationRelationToIsland
        TableToPlotBarRelationToIslandM1 <- data.frame(t(table(combinedResultsCorrelation$direction, 
                                          combinedResultsCorrelation$MethylationRelationToIsland))[ , - 3])
        if(ProduceImages == "Yes") {
        BarRelationToIsland <- TableToPlotBarRelationToIslandM1 %>%
                                ggplot2::ggplot(aes(fill = Var1, y = Freq, x = Var2)) +
                                ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
                                labs(x = "Relation To Island",
                                     y = "Percentage WRT Relation To Island",
                                     fill = "Relation To Island") +
                                scale_y_continuous(labels = scales::percent) +
                                scale_fill_manual(values = RelationToIsland) +
                                theme_bw(base_rect_size = 1) + theme(axis.text.y = element_text(colour = "black"),
                                                   axis.text.x = element_text(colour = "black", angle = 90),
                                                   panel.grid.major = element_blank(), 
                                                   panel.grid.minor = element_blank(),
                                                   axis.title.x = element_text(color = "black"),
                                                   axis.title.y = element_text(color = "black")) +
                                theme(aspect.ratio = 1, text = element_text(size = 15))
        ggplot2::ggsave(paste0(pathNow, "/img/39_BarRelationToIsland_Method1.", PNGorPDF),
                        width = 8.5, height = 8.5)
        }
      }
      
      


    
    
    
    
    
    
    
    
    
    # 2nd Method - match a selective probes out of multiple methylation probes that match to a gene
    # in methylation data and match that with RNAseq genes
    #   If majority of the probes have -ve FC, pick the probe with highest -ve FC
    #   If majority of the probes have +ve FC, pick the probe with highest +ve FC    
    #   If equal probes of -ve & +ve FC, mark that gene as 'noisy'
    #     containing multiple equal amounts of probes with -ve and +ve FC
    
    # For duplicated entries, find one value
    uniqueGenesMethylation <- unique(combinedResultsCorrelation$MethylationGeneName[which(combinedResultsCorrelation$colVector == "Sig")])
    # Empty matrix size of unique gene list to save data
    uniqueGenesMethylationDataFrame <- data.frame(matrix(nrow = length(uniqueGenesMethylation), 
                                                         ncol = ncol(methylationResults)))
    colnames(uniqueGenesMethylationDataFrame) <- colnames(methylationResults)
    
    # A loop that goes through each unique gene and see if there are multiple methylation probes
    # matching to that gene. If yes, pick the most representative probe. 
    for (n in 1:length(uniqueGenesMethylation)) {
      # Save gene entry data 
      saveGeneEntry <- methylationResults[which(methylationResults$GeneName == uniqueGenesMethylation[n]), ]
      # only save if P value is significant 
      if(length(which(saveGeneEntry$adj.P.Val < 0.05)) >= 1) {
        # select entries with significant p value
        
        saveGeneEntry <- saveGeneEntry[which(saveGeneEntry$adj.P.Val < 0.05),]
        
        # Obtain the gene, to see if multiple probes match to the gene name
        if (nrow(saveGeneEntry) > 1) { 
          # If only multiple probes match perform this analysis
          signOfEntries <- sign(saveGeneEntry$logFC) # Get the sign of the entries
          if(sum(signOfEntries) < 0) { 
            # If majority of signs are negative
            uniqueGenesMethylationDataFrame[n, ] <- saveGeneEntry[which(max(abs(saveGeneEntry$logFC[which(saveGeneEntry$logFC < 0)])) == abs(saveGeneEntry$logFC)), ]
          } else if(sum(signOfEntries) > 0) { 
            # If majority of signs are positive
            uniqueGenesMethylationDataFrame[n, ] <- saveGeneEntry[which(max(abs(saveGeneEntry$logFC[which(saveGeneEntry$logFC > 0)])) == abs(saveGeneEntry$logFC)), ]
          } else if(sum(signOfEntries) == 0) { 
            # If equal amount of negative and positive signs, mark it as a noisy entry
            uniqueGenesMethylationDataFrame[n, ] <- c("noisy", saveGeneEntry[1, c(2:11)])
          }
        } else {
          # If multiple probes don't match to the gene name, 
          # rather only one gene for one probe,  then save entry
          uniqueGenesMethylationDataFrame[n, ] <-  saveGeneEntry
        }
      }
      
    }
    
    cat("\n ###################################################################################")
    cat("\n Method 2: match only one representative probes of methylation with genes of RNAseq.")
    cat("\n The number of noisy genes are:", length(which(uniqueGenesMethylationDataFrame$logFC == "noisy")))
    
    # Remove noisy genes 
    uniqueGenesMethylationDataFrameNoNoisy <- uniqueGenesMethylationDataFrame[- which(uniqueGenesMethylationDataFrame$logFC == "noisy"), ]
    
    cat("\n Total methylation genes available for matching with RNAseq data after removing noisy genes are", 
        nrow(uniqueGenesMethylationDataFrameNoNoisy), ". Number of unqiue genes here is", 
        length(unique(uniqueGenesMethylationDataFrameNoNoisy$GeneName)), ".")
    cat("\n Total RNAseq protein-coding only genes available for matching with methylation is", nrow(RNAseqresults), ".
        Number of unqiue protein-coding genes here is", length(unique(RNAseqresults$GeneName)), ".")
    
    
    # Match gene names that correspond to methylation unique genes with RNAseq gene names
    compareRNAseqMethylationUniqueGenes <- match(uniqueGenesMethylationDataFrameNoNoisy$GeneName, 
                                                 RNAseqresults$GeneName)
    
    
    # Remove NA values from compareRNAseqMethylationGenes
    if (length(which(is.na(compareRNAseqMethylationUniqueGenes))) > 0) {
      compareRNAseqMethylationUniqueGenes <- compareRNAseqMethylationUniqueGenes[- which(is.na(compareRNAseqMethylationUniqueGenes))]
    }
    matchedRNAseqresultsMethylationUniqueGenes <- RNAseqresults[compareRNAseqMethylationUniqueGenes, ]
    # dim(matchedRNAseqresultsMethylationUniqueGenes) # 12937     9
    # Get only the matching entries from RNAseq data
    matchedMethylationResultsMethylationUniqueGenes <- uniqueGenesMethylationDataFrameNoNoisy[match(matchedRNAseqresultsMethylationUniqueGenes$GeneName, 
                                                                              uniqueGenesMethylationDataFrameNoNoisy$GeneName), ]
    # dim(matchedMethylationResultsMethylationUniqueGenes) #  12937    11
    
    # Printing the total number of probes/ genes involved for the user
    cat("\n The total number of matched genes between RNAseq and methylation: ", nrow(matchedRNAseqresultsMethylationUniqueGenes))
    
    
    
    # Vector to determine if both pvalues are significant or not among RNAseq and methylation
    colVectorMethylationUniqueGenes <- rep("NoSig", nrow(matchedMethylationResultsMethylationUniqueGenes))
    for(i in 1:nrow(matchedMethylationResultsMethylationUniqueGenes)) {
      if((matchedMethylationResultsMethylationUniqueGenes$adj.P.Val[i] < 0.05) && (matchedRNAseqresultsMethylationUniqueGenes$PValue[i] < 0.05)) {
        colVectorMethylationUniqueGenes[i] <- "Sig"
      }
    }
    colVectorMethylationUniqueGenes <- relevel(factor(colVectorMethylationUniqueGenes), ref = "Sig")
    cat("\n The number of significant and non significant entries, respectively: \n", table(colVectorMethylationUniqueGenes))
    

    # Combine all results with significance 
    combinedResultsMethylationUniqueGenes <- data.table(matchedMethylationResultsMethylationUniqueGenes, 
                                                        matchedRNAseqresultsMethylationUniqueGenes, 
                                                        colVectorMethylationUniqueGenes)
    colnames(combinedResultsMethylationUniqueGenes)[1:10] <- paste0("Methylation", colnames(combinedResultsMethylationUniqueGenes)[1:10])
    colnames(combinedResultsMethylationUniqueGenes)[11:18] <- paste0("RNAseq", colnames(combinedResultsMethylationUniqueGenes)[11:18])
    cat("\n Number of unique genes matched across methylation and RNAseq", length(unique(combinedResultsMethylationUniqueGenes$MethylationGeneName)), ".")
    
    
    if(ProduceImages == "Yes") {
      # plot of significant and non-siginficant entries 
      p2 <- ggplot2::ggplot(combinedResultsMethylationUniqueGenes, 
                            aes(x = as.numeric(MethylationlogFC), 
                                y = as.numeric(RNAseqlogFC), 
                                color = factor(colVectorMethylationUniqueGenes))) +
                            geom_point(size = 2) +
                            scale_y_continuous(name = "RNAseq logFC") +
                            labs(x = "Methylation logFC") +
                            theme_bw(base_rect_size = 1) + 
                            theme(text = element_text(size=20), 
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.text.x = element_text(colour = "black"),
                                  axis.text.y = element_text(colour = "black"),
                                  axis.title.x = element_text(color = "black"),
                                  axis.title.y = element_text(color = "black")) +
                            scale_color_manual(values = c("#b2182b", "#f0f0f0"),
                                               name = "Significance") +
                            theme(aspect.ratio = 1, 
                                  legend.position = "right", 
                                  panel.background = element_rect(colour = "black", size = 1.5))  
                                  # axis.title =  element_text(face = "bold"))
                          
      ggplot2::ggsave(paste0(pathNow, "/img/39_RNAseqVsMethFC_SigNonsig_Method2.", PNGorPDF),
                      width = 8.5, height = 8.5)
      
    }
    
    
    
    # Determine genes with inverse direction 
    inverseCorVectorMethylationUniqueGenes <- as.vector(combinedResultsMethylationUniqueGenes$colVectorMethylationUniqueGenes)
    for(i in 1:nrow(matchedMethylationResultsMethylationUniqueGenes)) {
      if(colVectorMethylationUniqueGenes[i] == "Sig") {
        if((matchedMethylationResultsMethylationUniqueGenes$logFC[i] > 0) && (matchedRNAseqresultsMethylationUniqueGenes$logFC[i] < 0)) {
          inverseCorVectorMethylationUniqueGenes[i] <- "PosMeth,NegExp"
        }
        if((matchedMethylationResultsMethylationUniqueGenes$logFC[i] < 0) && (matchedRNAseqresultsMethylationUniqueGenes$logFC[i] > 0)) {
          inverseCorVectorMethylationUniqueGenes[i] <- "NegMeth,PosExp"
        }
        if((matchedMethylationResultsMethylationUniqueGenes$logFC[i] < 0) && (matchedRNAseqresultsMethylationUniqueGenes$logFC[i] < 0)) {
          inverseCorVectorMethylationUniqueGenes[i] <- "NegMeth,NegExp"
        }
        if((matchedMethylationResultsMethylationUniqueGenes$logFC[i] > 0) && (matchedRNAseqresultsMethylationUniqueGenes$logFC[i] > 0)) {
          inverseCorVectorMethylationUniqueGenes[i] <- "PosMeth,PosExp"
        }
      }
    }
    
    # Combine all results with direction 
    combinedResultsCorrelationMethylationUniqueGenes <- data.table(combinedResultsMethylationUniqueGenes, 
                                                                   "direction" = inverseCorVectorMethylationUniqueGenes)
    
    
    cat("\n Summary of results:")
    print(table(combinedResultsCorrelationMethylationUniqueGenes$direction))
    cat("\n")
    
    if(ProduceImages == "Yes") {
      # Plot only significant entries with direction
      p3 <- ggplot2::ggplot(combinedResultsCorrelationMethylationUniqueGenes, 
                            aes(x = as.numeric(MethylationlogFC), 
                                y = as.numeric(RNAseqlogFC), 
                                color = factor(direction))) +
                            geom_point(size = 2) +
                            scale_y_continuous(name = "RNAseq logFC") +
                            labs(x = "Methylation logFC") +
                            theme_bw(base_rect_size = 1) + 
                            theme(text = element_text(size = 20), 
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.text.x = element_text(colour = "black"),
                                  axis.text.y = element_text(colour = "black"),
                                  axis.title.x = element_text(color = "black"),
                                  axis.title.y = element_text(color = "black"))+
                            scale_color_manual(values = c("#4292c6", "#de2d26", "#41ab5d", "#9970ab"),
                                               name = "Significance") +
                            theme(aspect.ratio = 1, 
                                  legend.position = "right", 
                                  panel.background = element_rect(colour = "black", size = 1.5)) 
                                  # axis.title =  element_text(face = "bold"))
        ggplot2::ggsave(paste0(pathNow, "/img/39_RNAseqVsMethFC_Direction_Method2.", PNGorPDF),
                                          width = 8.5, height = 8.5)
      
      
      
      # Plot stacked bar chart of distribution by UCSC Ref Genome or Relation to Island
      
      coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                          '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                          '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                          '#000075', '#808080')
      
      RelationToIsland <- c("Island" = "#a6cee3", "N_Shelf" = "#1f78b4",
                            "N_Shore" = "#fb9a99", "OpenSea" = "#fdbf6f", 
                            "S_Shelf" = "#ff7f00", "S_Shore" = "#cab2d6")
      
      GeneRegion <- c("1stExon" = "#8dd3c7", "3'UTR" = "#ffffb3", 
                      "5'UTR" = "#bebada", "Body" = "#fb8072",
                      "ExonBnd" = "#80b1d3", "NA" = "#fdb462", 
                      "TSS1500" = "#fccde5", "TSS200" = "#d9d9d9")
    }
      
      # if only one unique entry for MethylationUCSCRefGene
      if (length(unique(combinedResultsCorrelationMethylationUniqueGenes$MethylationUCSCRefGene)) == 1) {
        TableToPlotBarUCSCRefGeneGroupM2 <- data.frame(table(combinedResultsCorrelationMethylationUniqueGenes$direction, 
                                        combinedResultsCorrelationMethylationUniqueGenes$MethylationUCSCRefGene))
        
        if(ProduceImages == "Yes") {
        BarUCSCRefGeneGroup <- TableToPlotBarUCSCRefGeneGroupM2 %>%
                                ggplot2::ggplot(aes(fill = Var1, y = Freq, x = Var2)) +
                                ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
                                labs(x = "UCSC RefGene Group",
                                     y = "Percentage WRT Direction",
                                     fill = "Direction") +
                                scale_y_continuous(labels = scales::percent) +
                                scale_fill_manual(values = GeneRegion) +
                                theme_bw(base_rect_size = 1) + 
                                theme(axis.text.y = element_text(colour = "black"),
                                      axis.text.x = element_text(colour = "black", angle = 90),
                                      panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(),
                                      axis.title.x = element_text(color = "black"),
                                      axis.title.y = element_text(color = "black")) +
                                theme(aspect.ratio = 1, text = element_text(size = 15))
        ggplot2::ggsave(paste0(pathNow, "/img/39_BarUCSCRefGeneGroup_Method2.", PNGorPDF),
                        width = 8.5, height = 8.5)
        }
      } else { # if multiple entries for MethylationUCSCRefGene
        TableToPlotBarUCSCRefGeneGroupM2 <- data.frame(t(table(combinedResultsCorrelationMethylationUniqueGenes$direction, 
                                                               combinedResultsCorrelationMethylationUniqueGenes$MethylationUCSCRefGene)))
        if(ProduceImages == "Yes") {
        BarUCSCRefGeneGroup <- TableToPlotBarUCSCRefGeneGroupM2 %>%
                              ggplot2::ggplot(aes(fill = Var1, y = Freq, x = Var2)) +
                              ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
                              labs(x = "UCSC RefGene Group",
                                   y = "Percentage WRT UCSC RefGene Group",
                                   fill = "UCSC RefGene Group") +
                              scale_y_continuous(labels = scales::percent) +
                              scale_fill_manual(values = GeneRegion) +
                              theme_bw(base_rect_size = 1) + theme(
                                axis.text.y = element_text(colour = "black"),
                                axis.text.x = element_text(colour = "black", angle = 90),
                                panel.grid.major = element_blank(), 
                                                  panel.grid.minor = element_blank(),
                                axis.title.x = element_text(color = "black"),
                                axis.title.y = element_text(color = "black")) +
                              theme(aspect.ratio = 1, text = element_text(size = 15))
        ggplot2::ggsave(paste0(pathNow, "/img/39_BarUCSCRefGeneGroup_Method2.", PNGorPDF),
                        width = 8.5, height = 8.5)
        }
      }
      
      
      
      # if only one unique entry for MethylationRelationToIsland
      if (length(unique(combinedResultsCorrelationMethylationUniqueGenes$MethylationRelationToIsland)) == 1) {
        TableToPlotBarRelationToIslandM2 <- data.frame(table(combinedResultsCorrelationMethylationUniqueGenes$direction, 
                                                       combinedResultsCorrelationMethylationUniqueGenes$MethylationRelationToIsland))
        if(ProduceImages == "Yes") {
        BarRelationToIsland <- TableToPlotBarRelationToIslandM2 %>%
                                ggplot2::ggplot(aes(fill = Var1, y = Freq, x = Var2)) +
                                ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
                                labs(x = "Relation To Island",
                                     y = "Percentage WRT Direction",
                                     fill = "Direction") +
                                scale_y_continuous(labels = scales::percent) +
                                scale_fill_manual(values = RelationToIsland) +
                                theme_bw(base_rect_size = 1) + theme( 
                                  axis.text.y = element_text(colour = "black"),
                                  axis.text.x = element_text(colour = "black", angle = 90),
                                  panel.grid.major = element_blank(), 
                                                    panel.grid.minor = element_blank(),
                                  axis.title.x = element_text(color = "black"),
                                  axis.title.y = element_text(color = "black")) +
                                theme(aspect.ratio = 1, text = element_text(size = 15))
        ggplot2::ggsave(paste0(pathNow, "/img/39_BarRelationToIsland_Method2.", PNGorPDF),
                        width = 8.5, height = 8.5)
        }
        
      } else { # if multiple entries for MethylationRelationToIsland
        TableToPlotBarRelationToIslandM2 <- data.frame(t(table(combinedResultsCorrelationMethylationUniqueGenes$direction, 
                                          combinedResultsCorrelationMethylationUniqueGenes$MethylationRelationToIsland)))
        if(ProduceImages == "Yes") {
        BarRelationToIsland <- TableToPlotBarRelationToIslandM2 %>%
                                ggplot2::ggplot(aes(fill = Var1, y = Freq, x = Var2)) +
                                ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
                                labs(x = "Relation To Island",
                                     y = "Percentage WRT Relation To Island",
                                     fill = "Relation To Island") +
                                scale_y_continuous(labels = scales::percent) +
                                scale_fill_manual(values = RelationToIsland) +
                                theme_bw(base_rect_size = 1) + theme(axis.text.y = element_text(colour = "black"),
                                                                     axis.text.x = element_text(colour = "black", angle = 90),
                                                                     panel.grid.major = element_blank(), 
                                                                     panel.grid.minor = element_blank(),
                                                                     axis.title.x = element_text(color = "black"),
                                                                     axis.title.y = element_text(color = "black")) +
                                theme(aspect.ratio = 1, text = element_text(size = 15))
        ggplot2::ggsave(paste0(pathNow, "/img/39_BarRelationToIsland_Method2.", PNGorPDF),
                        width = 8.5, height = 8.5)
        }
      }
    
    
    
    
    
    
    
    # gProfiler analysis done based on method??
    gProfilerNotRun <- function() {
      # Running gprofiler analysis
      A <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,PosExp"), ]$RNAseqGeneName) 
      B <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,NegExp"), ]$RNAseqGeneName)
      C <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,NegExp"), ]$RNAseqGeneName)
      D <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName)
      E <- unfactor(combinedResults[which(inverseCorVector == "NoSig"), ]$RNAseqGeneName)
      
      gprofilerPosMethPosExp <- GProfilerAnalysis(GeneIDs = list(unique(A)),
                                                  Organism = "hsapiens",
                                                  OrderedQuery = TRUE,
                                                  PvalAlphaLevel = 0.01,
                                                  PositiveorNegFC = "positive", # with respect to expression
                                                  ConditionName = "PosMeth,PosExp",
                                                  ProduceImages = ProduceImages, 
                                                  PNGorPDF = "png")
      gprofilerPosMethNegExp <- GProfilerAnalysis(GeneIDs = list(unique(B)),
                                                  Organism = "hsapiens",
                                                  OrderedQuery = TRUE,
                                                  PvalAlphaLevel = 0.01,
                                                  PositiveorNegFC = "negative", # with respect to expression
                                                  ConditionName = "PosMeth,NegExp",
                                                  ProduceImages = ProduceImages, 
                                                  PNGorPDF = "png")
      gprofilerNegMethNegExp <- GProfilerAnalysis(GeneIDs = list(unique(C)),
                                                  Organism = "hsapiens",
                                                  OrderedQuery = TRUE,
                                                  PvalAlphaLevel = 0.01,
                                                  PositiveorNegFC = "negative", # with respect to expression
                                                  ConditionName = "NegMeth,NegExp",
                                                  ProduceImages = ProduceImages, 
                                                  PNGorPDF = "png")
      gprofilerNegMethPosExp <- GProfilerAnalysis(GeneIDs = list(unique(D)),
                                                  Organism = "hsapiens",
                                                  OrderedQuery = TRUE,
                                                  PvalAlphaLevel = 0.01,
                                                  PositiveorNegFC = "positive", # with respect to expression
                                                  ConditionName = "NegMeth,PosExp",
                                                  ProduceImages = ProduceImages, 
                                                  PNGorPDF = "png")
      gprofilerNoSig <- GProfilerAnalysis(GeneIDs = list(unique(E)),
                                          Organism = "hsapiens",
                                          OrderedQuery = TRUE,
                                          PvalAlphaLevel = 0.01,
                                          PositiveorNegFC = NA, # with respect to expression
                                          ConditionName = "NoSig",
                                          ProduceImages = ProduceImages, 
                                          PNGorPDF = "png")
      
    }
    
    
    
    
    
    
    
    RESULTS <- list(AllResultsWithDirectionMethod1 = combinedResultsCorrelation,
                    DirectionWRTUCSCRefGeneMethod1 = TableToPlotBarUCSCRefGeneGroupM1,
                    DirectionWRTRelationToIslandMethod1 = TableToPlotBarRelationToIslandM1,
                    AllResultsWithDirectionMethod2 = combinedResultsCorrelationMethylationUniqueGenes,
                    DirectionWRTUCSCRefGeneMethod2 = TableToPlotBarUCSCRefGeneGroupM2,
                    DirectionWRTRelationToIslandMethod2 = TableToPlotBarRelationToIslandM2  
                    
                    #gprofilerPosMethPosExp = gprofilerPosMethPosExp,
                    #gprofilerPosMethNegExp = gprofilerPosMethNegExp, 
                    #gprofilerNegMethNegExp = gprofilerNegMethNegExp,
                    #gprofilerNegMethPosExp = gprofilerNegMethPosExp, 
                    #gprofilerNoSig = gprofilerNoSig
                    )
    
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
    if(length(which(is.na(RNAseqresults$GeneName) == TRUE)) > 0) {
        RNAseqresults <- RNAseqresults[- which(is.na(RNAseqresults$GeneName) == TRUE), ]
    }
    
    
    
    
    # Match RNAseq genes with genes from DMRcate via comparing gene names
    compareRNAseqDMRcateGenes <- match(RNAseqresults$GeneName, saveDetails$GeneName)
    
    # Final matched DMRcate analysis
    if(length(which(is.na(compareRNAseqDMRcateGenes) == TRUE)) > 0) {
        matchedDMRcateResults <- saveDetails[compareRNAseqDMRcateGenes[! is.na(compareRNAseqDMRcateGenes)], ]
    }
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
    combinedResults$DMRcateMeanDiff <- round(as.numeric(combinedResults$DMRcateMeanDiff), 4)
    
    
    
    
    if(ProduceImages == "Yes") {
      # Plot of significant and non-siginficant entries, all values for DMRcate mean diff > 0.01
      p2 <- ggplot2::ggplot(combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ], 
                            aes(x = DMRcateMeanDiff, y = RNAseqlogFC, color = factor(colVectorDMRcate))) +
                            geom_point(size = 2) +
                            xlim(-0.12, 0.24) + 
                            scale_y_continuous(name = "RNAseq logFC") +
                            labs(x = "DMRcate Mean Difference")+
                            theme_bw(base_rect_size = 1) + 
                            theme(text = element_text(size = 20), 
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.text.x = element_text(color = "black"),
                                  axis.text.y = element_text(color = "black") ) +
                            scale_color_manual(values = c("#b2182b", "#f0f0f0"),
                                               name = "Significance")  +
                            theme(aspect.ratio = 1, 
                                  legend.position = "right", 
                                  panel.background = element_rect(colour = "black", size = 1.5),  
                                  axis.title =  element_text(face = "bold"))
                          
      
      ggplot2::ggsave(paste0(pathNow, "/img/39_RNAseqVsMethDMRcateFC_SigNonsig.", PNGorPDF),
                      width = 8.5, height = 8.5)
    }
    
    
    # Determine genes with inverse correlation 
    # All values for DMRcate mean diff > 0.01
    inverseCorVectorDMRcate <- rep("NoSig,DMRcateMeanDiff>0.01", 
                                   nrow(combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ]))
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
    cat("\n Breakdown of results:")
    print(table(combinedResultsCorrelation$direction))
    cat("\n")
    
    if(ProduceImages == "Yes") {
      # Plot only significant entries with direction
      
      p3 <- ggplot2::ggplot(combinedResultsCorrelation[- 
                            which(combinedResultsCorrelation$direction == "NoSig"), ], 
                            aes(x = DMRcateMeanDiff, y = RNAseqlogFC, color = factor(direction))) +
                            geom_point(size = 2) +
                            xlim(- 0.12, 0.24) + 
                            scale_y_continuous(name = "RNAseq logFC") +
                            labs(x = "DMRcate Mean Difference") +
                            theme_bw(base_rect_size = 1) + 
                            theme(text = element_text(size = 20), 
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.text.x = element_text(color = "black"),
                                axis.text.y = element_text(color = "black")) +
                            scale_color_manual(values = c("#4292c6", "#de2d26", "#41ab5d", "#9970ab"),
                                              name = "Significance") +
                            theme(aspect.ratio = 1, 
                                legend.position = "right", 
                                panel.background = element_rect(colour = "black", size = 1.5),  
                                axis.title =  element_text(face = "bold"))
      
      ggplot2::ggsave(paste0(pathNow, "/img/39_RNAseqVsMethDMRcateFC_Direction.", PNGorPDF),
                      width = 8.5, height = 8.5)
    }
    
    # Running gprofiler analysis
    A <- (combinedResultsCorrelation[which(combinedResultsCorrelation$direction == 
                                             "PosMeth,PosExp"), ]$RNAseqGeneName) 
    B <- (combinedResultsCorrelation[which(combinedResultsCorrelation$direction == 
                                             "PosMeth,NegExp"), ]$RNAseqGeneName)
    C <- (combinedResultsCorrelation[which(combinedResultsCorrelation$direction == 
                                             "NegMeth,NegExp"), ]$RNAseqGeneName)
    D <- (combinedResultsCorrelation[which(combinedResultsCorrelation$direction == 
                                             "NegMeth,PosExp"), ]$RNAseqGeneName)
    E <- (combinedResultsCorrelation[which(combinedResultsCorrelation$direction == 
                                             "NoSig"), ]$RNAseqGeneName)
    
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
# [END]
