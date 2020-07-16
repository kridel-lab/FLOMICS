# 10 July 2020
# Function: In take MutationCalls matrix and summarize genes with mutations as 0 (absent) or 1  
#           (present), after filtering for specified variant classification. User have the option 
#           to provide cluster labels or clinical category or none. If cluster labels are 
#           provided then association with mutation status will be performed using chi-squared 
#           test and Fisher's exact statistic. If ClinicalCategory is provded then association 
#           with mutation status will be performed based on clinical category. This may take
#           sometime to run, depending on the number of genes to be analyzed for mutations 
#           for each sample. 
# Author: Anjali Silva

# Input:
# RNAseqMutationCalls: A matrix of (SAMPLE_ID x Variants) x column listing SAMPLE_ID, Hugo_Symbol,
#                      Variant_Classification.
# VariantClassificationToRemove: A character vector listing type of variants to be removed prior to
#                                analyzing mutations. E.g., c("UTR3 .","UTR5 .", "downstream .", 
#                                "upstream .", "exonic unknown"). Default NA. 
# ClusterLabels: A matrix of size samples x 2, where the columns are "ID" and "Cluster". "ID" should 
#                list the unique SAMPLE IDs that are also used in RNAseqMutationCalls. "Cluster" 
#                should be a vector of integers specifying cluster. Default NA. If cluster labels 
#                are provided then association with mutation status will be performed using 
#                chi-squared test and Fisher's exact statistic. The "ID"s will be matched to be 
#                same length as number of SAMPLE IDs in RNAseqMutationCalls. Only one of ClusterLabels 
#                or ClinicalCategory could be provided. 
# ClinicalCategory: A matrix of size samples x 2, where the columns are "ID" and "ClinicalCategory".
#                   "ID" should list the unique SAMPLE IDs that are also used in RNAseqMutationCalls. 
#                   "ClinicalCategory" should be a vector specifying clinical category. Default NA. 
#                   If ClinicalCategory is provided then association with mutation status will be performed using 
#                   chi-squared test and Fisher's exact statistic. The "ID"s will be matched to be 
#                   same length as number of SAMPLE IDs in RNAseqMutationCalls. Only one of ClusterLabels 
#                   or ClinicalCategory could be provided. 
# ProduceImages: Produce images or not, options = "Yes" or "No"
# PNGorPDF: Output format of the image, options = "png" or "pdf"

# Output:
# MutationCalls01: A matrix of size samples x genes, containing 0 and 1. Here, 0 means absesnce of 
#                  mutation for the gene in column and 1 means presence of at least one mutation 
#                  for the gene in column. 
# AssociationMutationClustLabs: A matrix of size genes x categories listed below, associating mutations 
#                               with cluster labels. Will only be available if cluster labels are provided.
#                               Categories: Xsquared - Chi-squared test statistic from comparing mutation 
#                                                      status against cluster label; chisq.test(). 
#                                    XsquaredP.Value - Chi-squared test p-value; chisq.test(). 
#                                    XsquaredRes1-C1 - Pearson residuals, (observed - expected) / sqrt(expected) 
#                                                        for first cell in cross tabulation. 
#                                    XsquaredRes0-C1 - Pearson residuals, (observed - expected) / sqrt(expected) 
#                                                        for second cell in cross tabulation. 
#                                    XsquaredRes1-C2 - Pearson residuals, (observed - expected) / sqrt(expected) 
#                                                        for third cell in cross tabulation. 
#                                    XsquaredRes0-C2 - Pearson residuals, (observed - expected) / sqrt(expected) 
#                                                        for fourth cell in cross tabulation. 
#                                    fisherTOddsRatio - An estimate of the odds ratio from Fisher's exact test; fisher.test().
#                                    fisherTP.Value - P-value from Fisher's exact test; fisher.test().
#                                    MinOneMutation - Number of samples (patients) with at least one mutation for a given gene.
#                                    FrequencyMinOneMutation - Number of samples (patients) with at least one mutation for a 
#                                                                given gene, divided by the total number of samples (patients).
#                                    NumAllMutations - Some patients may have multiple mutations for the same gene. Here, 
#                                                        count of all such mutations are kept and summed for a given gene. If
#                                                        NumAllMutations for a given gene is > MinOneMutation, then one or more
#                                                        samples (patients) may have more than one mutation that given gene. 
# AssociationMutationClustLabsSigChiSquared: A matrix size genes x categories (see AssociationMutationClustLabs), but only 
#                                           for genes found to be significant for Chi-squared test p-value < 0.05. 
# AssociationMutationClustLabsSigFisher: A matrix size genes x categories (see AssociationMutationClustLabs), but only 
#                                           for genes found to be significant for Fisher exact test p-value < 0.05. 
# AssociationMutationClustLabsSigChiSquaredAndFisher: A matrix size genes x categories (see AssociationMutationClustLabs),   
#                                                    but for genes found to be significant for both Chi-squared test and 
#                                                    Fisher exact test p-value < 0.05. 
# AssociationMutationClustLabsPMpanel: A matrix size genes x categories, only for genes included in the PM (prespecified) 
#                                      panel.

# Visuals saved to img folder
# 38_RNAseqMutationCalls01SignficancePMPanelFiltered.* : -log2(odds ratio) vs -log10(P value) of genes from PM panel.


RNAseqToMutationCalls01 <- function(RNAseqMutationCalls, 
                                    VariantClassificationToRemove = NA, 
                                    ClusterLabels = NA,
                                    ClinicalCategory = NA, 
                                    ProduceImages = "Yes",
                                    PNGorPDF = "png") {

  
  # Checking user input
  if(! all(is.na(VariantClassificationToRemove))) {
    if(typeof(VariantClassificationToRemove) != "character") {
      stop("\n VariantClassificationToRemove should be NA or a character vector.")
    }
  }
  
  if(! all(is.na(ClusterLabels))  && ! all(is.na(ClinicalCategory))) {
    stop("\n Only one of ClusterLabels or ClinicalCategory could be provided.")
  }
  
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  RNAseqMutationCallsFiltered <- RNAseqMutationCalls
  
  # Removing VariantClassificationToRemove, if specified by user 
  if(! all(is.na(VariantClassificationToRemove))) {
    cat("\n Variants: '", VariantClassificationToRemove, "' will be removed.")
    for (i in 1:length(VariantClassificationToRemove)) {
      RNAseqMutationCallsFiltered <- RNAseqMutationCallsFiltered %>% 
        filter(Variant_Classification != VariantClassificationToRemove[i])
    }
  }
    
    # Creating a matrix to store 0 and 1 values 
    RNAseqMutationCallsFiltered01 <- matrix(0, nrow = length(unique(RNAseqMutationCallsFiltered$SAMPLE_ID)), 
                                            ncol = length(unique(RNAseqMutationCallsFiltered$Hugo_Symbol)))
    # Create a table to keep track of frequency of mutations
    FrequencyofRNAseqMutations <- RNAseqMutationCallsFiltered01 
    rownames(RNAseqMutationCallsFiltered01) <- rownames(FrequencyofRNAseqMutations) <- unique(RNAseqMutationCallsFiltered$SAMPLE_ID)
    colnames(RNAseqMutationCallsFiltered01) <- colnames(FrequencyofRNAseqMutations) <- unique(RNAseqMutationCallsFiltered$Hugo_Symbol)
    
    
    # In the 0 matrix, only write 1 if the muation is present 
    # for each sample
    cat("\n A total of ", length(unique(RNAseqMutationCallsFiltered$SAMPLE_ID)), "samples present to analyze. \n")
    if(length(unique(RNAseqMutationCallsFiltered$SAMPLE_ID)) > 50) {
      cat("\n This may take sometime....") 
    }
    for (i in 1:length(unique(RNAseqMutationCallsFiltered$SAMPLE_ID))) {
      cat("\n Running sample:", i)
      # for each sample, look at how many variants 
      for (j in 1:nrow(filter(RNAseqMutationCallsFiltered, SAMPLE_ID == unique(RNAseqMutationCallsFiltered$SAMPLE_ID)[i]))) {
        # cat("\n Running j:", j); for each variant of a given sample
        # locate where the gene for the variant is located at RNAseqMutationCallsFiltered01 matrix
        colName <- match(filter(RNAseqMutationCallsFiltered, SAMPLE_ID == unique(RNAseqMutationCallsFiltered$SAMPLE_ID)[i])$Hugo_Symbol[j], 
                         colnames(RNAseqMutationCallsFiltered01))
        if (RNAseqMutationCallsFiltered01[i, colName] != 0) {
          # cat("\n gene:", colnames(RNAseqMutationCallsFiltered01)[colName])
          # if RNAseqMutationCallsFiltered01 has already marked a mutation for this gene, keep count of that
          FrequencyofRNAseqMutations[i, colName] <- as.numeric(FrequencyofRNAseqMutations[i, colName]) + 1
          RNAseqMutationCallsFiltered01[i, colName]  <- 1
        } else if(RNAseqMutationCallsFiltered01[i, colName] == 0) {
          # if RNAseqMutationCallsFiltered01 has not marked a mutation for this gene, keep count of that
          RNAseqMutationCallsFiltered01[i, colName] <- FrequencyofRNAseqMutations[i, colName] <- 1
        }
      }
    }
    
  
    # If cluster labels are provided, then association can be tested
    if(! all(is.na(ClusterLabels))) {
      cat("\n Cluster labels have been provided.")
      ClusterLabelsMatched <- ClusterLabels$Cluster[match(rownames(RNAseqMutationCallsFiltered01), 
                                                          ClusterLabels$ID)]
      FrequencyofRNAseqMutationsClusterLabs <- data.frame(cbind(ClusterLabelsMatched, 
                                                                FrequencyofRNAseqMutations))
      RNAseqMutationCallsFiltered01ClusterLabs <- data.frame(cbind(ClusterLabelsMatched, 
                                                                   RNAseqMutationCallsFiltered01))
      colnames(RNAseqMutationCallsFiltered01ClusterLabs) <- c("Cluster", 
                                                              colnames(RNAseqMutationCallsFiltered01ClusterLabs)[- 1])
      colnames(FrequencyofRNAseqMutationsClusterLabs) <- c("Cluster", 
                                                           colnames(FrequencyofRNAseqMutationsClusterLabs)[- 1])
      
      # calculating for all patients
      RNAseqMutationCalls01ClusterLabsSignficance <- matrix(NA, 
                                                            ncol = 11, 
                                                            nrow = ncol(RNAseqMutationCallsFiltered01ClusterLabs))
      colnames(RNAseqMutationCalls01ClusterLabsSignficance) <- c("Xsquared", "XsquaredP.Value", 
                                                                 "XsquaredRes1-C1", "XsquaredRes0-C1", 
                                                                 "XsquaredRes1-C2", "XsquaredRes0-C2", 
                                                                 "fisherTOddsRatio", "fisherTP.Value",
                                                                 "MinOneMutation",
                                                                 "FrequencyMinOneMutation", 
                                                                 "NumAllMutations")
      rownames(RNAseqMutationCalls01ClusterLabsSignficance) <- colnames(RNAseqMutationCallsFiltered01ClusterLabs)
      # dim(RNAseqMutationCalls01ClusterLabsSignficance) #  7545    11
      for (k in 1:ncol(RNAseqMutationCallsFiltered01ClusterLabs)) { # for each gene/mutation
       #  cat("\nk is:", k)
        if(! all(RNAseqMutationCallsFiltered01ClusterLabs[, k] == 0, na.rm = T)) { # only if all entries are not zero
          Table <- table(MutationPresence = RNAseqMutationCallsFiltered01ClusterLabs[, k], 
                         Cluster = RNAseqMutationCallsFiltered01ClusterLabs$Cluster)
          if(k == 2) {
            # Note k =1 will be cluster against cluster and will be ignored
            # print table so user knows reference, which will be the first cell. 
            cat("\n The reference group will be the first cell \n")
            print(Table[c(2, 1), ])
          }
          chisqTest  <- stats::chisq.test(Table[c(2, 1), ])
          fisherTest <- stats::fisher.test(Table[c(2, 1), ])
          
          RNAseqMutationCalls01ClusterLabsSignficance[k, 1] <- chisqTest$statistic
          RNAseqMutationCalls01ClusterLabsSignficance[k, 2] <- chisqTest$p.value
          RNAseqMutationCalls01ClusterLabsSignficance[k, 3] <- chisqTest$residuals[1, 1]
          RNAseqMutationCalls01ClusterLabsSignficance[k, 4] <- chisqTest$residuals[2, 1]
          RNAseqMutationCalls01ClusterLabsSignficance[k, 5] <- chisqTest$residuals[1, 2]
          RNAseqMutationCalls01ClusterLabsSignficance[k, 6] <- chisqTest$residuals[2, 2]
          RNAseqMutationCalls01ClusterLabsSignficance[k, 7] <- fisherTest$estimate
          RNAseqMutationCalls01ClusterLabsSignficance[k, 8] <- fisherTest$p.value
          RNAseqMutationCalls01ClusterLabsSignficance[k, 9] <- length(which(RNAseqMutationCallsFiltered01ClusterLabs[, k] == 1))
          RNAseqMutationCalls01ClusterLabsSignficance[k, 10] <- (length(which(RNAseqMutationCallsFiltered01ClusterLabs[, k] == 1)) / 
                                                                   length(RNAseqMutationCallsFiltered01ClusterLabs[, k]))
          RNAseqMutationCalls01ClusterLabsSignficance[k, 11] <- (sum(FrequencyofRNAseqMutationsClusterLabs[, k]))
        }
      }
      # remove first row as it is for cluster against cluster results
      RNAseqMutationCalls01ClusterLabsSignficance <- RNAseqMutationCalls01ClusterLabsSignficance [- 1, ]
      
      
      # save results that are significant with respect to chi-squared or Fisher's Exact Test or both tests
      AssociationMutationClustLabsSigChiSquared <- RNAseqMutationCalls01ClusterLabsSignficance[which(RNAseqMutationCalls01ClusterLabsSignficance[, 2] < 0.05), ]
      AssociationMutationClustLabsSigFisher <- RNAseqMutationCalls01ClusterLabsSignficance[which(RNAseqMutationCalls01ClusterLabsSignficance[, 8] < 0.05), ]
      AssociationMutationClustLabsSigChiSquaredAndFisher <- RNAseqMutationCalls01ClusterLabsSignficance[
                                                            intersect(which(RNAseqMutationCalls01ClusterLabsSignficance[, 2] < 0.05), 
                                                            which(RNAseqMutationCalls01ClusterLabsSignficance[, 8] < 0.05)), ]

      
      # Genes from PM_FL_PANEL_4June2020
      PMPanelGenes <- c("ARID1A", "ATP6AP1", "ATP6V1B2", "B2M", "BCL7A", "BTG1", "BTG2", "BTK", "CARD11",   
                        "CCND3", "CD58", "CD79B", "CD83", "CDKN1B", "CDKN2A", "CREBBP", "CTSS", "DTX1",     
                        "EBF1", "EEF1A1", "EP300", "EZH2", "FAS", "FOXO1", "GNA13", "HIST1H1B", "HIST1H1C", 
                        "HIST1H1E", "HIST1H2AM", "HLA-DMB", "HVCN1", "IL4R", "IRF5", "IRF8", "ITPKB", "KLHL6",    
                        "KMT2D", "MAP2K1", "MAPK1", "MEF2B", "MEF2C", "MYC", "MYD88", "NM_020921", "NFKBIZ",   
                        "NOTCH1", "NOTCH2", "P2RY8", "PIM1", "POU2AF1", "PTEN", "RB1", "RHOA", "RRAGC",   
                        "RRAS", "S1PR2", "SGK1", "SMARCA4", "SOCS1", "SPEN", "STAT3", "STAT6", "TAF1",     
                        "TBL1XR1", "TCF3", "TNFAIP3", "TNFRSF14", "TP53", "U2AF2", "XBP1", "BCL6", 
                        "BCL2", "NFKBIZ")
      # select entries for only PM_FL_PANEL_4June2020
      matchPMwithRNAseqMutCalls <- match(PMPanelGenes, rownames(RNAseqMutationCalls01ClusterLabsSignficance))
      
      AssociationMutationClustLabsPMpanel <- RNAseqMutationCalls01ClusterLabsSignficance[matchPMwithRNAseqMutCalls[! is.na(matchPMwithRNAseqMutCalls)], ]
      
      if (ProduceImages == "Yes") {
        
        RNAseqMutationCalls01SignficancePMPanelFilteredPlottingT1 <- data.frame(
          AssociationMutationClinicalPMpanel, 
          NegLogOdds = -log2(AssociationMutationClustLabsPMpanel[, 7]),
          NegLogPval = -log10(AssociationMutationClustLabsPMpanel[, 8]))
        
        # Basic dot plot
        # plot P value against odds ratio for 60 matching genes of PM_FL_PANEL
        PlotPMPanel <- ggplot2::ggplot(data = RNAseqMutationCalls01SignficancePMPanelFilteredPlottingT1, 
                                       aes(x = NegLogOdds, 
                                           y = NegLogPval, 
                                           size = FrequencyMinOneMutation)) + 
                       ggplot2::geom_point() + 
                       scale_y_continuous(name = "-log10(P value)") +
                       # geom_text(aes(label = rownames(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlotting)), 
                       # hjust = 0, vjust = 0) + 
                       ggrepel::geom_label_repel(aes(label = rownames(RNAseqMutationCalls01SignficancePMPanelFilteredPlottingT1)),
                                                box.padding   = 0.35, 
                                                point.padding = 0.5,
                                                segment.color = 'grey50') +
                       labs(x = "-log2(odds ratio)") +
                       theme_bw() + 
                       theme(text = element_text(size=20), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())
        
        ggplot2::ggsave(paste0(pathNow, "/img/38_RNAseqMutationCalls01SignficancePMPanelFiltered.", PNGorPDF))
        grDevices::dev.off()
      }
      
      
      RESULTS <- list(MutationCalls01 = RNAseqMutationCallsFiltered01,
                     AssociationMutationClustLabs = RNAseqMutationCalls01ClusterLabsSignficance,
                     AssociationMutationClustLabsSigChiSquared = AssociationMutationClustLabsSigChiSquared,
                     AssociationMutationClustLabsSigFisher = AssociationMutationClustLabsSigFisher,
                     AssociationMutationClustLabsSigChiSquaredAndFisher = AssociationMutationClustLabsSigChiSquaredAndFisher,
                     AssociationMutationClustLabsPMpanel = AssociationMutationClustLabsPMpanel)
      

      
    } else if(! all(is.na(ClinicalCategory))) {
      cat("\n ClinicalCategory have been provided:", levels(factor(ClinicalCategory$ClinicalCategory)), ".")
      
      ClinicalCategoryMatched <- ClinicalCategory$ClinicalCategory[match(rownames(RNAseqMutationCallsFiltered01), 
                                                          ClinicalCategory$ID)]
      FrequencyofRNAseqMutationsClinicalCategory <- data.frame(cbind(as.character(ClinicalCategoryMatched), 
                                                                FrequencyofRNAseqMutations))
      FrequencyofRNAseqMutationsClinicalCategory <- FrequencyofRNAseqMutationsClinicalCategory[- which(is.na(FrequencyofRNAseqMutationsClinicalCategory$V1) == TRUE), ]
      
      RNAseqMutationCallsFiltered01ClinicalCategory <- data.frame(cbind(as.character(ClinicalCategoryMatched),
                                                                   RNAseqMutationCallsFiltered01))
      RNAseqMutationCallsFiltered01ClinicalCategory <- RNAseqMutationCallsFiltered01ClinicalCategory[- which(is.na(RNAseqMutationCallsFiltered01ClinicalCategory$V1) == TRUE), ]
      
      
      colnames(RNAseqMutationCallsFiltered01ClinicalCategory) <- c("ClinicalCategory", 
                                                                 colnames(RNAseqMutationCallsFiltered01ClinicalCategory)[- 1])
      colnames(RNAseqMutationCallsFiltered01ClinicalCategory) <- c("ClinicalCategory", 
                                                                 colnames(RNAseqMutationCallsFiltered01ClinicalCategory)[- 1])
      
      # Calculating for all patients
      RNAseqMutationCalls01ClinicalSignficance <- matrix(NA, 
                                                            ncol = 11, 
                                                            nrow = ncol(RNAseqMutationCallsFiltered01ClinicalCategory))
      colnames(RNAseqMutationCalls01ClinicalSignficance) <- c("Xsquared", "XsquaredP.Value", 
                                                                 paste0("XsquaredRes1-", levels(RNAseqMutationCallsFiltered01ClinicalCategory$ClinicalCategory)[1]), 
                                                                 paste0("XsquaredRes0-", levels(RNAseqMutationCallsFiltered01ClinicalCategory$ClinicalCategory)[1]),
                                                                 paste0("XsquaredRes1-", levels(RNAseqMutationCallsFiltered01ClinicalCategory$ClinicalCategory)[2]),
                                                                 paste0("XsquaredRes0-", levels(RNAseqMutationCallsFiltered01ClinicalCategory$ClinicalCategory)[2]),
                                                                 "fisherTOddsRatio", 
                                                                 "fisherTP.Value",
                                                                 "MinOneMutation",
                                                                 "FrequencyMinOneMutation", 
                                                                 "NumAllMutations")
      rownames(RNAseqMutationCalls01ClinicalSignficance) <- colnames(RNAseqMutationCallsFiltered01ClinicalCategory)
      # dim(RNAseqMutationCalls01ClinicalSignficance) #  7545    11
      for (k in 1:ncol(RNAseqMutationCallsFiltered01ClinicalCategory)) { # for each gene/mutation
         # cat("\nk is:", k)
        if(! all(RNAseqMutationCallsFiltered01ClinicalCategory[, k] == 0, na.rm = T)) { # only if all entries are not zero
          Table <- table(MutationPresence = RNAseqMutationCallsFiltered01ClinicalCategory[, k], 
                         ClinicalCategory = RNAseqMutationCallsFiltered01ClinicalCategory$ClinicalCategory)
          if(k == 2) {
            # Note k =1 will be cluster against cluster and will be ignored
            # print table so user knows reference, which will be the first cell. 
            cat("\n The reference group will be the first cell \n")
            print(Table[c(2, 1), ])
          }
          chisqTest  <- stats::chisq.test(Table[c(2, 1), ])
          fisherTest <- stats::fisher.test(Table[c(2, 1), ])
          
          RNAseqMutationCalls01ClinicalSignficance[k, 1] <- chisqTest$statistic
          RNAseqMutationCalls01ClinicalSignficance[k, 2] <- chisqTest$p.value
          RNAseqMutationCalls01ClinicalSignficance[k, 3] <- chisqTest$residuals[1, 1]
          RNAseqMutationCalls01ClinicalSignficance[k, 4] <- chisqTest$residuals[2, 1]
          RNAseqMutationCalls01ClinicalSignficance[k, 5] <- chisqTest$residuals[1, 2]
          RNAseqMutationCalls01ClinicalSignficance[k, 6] <- chisqTest$residuals[2, 2]
          RNAseqMutationCalls01ClinicalSignficance[k, 7] <- fisherTest$estimate
          RNAseqMutationCalls01ClinicalSignficance[k, 8] <- fisherTest$p.value
          RNAseqMutationCalls01ClinicalSignficance[k, 9] <- length(which(RNAseqMutationCallsFiltered01ClinicalCategory[, k] == 1))
          RNAseqMutationCalls01ClinicalSignficance[k, 10] <- (length(which(RNAseqMutationCallsFiltered01ClinicalCategory[, k] == 1)) / 
                                                                   length(RNAseqMutationCallsFiltered01ClinicalCategory[, k]))
          RNAseqMutationCalls01ClinicalSignficance[k, 11] <- (sum(as.numeric(FrequencyofRNAseqMutationsClinicalCategory[, k])))
        }
      }
      # remove first row as it is for cluster against cluster results
      RNAseqMutationCalls01ClinicalSignficance <- RNAseqMutationCalls01ClinicalSignficance [- 1, ]
      
      
      # save results that are significant with respect to chi-squared or Fisher's Exact Test or both tests
      AssociationMutationClinicalSigChiSquared <- RNAseqMutationCalls01ClinicalSignficance[which(RNAseqMutationCalls01ClinicalSignficance[, 2] < 0.05), ]
      AssociationMutationClinicalSigFisher <- RNAseqMutationCalls01ClinicalSignficance[which(RNAseqMutationCalls01ClinicalSignficance[, 8] < 0.05), ]
      AssociationMutationClinicalSigChiSquaredAndFisher <- RNAseqMutationCalls01ClinicalSignficance[
        intersect(which(RNAseqMutationCalls01ClinicalSignficance[, 2] < 0.05), 
                  which(RNAseqMutationCalls01ClinicalSignficance[, 8] < 0.05)), ]
      
      
      # Genes from PM_FL_PANEL_4June2020
      PMPanelGenes <- c("ARID1A", "ATP6AP1", "ATP6V1B2", "B2M", "BCL7A", "BTG1", "BTG2", "BTK", "CARD11",   
                        "CCND3", "CD58", "CD79B", "CD83", "CDKN1B", "CDKN2A", "CREBBP", "CTSS", "DTX1",     
                        "EBF1", "EEF1A1", "EP300", "EZH2", "FAS", "FOXO1", "GNA13", "HIST1H1B", "HIST1H1C", 
                        "HIST1H1E", "HIST1H2AM", "HLA-DMB", "HVCN1", "IL4R", "IRF5", "IRF8", "ITPKB", "KLHL6",    
                        "KMT2D", "MAP2K1", "MAPK1", "MEF2B", "MEF2C", "MYC", "MYD88", "NM_020921", "NFKBIZ",   
                        "NOTCH1", "NOTCH2", "P2RY8", "PIM1", "POU2AF1", "PTEN", "RB1", "RHOA", "RRAGC",   
                        "RRAS", "S1PR2", "SGK1", "SMARCA4", "SOCS1", "SPEN", "STAT3", "STAT6", "TAF1",     
                        "TBL1XR1", "TCF3", "TNFAIP3", "TNFRSF14", "TP53", "U2AF2", "XBP1", "BCL6", 
                        "BCL2", "NFKBIZ")
      # select entries for only PM_FL_PANEL_4June2020
      matchPMwithRNAseqMutCalls <- match(PMPanelGenes, rownames(RNAseqMutationCalls01ClinicalSignficance))
      
      AssociationMutationClinicalPMpanel <- RNAseqMutationCalls01ClinicalSignficance[matchPMwithRNAseqMutCalls[! is.na(matchPMwithRNAseqMutCalls)], ]
      
      if (ProduceImages == "Yes") {
      
        RNAseqMutationCalls01SignficancePMPanelFilteredPlottingT1 <- data.frame(
          AssociationMutationClinicalPMpanel, 
          NegLogOdds = -log2(AssociationMutationClustLabsPMpanel[, 7]),
          NegLogPval = -log10(AssociationMutationClustLabsPMpanel[, 8]))
        
        # Basic dot plot
        # plot P value against odds ratio for 60 matching genes of PM_FL_PANEL
        PlotPMPanel <- ggplot2::ggplot(data = RNAseqMutationCalls01SignficancePMPanelFilteredPlottingT1, 
                                       aes(x = NegLogOdds, 
                                           y = NegLogPval, 
                                           size = FrequencyMinOneMutation)) + 
                        ggplot2::geom_point() + 
                        scale_y_continuous(name = "-log10(P value)") +
                        # geom_text(aes(label = rownames(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlotting)), 
                        # hjust = 0, vjust = 0) + 
                        ggrepel::geom_label_repel(aes(label = rownames(RNAseqMutationCalls01SignficancePMPanelFilteredPlottingT1)),
                                                      box.padding   = 0.35, 
                                                      point.padding = 0.5,
                                                      segment.color = 'grey50') +
                                                  labs(x = "-log2(odds ratio)") +
                                                  theme_bw() + 
                                                  theme(text = element_text(size=20), 
                                                        panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank())
        
        ggplot2::ggsave(paste0(pathNow, "/img/38_RNAseqMutationCalls01SignficancePMPanelFiltered.", PNGorPDF))
        grDevices::dev.off()
      }
      
      
      
      
      RESULTS <- list(MutationCalls01 = RNAseqMutationCallsFiltered01,
                      AssociationMutationClinicalCategory = RNAseqMutationCalls01ClinicalSignficance,
                      AssociationMutationClinicalCategorySigChiSquared = AssociationMutationClinicalSigChiSquared,
                      AssociationMutationClinicalCategorySigFisher = AssociationMutationClinicalSigFisher,
                      AssociationMutationClinicalCategorySigChiSquaredAndFisher = AssociationMutationClinicalSigChiSquaredAndFisher,
                      AssociationMutationClinicalCategoryPMpanel = AssociationMutationClinicalPMpanel)
      
    } else { 
      # if cluster labels OR ClinicalCategory is not provided
      RESULTS <- list(MutationCalls01 = RNAseqMutationCallsFiltered01)
    }
  

  class(RESULTS) <- "RNAseqToMutationCalls01_ASilva"
  return(RESULTS)
}
  


