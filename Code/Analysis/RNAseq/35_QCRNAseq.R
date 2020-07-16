# 19 May 2020
# Function: Quality control of RNAseq count matrix using specified criteria. Return matrix
#           RNAseq counts matched to samples of BetaMatrix, removed of low quality samples, 
#           those with 0 library size. Also returns RNAseq counts with features (ENSEMBL IDs) 
#           removed. Finally returns RNAseq matrix of normalized counts after above steps  
#           (removed low quality samples and features). 
# Author: Anjali Silva


# Input:
# RNAseqCountMatrix: A matrix of counts (type "double") size features (e.g., genes) x patients, with 
#                    features in rows and patients as columns. The ncol(RNAseqCountMatrix)
#                    should match ncol(MvalueMatrix), ncol(BetaMatrix), nrow(ClinicalFile) and 
#                    nrow(QCMatrix).
# QCMatrix: A data frame of QC metrics, that has size patients x QC metric categories. The
#           nrow(QCMatrix) should match ncol(RNAseqCountMatrix).
# RNAseqAnnotationFile: A data frame with column name 'gene' giving the ENSEMBL IDs and column name 
#                      'chr' giving the choromosome to which the ENSEMBL ID belong to.  If provided
#                      this will be used to filter ENSEMBLIDs corresponding to autosomes (i.e., 
#                      remove those belonging to sex or mitochondrial chromosomes). Default: NA.
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns.
#             The RNAseqCountMatrix will be checked against samples present in BetaMatrix if provided.
# MvalueMatrix: A matrix of M values (type "double") size probes x patients, with probes in rows and
#               patients as columns. The ncol(MvalueMatrix) should match ncol(BetaMatrix).
# TumorPurity: A matrix of size samples x 1, where column titled "purity" lists the tumor purities. 
# ClinicalFile: A data frame of clinical categories, that has size patients x clinical categories. The
#               nrow(ClinicalFile) should match ncol(BetaMatrix). 
# SurvivalFile: A data frame of survival data,  with size patients x survival categories.
# RNAseqSampleCufoffUniqMapReadCount: An integer indicating the cut off value for uniquely mapped reads 
#                                      from QCMatrix. Samples with uniquely mapped reads below this cutoff 
#                                      will be removed. Default: samples with >30,000,000 will be retained.
# RNAseqSampleCufoffUniqMapReadPercentage: An integer indicating the cut off value for percentage of uniquely
#                                          mapped reads. Default: samples with >= 70 will be retained. 
# RNAseqSampleCufoffReadsMappedMultipleLoci: An integer indicating the cut off value for reads mapped to 
#                                            multiple loci. Default: < 20 will be retained. 
# RNAseqSampleCutoffRRNAcontam: An integer indicating the cut off value for rRNA contamination level from 
#                               QCMatrix. Samples with rRNA contamination above this cutoff will be removed. 
#                               Default: samples with > 40 will be removed. 
# RNAseqFeatureSelectionMethod: A character vector specifying the type of method used for feature
#                               selection for RNAseq data. The options include "SD" for standard
#                               deviation, "MAD" for median absolute deviation, "DifferentialExpression"
#                               for differential expression between FL and RLN, and "edgeR" for using 
#                               edgeR::filterByExpr() which filter out lowly expressed genes. Default
#                               method "edgeR". 
# RNAseqFeatureSelectionCutOff: A numeric value indicating the cutoff value to be used for "SD" or "MAD" 
#                               options in RNAseqFeatureSelectionMethod. Probes meeting a SD or MAD 
#                               value above the input cutoff will be returned. Only RNAseqFeatureSelectionCutOff
#                               or RNAseqFeatureSelectionNumberofProbes should be provided. 
# RNAseqFeatureSelectionNumberofProbes: An integer specifying the number of probes to be retained, if using 
#                                       "SD" or "MAD" options in RNAseqFeatureSelectionMethod. If this is 
#                                       provided, RNAseqFeatureSelectionCutOff will not be considered.
#                                       Only RNAseqFeatureSelectionNumberofProbes or RNAseqFeatureSelectionCutOff
#                                       should be provided. 
# RNAseqNormalizationMethod: A character vector specifying the type of normalization method to use. Options
#                            include "edgeR" for edgeR::cpm() and "StandardNormal", which is offered by 
#                            similarity network fusion (SNF) SNFtool::standardNormalization(). "edgeR" 
#                            could only be used if RNAseqFeatureSelectionMethod is also set to "edgeR". 
#                            Otherwise, StandardNormal will be used. Default method "edgeR".
# ImageName: A character vector specifying a unique name for the images to be produced. 
# PNGorPDF: Output format of the image, options = "png" or "pdf".




# Output
# RNASeqCountMatrixMatchedSampleFilteredNOFeatureFiltered: A matrix of raw counts, size transcripts (ENSEMBL IDs) x samples
#                                                        for RNAseq data, where only samples have been filtered.
#                                                        No features have been filtered.
# RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered: A matrix of raw counts, size transcripts (ENSEMBL IDs) x samples
#                                                        for RNAseq data, where transcripts have been filtered
#                                                        based on i) transcripts with 0 counts across all samples
#                                                        and ii) based on RNAseqFeatureSelectionMethod selected 
#                                                        by user. Samples are also filtered based on criteria 
#                                                        specified by user. 
# RNASeqCountMatrixMatchedSampleFilteredNOFeatureFilteredNoSexENSEMBLIDs: If RNAseqAnnotationFile is provided, then 
#                                                                         RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered
#                                                                         with ENSEMBL IDs corresponding to only autosomes 
#                                                                         (i.e., no sex or mitochondrial chromosomes) will be provided. 
# RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNormalized: A matrix of normalized values, size transcripts 
#                                                                 (ENSEMBL IDs) x samples for RNAseq data, where
#                                                                 transcripts have been filtered based on i) transcripts
#                                                                 with 0 counts across all samples and ii) 
#                                                                 based on RNAseqFeatureSelectionMethod selected 
#                                                                 by user. Samples are also filtered based on criteria 
#                                                                 specified by user. Normalized values are generated
#                                                                 based on method specified by user in 
#                                                                 RNAseqNormalizationMethod.
# RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized: If RNAseqAnnotationFile is provided, then 
#                                                                                 RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNormalized
#                                                                                 with ENSEMBL IDs corresponding to only autosomes 
#                                                                                 (i.e., no sex or mitochondrial chromosomes) 
#                                                                                 will be provided. 
# QCMatrixMatchedSampleFiltered: QCMatrix provided as an argument, updated according to selected sample size. Final
#                                output will be a matrix of samples x QC metrics. 
# BetaMatrixMatchedSampleFiltered: BetaMatrix provided as an argument, updated according to selected sample size. Final
#                                  output will be a matrix of probes x samples. 
# MvalueMatrixMatchedSampleFiltered: MvalueMatrix provided as an argument, updated according to selected sample size. Final
#                                  output will be a matrix of probes x samples. 
# ClinicalFileMatchedSampleFiltered: ClinicalFile provided as an argument, updated according to selected sample size. Final
#                                  output will be a matrix of samples x clinical categories. 
# SurvivalFileMatchedSampleFiltered: SurvivalFile provided as an argument, updated according to selected sample size. Final
#                                  output will be a matrix of samples x survival categories. 
# TumorPurityMatchedSampleFiltered: TumorPurity provided as an argument, updated according to selected sample size. Final
#                                   output will be a matrix of samples x survival categories.                               
# RNAseqAnnotationFileFeatureFilteredNoSexENSEMBLIDs: RNAseqAnnotationFile provided as an argument, updated only autosomes. 
#                                                     Correspond to ENSEMBL IDs in 
#                                                     RNASeqCountMatrixMatchedSampleFilteredNOFeatureFilteredNoSexENSEMBLIDs 
#                                                     and RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized.
                                                     

QCRNAseq <- function(RNAseqCountMatrix, 
                     QCMatrix, 
                     RNAseqAnnotationFile = NA, 
                     BetaMatrix = NA,  
                     MvalueMatrix = NA, 
                     TumorPurity = NA,
                     ClinicalFile = NA,
                     SurvivalFile,
                     RNAseqSampleCufoffUniqMapReadCount = 30000000, 
                     RNAseqSampleCufoffUniqMapReadPercentage = 70,
                     RNAseqSampleCufoffReadsMappedMultipleLoci = 20,
                     RNAseqSampleCutoffRRNAcontam = 40, 
                     RNAseqFeatureSelectionMethod = "edgeR",
                     RNAseqFeatureSelectionCutOff = NA,
                     RNAseqFeatureSelectionNumberofProbes = NA,
                     RNAseqNormalizationMethod = "edgeR",
                     ImageName = "35QCRNAseq",
                     PNGorPDF = "png",
                     ProduceImages = "Yes") {
  
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  library("SNFtool")
  library("edgeR")
  library("EnhancedVolcano")
  
  #  # Checking user input
  # Check if RNAseqCountMatrix and QCMatrix has the same samples, if not stop
  if (length(which((colnames(RNAseqCountMatrix) %in% QCMatrix$SAMPLE_ID) == FALSE)) != 0) {
    stop("\n RNAseqCountMatrix and QCMatrix should have the same samples.") }


  # Checking user input
  if(RNAseqFeatureSelectionMethod == "SD" || RNAseqFeatureSelectionMethod == "MAD") {
    if(is.na(RNAseqFeatureSelectionCutOff) && is.na(RNAseqFeatureSelectionNumberofProbes)) {
      stop("\n If RNAseqFeatureSelectionMethod is set to SD or MAD, need to provide
           a value for RNAseqFeatureSelectionCutOff or RNAseqFeatureSelectionNumberofProbes")
    }
  }
  
  # Match column names of RNAseqCountMatrix and QCMatrix 
  QCMatrixMatched <- QCMatrix[match(colnames(RNAseqCountMatrix), QCMatrix$SAMPLE_ID), ]
  
  
  # Check RNAseqCountMatrix against BetaMatrix, if samples don't match alter both 
  # RNAseqCountMatrix and QCMatrix
  if(class(BetaMatrix) != "logical") {
    if(ncol(RNAseqCountMatrix) != ncol(BetaMatrix)) {
      matchedColNames <- match(colnames(BetaMatrix), colnames(RNAseqCountMatrix))
      if(length(which(is.na(matchedColNames) == TRUE)) > 0) {
        RNAseqCountMatrix <- RNAseqCountMatrix[, matchedColNames[!is.na(matchedColNames)]]
        QCMatrixMatched <- QCMatrixMatched[matchedColNames[! is.na(matchedColNames)], ]
      } else {
        RNAseqCountMatrix <- RNAseqCountMatrix[, match(colnames(BetaMatrix), colnames(RNAseqCountMatrix))]
        QCMatrixMatched <- QCMatrixMatched[match(colnames(BetaMatrix), colnames(RNAseqCountMatrix)), ]
      }
    } 
  }

  
  # Letting the user know about number of observations
  cat("\n", c(nrow(QCMatrixMatched) - length(which(QCMatrixMatched$Uniquely.mapped > 
                                                     RNAseqSampleCufoffUniqMapReadCount))),
      "samples will be removed due to being below 'uniquely mapped count' cutoff of ", 
      RNAseqSampleCufoffUniqMapReadCount, ".")
  
  
  cat("\n", c(nrow(QCMatrixMatched) - length(which(QCMatrixMatched$Uniquely.mapped.reads.. >= 
                                                     RNAseqSampleCufoffUniqMapReadPercentage))),
      "samples will be removed due to being below 'uniquely mapped count percentage' cutoff of ", 
      RNAseqSampleCufoffUniqMapReadPercentage, "%.")
  
  
  cat("\n", c(nrow(QCMatrixMatched) - length(which(QCMatrixMatched$X..of.reads.mapped.to.multiple.loci < 
                                                     RNAseqSampleCufoffReadsMappedMultipleLoci))),
      "samples will be removed due to being above 'reads mapping to multiple loci percentage' cutoff of ", 
      RNAseqSampleCufoffReadsMappedMultipleLoci, ".")


  cat("\n", length(which(QCMatrixMatched$rrna_contam > RNAseqSampleCutoffRRNAcontam)) ,
      "samples will be removed due to being above 'rRNA contamination' cutoff of ", 
      RNAseqSampleCutoffRRNAcontam, ".")
  
  
  SamplesRemove <- unique(c(QCMatrixMatched$SAMPLE_ID[! (QCMatrixMatched$Uniquely.mapped > 
                                                           RNAseqSampleCufoffUniqMapReadCount)],
                          QCMatrixMatched$SAMPLE_ID[! (QCMatrixMatched$Uniquely.mapped.reads.. >= 
                                                         RNAseqSampleCufoffUniqMapReadPercentage)],
                          QCMatrixMatched$SAMPLE_ID[! (QCMatrixMatched$X..of.reads.mapped.to.multiple.loci < 
                                                         RNAseqSampleCufoffReadsMappedMultipleLoci)],
                          QCMatrixMatched$SAMPLE_ID[QCMatrixMatched$rrna_contam > 
                                                      RNAseqSampleCutoffRRNAcontam]))
  
  # length(SamplesRemove) # 88
  # table(substr(SamplesRemove, 4, 5))
  # 132-length(SamplesRemove)
  
  cat("\n", length(SamplesRemove), "samples will be removed all toegther. New total number of samples remaining is ", 
      ncol(RNAseqCountMatrix) - length(SamplesRemove), ".")

  
  
  SamplesRemove <- match(SamplesRemove, colnames(RNAseqCountMatrix))
  if(length(SamplesRemove) > 0) { # if samples to remove are present
    # Filtering samples
    RNASeqCountMatrixMatchedSampleFiltered <- RNAseqCountMatrix[, - c(SamplesRemove)]
    # dim(RNASeqCountMatrixMatchedSampleFiltered) # 57820    58
    
    QCMatrixMatchedSampleFiltered <- QCMatrixMatched[- c(SamplesRemove), ] 
    # dim(QCMatrixMatchedSampleFiltered) #  58 61
    
    if(class(BetaMatrix)[1] != "logical") {
      BetaMatrixMatchedSampleFiltered <- BetaMatrix[, match(colnames(RNASeqCountMatrixMatchedSampleFiltered), 
                                                            colnames(BetaMatrix))]
    }
    # dim(BetaMatrixMatchedSampleFiltered) # 595564     58
    
    if(class(MvalueMatrix)[1] != "logical") {
      MvalueMatrixMatchedSampleFiltered <- MvalueMatrix[, match(colnames(RNASeqCountMatrixMatchedSampleFiltered), 
                                                                colnames(MvalueMatrix))]
    }
    # dim(MvalueMatrixMatchedSampleFiltered) # 595564     58
    
    if(class(ClinicalFile)[1] != "logical") {
      ClinicalFileMatchedSampleFiltered <- ClinicalFile[match(colnames(RNASeqCountMatrixMatchedSampleFiltered), 
                                                              ClinicalFile$SAMPLE_ID), ]
    }
    # dim(ClinicalFileMatchedSampleFiltered) # 58 25
    
    if(class(SurvivalFile)[1] != "logical") {
      SurvivalFileMatchedSampleFiltered <- SurvivalFile[match(substr(colnames(RNASeqCountMatrixMatchedSampleFiltered), 1, 9), 
                                                              SurvivalFile$LY_FL_ID), ]
    }
    # dim(SurvivalFileMatchedSampleFiltered) # 58 23
    
    
    if(class(TumorPurity)[1] != "logical") {
      TumorPurityMatchedSampleFiltered <- as.matrix(TumorPurity[match(colnames(RNASeqCountMatrixMatchedSampleFiltered), rownames(TumorPurity)), ], ncol = 1)
      rownames(TumorPurityMatchedSampleFiltered) <- rownames(TumorPurity)[match(colnames(RNASeqCountMatrixMatchedSampleFiltered), rownames(TumorPurity))]
    }
  } else if(length(SamplesRemove) == 0){
    # Filtering samples
    RNASeqCountMatrixMatchedSampleFiltered <- RNAseqCountMatrix
    # dim(RNASeqCountMatrixMatchedSampleFiltered) # 57820    58
    
    QCMatrixMatchedSampleFiltered <- QCMatrixMatched
    # dim(QCMatrixMatchedSampleFiltered) #  58 61
    
    if(class(BetaMatrix)[1] != "logical") {
      BetaMatrixMatchedSampleFiltered <- BetaMatrix[, match(colnames(RNASeqCountMatrixMatchedSampleFiltered), 
                                                            colnames(BetaMatrix))]
    }
    # dim(BetaMatrixMatchedSampleFiltered) # 595564     58
    
    if(class(MvalueMatrix)[1] != "logical") {
      MvalueMatrixMatchedSampleFiltered <- MvalueMatrix[, match(colnames(RNASeqCountMatrixMatchedSampleFiltered), 
                                                                colnames(MvalueMatrix))]
    }
    # dim(MvalueMatrixMatchedSampleFiltered) # 595564     58
    
    if(class(ClinicalFile)[1] != "logical") {
      ClinicalFileMatchedSampleFiltered <- ClinicalFile[match(colnames(RNASeqCountMatrixMatchedSampleFiltered), 
                                                              ClinicalFile$SAMPLE_ID), ]
    }
    # dim(ClinicalFileMatchedSampleFiltered) # 58 25
    
    if(class(SurvivalFile)[1] != "logical") {
      SurvivalFileMatchedSampleFiltered <- SurvivalFile[match(substr(colnames(RNASeqCountMatrixMatchedSampleFiltered), 1, 9), 
                                                              SurvivalFile$LY_FL_ID), ]
    }
    # dim(SurvivalFileMatchedSampleFiltered) # 58 23
    
    
    if(class(TumorPurity)[1] != "logical") {
      TumorPurityMatchedSampleFiltered <- as.matrix(TumorPurity[match(colnames(RNASeqCountMatrixMatchedSampleFiltered), rownames(TumorPurity)), ], ncol = 1)
      rownames(TumorPurityMatchedSampleFiltered) <- rownames(TumorPurity)[match(colnames(RNASeqCountMatrixMatchedSampleFiltered), rownames(TumorPurity))]
    }
    
  }
  
  
  # RNAseq sample selection - remove samples with 0 library size
  if (length(which(colSums(RNASeqCountMatrixMatchedSampleFiltered) == 0)) > 0) {  
     
     cat("\n", length(which(colSums(RNASeqCountMatrixMatchedSampleFiltered) == 0)) ,
         "samples will be removed due to sample library size of 0 after feature filtering.
         The sample(s):", names(which(colSums(RNASeqCountMatrixMatchedSampleFiltered) == 0)))
     
    RNASeqCountMatrixMatchedSampleFiltered <- RNASeqCountMatrixMatchedSampleFiltered[, - which(colSums(RNASeqCountMatrixMatchedSampleFiltered) == 0)]
    # dim(RNASeqCountMatrixMatchedSampleFiltered) # 39778    
    QCMatrixMatchedSampleFiltered <- QCMatrixMatchedSampleFiltered[- which(colSums(RNASeqCountMatrixMatchedSampleFiltered) == 0), ]
    # dim(QCMatrixMatchedSampleFiltered) # 39778      
    BetaMatrixMatchedSampleFiltered <- BetaMatrixMatchedSampleFiltered[, - which(colSums(RNASeqCountMatrixMatchedSampleFiltered) == 0)]
    # dim(BetaMatrixMatchedSampleFiltered) # 595564     
    MvalueMatrixMatchedSampleFiltered <- MvalueMatrixMatchedSampleFiltered[, - which(colSums(RNASeqCountMatrixMatchedSampleFiltered) == 0)]
    # dim(MvalueMatrixMatchedSampleFiltered) # 595564     96
    ClinicalFileMatchedSampleFiltered <- ClinicalFileMatchedSampleFiltered[- which(colSums(RNASeqCountMatrixMatchedSampleFiltered) == 0), ]
    # dim(ClinicalFileMatchedSampleFiltered) # 58 25
    SurvivalFileMatchedSampleFiltered <- SurvivalFileMatchedSampleFiltered[- which(colSums(RNASeqCountMatrixMatchedSampleFiltered) == 0), ]
    # dim(SurvivalFileMatchedSampleFiltered) #  58 23
    TumorPurity <- TumorPurity[- which(colSums(RNASeqCountMatrixMatchedSampleFiltered) == 0)]
    # length(TumorPurity) # 96
  }
  # dim(RNASeqCountMatrixMatchedSampleFiltered) # 57820    58
  
  
  
  
  # save RNAseq matrix with no features saved
  RNASeqCountMatrixMatchedSampleFilteredNOFeatureFiltered <- RNASeqCountMatrixMatchedSampleFiltered
  
  # RNAseq feature selection
  # RNAseq sample selection - remove features with all zero counts
  # Filter probes with 0 counts across all samples
  if(length(which(rowSums(RNASeqCountMatrixMatchedSampleFiltered) == 0)) > 0) {
    RNASeqCountMatrixMatchedSampleFiltered <- RNASeqCountMatrixMatchedSampleFiltered[- which(rowSums(RNASeqCountMatrixMatchedSampleFiltered) == 0), ]
    # dim(RNASeqCountMatrixMatchedSampleFiltered) # 50911    58
  }

  
  
  set.seed(1234)
  if(RNAseqFeatureSelectionMethod == "SD") { # using standard deviation
    SDProbeFilterRNAseq <- SDeviation(BetaMatrix = RNASeqCountMatrixMatchedSampleFiltered, 
                                      MvalueMatrix = NA, 
                                      CutOff = RNAseqFeatureSelectionCutOff,
                                      NumberofProbes = RNAseqFeatureSelectionNumberofProbes)
    cat("\n", nrow(RNASeqCountMatrixMatchedSampleFiltered) - nrow(SDProbeFilterRNAseq$BetaMatrix_SD_Filtered), 
        "probes will be removed for RNAseq. New probe size is", nrow(SDProbeFilterRNAseq$BetaMatrix_SD_Filtered), ".")
    RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered <- SDProbeFilterRNAseq$BetaMatrix_SD_Filtered
    
  } else if (RNAseqFeatureSelectionMethod == "MAD") { # using median absolute deviance
    MADProbeFilterRNAseq <- MedianAbsoluteDeviation(BetaMatrix = RNASeqCountMatrixMatchedSampleFiltered, 
                                                    MvalueMatrix = NA, 
                                                    CutOff = RNAseqFeatureSelectionCutOff,
                                                    NumberofProbes = RNAseqFeatureSelectionNumberofProbes)
    # dim(MADProbeFilterRNAseq$BetaMatrix_SD_Filtered)
    cat("\n", nrow(RNASeqCountMatrixMatchedSampleFiltered) - nrow(MADProbeFilterRNAseq$BetaMatrix_MAD_Filtered), 
        "probes will be removed for RNAseq. New probe size is", nrow(MADProbeFilterRNAseq$BetaMatrix_MAD_Filtered), ".")
    RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered <- MADProbeFilterRNAseq$BetaMatrix_MAD_Filtered
    
  } else if (RNAseqFeatureSelectionMethod == "DifferentialExpression") {
    
    cat("\n Code is not set up for RNAseqFeatureSelectionMethod = DifferentialExpression.")
    
  } else if(RNAseqFeatureSelectionMethod == "edgeR") { # ClinicalFileMatchedSampleFiltered
    exprsRNAseq <- edgeR::DGEList(counts = RNASeqCountMatrixMatchedSampleFiltered, 
                                  group = factor(ClinicalFileMatchedSampleFiltered$TYPE))
    keepTranscripts <- edgeR::filterByExpr(y = exprsRNAseq$counts,
                                           group = factor(ClinicalFileMatchedSampleFiltered$TYPE))
    exprsRNAseq <- exprsRNAseq[keepTranscripts, ]
    exprsRNAseq$counts <- exprsRNAseq$counts
    RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered <- exprsRNAseq$counts
    cat("\n", nrow(RNASeqCountMatrixMatchedSampleFiltered) - nrow(exprsRNAseq$counts), 
        paste("ENSEMBL IDs will be removed for RNAseq based on RNAseq feature selection method:", 
              RNAseqFeatureSelectionMethod),
        "New ENSEMBL ID length is", nrow(exprsRNAseq$counts), ".")
    
  }


  # RNAseq Normalization
  if(RNAseqNormalizationMethod == "StandardNormalize" || is.na(RNAseqNormalizationMethod)) {
    # As per Dr. Bo Wang's response from 27 Oct 2018
    # "I would recommend to first take transformation of  (X = log2(1 + counts)), 
    # and then apply standardNormalization()"
    Data2RNAseqFilterededgeRNorm <- t(SNFtool::standardNormalization(log2(1 + Data2RNAseq_FilteredSamplesProbes)))
    # dim(Data2RNAseqFilterededgeRNorm) # 34749    58
  } else if (RNAseqNormalizationMethod == "edgeR") {
      if (RNAseqFeatureSelectionMethod == "edgeR") {
        # perform the TMM normalization and display the normalization factors
        Data2RNAseqFilterededgeRNormFactors <- edgeR::calcNormFactors(exprsRNAseq)
        # plotMDS(Data2RNAseqFilterededgeRNorm)
        # Generate matrix
        Data2RNAseqFilterededgeRNorm <- edgeR::cpm(Data2RNAseqFilterededgeRNormFactors, 
                                            normalized.lib.sizes = TRUE, log = TRUE)
        # dim(Data2RNAseqFilterededgeRNorm) # 34749    58
        # hist(exprs, main = "edgeR normalized data")
        
        # Another normalization based on transcripts per million 
        # https://gist.github.com/slowkow/c6ab0348747f86e2748b
        
      } else {
        cat("\n edgeR is not used for RNAseqFeatureSelectionMethod. 
        So edgeR cannot be used for RNAseqNormalizationMethod.
            Using SNFtool::standardNormalization() instead.")
        # As per Dr. Bo Wang's response from 27 Oct 2018
        # "I would recommend to first take transformation of  (X = log2(1 + counts)), 
        # and then apply standardNormalization()"
         Data2RNAseqFilterededgeRNorm <- t(SNFtool::standardNormalization(log2(1 + Data2RNAseq_FilteredSamplesProbes)))
         # dim(Data2RNAseqFilterededgeRNorm)
     }
  }
  
  
  
  # filtering for autosomes, if RNAseqAnnotationFile is provided
  if (all(is.na(RNAseqAnnotationFile)) == FALSE) {
    chrOfIDsTable <- data.frame(chromosome = RNAseqAnnotationFile$chr[match(rownames(RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), 
                                                                            RNAseqAnnotationFile$gene)],
                                name = RNAseqAnnotationFile$name[match(rownames(RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), 
                                                               RNAseqAnnotationFile$gene)],
                                type = RNAseqAnnotationFile$type[match(rownames(RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), 
                                                                      RNAseqAnnotationFile$gene)],
                                RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered)
    # remove all but autosomal chromosomes according to provided annotation files
    targetToKeep <- c(paste0("chr", 1:22), NA) # keep these
    chrOfIDsTableAutosomesOnly <- chrOfIDsTable %>%
                                  rownames_to_column('gene') %>%
                                  dplyr::filter(chromosome %in% targetToKeep) %>%
                                  column_to_rownames('gene')
    cat("\n", nrow(chrOfIDsTable) - nrow(chrOfIDsTableAutosomesOnly),
        "ENSEMBL IDs corresponding to sex (and mitochondrial) chromosomes have been removed.")
    
    # Saving results - if RNAseq annotation file is provided
    RESULTS <- list(RNASeqCountMatrixMatchedSampleFilteredNOFeatureFiltered = RNASeqCountMatrixMatchedSampleFilteredNOFeatureFiltered,
                    RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered = RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered,
                    RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDs = chrOfIDsTableAutosomesOnly[, c(4:ncol(chrOfIDsTableAutosomesOnly))],
                    RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations = chrOfIDsTableAutosomesOnly, 
                    RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNormalized = Data2RNAseqFilterededgeRNorm,
                    RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized = 
                      Data2RNAseqFilterededgeRNorm[match(rownames(chrOfIDsTableAutosomesOnly), rownames(Data2RNAseqFilterededgeRNorm)), ],
                    QCMatrixMatchedSampleFiltered = QCMatrixMatchedSampleFiltered,
                    BetaMatrixMatchedSampleFiltered = BetaMatrixMatchedSampleFiltered,
                    MvalueMatrixMatchedSampleFiltered = MvalueMatrixMatchedSampleFiltered,
                    ClinicalFileMatchedSampleFiltered = ClinicalFileMatchedSampleFiltered,
                    SurvivalFileMatchedSampleFiltered = SurvivalFileMatchedSampleFiltered,
                    TumorPurityMatchedSampleFiltered = TumorPurityMatchedSampleFiltered,
                    RNAseqAnnotationFileFeatureFiltered = chrOfIDsTableAutosomesOnly[, c(1:3)])
    
  } else {
    # Saving results - if NO RNAseq annotation file is provided
    RESULTS <- list(RNASeqCountMatrixMatchedSampleFilteredNOFeatureFiltered = RNASeqCountMatrixMatchedSampleFilteredNOFeatureFiltered,
                    RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered = RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered,
                    RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNormalized = Data2RNAseqFilterededgeRNorm,
                    QCMatrixMatchedSampleFiltered = QCMatrixMatchedSampleFiltered,
                    BetaMatrixMatchedSampleFiltered = BetaMatrixMatchedSampleFiltered,
                    MvalueMatrixMatchedSampleFiltered = MvalueMatrixMatchedSampleFiltered,
                    ClinicalFileMatchedSampleFiltered = ClinicalFileMatchedSampleFiltered,
                    SurvivalFileMatchedSampleFiltered = SurvivalFileMatchedSampleFiltered,
                    TumorPurityMatchedSampleFiltered = TumorPurityMatchedSampleFiltered)
  }
  


  class(RESULTS) <- "QCRNAseq_ASilva"
  return(RESULTS)
  
}
