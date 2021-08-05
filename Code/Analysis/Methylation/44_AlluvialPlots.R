# Updated 3 Aug 2021
# 4 January 2021
# Function: Generate Alluvial Plots
# Author: Anjali Silva


# Input:
# PNGorPDF: Output format of the image, options = "png" or "pdf". Default: png.
# ProduceImages: Produce images or not, options = "Yes" or "No". Default: Yes

# Output:
# FLResultsDataFrame: A data frame containing patients x categories, categories containing "epiCMIT", "Type",
#                    "Stage", "Sex", "BCL2Translocation", "Clusters", "EZH2Mut", "BCL2Mut", "KMT2DMut", "CREBBPMut",
#                    "EP300Mut", "epiCMIThyper", "epiCMIThypo", "indicator"


library(ggplot2)
library(tidyverse)
library(ggalluvial)
library(alluvial)

### Get data ####

# mutation data 
mutMergedT1Robert <- read.csv(paste0(TargetedDNAseqDirPath, "/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl_25Aug2020.csv"), row.names = 1)
dim(mutMergedT1Robert) # 55 138

# RNAseq data
RNAseqQC18June2020T3Samples132 <- QCRNAseq(RNAseqCountMatrix = RNAseqCountMatrixMatched, 
                                           QCMatrix = RNAseqQCFile, 
                                           RNAseqAnnotationFile = ENSMBLid, 
                                           BetaMatrix = BetaMatrix_T1,  
                                           MvalueMatrix = MvalueMatrix_T1, 
                                           TumorPurity = TumorPurity,
                                           ClinicalFile = ClinicalFile_T1,
                                           SurvivalFile = SurvivalFile,
                                           RNAseqSampleCufoffUniqMapReadCount = 10000000, 
                                           RNAseqSampleCufoffUniqMapReadPercentage = 0,
                                           RNAseqSampleCufoffReadsMappedMultipleLoci = 100,
                                           RNAseqSampleCutoffRRNAcontam = 100, 
                                           RNAseqFeatureSelectionMethod = "edgeR",
                                           RNAseqFeatureSelectionCutOff = NA,
                                           RNAseqFeatureSelectionNumberofProbes = NA,
                                           RNAseqNormalizationMethod = "edgeR",
                                           ImageName = "35_QCRNAseq_T3",
                                           PNGorPDF = "png",
                                           ProduceImages = "No")

dim(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized)
# 44219   132





# methylation data 
dim(BetaMatrix_T1[, c(11:165)]) # 595564    155

length(intersect(colnames(mutMergedT1Robert), colnames(BetaMatrix_T1[, c(11:165)]))) # 125
length(setdiff(colnames(BetaMatrix_T1[, c(11:165)]), colnames(mutMergedT1Robert))) # 30
samplesDiffMethDNA <- setdiff(colnames(mutMergedT1Robert), colnames(BetaMatrix_T1[, c(11:165)]))
length(samplesDiffMethDNA)
intersect(samplesDiffMethDNA, colnames(BetaMatrix_T1[, c(11:165)]))

ClinicalFile <- data.table::fread(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/Data_Clinical/sample_annotations_rcd6Nov2019.txt")
samplesDiffMethDNAClin <- match(samplesDiffMethDNA, ClinicalFile$SAMPLE_ID)

### Generate data frame ####

# empty data frame
threeDataType <- data.frame(Samples = c(colnames(BetaMatrix_T1), samplesDiffMethDNA),
                            Methylation = c(rep("Methylation", 170), rep("Not Included", 13)),
                            RNAseq = rep(NA, 183),
                            DNAseq = rep(NA, 183),
                            SNF = rep(NA, 183), 
                            Freq = rep(NA, 183),
                            Number = c(1:183),
                            Type = c(ClinicalFile_T1$TYPE, ClinicalFile$TYPE[samplesDiffMethDNAClin]),
                            Stage = c(ClinicalFile_T1$STAGE, ClinicalFile$STAGE[samplesDiffMethDNAClin]),
                            Institute = c(ClinicalFile_T1$INSTITUTION, ClinicalFile$INSTITUTION[samplesDiffMethDNAClin]))

for(num in 1:nrow(threeDataType)) {
  # check if methylation in RNAseq data
 RNAseq <- which(threeDataType[num, 1] == colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized))
 DNAseq <- which(threeDataType[num, 1] == colnames(mutMergedT1Robert))
 
 if(length(RNAseq) > 0) {
   threeDataType[num, 3] <- "RNAseq" 
 } else {
   threeDataType[num, 3] <- "Not Included"
 }
 
 if(length(DNAseq) > 0) {
   threeDataType[num, 4] <- "DNAseq"
 } else {
   threeDataType[num, 4] <- "Not Included"
 }

 
 if(length(which(threeDataType[num, ] == "Not Included")) == 1) {
   threeDataType[num, 5] <- "Not Included"
   threeDataType[num, 6] <- 2
 } else if(length(which(threeDataType[num, ] == "Not Included")) == 0) {
   threeDataType[num, 5] <- "SNF"
   threeDataType[num, 6] <- 3
 } else if(length(which(threeDataType[num, ] == "Not Included")) == 2) {
   threeDataType[num, 5] <- "Not Included"
   threeDataType[num, 6] <- 1
 } 
 
}

length(which(threeDataType$Freq == 3)) # 101
# replace NA values in Stage
threeDataType$Stage[which(is.na(threeDataType$Stage) == TRUE)] <- "Not Applicable"
# rename samples
threeDataType$Samples <- substr(threeDataType$Samples, 1, 2)

### Generate plot ####

alluvial::alluvial(threeDataType[, c(1:4, 8:10, 5)], freq = threeDataType$Freq,
         col = ifelse(threeDataType$SNF == "SNF", "blue", "grey"),
         border = ifelse(threeDataType$SNF == "SNF", "blue", "grey"),
         hide = threeDataType$Freq == 0,
         cex = 0.7,
         axis_labels = c("Samples (n=183)", "Methylation", "RNAseq", 
                         "DNAseq", "Type", 
                         "Stage", "Institute",
                         "SNF (n=101)"))


# [END]
