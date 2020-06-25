# Updated on 6 May 2019
# Reference: http://bioconductor.org/packages/devel/bioc/vignettes/missMethyl/inst/doc/missMethyl.html#testing-for-differential-methylation-using
# Author: Anjali Silva

########################################## #
# Source files
########################################## #
# ** Remember: If requesting analysis for "genes", need to provide summarized data matrices. ***
# ** Visualizations are typically done using Beta values and statistical calculations are done using M-values. ***
cat("\n Sourcing files")
source("0_Package_check.R") # This needs to be run separately as of May 9, 2019, if running on cluster
source("1_QCRemoveSamples.R")
source("2_LoadingMethylData.R")
source("3_DensityPlotMethylation.R")
source("4_SummarizationMvalues.R")
source("5_SummarizationBeta.R")
source("6_SummaryStatistics.R")
source("7_DifferentialMethylation.R")
source("8_DifferentialVariability.R")
source("9_DifferentiallyMethylatedRegions.R")
source("10_BoxPlotsMethylation.R") # previously 10_ViolinPlotsMethylation.R.R
source("11_HeatPlotMethylation.R")
source("12_EstimateCellCountsMethylation.R")
source("13_tSNEPlot.R")
source("14_Clustering_RPMM_InfiniumClust.R")
source("15_SurvivalAnalysis.R")
source("16_Huet23GeneModel.R")
#source("17_LinePlots.R") #* Under construction
source("18_StandardDeviation.R")
source("19_IsolateEntries.R")
source("20_GlmnetFeatureSelection.R")
source("21_ProportionVisualization.R") #* Under construction with bumphunter
source("22_MDSPlots.R")
source("23_Barplot.R")
source("24_Tumor_purity_check.R")
source("25_Variance.R")
source("26_Gprofiler.R")
source("27_CopyNumberAnalysis.R")
source("30_PieChart.R")
source("31_MeanSDPlot.R")
source("32_ClusterConfidencePlot.R")
source("33_SNFClustering.R")
source("34_MedianAbsoluteDeviation.R")
source("35_QCRNAseq.R")
source("36_DifferentialExpressionRNAseq.R")
source("37_ViolinPlotsRNAseq.R")
cat("\n Sourced files")

# Suggestions (to be done)
# MethylMix: Identifying methylation driven cancer genes
# https://www.bioconductor.org/packages/release/bioc/html/MethylMix.html
# Add bar plot

# Set pathways:
RNAseqDirPath <- "/Users/anjalisilva/Desktop/UHN/FLOMICS/RNAseq/ExtendedStudy2020/"
MethylationDirPath <- "/Users/anjalisilva/Desktop/UHN/FLOMICS/Methylation/Pipeline"
TargetedDNAseqDirPath <- "~/Desktop/UHN/FLOMICS/TargetedSequencing"
CodeDirPath <- "/Volumes/GoogleDrive/My Drive/UHN/FLOMICS/FLOMICS-Anjali/Methylation/Pipeline"

########################################## #
# Uploading needed datasets
########################################## #
# Annotation file sent by Alberto on 9 Nov 2018, should be in current directory
# library(readr)
# AnnotationFile <- read_csv(file="Ann_IlluminaHumanMethylationEPICanno.ilm10b2.hg19.csv")
# AnnotationFile <- data.table::fread(file = "Ann_IlluminaHumanMethylationEPICanno.ilm10b2.hg19.csv")
AnnotationFile <-  readRDS(file = paste0(MethylationDirPath, "/AnnotationFile_EPIC.rds"))
dim(AnnotationFile) # 866836     47
# write.csv(AnnotationFile, file = "Ann_IlluminaHumanMethylationEPICanno.ilm10b2.hg19.csv")
# AnnotationFile <- read.csv("Ann_IlluminaHumanMethylationEPICanno.ilm10b2.hg19.csv")

# SurvivalFile <- read_delim(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/Data_Survival/clinical_data_rcd12June2019.txt", delim = '\t')
# SurvivalFile <- read.delim(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/Data_Survival/clinical_data_rcd12June2019.txt")
# SurvivalFile <- data.table::fread(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/Data_Survival/clinical_data_rcd17Jan2020.txt")
# cat("\n Dimension of SurvivalFile is", dim(SurvivalFile)) # 496 23
SurvivalFile <- readRDS(file = paste0(MethylationDirPath, "/SurvivalFile.rds"))
dim(SurvivalFile) # 496  23


# From PipelineStep1
sheet <- readRDS(file = paste0(MethylationDirPath, "/1_MethylArraySheet_updSamples_Ordered_T1.rds"))
dim(sheet) #  170   8

# Reading clinical file, should be in current directory
# Altered location of clinical file on 10 Oct 2019 
# ClinicalFile2 <- read_delim(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/Data_Clinical/sample_annotations_rcd27Feb2019.txt", delim = '\t')
# ClinicalFile <- read.delim(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/Data_Clinical/sample_annotations_rcd27Feb2019.txt")
# ClinicalFile <- read.delim(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/Data_Clinical/sample_annotations_rcd16Oct2019.txt")
# ClinicalFile <- data.table::fread(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/Data_Clinical/sample_annotations_rcd16Oct2019.txt")
# The following file was generated from PipelineStep1
ClinicalFile_T1 <- readRDS(file = paste0(MethylationDirPath, "/1_ClinicalFile_updSamples_Ordered_T1.rds"))
dim(ClinicalFile_T1) # 170  25

BetaMatrix_T1 <- readRDS(file = paste0(MethylationDirPath, "/2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.rds"))
dim(BetaMatrix_T1) # 595564    170
range(BetaMatrix_T1) # 3.957731e-05 9.999547e-01
# write.csv(BetaMatrix_T1, "2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.csv")
# BetaMatrix_T12 <- as.matrix(read.csv("2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.csv", row.names = 1))
# dim(BetaMatrix_T12)
# range(BetaMatrix_T12) # 3.957731e-05 9.999547e-01


MvalueMatrix_T1 <- readRDS(file = paste0(MethylationDirPath, "/2_MvalueMatrix_updSamples_Ordered_T1_FilteredProbes.rds"))
dim(MvalueMatrix_T1) #   595564    170
range(MvalueMatrix_T1) # -7.93023  7.72350  
# write.csv(MvalueMatrix_T1, "2_MvalueMatrix_updSamples_Ordered_T1_FilteredProbes.csv", row.names = TRUE)
# MvalueMatrix_T12 <- data.table::fread("2_MvalueMatrix_updSamples_Ordered_T1_FilteredProbes.csv")
# dim(MvalueMatrix_T12[, -1]) # 595564    170
# range(MvalueMatrix_T12[, -1]) # -7.854857  7.662493

TumorPurity <- readRDS(file = paste0(MethylationDirPath, "/Purity_281probes_10Jan2020.rds"))
dim(TumorPurity) # 170   1
range(TumorPurity$purity) #  0.119587 0.816362
# write.csv(TumorPurity, "Purity_281probes_10Jan2020.csv")

library(InfiniumPurify)
data(iDMC) 
probes <- iDMC[["DLBC"]]
probes.true <- names(probes[probes == T])
beta.sel <- BetaMatrix_T1[row.names(BetaMatrix_T1) %in% probes.true,]
purity_OICR_withNormal <- getPurity(tumor.data = beta.sel, tumor.type = "DLBC")

# ClusterLabels2to4 <- readRDS(file = "InfiniumClustering2to4.rds")
# ClusterLabels2to4[[2]][[2]] # 2 cluster model 
# write.csv(data.frame(ID = colnames(MvalueMatrix_T1), Cluster = ClusterLabels2to4[[2]][[2]]), "InfiniumClustLabels2.csv")
InfiniumClustLabels <- readRDS(file = paste0(MethylationDirPath, "/InfiniumClustLabels2.RDS"))
dim(InfiniumClustLabels) # 170   2
range(InfiniumClustLabels$Cluster) # 1 2

EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020 <- readRDS(file = paste0(MethylationDirPath, "/EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020.rds"))
dim(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020) # 170 6
colnames(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020)
# "CD8T"  "CD4T"  "NK"    "Bcell" "Mono"  "Neu"
range(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020)
# -1.387779e-17  6.893340e-01
# write.csv(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020, file = "EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020.csv")



ENSMBLid <- read.csv(file = paste0(RNAseqDirPath, "/hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv"))
head(ENSMBLid$gene) #  ENSG00000000003 ENSG00000000005

# RNASeqCountMatrixMatched <- readRDS(file = paste0("RNASeqCountMatrix_Integrative.rds"))
# dim(RNASeqCountMatrixMatched) # 57820   132 No longer used
# As of 25 May 2020 don't use this
#RNAseqCountMatrix <- read.csv("2020-05-19_FL_136_samples_featureCounts_hg19_raw.csv", 
#                               row.names = 1)
# As of 18June2020  use this
# RNAseqCountMatrixMatched <- readRDS(file = paste0(RNAseqDirPath, "/2020-06-18STAR_quantmode_counts_matrix_FL_136_patients.RDS"))
# dim(RNAseqCountMatrixMatched) # 57820   132
# range(RNAseqCountMatrixMatched) # 0 9693321

RNAseqCountMatrixMatched <- data.frame(readRDS(file = paste0(RNAseqDirPath, "2020-06-18STAR_quantmode_counts_matrix_FL_132_patients.RDS")))
dim(RNAseqCountMatrixMatched) # 57820   132
range(RNAseqCountMatrixMatched) # 0 9693321
#RNAseqQCFile <- read.csv("FL_TGL_STAR_logQC_2020-05-13_summary_KI_ClusterContamAdded.csv",  
#                         row.names = 1) # removed
RNAseqQCFile <- readRDS(file = paste0(RNAseqDirPath, "/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.RDS"))
dim(RNAseqQCFile) # 132  60
range(RNAseqQCFile$Average.input.read.length) # 150 160


# RNAseqMutationCalls <- data.table::fread(file = paste0(RNAseqDirPath, "/2020-06-22_opossum_variant_FL_rna-seq_filtered.txt"))
# dim(RNAseqMutationCalls) # 150423     73 (previously # 99948    73)
# saveRDS(RNAseqMutationCalls, file = "2020-06-22_opossum_variant_FL_rna-seq_filtered.rds")
RNAseqMutationCalls <- readRDS(file = paste0(RNAseqDirPath, "/2020-06-22_opossum_variant_FL_rna-seq_filtered.rds"))
dim(RNAseqMutationCalls) # 150423     73
                               
                               
TargetSeqBC <- readRDS(file = paste0(TargetedDNAseqDirPath, "/BC_Cancer_capseq_data.rds"))
dim(TargetSeqBC) # 380  10
range(TargetSeqBC$Chromosome) # "1" "X"
BC_Cancer_capseq_data_with01ClusterLabs <- readRDS(file = paste0(TargetedDNAseqDirPath, "/BC_Cancer_capseq_data_with01ClusterLabs_AS_2June2020.rds"))
dim(BC_Cancer_capseq_data_with01ClusterLabs) # 31 67
range(BC_Cancer_capseq_data_with01ClusterLabs) # 0 2
################################## #






########################################## #
# Summarizing survival files
########################################## #
matchedIDs <- match(substr(colnames(BetaMatrix[, which(substr(colnames(BetaMatrix), 4, 5) == "FL")]), 1, 9), 
                    SurvivalFile$LY_FL_ID)
if (sum(is.na(matchedIDs) ) > 0) {
  # if NAs are present
  SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile[matchedIDs[- which(is.na(matchedIDs) == TRUE)], ]
} else {
  SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile[matchedIDs, ]
}
dim(SurvivalFile_OrderedbyBetaMatrixPatients) # 155  23

# AGE
median(SurvivalFile_OrderedbyBetaMatrixPatients$AGE_AT_DIAGNOSIS) # 58
IQR(SurvivalFile_OrderedbyBetaMatrixPatients$AGE_AT_DIAGNOSIS) # 17.5
range(SurvivalFile_OrderedbyBetaMatrixPatients$AGE_AT_DIAGNOSIS) # 20 86
quantile(SurvivalFile_OrderedbyBetaMatrixPatients$AGE_AT_DIAGNOSIS)
#  0%  25%  50%  75% 100% 
# 20.0 49.5 58.0 67.0 86.0 

# Sex
(table(SurvivalFile_OrderedbyBetaMatrixPatients$SEX)/155) * 100
# F  M 
# 80 75 

# Translocation
(table(SurvivalFile_OrderedbyBetaMatrixPatients$T_14_18)/155)*100
length(SurvivalFile_OrderedbyBetaMatrixPatients$T_14_18)

# FLIPI binary
(table(SurvivalFile_OrderedbyBetaMatrixPatients$FLIPI_BINARY)/155)*100
length(SurvivalFile_OrderedbyBetaMatrixPatients$T_14_18)

# FLIPI binary
(table(ClinicalFile_T1$STAGE)/155)*100
length(SurvivalFile_OrderedbyBetaMatrixPatients$T_14_18)

########################################## #
# Running functions
########################################## #
# ClinicalFile16Oct2019 <- data.table::fread(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/Data_Clinical/sample_annotations_rcd16Oct2019.txt")
#  PipelineStep1 <- QCRemoveSamples(ClinicalFile = as.data.frame(ClinicalFile_T1))
#  saveRDS(PipelineStep1$ClinicalFile_updSamples_170_Ordered_T1, file = "1_ClinicalFile_updSamples_Ordered_T1.rds")
#  saveRDS(PipelineStep1$ClinicalFile_updSamples_170_Ordered_T2, file = "1_ClinicalFile_updSamples_Ordered_T2.rds")
#  saveRDS(PipelineStep1$Sheet_updSamples_170_Ordered_T1, file = "1_MethylArraySheet_updSamples_Ordered_T1.rds")
#  saveRDS(PipelineStep1$Sheet_updSamples_170_Ordered_T2, file = "1_MethylArraySheet_updSamples_Ordered_T2.rds")

#  write.csv(PipelineStep1$ClinicalFile_updSamples_170_Ordered_T1, file = "1_ClinicalFile_updSamples_Ordered_T1.csv")
#  write.csv(PipelineStep1$ClinicalFile_updSamples_170_Ordered_T2, file = "1_ClinicalFile_updSamples_Ordered_T2.csv")
#  write.csv(PipelineStep1$Sheet_updSamples_170_Ordered_T1, file = "1_MethylArraySheet_updSamples_Ordered_T1.csv")
#  write.csv(PipelineStep1$Sheet_updSamples_170_Ordered_T2, file = "1_MethylArraySheet_updSamples_Ordered_T2.csv")


  PipelineStep2 <- LoadingMethylData(Sheet = sheet, 
                                     ClinicalFile = ClinicalFile_T1, 
                                     AnnotationFile = AnnotationFile, 
                                     FigureGenerate = "No", 
                                     TableGenerate = "No", 
                                     PNGorPDF = "png")
# saveRDS(PipelineStep2$RGChannelSet, file = "2_RGChannelSet_25June2020.rds")
# saveRDS(PipelineStep2, file = "2_PipelineStep2_Output_25June2020.rds")
# saveRDS(PipelineStep2$BetaMatrix, file = "2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.rds")
# write.csv(PipelineStep2$BetaMatrix, file = "2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.csv") # 719990 170
# saveRDS(Output_Data_2$MvalueMatrix, file = "2_MvalueMatrix_updSamples_Ordered_T1_FilteredProbes.rds")
# save.image("Output_Data_2_23Jan2020.RData")
# table(Output_Data_2$PredictedSex[which(ClinicalFile_T1$TYPE == "FL")], ClinicalFile_T1$SEX[which(ClinicalFile_T1$TYPE == "FL")]) 
#   F  M
# F 73  0
# M  7 75

# which(Output_Data_2$PredictedSex[which(ClinicalFile_T1$TYPE == "FL")] != ClinicalFile_T1$SEX[which(ClinicalFile_T1$TYPE == "FL")])
# 10 21 42 47 48 50 53
# Output_Data_2$PredictedSex[which(ClinicalFile_T1$TYPE == "FL")][c(10, 21, 42, 47, 48, 50, 53)]
# ClinicalFile_T1$SAMPLE_ID[which(ClinicalFile_T1$TYPE == "FL")[c(10, 21, 42, 47, 48, 50, 53)]]
# ClinicalFile_T1$TYPE_BIOPSY[which(ClinicalFile_T1$TYPE == "FL")[c(10, 21, 42, 47, 48, 50, 53)]]
# ClinicalFile_T1$INSTITUTION[which(ClinicalFile_T1$TYPE == "FL")[c(10, 21, 42, 47, 48, 50, 53)]]
# ClinicalFile_T1$STAGE[which(ClinicalFile_T1$TYPE == "FL")[c(10, 21, 42, 47, 48, 50, 53)]]


# BetaMatrix_T1 <- readRDS(file = "2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.rds")
# dim(BetaMatrix_T1) # 595564    170
# MvalueMatrix_T1 <- readRDS(file = "2_MvalueMatrix_updSamples_Ordered_T1_FilteredProbes.rds")
# dim(MvalueMatrix_T1) # 595564    170
# TumorPurity <- readRDS(file = "Purity_281probes_10Jan2020.rds")
# dim(TumorPurity) # 170   1
# RNASeqCountMatrixMatched <- readRDS(file = paste0("RNASeqCountMatrix_Integrative.rds"))
# dim(RNASeqCountMatrixMatched) # 57820   132


# Output_Visuals_3_TYPE <- MethylationDensityPlot(AnnotationCategoryToVisualize="Relation_to_Island", PlotWithinAnnotationCategories="Yes", ClinicalCategoryToVisualize="TYPE", BetaMatrix=Output_Remove_2$BetaMatrix_updSamples, AnnotationFile=AnnotationFile, ClinicalFile=Output_Remove_2$ClinicalFile_updSamples, SampleSheet=Output_Data_1$SampleSheet, FigureGenerate="Yes",  ImageName = "BetaMatrix_updSamples_rmdSexChsomes_TYPE", PNGorPDF="png")
# Output_Visuals_3_TYPE2 <- MethylationDensityPlot(AnnotationCategoryToVisualize="None", PlotWithinAnnotationCategories="No", ClinicalCategoryToVisualize="TYPE", BetaMatrix=Output_Remove_2$BetaMatrix_updSamples, AnnotationFile=AnnotationFile, ClinicalFile=Output_Remove_2$ClinicalFile_updSamples, SampleSheet=Output_Data_1$SampleSheet, FigureGenerate="Yes", PNGorPDF="png")
# Output_Visuals_3_STAGE <- MethylationDensityPlot(AnnotationCategoryToVisualize="Relation_to_Island", PlotWithinAnnotationCategories="Yes", ClinicalCategoryToVisualize="STAGE", BetaMatrix=Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, AnnotationFile=AnnotationFile, ClinicalFile=Output_Remove_2$ClinicalFile_updSamples, SampleSheet=Output_Data_1$SampleSheet, FigureGenerate="Yes", PNGorPDF="png")
# Output_Visuals_3_STAGE2 <- MethylationDensityPlot(AnnotationCategoryToVisualize="None", PlotWithinAnnotationCategories="No", ClinicalCategoryToVisualize="STAGE", BetaMatrix=Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, AnnotationFile=AnnotationFile, ClinicalFile=Output_Remove_2$ClinicalFile_updSamples, SampleSheet=Output_Data_1$SampleSheet, FigureGenerate="Yes", PNGorPDF="png")
# *** Not done  Output_SummarizedMvalues_4 <- SummarizationMvalues(MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ClinicalFile=Output_Remove_2$ClinicalFile_updSamples)
# *** Not done  Output_SummarizedBetavalues_5 <- SummarizationBeta(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples, ClinicalFile=Output_Remove_2$ClinicalFile_updSamples)
# Output_Summary_6 <- SummaryStatistics(BetaMatrixProbes=Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, BetaMatrixGenes=NA, ClinicalCategoryToVisualize=NA, ClinicalFile=NA, MvalueMatrixProbes=Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, MvalueMatrixGenes=NA, ProduceImages="Yes", PNGorPDF="png")
# Output_Differential_probes_TYPE_7 <- DifferentialMethylation(ClinicalFile=Output_Remove_2$ClinicalFile_updSamples, MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ProbeOrGene = "Probe", ContrastColumnName="TYPE", ProduceImages="No", PNGorPDF="png")
# *** Not done  Output_Differential_genes_TYPE_7  <- DifferentialMethylation(ClinicalFile=Output_Remove_2$ClinicalFile_updSamples, MvalueMatrix=Output_SummarizedMvalues_4$SummarizedMvalueMatrix, ProbeOrGene = "Gene", ContrastColumnName="TYPE", ProduceImages="Yes", PNGorPDF="png")
# Output_Differential_probes_STAGE_7 <- DifferentialMethylation(ClinicalFile = Output_Remove_2_ClinicalFile, 
#                                                               MvalueMatrix = Output_Remove_2_MvalueMatrix, 
#                                                               BetaMatrix = Output_Remove_2_BetaMatrix,
#                                                               ProbeOrGene = "Probe", 
#                                                               ContrastColumnName = "STAGE", 
#                                                               RGChannelSet = Output_Data_1$RGChannelSet, 
#                                                               SampleSheet = Output_Data_1$SampleSheet, 
#                                                               ProduceImages = "No", PNGorPDF = "png")
# *** Not done  Output_Differential_genes_STAGE_7 <- DifferentialMethylation(ClinicalFile=Output_Remove_2$ClinicalFile_updSamples, MvalueMatrix=Output_SummarizedMvalues_4$SummarizedMvalueMatrix, ProbeOrGene = "Gene", ContrastColumnName="STAGE", PNGorPDF="png")

# Output_Differential_Variability_probes_STAGE_8 <-DifferentialVariability(ClinicalFile=Output_Remove_2$ClinicalFile_updSamples, MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ProbeOrGene="Probe", ContrastColumnName="STAGE", PNGorPDF="png")
# *** Not done Output_Differential_Variability_genes_STAGE_8 <-DifferentialVariability(ClinicalFile=Output_Remove_2$ClinicalFile_updSamples, MvalueMatrix=Output_SummarizedMvalues_4$SummarizedMvalueMatrix, ProbeOrGene="Gene", ContrastColumnName="STAGE", PNGorPDF="png")
#  Output_Visuals_Relation_to_Island_STAGE_2 <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, AnnotationFile = AnnotationFile, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, ClinicalCategory = "STAGE", PlotWithinCategories="Yes", FigureGenerate="Yes", PNGorPDF="png")
#  Output_Visuals_Relation_to_Island_TYPE_21 <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, AnnotationFile = AnnotationFile, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, ClinicalCategory = "TYPE", PlotWithinCategories="No", FigureGenerate="Yes", PNGorPDF="png")
# *** Not done Output_Visuals_Chromosome_21 <- ProportionVisualization(CategoryToVisualize = "chr", BetaMatrix = Output_Remove_2$BetaMatrix, AnnotationFile = AnnotationFile, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, ClinicalCategory = "TYPE", PlotWithinCategories="Yes", FigureGenerate = "Yes", PNGorPDF="png")
# *** Not done Output_Visuals_Differentially_Methylated_Regions_9 <- ProportionVisualization(CategoryToVisualize = "DMR", BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, AnnotationFile = AnnotationFile, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, ClinicalCategory = "TYPE", PlotWithinCategories="Yes", FigureGenerate="Yes", PNGorPDF="png")
#  Output_Visuals_Type_10 <- ViolinPlot(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, CategoryToVisualize="TYPE", PNGorPDF="png")
#  Output_Visuals_Stage_10 <- ViolinPlot(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, CategoryToVisualize="STAGE", PNGorPDF="png")
#  Output_Visuals_Type_11 <- HeatPlot(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, CategoryToVisualize = "TYPE", PNGorPDF="png")
#  Output_Visuals_Stage_11 <- HeatPlot(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, CategoryToVisualize = "STAGE", PNGorPDF="png")
# *** No done  Output_Visuals_Stage_12 <- BoxPlot(BetaMatrix = Output_SummarizedBetavalues_5$SummarizedBetaMatrix, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, CategoryToVisualize = "STAGE", GeneOrProbeName="APOC1", PNGorPDF="png")
# Output_tSNE_13 <- tSNEPlotGeneration(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, PerplexityParameter=6, PNGorPDF="png")
# Output_Clustering_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = Output_Differential_probes_TYPE_7$MultipleComparisonCommonProbes$ListPValues, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="DiffMethylation" )
# *** No done Output_SurvivalAnalysis_15 <- SurvivalAnalysis ()
# *** No done Output_Huet23GeneModel_16 <- Huet23GeneModel(BetaMatrixSummarized = Output_SummarizedBetavalues_5$SummarizedBetaMatrix, BetaMatrixNotSummarized = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrixSummarized = Output_SummarizedMvalues_4$SummarizedMvalueMatrix, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png")
# *** No done Output_LinePlots_17 <- LinePlot(dataset = Output_Clustering_14$BetaMatrixProbesUsedClustering, ClusterMembershipVector = Output_Clustering_14$RPMM$Labels, name = "RPMM")

#*** No done Output_SDeviation_0.3_18 <- SDeviation(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples, standard_deviation=0.3)
#*** No done Output_Clustering_SD_0.3_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_SDeviation_0.3_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.3SD")
#*** No done cat("\n Running SD of 0.25")
#*** No done Output_SDeviation_0.25_18 <- SDeviation(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples, standard_deviation=0.25)
#*** No done cat("\n Clustering on SD of 0.25")
#*** No done Output_Clustering_SD_0.25_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_SDeviation_0.25_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.25SD")






# cat("\n Running SD of 0.25 with no sex xsomes of all samples = 607 probes")
#  Output_SDeviation_0.25SD_noSexChsomes_18 <- SDeviation(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, standard_deviation=0.25)
# cat("\n Clustering on SD of 0.25 with no sex xsomes")
#  Output_Clustering_SD_0.25_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ListofProbes = rownames(Output_SDeviation_0.25SD_noSexChsomes_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="No", PNGorPDF="png", ImageName="0.25SDNoSexChsomes")
# Saving labels
# Probes607_3ClusterModel_RPMM_Labels <- cbind(colnames(Output_Clustering_SD_0.25_Output_14$BetaMatrixProbesUsedClustering), Output_Clustering_SD_0.25_Output_14$RPMM$RPMMoutputLabels)
# write.csv(Probes607_3ClusterModel_RPMM_Labels, file="607Probes_3ClusterModel_RPMM_Labels.csv")
# Output_tSNE_SDeviation_0.25SD_noSexChsomes_Allsamples_Allprobes_13 <- tSNEPlotGeneration(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, PerplexityParameter=6, ClusterLabels=Output_Clustering_SD_0.25_Output_14$RPMM$RPMMoutputLabels, FigureGenerate="Yes", PNGorPDF="png", ImageName="NoSexChsomes_Allsamples_AllProbes")
# Output_tSNE_SDeviation_0.25SD_noSexChsomes_Allsamples_13 <- tSNEPlotGeneration(BetaMatrix=Output_SDeviation_0.25SD_noSexChsomes_18$BetaMatrix_SD_Filtered, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, PerplexityParameter=6, ClusterLabels=Output_Clustering_SD_0.25_Output_14$RPMM$RPMMoutputLabels, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.25SD_NoSexChsomes_Allsamples")
# Output_Survival_SDeviation_0.25SD_noSexChsomes_Allsamples_15 <- SurvivalAnalysis(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, ClinicalFile=Output_Remove_2$ClinicalFile_updSamples, SurvivalFile=SurvivalFile, ClusterLabels=Output_Clustering_SD_0.25_Output_14$RPMM$RPMMoutputLabels, FigureGenerate="Yes", PNGorPDF="png")
# Output_SummarizationBeta_SD_0.25_Output_5 <- SummarizationBeta(BetaMatrix = Output_SDeviation_0.25SD_noSexChsomes_18$BetaMatrix_SD_Filtered, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples)
# Output_Density_Clustering_SD_0.25_Allsamples_3 <- MethylationDensityPlot(ClusterLabels = Output_Clustering_SD_0.25_Output_14$RPMM$RPMMoutputLabels, BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, FigureGenerate = "Yes", PNGorPDF = "png")
# Output_ProportionVisualization_Clustering_SD_0.25_Allsamples_9 <- ProportionVisualization(CategoryToVisualize=NA, ClusterLabels = Output_Clustering_SD_0.25_Output_14$RPMM$RPMMoutputLabels, BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, AnnotationFile = NA, ClinicalFile = NA, ClinicalCategory = NA, PlotWithinCategories = NA, FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "3RPMMmodel")




#*** No done Output_SDeviation_0.3_18 <- SDeviation(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, standard_deviation=0.3)
#*** No done Output_Clustering_SD_0.3_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ListofProbes = rownames(Output_SDeviation_0.3_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.3SD")
#*** No done Output_SDeviation_0.4_18 <- SDeviation(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples, standard_deviation=0.4)
#*** No done Output_Clustering_SD_0.4_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_SDeviation_0.4_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.4SD")

#*** No done Output_SDeviation_18 <- SDeviation(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples, standard_deviation=0.2)
#*** No done Output_Clustering_SD_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_SDeviation_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="pdf", ImageName="0.2SD")
#*** No done Output_SDeviation_0.3_18 <- SDeviation(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples, standard_deviation=0.3)
#*** No done Output_Clustering_SD_0.3_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_SDeviation_0.3_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.3SD")



#*** No done Output_IsolateEntries_RelationToIsland_19 <- IsolateEntries(CategoryToVisualize = "Relation_to_Island", BetaMatrix = Output_Remove_2$BetaMatrix, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, AnnotationFile = AnnotationFile, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, ProvideWithinCategories="Yes")
#*** No done Output_Clustering_Island_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_IsolateEntries_RelationToIsland_19$Matrix[[1]]$BetaMatrix_AnnotationCategory)[c(1:100)], ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="Island")
#*** No done Output_Clustering_N_Shore_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_IsolateEntries_RelationToIsland_19$Matrix[[2]]$BetaMatrix_AnnotationCategory)[c(1:100)], ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="N_Shore")
#*** No done Output_Clustering_OpenSea_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_IsolateEntries_RelationToIsland_19$Matrix[[3]]$BetaMatrix_AnnotationCategory)[c(1:100)], ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="OpenSea")
#*** No done Output_Clustering_N_Shelf_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_IsolateEntries_RelationToIsland_19$Matrix[[4]]$BetaMatrix_AnnotationCategory)[c(1:100)], ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="N_Shelf")
#*** No done Output_Clustering_S_Shelf_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_IsolateEntries_RelationToIsland_19$Matrix[[5]]$BetaMatrix_AnnotationCategory)[c(1:100)], ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="S_Shelf")
#*** No done Output_Clustering_S_Shore_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_IsolateEntries_RelationToIsland_19$Matrix[[6]]$BetaMatrix_AnnotationCategory)[c(1:100)], ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="S_Shore")


# Output_Glmnet_20 <- Glmnet_FeatureSelection(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, FigureGenerate="Yes", PNGorPDF="png")
# dim(Output_Glmnet_20$BetaMatrix_Lambda.min) # 29 x 171
# Output_Clustering_Glmnet_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_Glmnet_20$BetaMatrix_Lambda.min), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="Glmnet")


## Isolate Probes that belong to Island and has an SD of 0.2
# Originally 117461 probes
# dim(Output_IsolateEntries_RelationToIsland_19$Matrix[[1]]$BetaMatrix_AnnotationCategory)
# Brought down to 4871 probes

# cat("\n Running SD0.25 Island")
# Output_SDeviation_0.25_Island_18 <- SDeviation(BetaMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[1]]$BetaMatrix_AnnotationCategory, MvalueMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[1]]$MvalueMatrix_AnnotationCategory, standard_deviation=0.25)
#  cat("\n Running Clustering after SD0.25 Island")
#  Output_Clustering_SD_0.25_Island_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_SDeviation_0.25_Island_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.25SD_Island")

# Look into how to isolate for islands??
############3 EDITTTTT #########
# cat("\n EDITTTTT")
# cat("\n Running SD of 0.25 Island with no sex xsomes")
# Output_SDeviation_0.25SD_Island_noSexChsomes_18 <- SDeviation(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, standard_deviation=0.25)
# cat("\n Clustering on SD of 0.25 Island  with no sex xsomes")
# Output_Clustering_SD_0.25_Island_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ListofProbes = rownames(Output_SDeviation_0.25SD_noSexChsomes_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.25SDNoSexChsomes")


# Output_SDeviation_0.2_Island_18 <- SDeviation(BetaMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[1]]$BetaMatrix_AnnotationCategory, MvalueMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[1]]$MvalueMatrix_AnnotationCategory, standard_deviation=0.2)
# Output_Clustering_SD_0.2_Island_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_SDeviation_0.2_Island_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.2SD_Island")


## Isolate Probes that belong to Island and has an SD of 0.3
# Originally 117461 probes
# Brought down to 108 probes
# Output_SDeviation_0.3_Island_18 <- SDeviation(BetaMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[1]]$BetaMatrix_AnnotationCategory, MvalueMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[1]]$MvalueMatrix_AnnotationCategory, standard_deviation=0.3)
# Output_Clustering_SD_0.3_Island_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_SDeviation_0.3_Island_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.3SD_Island")

## Isolate Probes that belong to N_Shore and has an SD of 0.2
# Originally 54361 probes
# dim(Output_IsolateEntries_RelationToIsland_19$Matrix[[2]]$BetaMatrix_AnnotationCategory)
# Brought down to 171 probes

# Output_SDeviation_0.25_N_Shore_18 <- SDeviation(BetaMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[2]]$BetaMatrix_AnnotationCategory, MvalueMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[2]]$MvalueMatrix_AnnotationCategory, standard_deviation=0.25)
# Output_Clustering_SD_0.25_N_Shore_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_SDeviation_0.25_N_Shore_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.25N_Shore")

# Output_SDeviation_0.2_N_Shore_18 <- SDeviation(BetaMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[2]]$BetaMatrix_AnnotationCategory, MvalueMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[2]]$MvalueMatrix_AnnotationCategory, standard_deviation=0.2)
# Output_Clustering_SD_0.2_N_Shore_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_SDeviation_0.2_N_Shore_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.2N_Shore")


## Isolate Probes that belong to N_Shore and has an SD of 0.3
# Originally 54361 probes
# Brought down to 171 probes
# Output_SDeviation_0.3_N_Shore_18 <- SDeviation(BetaMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[2]]$BetaMatrix_AnnotationCategory, MvalueMatrix=Output_IsolateEntries_RelationToIsland_19$Matrix[[2]]$MvalueMatrix_AnnotationCategory, standard_deviation=0.3)
# Output_Clustering_SD_0.3_N_Shore_Output_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_SDeviation_0.3_N_Shore_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.3N_Shore")

# Output_Glmnet_20 <- Glmnet_FeatureSelection(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, FigureGenerate="Yes", PNGorPDF="png")
# Output_Clustering_Glmnet_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ListofProbes = rownames(Output_Glmnet_20$BetaMatrix_Lambda.min), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, FigureGenerate="Yes", PNGorPDF="png", ImageName="Glmnet")

# Diff_MethylatedRegions(Method = "DMRcate", BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples, ContrastColumnName="TYPE", ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, AnnotationFile = AnnotationFile, ProduceImages="Yes", PNGorPDF="png")

# MDSplots(ClinicalCategoryToVisualize="TYPE", BetaMatrix=Output_Remove_2$BetaMatrix_updSamples, AnnotationFile=AnnotationFile, ClinicalFile=Output_Remove_2$ClinicalFile_updSamples, SampleSheet=Output_Data_1$SampleSheet, FigureGenerate="Yes", PNGorPDF="png")

# cat("\n Running Diff_MethylatedRegions")
#Diff_MethylatedRegions_0.25SD <- Diff_MethylatedRegions(Method = "DMRcate", BetaMatrix = Output_SDeviation_0.25$BetaMatrix_SD_Filtered, MvalueMatrix = Output_SDeviation_0.25$MvalueMatrix_SD_Filtered, ContrastColumnName="TYPE", ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, AnnotationFile = AnnotationFile, ProduceImages="No", DMR=1, PNGorPDF="pdf")

#Diff_MethylatedRegions_Type_All <- Diff_MethylatedRegions(Method = "DMRcate", BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ContrastColumnName="TYPE", ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, AnnotationFile = AnnotationFile, ProduceImages="Yes", DMR=1, PNGorPDF="png")

#  Diff_MethylatedRegions_STAGE_rmdSexChsomes <- Diff_MethylatedRegions(Method = "DMRcate", BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ContrastColumnName="STAGE", ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, AnnotationFile = AnnotationFile, ProduceImages="No", DMR=1, PNGorPDF="png")
#  Diff_MethylatedRegions_TYPE_rmdSexChsomes <- Diff_MethylatedRegions(Method = "DMRcate", BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ContrastColumnName="TYPE", ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, AnnotationFile = AnnotationFile, ProduceImages="No", DMR=1, PNGorPDF="png")


# Output_Visuals_Regulatory_Feature_Group_9 <- ProportionVisualization(CategoryToVisualize = "Regulatory_Feature_Group", BetaMatrix = Output_Remove_2$BetaMatrix, AnnotationFile = AnnotationFile, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, ClinicalCategory = "TYPE", PlotWithinCategories="Yes", FigureGenerate = "Yes", PNGorPDF="png")
# Output_Visuals_Relation_to_Island_TYPE_9 <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", BetaMatrix = Output_Remove_2$BetaMatrix_updSamples, AnnotationFile = AnnotationFile, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, ClinicalCategory = "TYPE", PlotWithinCategories="No", FigureGenerate="Yes", PNGorPDF="png")

# Diff_MethylatedRegions_CLUSTER_AllProbes <- Diff_MethylatedRegions(Method = "DMRcate", BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes, MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes, ContrastColumnName = "CLUSTER", ClusterLabels = Output_Clustering_SD_0.25_Output_14$RPMM$RPMMoutputLabels, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, AnnotationFile = AnnotationFile, ProduceImages="Yes", DMR=1, PNGorPDF="png")
# table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast <- cbind(Diff_MethylatedRegions_CLUSTER_AllProbes$DifferentiallyMethylatedRegions[[1]]$results, Diff_MethylatedRegions_CLUSTER_AllProbes$GRangesObject[[1]])
# table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast <- cbind(Diff_MethylatedRegions_CLUSTER_AllProbes$DifferentiallyMethylatedRegions[[2]]$results, Diff_MethylatedRegions_CLUSTER_AllProbes$GRangesObject[[2]])
# table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast <- cbind(Diff_MethylatedRegions_CLUSTER_AllProbes$DifferentiallyMethylatedRegions[[3]]$results, Diff_MethylatedRegions_CLUSTER_AllProbes$GRangesObject[[3]])
#write.csv(table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast, file="Table_Diff_MethylatedRegions_CLUSTER_AllProbes_C2vsC3.csv")

# Function extract gene regions list from Diff_MethylatedRegions
ExtractingElements <- function(x, FILE) {
  
  example <- unlist(strsplit(as.character(FILE$overlapping.promoters[x]), "[,]"))
  unique_elements <- unique(trimws(sub("\\-.*", "", example)))
  
  GettingTable <- function(i, n, numb_unique_elements){
    element_matrix <- matrix(ncol=2, nrow=numb_unique_elements)
    element_matrix[i,1] <- unique_elements[i]
    element_matrix[i,2] <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast$Stouffer[n]
    return(element_matrix)
  }
  
  element_matrix <- GettingTable(i=c(1:length(unique_elements)), n=x, numb_unique_elements=length(unique_elements))
  return(element_matrix)
}

# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast <- sapply(c(1:500), function(i) ExtractingElements(i, FILE=table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast))
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast, quote=FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast, file="GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast.csv")
# # Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast <- sapply(c(1:500), function(i) ExtractingElements(i, FILE=table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast))
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast, quote=FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast, file="GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast.csv")
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast <- sapply(c(1:500), function(i) ExtractingElements(i, FILE=table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast))
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast, quote=FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast, file="GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast.csv")







# cat("\n Running SD of 0.25 with no sex xsomes, FL samples only =  585 probes")
# Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) == "FL")], standard_deviation = 0.25)
# Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], ListofProbes = rownames(Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL"),], FigureGenerate="No", PNGorPDF="png", ImageName="0.25SDNoSexChsomes_FLonly")
# Saving labels
# Probes585_2ClusterModel_RPMM_Labels<-cbind(colnames(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$BetaMatrixProbesUsedClustering), Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels)
# write.csv(Probes585_2ClusterModel_RPMM_Labels, file="585Probes_2ClusterModel_RPMM_Labels.csv")
# Output_tSNE_SDeviation_0.25SD_noSexChsomes_FLonlysamples_13 <- tSNEPlotGeneration(BetaMatrix=Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL"),], PerplexityParameter=6, ClusterLabels=Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels, FigureGenerate="Yes", PNGorPDF="png", ImageName="0.25SD_NoSexChsomes_FLonly")
# Output_tSNE_SDeviation_0.25SD_noSexChsomes_FLonlysamples_Allprobes_13 <- tSNEPlotGeneration(BetaMatrix=Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL"),], PerplexityParameter=6, ClusterLabels=Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels, FigureGenerate="Yes", PNGorPDF="png", ImageName="_NoSexChsomes_FLonly_AllProbes")
# Output_SummarizationBeta_SDeviation_0.25SD_noSexChsomes_FLonlysamples_5 <- SummarizationBeta(Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples)
# Output_SummarizationBeta_SDeviation_0.25SD_noSexChsomes_FLonlysamples_15 <- SurvivalAnalysis(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL")],
#                  ClinicalFile = Output_Remove_2$ClinicalFile_updSamples,
#                  ClusterLabels = Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels,
#                  SurvivalFile = SurvivalFile, FigureGenerate = "No",
#                  PNGorPDF = "png")
# Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_3 <- MethylationDensityPlot(ClusterLabels = Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels, BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], FigureGenerate = "Yes", ImageName = "2RPMMmodel", PNGorPDF = "png")

# Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM <- Diff_MethylatedRegions(Method = "DMRcate", BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], ContrastColumnName = "CLUSTER", ClusterLabels = Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, AnnotationFile = AnnotationFile, ProduceImages="No", DMR=1, PNGorPDF="png")
# table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_2RPMM <- cbind(Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM$DifferentiallyMethylatedRegions[[1]]$results, Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM$GRangesObject[[1]])
# write.csv(table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_2RPMM, file="Table_Diff_MethylatedRegions_CLUSTER_AllProbes_C1vsC2_2RPMM.csv")

# Function extract gene regions list
ExtractingElements <-function(x, FILE) {
  
  example <- unlist(strsplit(as.character(FILE$overlapping.promoters[x]), "[,]"))
  unique_elements <- unique(trimws(sub("\\-.*", "", example)))
  
  GettingTable <- function(i, n, numb_unique_elements){
    element_matrix <- matrix(ncol=2, nrow=numb_unique_elements)
    element_matrix[i,1] <- unique_elements[i]
    element_matrix[i,2] <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_2RPMM$Stouffer[n]
    return(element_matrix)
  }
  
  element_matrix <- GettingTable(i=c(1:length(unique_elements)), n=x, numb_unique_elements=length(unique_elements))
  return(element_matrix)
}
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_2RPMM <- sapply(c(1:500), function(i) ExtractingElements(i, FILE=table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_2RPMM))
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_2RPMM <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_2RPMM , quote=FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_2RPMM, file="GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_2RPMM .csv")

# Preping for proportion visualization for only 585 probes from 2 RPMM model
# matchedIDs <- match(ListofProbes, rownames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes))
# BetaMatrix_OrderListofProbes <- Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[matchedIDs, which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")]
# Output_ProportionVisualization_Clustering_SD_0.25_FLOnly_21_2RPMM <- ProportionVisualization(CategoryToVisualize=NA, ClusterLabels = Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels, BetaMatrix = BetaMatrix_OrderListofProbes, AnnotationFile=NA, ClinicalFile=NA, ClinicalCategory=NA, PlotWithinCategories=NA, FigureGenerate="Yes", PNGorPDF="png", ImageName="2RPMMmodel")

# Preping for proportion visualization for only ALL probes from 2 RPMM model
# Output_ProportionVisualization_Clustering_SD_0.25_FLOnly_Allprobes_21_2RPMM <- ProportionVisualization(CategoryToVisualize=NA, ClusterLabels = Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels, BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[, which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], AnnotationFile=NA, ClinicalFile=NA, ClinicalCategory=NA, PlotWithinCategories=NA, FigureGenerate="Yes", PNGorPDF="png", ImageName="2RPMMmodel_Allprobes")


# ########################## Looking at overall methylation in Cluster 1 and Cluster 2 of 2 RPMM model (585 probes)
# BetaMatrix_2RPMMmodel <- Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL")]
# MvalueMatrix_2RPMMmodel <- Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL")]
# ListofProbes_2RPMMmodel = rownames(Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered)
# From DifferentialExpressionResults select top probes (arranged by p-value) with p values less than 0.01
# matchedIDs <- match(ListofProbes, rownames(BetaMatrix_2RPMMmodel))
# BetaMatrix_OrderListofProbes_2RPMMmodel <- BetaMatrix_2RPMMmodel[matchedIDs,]

# dim(BetaMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])
# mean(BetaMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])
# median(BetaMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])
# range(BetaMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])

# dim(BetaMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])
# mean(BetaMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])
# median(BetaMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])
# range(BetaMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])

# matchedIDs <- match(ListofProbes, rownames(MvalueMatrix_2RPMMmodel))
# MvalueMatrix_OrderListofProbes_2RPMMmodel <- MvalueMatrix_2RPMMmodel[matchedIDs,]

# dim(MvalueMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])
# mean(MvalueMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])
# median(MvalueMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])
# range(MvalueMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])

# dim(MvalueMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])
# mean(MvalueMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])
# median(MvalueMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])
# range(MvalueMatrix_OrderListofProbes_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])

# ########################## Looking at overall methylation in Cluster 1 and Cluster 2 of 2 RPMM model (all probes)
# mean(BetaMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])
# median(BetaMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])
# range(BetaMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])

# mean(BetaMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])
# median(BetaMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])
# range(BetaMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])

# mean(MvalueMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])
# median(MvalueMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])
# range(MvalueMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 1)])

# mean(MvalueMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])
# median(MvalueMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])
# range(MvalueMatrix_2RPMMmodel[,which(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels == 2)])



# Comparing 3 and 2 RPMM clusters with 2 cluster model from SNF
#  SNFClusters <- read.csv(file = "SNF_cluster_labels.csv")
#  RPMM2model <- read.csv(file = "585Probes_2ClusterModel_RPMM_Labels.csv")
#  RPMM3model <- read.csv(file = "607Probes_3ClusterModel_RPMM_Labels.csv")

# Match samples from SNF with 2 RPMM model
#  match_RPMM2_bySNF <- match(SNFClusters$patients, substr(RPMM2model$V1, 1, 9))
# table(SNFClusters$cluster,RPMM2model$V2[match_RPMM2_bySNF])
# chisq.test(table(SNFClusters$cluster, RPMM2model$V2[match_RPMM2_bySNF]))

#  Output_tSNE_13_SNFClusters_Allprobes <- tSNEPlotGeneration(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,match_RPMM2_bySNF], PerplexityParameter = 6, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[match_RPMM2_bySNF,], ClusterLabels=SNFClusters$cluster, Filter=NA, FigureGenerate="Yes", PNGorPDF="png", ImageName="30FLsamples_SNFclustering")


# Match samples from SNF with 3 RPMM model
# match_RPMM3_bySNF <- match(SNFClusters$patients, substr(RPMM3model$V1, 1, 9))
# table(SNFClusters$cluster, RPMM3model$V2[match_RPMM3_bySNF])

# Trying patients from SNF in RPMM
# cat("\n Running SD of 0.25 with no sex xsomes for 30 samples used in SNF")
# match_RPMM2_bySNF <- match(SNFClusters$patients, substr(RPMM2model$V1, 1, 9))
# Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_SNF_18 <- SDeviation(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[ , match_RPMM2_bySNF], MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[ , match_RPMM2_bySNF], standard_deviation=0.25)
# Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_SNF_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[ , match_RPMM2_bySNF], MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[ , match_RPMM2_bySNF], ListofProbes = rownames(Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_SNF_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[match_RPMM2_bySNF, ], FigureGenerate="No", PNGorPDF="png", ImageName="0.25SDNoSexChsomes_FLonly_SNF30samples")
# Compare RPMM results of 30 samples with SNF original results
# table(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_SNF_14$RPMM$RPMMoutputLabels, SNFClusters$cluster)
# Compare RPMM results of 30 samples with original RPMM 2 cluster results (30 cases separated)
# table(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_SNF_14$RPMM$RPMMoutputLabels, RPMM2model$V2[match_RPMM2_bySNF])
#  Output_tSNE_13_SNFClusters_Allprobes <- tSNEPlotGeneration(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,match_RPMM2_bySNF], PerplexityParameter = 6, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[match_RPMM2_bySNF,], ClusterLabels = Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_SNF_14$RPMM$RPMMoutputLabels, Filter=NA, FigureGenerate="Yes", PNGorPDF="png", ImageName="30FLsamples_RPMMClustering_Allprobes")
#  Output_tSNE_13_SNFClusters_0.25SD <- tSNEPlotGeneration(BetaMatrix = Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_SNF_18$BetaMatrix_SD_Filtered[,match_RPMM2_bySNF], PerplexityParameter = 6, ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[match_RPMM2_bySNF,], ClusterLabels = Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_SNF_14$RPMM$RPMMoutputLabels, Filter=NA, FigureGenerate="Yes", PNGorPDF="png", ImageName="30FLsamples_RPMMClustering")


# cat("\n Running all probes with no sex xsomes for 30 samples used in SNF") = don't work; error "vector memory exhausted (limit reached?)"
# Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_Allprobes_SNF_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[ , match_RPMM2_bySNF], MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[ , match_RPMM2_bySNF], ListofProbes = rownames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes)[ , ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[match_RPMM2_bySNF, ], FigureGenerate="Yes", PNGorPDF="png", ImageName="Halfprobes_NoSexChsomes_FLonly_SNF30samples")

# Using 585 probes from 2 RPMM model, run SNF 30 samples
# Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_2RPMM585probes_SNF30samples_14 <- Clustering(BetaMatrix = Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered[ , match_RPMM2_bySNF], MvalueMatrix = Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_18$MvalueMatrix_SD_Filtered[ , match_RPMM2_bySNF], ListofProbes = rownames(Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_18$MvalueMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[match_RPMM2_bySNF, ], FigureGenerate="Yes", PNGorPDF="png", ImageName="585probes_NoSexChsomes_FLonly_SNF30samples")






# 26 July 2019
# During meeting, asked to try different cutoffs of SD/ probe selection
# Trying 607 probes from 3RPMM model with FL only samples
# Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_607probes_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], ListofProbes = rownames(Output_SDeviation_0.25SD_noSexChsomes_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL"),], FigureGenerate="No", PNGorPDF="png", ImageName="0.25SDNoSexChsomes_FLonly_with607ProbesFromRPMM")
# CROSS tabulation between 607 probes and 585 probes, both for 2 RPMM model
# table(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_607probes_14$RPMM$RPMMoutputLabels, Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels)
#   1  2
#1 72  0
#2  2 82


# cat("\n Running SD of 0.2 with no sex xsomes, FL samples only =  5394 probes")
# Output_SDeviation_0.2SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], MvalueMatrix=Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], standard_deviation=0.2)
# Output_Clustering_0.2SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], ListofProbes = rownames(Output_SDeviation_0.2SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes),4,5) =="FL"),], FigureGenerate="No", PNGorPDF="png", ImageName="0.2SDNoSexChsomes_FLonly")
# table(Output_Clustering_0.2SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels, Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels)



# 6 August 2019
# Robert suggested to cluster 10 - 50 % of dataset by relaxing SD
# BetaMatrix_updSamples_rmdSexChsomes <- Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes
# MvalueMatrix_updSamples_rmdSexChsomes <- Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes
# ClinicalFile_updSamples <- Output_Remove_2$ClinicalFile_updSamples
# write.csv(BetaMatrix_updSamples_rmdSexChsomes, file = "BetaMatrix_updSamples_rmdSexChsomes.csv")
# write.csv(MvalueMatrix_updSamples_rmdSexChsomes, file = "MvalueMatrix_updSamples_rmdSexChsomes.csv")
# write.csv(ClinicalFile_updSamples, file = "ClinicalFile_updSamples.csv")

# cat("\n Reading files")
# BetaMatrix_updSamples_rmdSexChsomes <- read.csv(file="BetaMatrix_updSamples_rmdSexChsomes.csv", row.names = 1)
# MvalueMatrix_updSamples_rmdSexChsomes <- read.csv(file="MvalueMatrix_updSamples_rmdSexChsomes.csv")
# ClinicalFile_updSamples <- read.csv(file = "ClinicalFile_updSamples.csv")

# cat("\n Running SD of 0.3 with no sex xsomes, FL samples only =  probes")
# Output_SDeviation_0.3SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
#                                                                     MvalueMatrix = NA, 
#                                                                     standard_deviation = 0.3)
# Output_Clustering_0.3SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]), 
#                                                                      MvalueMatrix = NA, 
#                                                                      ListofProbes = rownames(Output_SDeviation_0.3SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                      ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL"), ], 
#                                                                      FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.3SDNoSexChsomes_FLonly")
# write.csv(cbind(colnames(Output_Clustering_0.3SD_noSexChsomes_FLonlysamples_14$BetaMatrixProbesUsedClustering), Output_Clustering_0.3SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels), file = "70_point3SD_Probes_4ClusterModel_RPMM_Labels.csv")
# Using 70 probes obtained from FL only samples across ALL samples
# Output_Clustering_0.3SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                    MvalueMatrix = NA, 
#                                                                    ListofProbes = rownames(Output_SDeviation_0.3SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                   ClinicalFile = ClinicalFile_updSamples, 
#                                                                    FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.3SD(fromFL)_NoSexChsomes_ALLsamples")


# cat("\n Running SD of 0.29 with no sex xsomes, FL samples only = 127  probes")
#  Output_SDeviation_0.29SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
#                                                                      MvalueMatrix = NA, 
#                                                                      standard_deviation = 0.29)
# Output_Clustering_0.29SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]), 
#                                                                      MvalueMatrix = NA, 
#                                                                      ListofProbes = rownames(Output_SDeviation_0.29SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                      ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL"), ], 
#                                                                      FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.3SDNoSexChsomes_FLonly")
# write.csv(cbind(colnames(Output_Clustering_0.29SD_noSexChsomes_FLonlysamples_14$BetaMatrixProbesUsedClustering), Output_Clustering_0.29SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels), file = "127 _point29SD_Probes_XClusterModel_RPMM_Labels.csv")
# Using  probes obtained from FL only samples across ALL samples
# Output_Clustering_0.29SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                   MvalueMatrix = NA, 
#                                                                   ListofProbes = rownames(Output_SDeviation_0.29SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                   ClinicalFile = ClinicalFile_updSamples, 
#                                                                   FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.29SD(fromFL)_NoSexChsomes_ALLsamples")


# cat("\n Running SD of 0.28 with no sex xsomes, FL samples only =  185 probes")
#   Output_SDeviation_0.28SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
#                                                                       MvalueMatrix = NA, 
#                                                                       standard_deviation = 0.28)
# Output_Clustering_0.28SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]), 
#                                                                      MvalueMatrix = NA, 
#                                                                      ListofProbes = rownames(Output_SDeviation_0.28SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                      ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL"), ], 
#                                                                      FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.3SDNoSexChsomes_FLonly")
# write.csv(cbind(colnames(Output_Clustering_0.28SD_noSexChsomes_FLonlysamples_14$BetaMatrixProbesUsedClustering), Output_Clustering_0.29SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels), file = "127 _point29SD_Probes_XClusterModel_RPMM_Labels.csv")
# Using  probes obtained from FL only samples across ALL samples
# Output_Clustering_0.28SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                   MvalueMatrix = NA, 
#                                                                   ListofProbes = rownames(Output_SDeviation_0.28SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                   ClinicalFile = ClinicalFile_updSamples, 
#                                                                   FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.28SD(fromFL)_NoSexChsomes_ALLsamples")


# cat("\n Running SD of 0.27 with no sex xsomes, FL samples only = 265 probes")
# Output_SDeviation_0.27SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
#                                                                        MvalueMatrix = NA, 
#                                                                        standard_deviation = 0.27)
# Output_Clustering_0.27SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]), 
#                                                                      MvalueMatrix = NA, 
#                                                                      ListofProbes = rownames(Output_SDeviation_0.27SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                      ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL"), ],   #                                                                      FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.27SDNoSexChsomes_FLonly")
# write.csv(cbind(colnames(Output_Clustering_0.27SD_noSexChsomes_FLonlysamples_14$BetaMatrixProbesUsedClustering), Output_Clustering_0.27SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels), file = "127 _point29SD_Probes_XClusterModel_RPMM_Labels.csv")
# Using  probes obtained from FL only samples across ALL samples
# Output_Clustering_0.27SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                   MvalueMatrix = NA, 
#                                                                   ListofProbes = rownames(Output_SDeviation_0.27SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                   ClinicalFile = ClinicalFile_updSamples, 
#                                                                   FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.27SD(fromFL)_NoSexChsomes_ALLsamples")



# cat("\n Running SD of 0.26 with no sex xsomes, FL samples only = 265 probes")
# Output_SDeviation_0.26SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
#                                                                      MvalueMatrix = NA, 
#                                                                      standard_deviation = 0.27)
# Output_Clustering_0.26SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]), 
#                                                                      MvalueMatrix = NA, 
#                                                                      ListofProbes = rownames(Output_SDeviation_0.26SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                      ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL"), ],   #                                                                      FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.27SDNoSexChsomes_FLonly")
# write.csv(cbind(colnames(Output_Clustering_0.26SD_noSexChsomes_FLonlysamples_14$BetaMatrixProbesUsedClustering), Output_Clustering_0.26SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels), file = "127 _point29SD_Probes_XClusterModel_RPMM_Labels.csv")
# Using  probes obtained from FL only samples across ALL samples
# Output_Clustering_0.26SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                   MvalueMatrix = NA, 
#                                                                   ListofProbes = rownames(Output_SDeviation_0.26SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                   ClinicalFile = ClinicalFile_updSamples, 
#                                                                   FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.26SD(fromFL)_NoSexChsomes_ALLsamples")



# cat("\n Running SD of 0.25 with no sex xsomes, FL samples only = 585 probes")
Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
                                                                     MvalueMatrix = NA, 
                                                                     standard_deviation = 0.25)
# Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]), 
#                                                                       MvalueMatrix = NA, 
#                                                                       ListofProbes = rownames(Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                      ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes),4,5) == "FL"), ], 
#                                                                       FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.25SDNoSexChsomes_FLonly")
# Using 585 probes obtained from FL only samples across ALL samples
Output_Clustering_0.25SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
                                                                  MvalueMatrix = NA, 
                                                                  ListofProbes = rownames(Output_SDeviation_0.25SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
                                                                  ClinicalFile = ClinicalFile_updSamples, 
                                                                  FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.25SD(fromFL)_NoSexChsomes_ALLsamples")
# Using 607 probes obtained from ALL samples
# Output_SDeviation_0.25SD_noSexChsomes_Allsamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes, 
#                                                                       MvalueMatrix = NA, 
#                                                                       standard_deviation = 0.25)
# Output_Clustering_0.25SD_noSexChsomes_AllSamples_14_607probes <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                             MvalueMatrix = NA, 
#                                                                             ListofProbes = rownames(Output_SDeviation_0.25SD_noSexChsomes_Allsamples_18$BetaMatrix_SD_Filtered), 
#                                                                             ClinicalFile = ClinicalFile_updSamples, 
#                                                                             FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.25SD_NoSexChsomes_ALLsamples")



# cat("\n Running SD of 0.24 with no sex xsomes, FL samples only = 897 probes")
# Output_SDeviation_0.24SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
#                                                                      MvalueMatrix = NA, 
#                                                                      standard_deviation = 0.24)
# Output_Clustering_0.24SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]), 
#                                                                      MvalueMatrix = NA, 
#                                                                      ListofProbes = rownames(Output_SDeviation_0.24SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                      ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes),4,5) == "FL"), ], 
#                                                                      FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.24SDNoSexChsomes_FLonly")
# Using 897 probes obtained from FL only samples across ALL samples
# Output_Clustering_0.24SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                  MvalueMatrix = NA, 
#                                                                  ListofProbes = rownames(Output_SDeviation_0.24SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                  ClinicalFile = ClinicalFile_updSamples, 
#                                                                  FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.24SD(fromFL)_NoSexChsomes_ALLsamples")
# write.csv(cbind(colnames(Output_Clustering_0.24SD_noSexChsomes_AllSamples_14$BetaMatrixProbesUsedClustering), Output_Clustering_0.24SD_noSexChsomes_AllSamples_14$RPMM$RPMMoutputLabels), file = "897_point24SD_Probes_3ClusterModel_RPMM_Allsamples_Labels.csv")


# cat("\n Running SD of 0.23 with no sex xsomes, FL samples only = 1381 probes")
# Output_SDeviation_0.23SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
#                                                                      MvalueMatrix = NA, 
#                                                                      standard_deviation = 0.23)
# Output_Clustering_0.23SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]), 
#                                                                      MvalueMatrix = NA, 
#                                                                      ListofProbes = rownames(Output_SDeviation_0.23SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                      ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes),4,5) == "FL"), ], 
#                                                                      FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.23SDNoSexChsomes_FLonly")
# Using 1381 probes obtained from FL only samples across ALL samples
# Output_Clustering_0.23SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                 MvalueMatrix = NA, 
#                                                                 ListofProbes = rownames(Output_SDeviation_0.23SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                 ClinicalFile = ClinicalFile_updSamples, 
#                                                                 FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.23SD(fromFL)_NoSexChsomes_ALLsamples")
# write.csv(cbind(colnames(Output_Clustering_0.23SD_noSexChsomes_AllSamples_14$BetaMatrixProbesUsedClustering), Output_Clustering_0.23SD_noSexChsomes_AllSamples_14$RPMM$RPMMoutputLabels), file = "1381_point23SD_Probes_3ClusterModel_RPMM_Allsamples_Labels.csv")




# cat("\n Running SD of 0.2 with no sex xsomes, FL samples only = 5622 probes")
Output_SDeviation_0.2SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_T1[, c(11:165)], 
                                                                    MvalueMatrix = NA, 
                                                                    standard_deviation = 0.2)
Output_Clustering_0.2SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = BetaMatrix_T1[, c(11:165)], 
                                                                    MvalueMatrix = NA, 
                                                                    ListofProbes = rownames(Output_SDeviation_0.2SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
                                                                    ClinicalFile = ClinicalFile_T1[c(11:165), ], 
                                                                    FigureGenerate = "Yes", 
                                                                    PNGorPDF = "png", 
                                                                    ImageName = "0.2SDNoSexChsomes_FLonly")
# write.csv(cbind(colnames(Output_Clustering_0.2SD_noSexChsomes_FLonlysamples_14$BetaMatrixProbesUsedClustering), Output_Clustering_0.2SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels), file = "5394_point2SD_Probes_4ClusterModel_RPMM_Labels.csv")
# Using 5394 probes obtained from FL only samples across ALL samples
Output_Clustering_0.2SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
                                                                 MvalueMatrix = NA, 
                                                                 ListofProbes = rownames(Output_SDeviation_0.2SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
                                                                 ClinicalFile = ClinicalFile_updSamples, 
                                                                 FigureGenerate = "No", PNGorPDF = "png", ImageName = "0.2SD(fromFL)_NoSexChsomes_ALLsamples")


# cat("\n Running SD of 0.19 with no sex xsomes, FL samples only = 8472 probes")
#  Output_SDeviation_0.19SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], 
#                                                                       MvalueMatrix = NA, 
#                                                                      standard_deviation = 0.19)
# Output_Clustering_0.19SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")]), 
#                                                                      MvalueMatrix = NA, 
#                                                                      ListofProbes = rownames(Output_SDeviation_0.19SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                      ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL"),], 
#                                                                      FigureGenerate="No", PNGorPDF="png", ImageName="0.19SDNoSexChsomes_FLonly")
# Using 8472 probes obtained from FL only samples across ALL samples
#  Output_Clustering_0.19SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                    MvalueMatrix = NA, 
#                                                                    ListofProbes = rownames(Output_SDeviation_0.19SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                    ClinicalFile = ClinicalFile_updSamples, 
#                                                                    FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.19SD(fromFL)_NoSexChsomes_ALLsamples")




#  cat("\n Running SD of 0.18 with no sex xsomes, FL samples only = 13615 probes")
#  Output_SDeviation_0.18SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], 
#                                                                      MvalueMatrix = NA, 
#                                                                      standard_deviation = 0.18)
# Output_Clustering_0.18SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes),4,5) =="FL")], 
#                                                                       MvalueMatrix = NA, 
#                                                                       ListofProbes = rownames(Output_SDeviation_0.18SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                       ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes),4,5) == "FL"),], 
#                                                                       FigureGenerate="Yes", PNGorPDF="png", ImageName="0.18SDNoSexChsomes_FLonly")
# Using 13615 probes obtained from FL only samples across ALL samples
# Output_Clustering_0.18SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                   MvalueMatrix = NA, 
#                                                                   ListofProbes = rownames(Output_SDeviation_0.18SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                   ClinicalFile = ClinicalFile_updSamples, 
#                                                                   FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.18SD(fromFL)_NoSexChsomes_ALLsamples")





#  cat("\n Running SD of 0.17 with no sex xsomes, FL samples only = 22839 probes")
# Output_SDeviation_0.17SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")],
#                                                                       MvalueMatrix = NA,
#                                                                       standard_deviation = 0.17)
# Output_Clustering_0.17SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]),
#                                                                      MvalueMatrix = NA,
#                                                                      ListofProbes = rownames(Output_SDeviation_0.17SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered),
#                                                                     ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) =="FL"),],
#                                                                     FigureGenerate = "Yes",
#                                                                     PNGorPDF = "png",
#                                                                     ImageName = "0.17SDNoSexChsomes_FLonly")
# Using 22839 probes obtained from FL only samples across ALL samples
# Output_Clustering_0.17SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                    MvalueMatrix = NA, 
#                                                                    ListofProbes = rownames(Output_SDeviation_0.17SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                    ClinicalFile = ClinicalFile_updSamples, 
#                                                                    FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.17SD(fromFL)_NoSexChsomes_ALLsamples")



# cat("\n Running SD of 0.16 with no sex xsomes, FL samples only = 39260 probes")
# Output_SDeviation_0.16SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")],
#                                                                     MvalueMatrix = NA,
#                                                                     standard_deviation = 0.16)
# Output_Clustering_0.16SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]),
#                                                                    MvalueMatrix = NA,
#                                                                    ListofProbes = rownames(Output_SDeviation_0.16SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered),
#                                                                    ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) =="FL"),],
#                                                                    FigureGenerate = "Yes",
#                                                                    PNGorPDF = "png",
#                                                                    ImageName = "0.16SDNoSexChsomes_FLonly")
# Using 39260 probes obtained from FL only samples across ALL samples
# Output_Clustering_0.16SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                   MvalueMatrix = NA, 
#                                                                   ListofProbes = rownames(Output_SDeviation_0.16SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                   ClinicalFile = ClinicalFile_updSamples, 
#                                                                   FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.16SD(fromFL)_NoSexChsomes_ALLsamples")



# cat("\n Running SD of 0.15 with no sex xsomes, FL samples only = 64613 probes")
#  Output_SDeviation_0.15SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")],
#                                                                        MvalueMatrix = NA,
#                                                                        standard_deviation = 0.15)
# Output_Clustering_0.15SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]),
#                                                                      MvalueMatrix = NA,
#                                                                      ListofProbes = rownames(Output_SDeviation_0.15SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered),
#                                                                      ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) =="FL"),],
#                                                                      FigureGenerate = "Yes",
#                                                                      PNGorPDF = "png",
#                                                                      ImageName = "0.15SDNoSexChsomes_FLonly")
# Using 64613 probes obtained from FL only samples across ALL samples
# Output_Clustering_0.15SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                   MvalueMatrix = NA, 
#                                                                   ListofProbes = rownames(Output_SDeviation_0.15SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                   ClinicalFile = ClinicalFile_updSamples, 
#                                                                   FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.15SD(fromFL)_NoSexChsomes_ALLsamples")



# cat("\n Running SD of 0.14 with no sex xsomes, FL samples only =  probes")
# Output_SDeviation_0.14SD_noSexChsomes_FLonlysamples_18 <- SDeviation(BetaMatrix = BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")],
#                                                                      MvalueMatrix = NA,
#                                                                      standard_deviation = 0.14)
# Output_Clustering_0.14SD_noSexChsomes_FLonlysamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes[,which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")]),
#                                                                     MvalueMatrix = NA,
#                                                                     ListofProbes = rownames(Output_SDeviation_0.14SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered),
#                                                                     ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) =="FL"),],
#                                                                     FigureGenerate = "Yes",
#                                                                     PNGorPDF = "png",
#                                                                     ImageName = "0.14SDNoSexChsomes_FLonly")
# Using 98523 probes obtained from FL only samples across ALL samples
# Output_Clustering_0.14SD_noSexChsomes_AllSamples_14 <- Clustering(BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), 
#                                                                   MvalueMatrix = NA, 
#                                                                   ListofProbes = rownames(Output_SDeviation_0.14SD_noSexChsomes_FLonlysamples_18$BetaMatrix_SD_Filtered), 
#                                                                   ClinicalFile = ClinicalFile_updSamples, 
#                                                                   FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.14SD(fromFL)_NoSexChsomes_ALLsamples")


# Comparing results with 2 RPMM model 
#RPMM_2RPMMmodel <- read.csv(file = "585Probes_2ClusterModel_RPMM_Labels.csv", row.names = 1)

# RPMM_0.3RPMMmodel <- read.csv(file = "70_point3SD_Probes_4ClusterModel_RPMM_Labels.csv", row.names = 1)
# table(RPMM_2RPMMmodel$V2, RPMM_0.3RPMMmodel$V2)
# print(chisq.test(table(RPMM_2RPMMmodel$V2, RPMM_0.3RPMMmodel$V2))) # X-squared = 21.02, df = 1, p-value = 4.546e-06
# adjustedRandIndex(RPMM_2RPMMmodel$V2, RPMM_0.3RPMMmodel$V2) #  0.1066499
# ((sum(diag(table(RPMM_2RPMMmodel$V2, RPMM_0.3RPMMmodel$V2)))) / 156) * 100 # 66.7
# ClinicalFile = ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) =="FL"),]
# adjustedRandIndex(RPMM_0.3RPMMmodel$V2, ClinicalFile$STAGE)
# table(RPMM_0.3RPMMmodel$V2, ClinicalFile$STAGE)


#  RPMM_0.25RPMMmodel_all <- read.csv(file = "607_point25_Probes_3ClusterModel_RPMM_Allsamples_Labels.csv")
#  RPMM_0.25RPMMmodel_all_585probes <- read.csv(file = "585_point25SD_Probes_3ClusterModel_RPMM_Allsamples_Labels.csv")


# RPMM_0.2RPMMmodel <- read.csv(file = "5394_point2SD_Probes_4ClusterModel_RPMM_Labels.csv", row.names = 1)
# table(RPMM_2RPMMmodel$V2, RPMM_0.2RPMMmodel$V2)
# print(chisq.test(table(RPMM_2RPMMmodel$V2, RPMM_0.2RPMMmodel$V2))) # X-squared = 110.39, df = 3, p-value < 2.2e-16
# adjustedRandIndex(RPMM_2RPMMmodel$V2, RPMM_0.2RPMMmodel$V2) #  0.356884
#((sum(diag(table(RPMM_2RPMMmodel$V2, RPMM_0.2RPMMmodel$V2))))/ 156 ) * 100 # 19.9
# RPMM_0.2RPMMmodel_all <- as.matrix(read.csv(file = "5394_point2SD_Probes_4ClusterModel_RPMM_Allsamples_Labels.csv", row.names = 1))
# table(as.numeric(RPMM_0.2RPMMmodel_all[match(ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 1], rownames(RPMM_0.2RPMMmodel_all)), ]), ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 16])
# Output_Clustering_0.2SD_noSexChsomes_FLonlysamples_14 <- read.csv(file = "5394_point2SD_Probes_4ClusterModel_RPMM_FLonly_Labels.csv", row.names = 1)



# RPMM_0.19RPMMmodel <- read.csv(file = "8472_point19SD_Probes_4ClusterModel_RPMM_Labels.csv", row.names = 1)
# table(RPMM_2RPMMmodel$V2, RPMM_0.19RPMMmodel$V2)
# print(chisq.test(table(RPMM_2RPMMmodel$V2, RPMM_0.19RPMMmodel$V2))) # X-squared = 107.73, df = 3, p-value < 2.2e-16
# adjustedRandIndex(RPMM_2RPMMmodel$V2, RPMM_0.19RPMMmodel$V2) #  0.3406146
# ((sum(diag(table(RPMM_2RPMMmodel$V2, RPMM_0.19RPMMmodel$V2)))) / 156) * 100 # 23.08
# RPMM_0.19RPMMmodel_all <- as.matrix(read.csv(file = "8472_point19SD_Probes_4ClusterModel_RPMM_Allsamples_Labels.csv", row.names = 1))
# table(as.numeric(RPMM_0.19RPMMmodel_all[match(ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 1], rownames(RPMM_0.19RPMMmodel_all)), ]), ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 16])


# RPMM_0.18RPMMmodel <- read.csv(file = "13615_point18SD_Probes_4ClusterModel_RPMM_Labels.csv", row.names = 1)
# table(RPMM_2RPMMmodel$V2, RPMM_0.18RPMMmodel$V2)
# print(chisq.test(table(RPMM_2RPMMmodel$V2, RPMM_0.18RPMMmodel$V2))) # X-squared = 93.464, df = 3, p-value < 2.2e-16
# adjustedRandIndex(RPMM_2RPMMmodel$V2, RPMM_0.18RPMMmodel$V2) #  0.2879752
# ((sum(diag(table(RPMM_2RPMMmodel$V2, RPMM_0.18RPMMmodel$V2)))) / 156) * 100 # 21.79
# RPMM_0.18RPMMmodel_all <- as.matrix(read.csv(file = "13615_point18SD_Probes_4ClusterModel_RPMM_Allsamples_Labels.csv", row.names = 1))
# table(as.numeric(RPMM_0.18RPMMmodel_all[match(ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 1], rownames(RPMM_0.18RPMMmodel_all)), ]), ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 16])




# RPMM_0.17RPMMmodel <- read.csv(file = "22839_point17SD_Probes_4ClusterModel_RPMM_Labels.csv", row.names = 1)
# table(RPMM_2RPMMmodel$V2, RPMM_0.17RPMMmodel$V2)
# print(chisq.test(table(RPMM_2RPMMmodel$V2, RPMM_0.17RPMMmodel$V2))) # X-squared = 87.191, df = 3, p-value < 2.2e-16
# adjustedRandIndex(RPMM_2RPMMmodel$V2, RPMM_0.17RPMMmodel$V2) # 0.2504224
#((sum(diag(table(RPMM_2RPMMmodel$V2, RPMM_0.17RPMMmodel$V2)))) / 156) * 100 # 26.28
# RPMM_0.17RPMMmodel_all <- as.matrix(read.csv(file = "22839_point17SD_Probes_5ClusterModel_RPMM_Allsamples_Labels.csv", row.names = 1))
# table(as.numeric(RPMM_0.17RPMMmodel_all[match(ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 1], rownames(RPMM_0.17RPMMmodel_all)), ]), ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 16])



# RPMM_0.16RPMMmodel <- read.csv(file = "39260_point16SD_Probes_4ClusterModel_RPMM_Labels.csv", row.names = 1)
# table(RPMM_2RPMMmodel$V2, RPMM_0.16RPMMmodel$V2)
# print(chisq.test(table(RPMM_2RPMMmodel$V2, RPMM_0.16RPMMmodel$V2))) # X-squared = 63.393, df = 3, p-value = 1.107e-13
# adjustedRandIndex(RPMM_2RPMMmodel$V2, RPMM_0.16RPMMmodel$V2) # 0.1577759
# ((sum(diag(table(RPMM_2RPMMmodel$V2, RPMM_0.16RPMMmodel$V2)))) / 156) * 100 # 30.77
# RPMM_0.16RPMMmodel_all <- as.matrix(read.csv(file = "39260_point16SD_Probes_5ClusterModel_RPMM_Allsamples_Labels.csv", row.names = 1))
# table(as.numeric(RPMM_0.16RPMMmodel_all[match(ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 1], rownames(RPMM_0.16RPMMmodel_all)), ]), ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 16])



# RPMM_0.15RPMMmodel <- read.csv(file = "64613_point15SD_Probes_4ClusterModel_RPMM_Labels.csv", row.names = 1)
# table(RPMM_2RPMMmodel$V2, RPMM_0.15RPMMmodel$V2)
# print(chisq.test(table(RPMM_2RPMMmodel$V2, RPMM_0.15RPMMmodel$V2))) # X-squared = 60.031, df = 3, p-value = 5.791e-13
# adjustedRandIndex(RPMM_2RPMMmodel$V2, RPMM_0.15RPMMmodel$V2) # 0.1470181
#((sum(diag(table(RPMM_2RPMMmodel$V2, RPMM_0.15RPMMmodel$V2)))) / 156) * 100 # 39.10
# RPMM_0.15RPMMmodel_all <- as.matrix(read.csv(file = "64613_point15SD_Probes_5ClusterModel_RPMM_Allsamples_Labels.csv", row.names = 1))
# table(as.numeric(RPMM_0.15RPMMmodel_all[match(ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 1], rownames(RPMM_0.15RPMMmodel_all)), ]), ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 16])


# RPMM_0.14RPMMmodel <- read.csv(file = "98523_point14SD_Probes_3ClusterModel_RPMM_Labels.csv", row.names = 1)
# table(RPMM_2RPMMmodel$V2, RPMM_0.14RPMMmodel$V2)
# print(chisq.test(table(RPMM_2RPMMmodel$V2, RPMM_0.14RPMMmodel$V2))) # X-squared = 60.031, df = 3, p-value = 5.791e-13
# adjustedRandIndex(RPMM_2RPMMmodel$V2, RPMM_0.14RPMMmodel$V2) # 0.1470181
# ((sum(diag(table(RPMM_2RPMMmodel$V2, RPMM_0.14RPMMmodel$V2)))) / 156) * 100  # 48.72
# RPMM_0.14RPMMmodel_all <- as.matrix(read.csv(file = "98523_point14SD_Probes_4ClusterModel_RPMM_Allsamples_Labels.csv", row.names = 1))
# table(as.numeric(RPMM_0.14RPMMmodel_all[match(ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 1], rownames(RPMM_0.14RPMMmodel_all)), ]), ClinicalFile[which(! is.na(ClinicalFile$COO) == TRUE), 16])


# Output_ProportionVisualization_Clustering_SD_0.25_Allprobes_3RPMM <- Varaince(CategoryToVisualize=NA, ClusterLabels = Output_Clustering_0.25SD_noSexChsomes_AllSamples_14$RPMM$RPMMoutputLabels, BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), AnnotationFile=NA, ClinicalFile=NA, ClinicalCategory=NA, PlotWithinCategories=NA, FigureGenerate="Yes", PNGorPDF="png", ImageName="3RPMMmodel_0.25SD_Allprobes")
Output_ProportionVisualization_Clustering_SD_0.25_Allprobes_3RPMM <- Varaince(CategoryToVisualize=NA, ClusterLabels = Output_Clustering_0.25SD_noSexChsomes_AllSamples_14$RPMM$RPMMoutputLabels, BetaMatrix = as.matrix(BetaMatrix_updSamples_rmdSexChsomes), AnnotationFile=NA, ClinicalFile=NA, ClinicalCategory=NA, PlotWithinCategories=NA, FigureGenerate="Yes", PNGorPDF="png", ImageName="3RPMMmodel_0.25SD_Allprobes")



# Looking into copy number data - 12 Sept 2019
# copy_number_data <- read.delim(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/CopyNumberAnalysis/all_lesions.txt", row.names = 1)
copy_number_data_conf90 <- read.delim(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/CopyNumberAnalysis/all_lesions.conf_90.txt", row.names = 1)
dim(copy_number_data_conf90) # 146 181
dim(copy_number_data_conf90[,c(9:180)]) # 146 172
colnames_copy_number_data_conf90 <- colnames(copy_number_data_conf90)

# correct the column names
colnames_KRI <- which(substr(colnames(copy_number_data_conf90), 1, 4) == "KRI_")
colnames_copy_number_data_conf90[colnames_KRI] <- 
  substr(colnames(copy_number_data_conf_90[, colnames_KRI]), 5, 16)

colnames_copy_number_data_conf90[which(substr(colnames_copy_number_data_conf90, 11, 12) == "")] <-
  paste0(colnames_copy_number_data_conf90[which(substr(colnames_copy_number_data_conf90, 11, 12) == "")], "_T1")

# assign the correct column names
colnames(copy_number_data_conf90) <- colnames_copy_number_data_conf90
# remove "LY_FL_159_T1"
copy_number_data_conf90 <- copy_number_data_conf90[, - which(colnames(copy_number_data_conf90) == "LY_FL_159_T1")]
dim(copy_number_data_conf90[, 9:179]) # 146 171
copy_number_data_conf90 <- copy_number_data_conf90[, 9:179]
colnames(copy_number_data_conf90)[which(colnames(copy_number_data_conf90) == "LY_FL_159_T1_rep")] <- "LY_FL_159_T1"
dim(copy_number_data_conf90) # 146 171

# Compare sample names with that of processed beta matrix
MatchingNames <- match(colnames(Output_Remove_2_BetaMatrix), colnames(copy_number_data_conf90))
colnames(copy_number_data_conf90[, MatchingNames[! is.na(MatchingNames)]])
dim(copy_number_data_conf90[, MatchingNames[! is.na(MatchingNames)]]) # 146 155
copy_number_data_conf90 <- copy_number_data_conf90[, MatchingNames[! is.na(MatchingNames)]]
dim(copy_number_data_conf90) # 146 155

# 4RPMM model
# ClusterLabelMatrix_4RPMM <- read.csv(file ="/Volumes/GoogleDrive/My Drive/UHN/FLOMICS/FLOMICS/Methylation/Pipeline/DataFiles/5394_point2SD_Probes_4ClusterModel_RPMM_FLonly_Labels.csv", row.names = 1)
# dim(ClusterLabelMatrix_4RPMM) # 156   1
# Select copy number data for only the 156 patients

FourRPMMmodel_CopyNumber <- CopyNumberAnalysis(CopyNumberData = copy_number_data, ClusterLabels = ClusterLabelMatrix_4RPMM)



# 3RPMM model
# ClusterLabelMatrix_3RPMM <- matrix(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels)
# rownames(ClusterLabelMatrix_3RPMM) <- colnames(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$BetaMatrixProbesUsedClustering)
# ThreeRPMMmodel_CopyNumber <- CopyNumberAnalysis(CopyNumberData = copy_number_data, ClusterLabels = ClusterLabelMatrix_3RPMM)



# 17 Sept 2019
# Decided to choose 0.2SD cutoff and run DMRCate analysis (with 6 contrasts)

# Output_SummarizationBeta_SDeviation_0.2SD_noSexChsomes_FLonlysamples_15 <- SurvivalAnalysis(BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[, which(substr(colnames(Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")],
#                    ClinicalFile = Output_Remove_2$ClinicalFile_updSamples,
#                    ClusterLabels = Output_Clustering_0.2SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels,
#                    SurvivalFile = SurvivalFile, 
#                    FigureGenerate = "No",
#                    PNGorPDF = "png")


# Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD <- Diff_MethylatedRegions(Method = "DMRcate", 
#                                                                                BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[, which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
#                                                                                MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[, which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
#                                                                                ContrastColumnName = "CLUSTER", 
#                                                                                ClusterLabels = Output_Clustering_0.2SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels, 
#                                                                                ClinicalFile = Output_Remove_2$ClinicalFile_updSamples, 
#                                                                                AnnotationFile = AnnotationFile, 
#                                                                                ProduceImages = "No", DMR = 1, PNGorPDF = "png", ExpressionFile = NA)

# Filter out probes 2 nucleotides or closer to a SNP that have a minor allele frequency greater than 0.05
# After filtering, 35897 probes were removed from the original 559110 probes. 
# Now there are 523213 probes.
# Contrast considered is: one-two 
# Your contrast returned 36353 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

# Contrast considered is: one-three 
# Your contrast returned 151197 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

# Contrast considered is: one-four 
# Your contrast returned 217048 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

# Contrast considered is: two-three 
# Your contrast returned 44316 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

# Contrast considered is: two-four 
# Your contrast returned 173280 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

# Contrast considered is: three-four 
# Your contrast returned 65864 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

# Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD_Betaplots <- MethylationDensityPlot(ClusterLabels = Output_Clustering_0.2SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels, BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[, which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes), 4, 5) =="FL")], FigureGenerate = "Yes", ImageName = "0.2SD_FLonly_2RPMM", PNGorPDF = "png")

#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM <- cbind(Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$DifferentiallyMethylatedRegions[[1]]$results, Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$GRangesObject[[1]])
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_positive <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM$maxbetafc > 0, ]
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_negative <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM$maxbetafc < 0, ]

#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM <- cbind(Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$DifferentiallyMethylatedRegions[[2]]$results, Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$GRangesObject[[2]])
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM_positive <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM$maxbetafc > 0, ]
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM_negative <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM$maxbetafc < 0, ]

#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM <- cbind(Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$DifferentiallyMethylatedRegions[[3]]$results, Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$GRangesObject[[3]])
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM_positive <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM$maxbetafc > 0, ]
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM_negative <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM$maxbetafc < 0, ]

#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM <- cbind(Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$DifferentiallyMethylatedRegions[[4]]$results, Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$GRangesObject[[4]])
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM_positive <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM$maxbetafc > 0, ]
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM_negative <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM$maxbetafc < 0, ]

#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM <- cbind(Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$DifferentiallyMethylatedRegions[[5]]$results, Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$GRangesObject[[5]])
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM_positive <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM$maxbetafc > 0, ]
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM_negative <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM$maxbetafc < 0, ]

#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM <- cbind(Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$DifferentiallyMethylatedRegions[[6]]$results, Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$GRangesObject[[6]])
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM_positive <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM$maxbetafc > 0, ]
#  table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM_negative <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM[table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM$maxbetafc < 0, ]


# write.csv(table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM, file="Table_Diff_MethylatedRegions_CLUSTER_AllProbes_C1vsC2_4RPMM.csv")
# write.csv(table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM, file="Table_Diff_MethylatedRegions_CLUSTER_AllProbes_C1vsC3_4RPMM.csv")
# write.csv(table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM, file="Table_Diff_MethylatedRegions_CLUSTER_AllProbes_C1vsC4_4RPMM.csv")
# write.csv(table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM, file="Table_Diff_MethylatedRegions_CLUSTER_AllProbes_C2vsC3_4RPMM.csv")
# write.csv(table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM, file="Table_Diff_MethylatedRegions_CLUSTER_AllProbes_C2vsC4_4RPMM.csv")
# write.csv(table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM, file="Table_Diff_MethylatedRegions_CLUSTER_AllProbes_C3vsC4_4RPMM.csv")


# Function extract gene regions list
ExtractingElements <- function(x, FILE) {
  
  example <- unlist(strsplit(as.character(FILE$overlapping.promoters[x]), "[,]"))
  unique_elements <- unique(trimws(sub("\\-.*", "", example)))
  
  GettingTable <- function(i, n, numb_unique_elements) {
    element_matrix <- matrix(ncol = 2, nrow = numb_unique_elements)
    element_matrix[i, 1] <- unique_elements[i]
    element_matrix[i, 2] <- FILE$Stouffer[n]
    return(element_matrix)
  }
  
  element_matrix <- GettingTable(i = c(1:length(unique_elements)), n = x, numb_unique_elements = length(unique_elements))
  return(element_matrix)
}

# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM <- sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM))
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM, quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM, file="GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM.csv")
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_positive <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_positive)), quote = FALSE)
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_negative <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_negative)), quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_positive, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_positive.csv")
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_negative, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_negative.csv")

# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM <- sapply(c(1:500), function(i) ExtractingElements(i, FILE=table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM))
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM, quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM, file="GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM.csv")
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM_positive <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM_positive)), quote = FALSE)
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM_negative <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM_negative)), quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM_positive, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM_positive.csv")
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM_negative, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_2ndContrast_4RPMM_negative.csv")


# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM <- sapply(c(1:500), function(i) ExtractingElements(i, FILE=table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM))
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM, quote=FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM, file="GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM.csv")
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM_positive <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM_positive)), quote = FALSE)
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM_negative <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM_negative)), quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM_positive, file="GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM_positive.csv")
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM_negative, file="GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_3rdContrast_4RPMM_negative.csv")


# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM <- sapply(c(1:500), function(i) ExtractingElements(i, FILE=table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM))
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM, quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM.csv")
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM_positive <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM_positive)), quote = FALSE)
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM_negative <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM_negative)), quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM_positive, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM_positive.csv")
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM_negative, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_4thContrast_4RPMM_negative.csv")


# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM <- sapply(c(1:500), function(i) ExtractingElements(i, FILE=table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM))
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM, quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM.csv")
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM_positive <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM_positive)), quote = FALSE)
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM_negative <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM_negative)), quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM_positive, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM_positive.csv")
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM_negative, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_5thContrast_4RPMM_negative.csv")


# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM <- sapply(c(1:500), function(i) ExtractingElements(i, FILE=table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM))
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM, quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM.csv")
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM_positive <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM_positive)), quote = FALSE)
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM_negative <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM_negative)), quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM_positive, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM_positive.csv")
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM_negative, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_6thContrast_4RPMM_negative.csv")

# Getting the path of the file, which should contain a folder called "img"
# pathNow <- getwd() 

# GProfiler_input <- read.csv(paste0(pathNow,"/DataFiles/GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM.csv"), row.names = 1)
# geneIDs_AllProbes_1stContrast_4RPMM <- as.vector(GProfiler_input$V1[! is.na(GProfiler_input$V1)])
# gprofiler(genes, organism = "hsapiens", max_p_value = 0.01,
# min_set_size = 5, max_set_size = 250, min_isect_size = 2, src_filter=c("GO:BP", "REAC"))
# gProfileR_AllProbes_1stContrast_4RPMM <- gProfileR::gprofiler(query = geneIDs_AllProbes_1stContrast_4RPMM, ordered_query = TRUE)


# geneIDs_AllProbes_1stContrast_4RPMM_positive <- as.vector(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_positive[, 1][!is.na(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_4RPMM_positive[, 1])])
# gProfileR_AllProbes_1stContrast_4RPMM_positive <- gProfileR::gprofiler(query = geneIDs_AllProbes_1stContrast_4RPMM_positive, ordered_query = TRUE, png_fn = TRUE, include_graph = TRUE)
# gProfileR_AllProbes_1stContrast_4RPMM_positive$intersection


# DMRcate analysis between advanced-stage and limited-stage patients; 8 October 2019
# Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD <- Diff_MethylatedRegions(Method = "DMRcate", 
#                                                                                 BetaMatrix = Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes[, which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
#                                                                                 MvalueMatrix = Output_Remove_2$MvalueMatrix_updSamples_rmdSexChsomes[, which(substr(colnames(Output_Remove_2$BetaMatrix_updSamples_rmdSexChsomes), 4, 5) == "FL")], 
#                                                                                 ContrastColumnName = "STAGE", 
#                                                                                ClusterLabels = NA, 
#                                                                                ClinicalFile = Output_Remove_2$ClinicalFile_updSamples[which(substr(colnames(BetaMatrix_updSamples_rmdSexChsomes), 4, 5) =="FL"),], 
#                                                                                 AnnotationFile = AnnotationFile, 
#                                                                                 ProduceImages = "No", DMR = 1, PNGorPDF = "png", ExpressionFile = NA)
# Your contrast returned 223545 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

# table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE <- cbind(Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$DifferentiallyMethylatedRegions[[1]]$results, Diff_MethylatedRegions_CLUSTER_AllProbes_2RPMM_0.2SD$GRangesObject[[1]])
# table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE_positive <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE[table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE$maxbetafc > 0, ]
# table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE_negative <- table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE[table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE$maxbetafc < 0, ]

# Function extract gene regions list
ExtractingElements <- function(x, FILE) {
  
  example <- unlist(strsplit(as.character(FILE$overlapping.promoters[x]), "[,]"))
  unique_elements <- unique(trimws(sub("\\-.*", "", example)))
  
  GettingTable <- function(i, n, numb_unique_elements) {
    element_matrix <- matrix(ncol = 2, nrow = numb_unique_elements)
    element_matrix[i, 1] <- unique_elements[i]
    element_matrix[i, 2] <- FILE$Stouffer[n]
    return(element_matrix)
  }
  
  element_matrix <- GettingTable(i = c(1:length(unique_elements)), n = x, numb_unique_elements = length(unique_elements))
  return(element_matrix)
}

# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE_positive <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE_positive)), quote = FALSE)
# Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE_negative <- do.call(rbind, sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE_negative)), quote = FALSE)
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE_positive, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE_positive.csv")
# write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE_negative, file = "GProfiler_Input_overlapping_promoters_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_STAGE_negative.csv")


# compare IDs - 10 October 2019
venndiagram_function <- function() {
  
  library(eulerr)
  TargetedSeq_Mehran <- read.csv("/Users/anjalisilva/Desktop/UHN/FLOMICS/TargetedSequencing/Sample_sumbission_Track_no._rcd7Oct2019/tra52892_Kridel92sample_15_failed_QC_190624TZ.csv")
  TargetedSeq2_Mehran <- read.csv("/Users/anjalisilva/Desktop/UHN/FLOMICS/TargetedSequencing/Sample_sumbission_Track_no._rcd7Oct2019/tra52893_Kridel41sample_5_failed_QC_190624TZ.csv")
  TargetedSeq_Mehran_ID <- as.character(TargetedSeq_Mehran$X.7)[2:93]
  TargetedSeq2_Mehran_ID <- as.character(TargetedSeq2_Mehran$External_ID)
  TargetedSeq_Mehran_ID <- unique(c(TargetedSeq_Mehran_ID, TargetedSeq2_Mehran_ID))
  
  TargetedSeq_BCCA <- read.csv("/Users/anjalisilva/Desktop/UHN/FLOMICS/TargetedSequencing/Results_BCCA_Sept2019/IX7722_CE010ANXX_3_gsc_library.summary.csv")
  TargetedSeq_BCCA2 <- read.csv("/Users/anjalisilva/Desktop/UHN/FLOMICS/TargetedSequencing/Results_BCCA_Sept2019/IX7723_CDYGFANXX_4_gsc_library.summary.csv")
  TargetedSeq_BCCA3 <- read.csv("/Users/anjalisilva/Desktop/UHN/FLOMICS/TargetedSequencing/Results_BCCA_Sept2019/IX7724_CDYGFANXX_5_gsc_library.summary.csv")
  TargetedSeq_BCCA_ID <- unique(c(as.character(TargetedSeq_BCCA$X.13)[-c(1:21,63:66)], as.character(TargetedSeq_BCCA2$X.13)[-c(1:21, 67:70)], as.character(TargetedSeq_BCCA3$X.13)[-c(1:21, 67:70)]))
  # TargetedSeq_BCCA_ID <- TargetedSeq_BCCA_ID[-which(substr(TargetedSeq_BCCA_ID, 11, 12) == "T2")]
  #  setdiff(TargetedSeq_BCCA_ID, colnames(BetaMatrix_updSamples_rmdSexChsomes_Ordered_T1))
  # Reduce(intersect, list(TargetedSeq_BCCA_ID,colnames(BetaMatrix_updSamples_rmdSexChsomes_Ordered_T1)))
  
  BetaMatrix_updSamples_rmdSexChsomes_Ordered_T1 <- readRDS(file = '2_BetaMatrix_updSamples_rmdSexChsomes_Ordered_T1.rds')
  ClinicalFile_updSamples_rmdSexChsomes_Ordered_T1 <- readRDS(file = '2_ClinicalFile_updSamples_rmdSexChsomes_Ordered_T1.rds')
  
  # CopyNumber_IDs <- colnames(read.delim(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/CopyNumberAnalysis/all_lesions.txt", row.names = 1))
  CopyNumber_IDs <- colnames(read.delim(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/CopyNumberAnalysis/all_lesions.txt", row.names = 1))
  
  library(eulerr)
  vd <- euler(c("A" = length(TargetedSeq_BCCA_ID), 
                "B" = length(colnames(BetaMatrix_updSamples_rmdSexChsomes_Ordered_T1)), 
                "C" = length(ClinicalFile_updSamples_rmdSexChsomes_Ordered_T1$SAMPLE_ID), 
                "D" = length(CopyNumber_IDs),
                "A&B" = length(Reduce(intersect, list(TargetedSeq_BCCA_ID, colnames(BetaMatrix_updSamples_rmdSexChsomes_Ordered_T1)))), 
                "A&C" = length(Reduce(intersect, list(TargetedSeq_BCCA_ID, ClinicalFile_updSamples_rmdSexChsomes_Ordered_T1$SAMPLE_ID))),
                "A&D" = length(Reduce(intersect, list(TargetedSeq_Mehran_ID, CopyNumber_IDs))), 
                "B&C" = length(Reduce(intersect, list(colnames(BetaMatrix_updSamples_rmdSexChsomes_Ordered_T1), ClinicalFile_updSamples_rmdSexChsomes_Ordered_T1$SAMPLE_ID))),
                "B&D" = length(Reduce(intersect, list(TargetedSeq_BCCA_ID, CopyNumber_IDs))),
                "C&D" = length(Reduce(intersect, list(Methylation_ID, CopyNumber_IDs))),
                "A&B&C" = length(Reduce(intersect, list(TargetedSeq_Mehran_ID, TargetedSeq_BCCA_ID, Methylation_ID))),
                "A&B&D" = length(Reduce(intersect, list(TargetedSeq_Mehran_ID, TargetedSeq_BCCA_ID, Methylation_ID))), 
                "B&C&D" = length(Reduce(intersect, list(TargetedSeq_BCCA_ID, Methylation_ID, CopyNumber_IDs))),
                "A&C&D" = length(Reduce(intersect, list(TargetedSeq_Mehran_ID, Methylation_ID, CopyNumber_IDs))),
                "A&B&C&D" = length(Reduce(intersect, list(TargetedSeq_Mehran_ID, TargetedSeq_BCCA_ID, Methylation_ID, CopyNumber_IDs)))))
  
  plot(vd, key = TRUE, counts = TRUE, quantities = TRUE, labels = c("TargetedSeq_Submitted", "TargetedSeq_BCCA", "Methylation_EPIC", "CopyNumber"))
  
  
  vd2_methylation_used <- euler(c("A" = length(TargetedSeq_BCCA_ID), 
                                  "B" = length(colnames(BetaMatrix_updSamples_rmdSexChsomes_Ordered_T1)), 
                                  "C" = length(ClinicalFile_updSamples_rmdSexChsomes_Ordered_T1$SAMPLE_ID), 
                                  "A&B" = length(Reduce(intersect, list(TargetedSeq_BCCA_ID, colnames(BetaMatrix_updSamples_rmdSexChsomes_Ordered_T1)))), 
                                  "A&C" = length(Reduce(intersect, list(TargetedSeq_BCCA_ID, ClinicalFile_updSamples_rmdSexChsomes_Ordered_T1$SAMPLE_ID))),
                                  "B&C" = length(Reduce(intersect, list(colnames(BetaMatrix_updSamples_rmdSexChsomes_Ordered_T1), ClinicalFile_updSamples_rmdSexChsomes_Ordered_T1$SAMPLE_ID))),
                                  "A&B&C" = length(Reduce(intersect, list(TargetedSeq_BCCA_ID, colnames(BetaMatrix_updSamples_rmdSexChsomes_Ordered_T1), ClinicalFile_updSamples_rmdSexChsomes_Ordered_T1$SAMPLE_ID)))))
  
  plot(vd2_methylation_used, key = TRUE, counts = TRUE, quantities = TRUE, labels = c("TargetedSeq_BCCA", "Methylation_EPIC", "ClinicalFile"))
  
}


# 15 Oct 2019
# Going back to performing differential methylation and removing unwanted variation

Output_Differential_probes_STAGE_7 <- DifferentialMethylation(ClinicalFile = Output_Remove_2_ClinicalFile, 
                                                              MvalueMatrix = Output_Remove_2_MvalueMatrix, 
                                                              BetaMatrix = Output_Remove_2_BetaMatrix,
                                                              ProbeOrGene = "Probe", 
                                                              ContrastColumnName = "STAGE", 
                                                              RGChannelSet = Output_Data_1$RGChannelSet, 
                                                              SampleSheet = Output_Data_1$SampleSheet, 
                                                              ProduceImages = "No", PNGorPDF = "png")

Comparison_RUVvsRegular <- Reduce(intersect, list(Output_Differential_probes_STAGE_7$SignificantProbes, rownames(MvalueMatrix_rfit4)))


# 2 Dec 2019
# Trying deconvolution via method used by Queiros et al., 2016
# Email received 29 Nov 2019

Y = B1 %*% T1 + B2 %*% T2
Y <- BetaMatrix_T1
# B1 <- # Unknown # adjustedBetaMatrix_T1
T1 <- purity_OICR_withNormal

# EstimatedCounts_170Samples <- readRDS(file = "EstimatedCellCOunts_FlowSortedBloodEPIC_13Dec2019.rds")
B2 <- 
  T2 <- EstimatedCounts_170Samples # patients vs cell count

# Sent by Marti
B <- Y
# Construct B2O2 term
B2O2 <- B2%*%t(O2)
# Proportion of pure cell types in mixture(cancer) samples
O1
# B1 term (in silico purified betas from Y)
B1 <- t(t(B-B2O2)/O1)


#13 Dec 2019
Output_Visuals_Relation_to_Island_STAGE_2 <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", BetaMatrix = BetaMatrix_T1, AnnotationFile = AnnotationFile, ClinicalFile = ClinicalFile_T1, ClinicalCategory = "STAGE", PlotWithinCategories="Yes", FigureGenerate="Yes", PNGorPDF="png", ImageName = "Relation_to_Island_STAGE_30Jan2020")
Output_Visuals_Relation_to_Island_TYPE_2 <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", BetaMatrix = BetaMatrix_T1, AnnotationFile = AnnotationFile, ClinicalFile = ClinicalFile_T1, ClinicalCategory = "TYPE", PlotWithinCategories="Yes", FigureGenerate="Yes", PNGorPDF="png", ImageName = "Relation_to_Island_TYPE_30Jan2020")
Density_RelationToIsland_TYPE <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", ClusterLabels = NA, PlotWithinAnnotationCategories = "Yes", ClinicalCategoryToVisualize = "TYPE", BetaMatrix = BetaMatrix_T1, AnnotationFile = AnnotationFile, ClinicalFile = ClinicalFile_T1, SampleSheet = sheet, FigureGenerate = "Yes", ImageName = NA, PNGorPDF = "png") 
Density_RelationToIsland_TYPE <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", ClusterLabels = NA, PlotWithinAnnotationCategories = "No", ClinicalCategoryToVisualize = "TYPE", BetaMatrix = BetaMatrix_T1, AnnotationFile = AnnotationFile, ClinicalFile = ClinicalFile_T1, SampleSheet = sheet, FigureGenerate = "Yes", ImageName = "3_DensityBeta_Relation_to_Island_NoWithinCategories", PNGorPDF = "png") 
Density_RelationToIsland_STAGE <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", ClusterLabels = NA, PlotWithinAnnotationCategories = "Yes", ClinicalCategoryToVisualize = "STAGE", BetaMatrix = BetaMatrix_T1, AnnotationFile = AnnotationFile, ClinicalFile = ClinicalFile_T1, SampleSheet = sheet, FigureGenerate = "Yes", ImageName = "3_DensityBeta_Relation_to_Island_NoWithinCategories", PNGorPDF = "png") 
Density_RelationToIsland_STAGE <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", ClusterLabels = NA, PlotWithinAnnotationCategories = "No", ClinicalCategoryToVisualize = "STAGE", BetaMatrix = BetaMatrix_T1, AnnotationFile = AnnotationFile, ClinicalFile = ClinicalFile_T1, SampleSheet = sheet, FigureGenerate = "Yes", ImageName = "3_DensityBeta_Relation_to_Island_NoWithinCategories", PNGorPDF = "png") 
ViolinPlot(BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, ClinicalFile = ClinicalFile_T1, CategoryToVisualize = "TYPE", PNGorPDF = "png")
ViolinPlot(BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, ClinicalFile = ClinicalFile_T1, CategoryToVisualize = "STAGE", PNGorPDF = "png")

cat("\n Running SD of 0.2 with no sex xsomes, no SNPs, of FL only samples =  5408  probes")
Output_SDeviation_0.2SD <- SDeviation(BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, CutOff = 0.2)
Output_SDeviation_10000 <- SDeviation(BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, NumberofProbes = 10000)

Output_Clustering_SD_0.2SD <- Clustering(BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, 
                                         ListofProbes = rownames(Output_SDeviation_10000$BetaMatrix_SD_Filtered), 
                                         ClinicalFile = ClinicalFile_T1, 
                                         FigureGenerate="No", 
                                         PNGorPDF="png", ImageName="10000_FLonly")
# FL only samples
Output_SDeviation_10000_FLonly <- SDeviation(BetaMatrix = BetaMatrix_T1[, 11:165], MvalueMatrix = MvalueMatrix_T1[, 11:165], NumberofProbes = 10000)
Output_SDeviation_10000_FLonly <- Clustering(BetaMatrix = BetaMatrix_T1[, 11:165], MvalueMatrix = MvalueMatrix_T1[, 11:165], ListofProbes = rownames(Output_SDeviation_10000$BetaMatrix_SD_Filtered), ClinicalFile = ClinicalFile_T1[11:165, ], FigureGenerate="No", PNGorPDF="png", ImageName="0.2SD_FLonly")



mean(MvalueMatrix_T1[,which(ClinicalFile_T1$TYPE == "FL")]) # 0.20
median(MvalueMatrix_T1[,which(ClinicalFile_T1$TYPE == "FL")]) # 0.90
range(MvalueMatrix_T1[,which(ClinicalFile_T1$TYPE == "FL")]) # -7.40  7.40


# 23 Jan 2020
# Surrogate variable analysis to batch effects

# 7.4 Batch effects correction with SVA
# https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#batch-effects-correction-with-sva
# SVA vignette
# http://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf

exampleSVA <- function() {
  library(sva)
  mval <- getM(GRset)[1:5000,]
  pheno <- pData(GRset)
  mod <- model.matrix(~as.factor(status), data=pheno) # status is normal or cancer
  mod0 <- model.matrix(~1, data=pheno)
  sva.results <- sva::sva(mval, mod, mod0)
}

ApplyingSVA <- function() {
  # from https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#batch-effects-correction-with-sva
  # from ? sva documentation
  # # http://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf
  # It is also possible to use the sva function with the ComBat function to
  # remove both known batch effects and other potential latent sources of variation.
  
  # identifying and estimating surrogate variables for unknown sources of variation in 
  # highthroughput experiments and
  # (1) adjustment variables = age of the patients, the sex of
  # the patients, and a variable like the date the arrays were processed.
  # (2) variables of interest = indicator of cancer versus control 
  
  ##### INSTITUTE = adjustment variable, TYPE = variable of interest ####
  # With TYPE
  library(sva)
  pheno <- Biobase::pData(Output_Data_2$GenomicMethylSet)
  # Full model - including both the adjustment variables and the variable of 
  # interest (cancer status).
  # In our case cancer type + insititution
  # ** Note: ClinicalFile_T1$STAGE could not be tested as only 155 values. 
  mod_TYPE_INSTITUTE <- stats::model.matrix(~(as.factor(ClinicalFile_T1$TYPE) + as.factor(ClinicalFile_T1$INSTITUTION)), data = pheno)
  dim(mod_TYPE_INSTITUTE) # 170   5
  # Null model contains only the adjustment variables.
  # In our case insititution
  mod0_TYPE_INSTITUTE <- stats::model.matrix(~as.factor(ClinicalFile_T1$INSTITUTION), data = pheno)
  dim(mod0_TYPE_INSTITUTE) # 170   5
  svaResults_TYPE_INSTITUTE <- sva::sva(MvalueMatrix_T1, mod_TYPE_INSTITUTE, mod0_TYPE_INSTITUTE)
  # Number of significant surrogate variables is:  25 
  # Iteration (out of 5 ):1  2  3  4  5  
  names(svaResults_TYPE_INSTITUTE) # "sv"        "pprob.gam" "pprob.b"   "n.sv"  
  # pprob.b is the posterior probability that each gene is associated
  # with the variables of interest
  head(svaResults_TYPE_INSTITUTE$pprob.b)
  length(svaResults_TYPE_INSTITUTE$pprob.b) # 595564
  # n.sv The number of significant surrogate variables
  svaResults_TYPE_INSTITUTE$n.sv #25
  # pprob.gam is the posterior probability that each gene is associated with one
  # or more latent variables
  svaResults_TYPE_INSTITUTE$pprob.gam
  length(svaResults_TYPE_INSTITUTE$pprob.gam) # 595564
  # sv is a matrix whose columns correspond to the estimated surrogate variables
  dim(svaResults_TYPE_INSTITUTE$sv) # 170  25
  
  # Adjusting for surrogate variables using the f.pvalue function
  pValues_TYPE_INSTITUTE = sva::f.pvalue(MvalueMatrix_T1, mod_TYPE_INSTITUTE, mod0_TYPE_INSTITUTE)
  qValues_TYPE_INSTITUTE = stats::p.adjust(pValues_TYPE_INSTITUTE, method = "BH")
  # length(which(qValues < 0.05)) # 245660
  # (245660/ 595564)*100 = 41.2%
  # Note that nearly 41% of the genes are strongly differentially expressed at an
  # FDR of less than 5% between groups. 
  
  # perform the same analysis, but adjusting for surrogate variables. 
  # The first step is to include the surrogate variables in both the null and 
  # full models. The reason is that we want to adjust for the surrogate 
  # variables, so we treat them as adjustment variables that must
  # be included in both models.
  modSv_TYPE_INSTITUTE = cbind(mod, svaResults_TYPE_INSTITUTE$sv)
  mod0Sv_TYPE_INSTITUTE = cbind(mod0, svaResults_TYPE_INSTITUTE$sv)
  pValuesSv_TYPE_INSTITUTE = sva::f.pvalue(MvalueMatrix_T1, modSv_TYPE_INSTITUTE, mod0Sv_TYPE_INSTITUTE)
  qValuesSv_TYPE_INSTITUTE = p.adjust(pValuesSv_TYPE_INSTITUTE, method = "BH")
  # length(which(qValuesSv_TYPE_INSTITUTE < 0.05)) # 34528
  # (34528/ 595564)*100 = 5.79753 %
  
  # calculating purity
  library(InfiniumPurify)
  data(iDMC) 
  probes <- iDMC[["DLBC"]]
  probes.true <- names(probes[probes == T])
  beta.sel <- BetaMatrix_T1[row.names(BetaMatrix_T1) %in% probes.true,]
  purity_OICR_withNormal <- getPurity(tumor.data = beta.sel, tumor.type = "DLBC")
  
  
  set.seed(1234)
  InfiniumClust_TYPE_INSTITUTE <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_T1[which(qValuesSv_TYPE_INSTITUTE < 0.05), ], 
                                                                purity = purity_OICR_withNormal, 
                                                                K = 2, maxiter = 100, tol = 0.0001)
  table(mclust::map(InfiniumClust_TYPE_INSTITUTE$Z))
  # 1  2 
  # 91 79 
  
  table(mclust::map(InfiniumClust_TYPE_INSTITUTE$Z), ClinicalFile_T1$STAGE)
  # ADVANCED LIMITED
  # 1  56      25      
  # 2  24      50     
  chisq.test(table(mclust::map(InfiniumClust_TYPE_INSTITUTE$Z), ClinicalFile_T1$STAGE))
  # X-squared = 19.416, df = 1, p-value = 1.051e-05
  
  
  table(mclust::map(InfiniumClust_TYPE_INSTITUTE$Z), ClinicalFile_T1$TRANSLOC_14_18)
  #  0  1
  # 1 16 32
  # 2 15  7
  
  chisq.test(table(mclust::map(InfiniumClust_TYPE_INSTITUTE$Z), ClinicalFile_T1$TRANSLOC_14_18))
  # X-squared = 6.0799, df = 1, p-value = 0.01367
  
  table(mclust::map(InfiniumClust_TYPE_INSTITUTE$Z), ClinicalFile_T1$TYPE)
  #    DLBCL FL RLN
  # 1   10 81   0  
  # 2   0 74   5   
  chisq.test(table(mclust::map(InfiniumClust_TYPE_INSTITUTE$Z), ClinicalFile_T1$TYPE))
  # X-squared = 14.542, df = 2, p-value = 0.0006956
  
  InfiniumClust_TYPE_INSTITUTE_WithPurity <- data.frame(cluster = mclust::map(InfiniumClust_TYPE_INSTITUTE$Z), 
                                                        purity = purity_OICR_withNormal,
                                                        stage = ClinicalFile_T1$STAGE,
                                                        translocation = ClinicalFile_T1$TRANSLOC_14_18,
                                                        type = ClinicalFile_T1$TYPE)
  
  my_comparisons <- list(c("1", "2"))
  p2_TYPE_INSTITUTE <- ggpubr::ggviolin(InfiniumClust_TYPE_INSTITUTE_WithPurity, 
                                        x = "cluster", y = "purity", fill = "cluster",
                                        add = "boxplot", ylab=" Purity") +
    ggtitle("All samples") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  p1_TYPE_INSTITUTE <- ggplot2::ggplot(InfiniumClust_TYPE_INSTITUTE_WithPurity, 
                                       aes(x = factor(cluster), y = purity, fill=stage)) +
    geom_boxplot(size = 1, outlier.size = 3) +
    ggtitle("All samples") +
    scale_y_continuous(name = "Purity") +
    # geom_jitter(aes(color = factor(stage)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  
  
  p3_TYPE_INSTITUTE <- ggplot2::ggplot(InfiniumClust_TYPE_INSTITUTE_WithPurity, 
                                       aes(x = factor(cluster), y = purity, fill=factor(translocation) )) +
    geom_boxplot(size = 1, outlier.size = 3) +
    ggtitle("All samples") +
    scale_y_continuous(name = "Purity") +
    # geom_jitter(aes(color = factor(stage)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    labs(fill = "translocation") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  # Outcome analysis
  SurvivalFile_All <- SurvivalFile[match(substr(colnames(BetaMatrix_T1), 1, 9), SurvivalFile$LY_FL_ID), ]
  BetaMatrix = BetaMatrix_T1[which(qValuesSv < 0.05), ]
  ClinicalFile = ClinicalFile_T1
  SurvivalFile = SurvivalFile_All
  ClusterLabels = InfiniumClust_TYPE_INSTITUTE_WithPurity$cluster
  # Then see 15_SurivivalAnalysis.R
  
  
  
  ##### INSTITUTE = adjustment variable, STAGE = variable of interest ####
  mod_STAGE_INSTITUTE <- stats::model.matrix(~(as.factor(ClinicalFile_T1$STAGE[- which(is.na(ClinicalFile_T1$STAGE) == TRUE)]) + as.factor(ClinicalFile_T1$INSTITUTION[- which(is.na(ClinicalFile_T1$STAGE) == TRUE)])), data = pheno[- which(is.na(ClinicalFile_T1$STAGE) == TRUE), ])
  dim(mod_STAGE_INSTITUTE) # 155   4
  # Null model contains only the adjustment variables.
  # In our case insititution
  mod0_STAGE_INSTITUTE <- stats::model.matrix(~as.factor(ClinicalFile_T1$INSTITUTION[- which(is.na(ClinicalFile_T1$STAGE) == TRUE)]), data = pheno[- which(is.na(ClinicalFile_T1$STAGE) == TRUE), ])
  dim(mod0_STAGE_INSTITUTE) # 155   3
  svaResults_STAGE_INSTITUTE <- sva::sva(MvalueMatrix_T1[, - which(is.na(ClinicalFile_T1$STAGE) == TRUE)], mod_STAGE_INSTITUTE, mod0_STAGE_INSTITUTE)
  # Number of significant surrogate variables is:  23 
  # Iteration (out of 5 ):1  2  3  4  5
  
  # perform the same analysis, but adjusting for surrogate variables. 
  modSv_STAGE_INSTITUTE = cbind(mod_STAGE_INSTITUTE, svaResults_STAGE_INSTITUTE$sv)
  mod0Sv_STAGE_INSTITUTE = cbind(mod0_STAGE, svaResults_STAGE_INSTITUTE$sv)
  pValuesSv_STAGE_INSTITUTE = sva::f.pvalue(MvalueMatrix_T1[, - which(is.na(ClinicalFile_T1$STAGE) == TRUE)], modSv_STAGE_INSTITUTE, mod0Sv_STAGE_INSTITUTE)
  qValuesSv_STAGE_INSTITUTE = p.adjust(pValuesSv_STAGE_INSTITUTE, method = "BH")
  # length(which(qValuesSv_STAGE < 0.05)) # 34528
  # (34528/ 595564)*100 = 5.79753 %
  
  
  # calculating purity
  library(InfiniumPurify)
  data(iDMC) 
  probes <- iDMC[["DLBC"]]
  probes.true <- names(probes[probes == T])
  beta.sel <- BetaMatrix_T1[row.names(BetaMatrix_T1) %in% probes.true,]
  purity_OICR_withNormal <- getPurity(tumor.data = beta.sel, tumor.type = "DLBC")
  
  
  set.seed(1234)
  InfiniumClust_STAGE_INSTITUTE <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_T1[which(qValuesSv_STAGE_INSTITUTE < 0.05), ], 
                                                                 purity = purity_OICR_withNormal, 
                                                                 K = 2, maxiter = 100, tol = 0.0001)
  table(mclust::map(InfiniumClust_STAGE_INSTITUTE$Z))
  # 1  2 
  # 81 89 
  
  table(mclust::map(InfiniumClust_STAGE_INSTITUTE$Z), ClinicalFile_T1$STAGE)
  # ADVANCED LIMITED
  # 1  61      13   
  # 2  19      62  
  chisq.test(table(mclust::map(InfiniumClust_STAGE_INSTITUTE$Z), ClinicalFile_T1$STAGE))
  # X-squared = 51.521, df = 1, p-value = 7.082e-13
  
  table(mclust::map(InfiniumClust_STAGE_INSTITUTE$Z), ClinicalFile_T1$TRANSLOC_14_18)
  #  0  1
  # 1 10 23
  # 2 21 16
  
  chisq.test(table(mclust::map(InfiniumClust_STAGE_INSTITUTE$Z), ClinicalFile_T1$TRANSLOC_14_18))
  # X-squared = 3.9332, df = 1, p-value = 0.04734
  
  
  InfiniumClust_STAGE_INSTITUTE_WithPurity <- data.frame(cluster = mclust::map(InfiniumClust_STAGE_INSTITUTE$Z), 
                                                         purity = purity_OICR_withNormal,
                                                         stage = ClinicalFile_T1$STAGE,
                                                         translocation = ClinicalFile_T1$TRANSLOC_14_18,
                                                         type = ClinicalFile_T1$TYPE)
  
  my_comparisons <- list(c("1", "2"))
  p2_STAGE_INSTITUTE <- ggpubr::ggviolin(InfiniumClust_STAGE_INSTITUTE_WithPurity, 
                                         x = "cluster", y = "purity", fill = "cluster",
                                         add = "boxplot", ylab=" Purity") +
    ggtitle("FL samples") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  p1_STAGE_INSTITUTE <- ggplot2::ggplot(InfiniumClust_STAGE_INSTITUTE_WithPurity, 
                                        aes(x = factor(cluster), y = purity, fill=stage)) +
    geom_boxplot(size = 1, outlier.size = 3) +
    ggtitle("FL samples") +
    scale_y_continuous(name = "Purity") +
    # geom_jitter(aes(color = factor(stage)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  
  p3_STAGE_INSTITUTE <- ggplot2::ggplot(InfiniumClust_STAGE_INSTITUTE_WithPurity, 
                                        aes(x = factor(cluster), y = purity, fill=factor(translocation) )) +
    geom_boxplot(size = 1, outlier.size = 3) +
    ggtitle("FL samples") +
    scale_y_continuous(name = "Purity") +
    # geom_jitter(aes(color = factor(stage)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    labs(fill = "translocation") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  # Outcome analysis
  SurvivalFile_All <- SurvivalFile[match(substr(colnames(BetaMatrix_T1), 1, 9), SurvivalFile$LY_FL_ID), ]
  BetaMatrix = BetaMatrix_T1[which(qValuesSv_STAGE_INSTITUTE < 0.05), ]
  ClinicalFile = ClinicalFile_T1
  SurvivalFile = SurvivalFile_All
  ClusterLabels = InfiniumClust_TYPE_INSTITUTE_WithPurity$cluster 
  # Then see 15_SurivivalAnalysis.R
  
  
  ##### DATE = adjustment variable, TYPE = variable of interest ####
  # The 30 samples were received early
  # The 140 samples were received later
  
  Patient_RNAseq_TDNAseq_Methylation.xlsx
  
  # with TYPE
  pheno <- Biobase::pData(Output_Data_2$GenomicMethylSet)
  # Full model - including both the adjustment variables and the variable of 
  # interest (cancer status).
  # In our case cancer type + DATE
  mod_TYPE_DATE <- stats::model.matrix(~(as.factor(ClinicalFile_T1$TYPE) + as.factor(DATE)), data = pheno)
  dim(mod_TYPE_DATE) # 170   5
  # Null model contains only the adjustment variables.
  # In our case DATE
  mod0_TYPE_DATE <- stats::model.matrix(~as.factor(DATE), data = pheno)
  dim(mod0_TYPE_DATE) # 170   5
  svaResults_TYPE_DATE <- sva::sva(MvalueMatrix_T1, mod_TYPE_DATE, mod0_TYPE_DATE)
  
  modSv_TYPE_DATE = cbind(mod_TYPE_DATE, svaResults_TYPE_DATE$sv)
  mod0Sv_TYPE_DATE = cbind(mod0_TYPE_DATE, svaResults_TYPE_DATE$sv)
  pValuesSv_TYPE_DATE = sva::f.pvalue(MvalueMatrix_T1, modSv_TYPE_DATE, mod0Sv_TYPE_DATE)
  qValuesSv_TYPE_DATE = p.adjust(pValuesSv_TYPE_DATE, method = "BH")
  # length(which(qValuesSv_TYPE_DATE < 0.05)) # 
  # (x / 595564)*100 =  %
  
  
  
  ##### DATE = adjustment variable, STAGE = variable of interest ####
  mod_STAGE_DATE <- stats::model.matrix(~(as.factor(ClinicalFile_T1$STAGE[- which(is.na(ClinicalFile_T1$STAGE) == TRUE)]) + as.factor(DATE[- which(is.na(ClinicalFile_T1$STAGE) == TRUE)])), data = pheno[- which(is.na(ClinicalFile_T1$STAGE) == TRUE), ])
  dim(mod_STAGE_DATE) # 155   4
  # Null model contains only the adjustment variables.
  # In our case insititution
  mod0_STAGE_DATE <- stats::model.matrix(~as.factor(DATE[- which(is.na(ClinicalFile_T1$STAGE) == TRUE)]), data = pheno[- which(is.na(ClinicalFile_T1$STAGE) == TRUE), ])
  dim(mod0_STAGE_DATE) # 155   3
  svaResults_STAGE_DATE <- sva::sva(MvalueMatrix_T1[, - which(is.na(ClinicalFile_T1$STAGE) == TRUE)], mod_STAGE_DATE, mod0_STAGE_DATE)
  
  # Adjusting for surrogate variables using the f.pvalue function
  pValues_STAGE_DATE = sva::f.pvalue(MvalueMatrix_T1[, - which(is.na(ClinicalFile_T1$STAGE) == TRUE)], mod_STAGE_DATE, mod0_STAGE_DATE)
  qValues_STAGE_DATE = stats::p.adjust(pValues_STAGE_DATE, method = "BH")
  # length(which(qValues_STAGE_DATE < 0.05)) # 
  # (X / 595564)*100 = 29.08957%
  
  # perform the same analysis, but adjusting for surrogate variables. 
  modSv_STAGE_DATE = cbind(mod_STAGE_DATE, svaResults_STAG_DATEE$sv)
  mod0Sv_STAGE_DATE = cbind(mod0_STAGE_DATE, svaResults_STAGE_DATE$sv)
  pValuesSv_STAGE_DATE = sva::f.pvalue(MvalueMatrix_T1[, - which(is.na(ClinicalFile_T1$STAGE) == TRUE)], modSv_STAGE_DATE, mod0Sv_STAGE_DATE)
  qValuesSv_STAGE_DATE = p.adjust(pValuesSv_STAGE_DATE, method = "BH")
  # length(which(qValuesSv_STAGE_DATE < 0.05)) # 
  # (X / 595564)*100 = 5.79753 %
  
  
  ##### DATE + INSTITUTE = adjustment variable, TYPE = variable of interest ####
  
  ##### DATE + INSTITUTE = adjustment variable, STAGE = variable of interest ####
}



# Once the surrogate variables are computed, one can include them in the downstream
# analysis to adjust for unknown unwanted variation. See the sva vignette for a
# more comprehensive use of sva.


# 31 Jan 2020
BatchInfo <- read.csv("Patient_RNAseq_TDNAseq_Methylation.csv", header = TRUE)
BatchInfo$RNAseq[which(BatchInfo$RNAseq_First_Time == TRUE)] # get names from 1st batch)
matchClinicalBatch <- match(BatchInfo$RNAseq[which(BatchInfo$RNAseq_First_Time == TRUE)], ClinicalFile_T1$SAMPLE_ID)
ClinicalFile_T1 <- readRDS(file = "1_ClinicalFile_updSamples_Ordered_T1.rds")
ClinicalFile_T1$SAMPLE_ID[matchClinicalBatch]
BatchVector <- rep("Second", 170)
BatchVector[matchClinicalBatch] <- "First"
ClinicalFile_T1 <- ClinicalFile_T1 %>% tibble::add_column(BATCH = factor(BatchVector))


# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
BetaMatrixPCA <- stats::prcomp(t(BetaMatrix_T1), scale = TRUE)
summary(BetaMatrixPCA) #  PC1 explains 14.8% of the total variance
# 170 Components in total   
devtools::install_github("vqv/ggbiplot")
library(ggbiplot)
# biplot: shows both the position of each sample in terms of PC1 and PC2
ggbiplot::ggbiplot(BetaMatrixPCA)


factoextra::fviz_eig(BetaMatrixPCA) 
# first dimension account for 88% of variability
# Add labels
fviz_eig(BetaMatrixPCA, addlabels=TRUE, hjust = -0.3) + ylim(0, 95)

colour_ind <- c("#00AFBB", "#E7B800", "#FC4E07")
# Visualize Principal Component Analysis- Graph of individuals
fviz_pca_ind(X = BetaMatrixPCA,
             col.ind = "y", # Color by the quality of representation
             gradient.cols = colour_ind,
             repel = TRUE ,    # Avoid text overlapping
             habillage ="none"
)


fviz_pca_var(BetaMatrixPCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Extract all the results (coordinates, squared cosine, contributions) for 
# the active individuals/variables from Principal Component Analysis (PCA) outputs.
BetaMatrixPCAind <- factoextra::get_pca_ind(BetaMatrixPCA)
BetaMatrixPCAind <- BetaMatrixPCAind$coord
BetaMatrixPCAind <- data.frame(BetaMatrixPCAind)
BetaMatrixPCAind$SAMPLE_ID <- row.names(BetaMatrixPCAind)
# Retain only the first two dimensions
BetaMatrixPCAind <- BetaMatrixPCAind %>%
  select(SAMPLE_ID, BetaMatrixPCAindDim1 = Dim.1, BetaMatrixPCAindDim2 = Dim.2)

BetaMatrixPCAind <- BetaMatrixPCAind %>%
  left_join(BetaMatrixPCAind) %>%
  left_join(ClinicalFile_T1[,c("SAMPLE_ID", "TYPE", "INSTITUTION", "BATCH")])

Beta_plots_to_save <- list()

Beta_plots_to_save[[1]] <- BetaMatrixPCAind %>%
  ggplot(aes(BetaMatrixPCAindDim1, BetaMatrixPCAindDim2, col = TYPE)) +
  labs(x = "Dim 1 (14.8%)", y = "Dim 2 (12.7%)") +
  geom_point(size = 2) +
  # geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  # stat_cor() +
  ggtitle("PCA plot by TYPE based on Beta values") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

Beta_plots_to_save[[2]] <- BetaMatrixPCAind %>%
  ggplot(aes(BetaMatrixPCAindDim1, BetaMatrixPCAindDim2, col = TYPE)) +
  labs(x = "Dim 1 (14.8%)", y = "Dim 2 (12.7%)") +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  stat_cor() +
  ggtitle("PCA plot by TYPE based on Beta values") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

Beta_plots_to_save[[3]] <- BetaMatrixPCAind %>%
  ggplot(aes(BetaMatrixPCAindDim1, BetaMatrixPCAindDim2, col = INSTITUTION)) +
  labs(x = "Dim 1 (14.8%)", y = "Dim 2 (12.7%)") +
  geom_point(size = 2) +
  #geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  #stat_cor() +
  ggtitle("PCA plot by INSTITUTION based on Beta values") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


Beta_plots_to_save[[4]] <- BetaMatrixPCAind %>%
  ggplot(aes(BetaMatrixPCAindDim1, BetaMatrixPCAindDim2, col = INSTITUTION)) +
  labs(x = "Dim 1 (14.8%)", y = "Dim 2 (12.7%)") +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  stat_cor() +
  ggtitle("PCA plot by INSTITUTION based on Beta values") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


Beta_plots_to_save[[5]] <- BetaMatrixPCAind %>%
  ggplot(aes(BetaMatrixPCAindDim1, BetaMatrixPCAindDim2, col = BATCH),
         main = "PCA plot") +
  labs(x = "Dim 1 (14.8%)", y = "Dim 2 (12.7%)") +
  geom_point(size = 2) +
  #geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  #stat_cor() +
  ggtitle("PCA plot by BATCH based on Beta values") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


Beta_plots_to_save[[6]] <- BetaMatrixPCAind %>%
  ggplot(aes(BetaMatrixPCAindDim1, BetaMatrixPCAindDim2, col = BATCH),
         main = "PCA plot") +
  labs(x = "Dim 1 (14.8%)", y = "Dim 2 (12.7%)") +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  stat_cor() +
  ggtitle("PCA plot by BATCH based on Beta values") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


# PCA analysis based on Mvalues
MvalueMatrixPCA <- stats::prcomp(t(MvalueMatrix_T1), scale = TRUE)
summary(MvalueMatrixPCA) #  PC1 explains 17.4% of the total variance
# PC2 explains 12.6% of total variance
fviz_eig(MvalueMatrixPCA, addlabels=TRUE, hjust = -0.3) + ylim(0, 95)


MvalueMatrixPCAind <- factoextra::get_pca_ind(MvalueMatrixPCA)
MvalueMatrixPCAind <- MvalueMatrixPCAind$coord
MvalueMatrixPCAind <- data.frame(MvalueMatrixPCAind)
MvalueMatrixPCAind$SAMPLE_ID <- row.names(MvalueMatrixPCAind)
# Retain only the first two dimensions
MvalueMatrixPCAind <- MvalueMatrixPCAind %>%
  select(SAMPLE_ID, MvalueMatrixPCAindDim1 = Dim.1, MvalueMatrixPCAindDim2 = Dim.2)

MvalueMatrixPCAind <- MvalueMatrixPCAind %>%
  left_join(MvalueMatrixPCAind) %>%
  left_join(ClinicalFile_T1[,c("SAMPLE_ID", "TYPE", "INSTITUTION", "BATCH")])

Mvalue_plots_to_save <- list()

Mvalue_plots_to_save[[1]] <- MvalueMatrixPCAind %>%
  ggplot(aes(MvalueMatrixPCAindDim1, MvalueMatrixPCAindDim2, col = TYPE)) +
  labs(x = "Dim 1 (17.4%)", y = "Dim 2 (12.6%)") +
  geom_point(size = 2) +
  # geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  # stat_cor() +
  ggtitle("PCA plot by TYPE based on Mvalues") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

Mvalue_plots_to_save[[2]] <- MvalueMatrixPCAind %>%
  ggplot(aes(MvalueMatrixPCAindDim1, MvalueMatrixPCAindDim2, col = TYPE)) +
  labs(x = "Dim 1 (17.4%)", y = "Dim 2 (12.6%)") +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  stat_cor() +
  ggtitle("PCA plot by TYPE based on Mvalues") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

Mvalue_plots_to_save[[3]] <- MvalueMatrixPCAind %>%
  ggplot(aes(MvalueMatrixPCAindDim1, MvalueMatrixPCAindDim2, col = INSTITUTION)) +
  labs(x = "Dim 1 (17.4%)", y = "Dim 2 (12.6%)") +
  geom_point(size = 2) +
  #geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  #stat_cor() +
  ggtitle("PCA plot by INSTITUTION based on Mvalues") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


Mvalue_plots_to_save[[4]] <- MvalueMatrixPCAind %>%
  ggplot(aes(MvalueMatrixPCAindDim1, MvalueMatrixPCAindDim2, col = INSTITUTION)) +
  labs(x = "Dim 1 (17.4%)", y = "Dim 2 (12.6%)") +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  stat_cor() +
  ggtitle("PCA plot by INSTITUTION based on Mvalues") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


Mvalue_plots_to_save[[5]] <- MvalueMatrixPCAind %>%
  ggplot(aes(MvalueMatrixPCAindDim1, MvalueMatrixPCAindDim2, col = BATCH)) +
  labs(x = "Dim 1 (17.4%)", y = "Dim 2 (12.6%)") +
  geom_point(size = 2) +
  # geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  # stat_cor() +
  ggtitle("PCA plot by BATCH based on Mvalues") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

Mvalue_plots_to_save[[6]] <- MvalueMatrixPCAind %>%
  ggplot(aes(MvalueMatrixPCAindDim1, MvalueMatrixPCAindDim2, col = BATCH)) +
  labs(x = "Dim 1 (17.4%)", y = "Dim 2 (12.6%)") +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  stat_cor() +
  ggtitle("PCA plot by BATCH based on Mvalues") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# TYPE
# Differentially methylated regions for FL vs RLN (across ALL probes)
Diff_MethylatedRegions_FLvRLN <- Diff_MethylatedRegions(Method = "DMRcate", 
                                                        BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) != "DL")], 
                                                        MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(MvalueMatrix_T1), 4, 5) != "DL")], 
                                                        ContrastColumnName = "TYPE", 
                                                        ClusterLabels = NA, 
                                                        ClinicalFile = ClinicalFile_T1[which(substr(colnames(MvalueMatrix_T1), 4, 5) != "DL"), ], 
                                                        AnnotationFile = AnnotationFile, 
                                                        ProduceImages = "No", DMR = 1, PNGorPDF = "png", ExpressionFile = NA)


# Differentially methylated regions for FL vs DLBCL (across ALL probes)
Diff_MethylatedRegions_FLvRLN <- Diff_MethylatedRegions(Method = "DMRcate", 
                                                        BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) != "RL")], 
                                                        MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(MvalueMatrix_T1), 4, 5) != "RL")], 
                                                        ContrastColumnName = "TYPE", 
                                                        ClusterLabels = NA, 
                                                        ClinicalFile = ClinicalFile_T1[which(substr(colnames(MvalueMatrix_T1), 4, 5) != "RL"), ], 
                                                        AnnotationFile = AnnotationFile, 
                                                        ProduceImages = "Yes", DMR = 1, PNGorPDF = "png", ExpressionFile = NA)

# Differentially methylated regions for RLN vs DLBCL (across ALL probes)
Diff_MethylatedRegions_FLvRLN <- Diff_MethylatedRegions(Method = "DMRcate", 
                                                        BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) != "FL")], 
                                                        MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(MvalueMatrix_T1), 4, 5) != "FL")], 
                                                        ContrastColumnName = "TYPE", 
                                                        ClusterLabels = NA, 
                                                        ClinicalFile = ClinicalFile_T1[which(substr(colnames(MvalueMatrix_T1), 4, 5) != "FL"), ], 
                                                        AnnotationFile = AnnotationFile, 
                                                        ProduceImages = "Yes", DMR = 1, PNGorPDF = "png", ExpressionFile = NA)


Output_Differential_probes_TYPE <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1, 
                                                           MvalueMatrix = MvalueMatrix_T1, 
                                                           BetaMatrix = BetaMatrix_T1,
                                                           ProbeOrGene = "Probe", 
                                                           ContrastColumnName = "TYPE", 
                                                           RGChannelSet = Output_Data_1$RGChannelSet, 
                                                           SampleSheet = Output_Data_1$SampleSheet, 
                                                           ProduceImages = "No", PNGorPDF = "png")


# STAGE
# Differentially methylated regions for FL (across ALL probes)
Diff_MethylatedRegions_LimAdv <- Diff_MethylatedRegions(Method = "DMRcate", 
                                                        BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                                        MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(MvalueMatrix_T1), 4, 5) == "FL")], 
                                                        ContrastColumnName = "STAGE", 
                                                        ClusterLabels = NA, 
                                                        ClinicalFile = ClinicalFile_T1[which(substr(colnames(MvalueMatrix_T1), 4, 5) == "FL"), ], 
                                                        AnnotationFile = AnnotationFile, 
                                                        ProduceImages = "No", DMR = 1, PNGorPDF = "png", ExpressionFile = NA)


Output_Differential_probes_LimAdv <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1[which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL"), ], 
                                                             MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                                             BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")],
                                                             ProbeOrGene = "Probe", 
                                                             ContrastColumnName = "STAGE", 
                                                             RGChannelSet = Output_Data_1$RGChannelSet, 
                                                             SampleSheet = Output_Data_1$SampleSheet, 
                                                             ProduceImages = "No", PNGorPDF = "png")

# islands
ContrastColumnName = "STAGE"
BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")]
MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")]
ClinicalFile = ClinicalFile_T1[which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL"), ]
# Matching rownames in methylation file with corresponding location in annotation file
match_ids_methylation_beta <- match(rownames(BetaMatrix), AnnotationFile$V1)

# Bringing annotation file in the same order as BetaMatrix
AnnotationFile2 <- data.frame(AnnotationFile[match_ids_methylation_beta, ])


# Identify column number matching the CategoryToVisualize
columnNumber = as.numeric(which(colnames(AnnotationFile2) == "Relation_to_Island"))

# Only take in the non empty entries of the selected column
non_empty_entries <- which(AnnotationFile2[, columnNumber] != "")

# rows corresponding to the current AnnotationCategory, e.g. Island from the list Island  N_Shore OpenSea N_Shelf S_Shelf S_Shore
correspondingRows <- which(AnnotationFile2[non_empty_entries, columnNumber] == "S_Shore")

MvalueMatrix = MvalueMatrix_T1[correspondingRows, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")]
BetaMatrix = BetaMatrix_T1[correspondingRows, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")]
# now go to 7_DifferentialMethylation.R


# Infinium Clustering
Output_SDeviation_0.3_18 <- SDeviation(BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, standard_deviation = 0.3)
Output_SDeviation_0.25_18 <- SDeviation(BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, standard_deviation = 0.25)
Output_SDeviation_0.23_18 <- SDeviation(BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, standard_deviation = 0.23)
Output_SDeviation_0.20_18 <- SDeviation(BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, standard_deviation = 0.2)
Output_SDeviation_0.18_18 <- SDeviation(BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, standard_deviation = 0.18)
Output_SDeviation_0.16_18 <- SDeviation(BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, standard_deviation = 0.16)

library(InfiniumPurify)
tumor.type <- "DLBC"
data(iDMC) 
probes <- iDMC[[tumor.type]]
probes.true <- names(probes[probes == T])
beta.sel <- BetaMatrix_T1[row.names(BetaMatrix_T1) %in% probes.true,]
# data(beta.emp) # An example data set for InfiniumClust and InfiniumPurify.
# purity_OICR_withNormal <- data.frame(purity = getPurity(tumor.data = BetaMatrix_T1, tumor.type = tumor.type))
purity_OICR_withNormal <- InfiniumPurify::getPurity(tumor.data = beta.sel, tumor.type = tumor.type)

# All samples 
set.seed(1234)
Output_SDeviation_0.25_18_clust <- Clustering(TypeofClustering = "InfiniumClust", Purity = purity_OICR_withNormal, 
                                              BetaMatrix = BetaMatrix_T1, MvalueMatrix = MvalueMatrix_T1, 
                                              ListofProbes = rownames(Output_SDeviation_0.25_18$BetaMatrix_SD_Filtered), 
                                              ClinicalFile = ClinicalFile_T1, FigureGenerate = "No", 
                                              PNGorPDF = "png", ImageName = "InfiniumClust_point3SD") 



# FL only samples 
Output_FLonly_SDeviation_0.3_18 <- SDeviation(BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                              MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                              standard_deviation = 0.3)
Output_FLonly_SDeviation_0.25_18 <- SDeviation(BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                               MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                               standard_deviation = 0.25)
Output_FLonly_SDeviation_0.23_18 <- SDeviation(BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                               MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                               standard_deviation = 0.23)
Output_FLonly_SDeviation_0.20_18 <- SDeviation(BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                               MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                               standard_deviation = 0.2)
Output_FLonly_SDeviation_0.18_18 <- SDeviation(BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")],
                                               MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")],
                                               standard_deviation = 0.18)
Output_FLonly_SDeviation_0.16_18 <- SDeviation(BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                               MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                               standard_deviation = 0.16)


# clustering FL only samp,es 
Output_FLonly_SDeviation_0.25_18_clust <- Clustering(TypeofClustering = "InfiniumClust", Purity = purity_OICR_withNormal[c(11:165)], 
                                                     BetaMatrix = BetaMatrix_T1[, c(11:165)], MvalueMatrix = MvalueMatrix_T1[, c(11:165)], 
                                                     ListofProbes = rownames(Output_FLonly_SDeviation_0.25_18$BetaMatrix_SD_Filtered), 
                                                     ClinicalFile = ClinicalFile_T1[c(11:165), ], FigureGenerate = "No", 
                                                     PNGorPDF = "png", ImageName = "None") 
set.seed(1234)
Output_FLonly_SDeviation_0.23_18_clust <- Clustering(TypeofClustering = "InfiniumClust", Purity = purity_OICR_withNormal[c(11:165)], 
                                                     BetaMatrix = BetaMatrix_T1[, c(11:165)], MvalueMatrix = MvalueMatrix_T1[, c(11:165)], 
                                                     ListofProbes = rownames(Output_FLonly_SDeviation_0.23_18$BetaMatrix_SD_Filtered), 
                                                     ClinicalFile = ClinicalFile_T1[c(11:165), ], FigureGenerate = "No", 
                                                     PNGorPDF = "png", ImageName = "None") 
set.seed(1234)
Output_FLonly_SDeviation_0.20_18_clust <- Clustering(TypeofClustering = "InfiniumClust", Purity = purity_OICR_withNormal[c(11:165)], 
                                                     BetaMatrix = BetaMatrix_T1[, c(11:165)], MvalueMatrix = MvalueMatrix_T1[, c(11:165)], 
                                                     ListofProbes = rownames(Output_FLonly_SDeviation_0.20_18$BetaMatrix_SD_Filtered), 
                                                     ClinicalFile = ClinicalFile_T1[c(11:165), ], FigureGenerate = "No", 
                                                     PNGorPDF = "png", ImageName = "None") 
set.seed(1234)
Output_FLonly_SDeviation_0.18_18_clust <- Clustering(TypeofClustering = "InfiniumClust", Purity = purity_OICR_withNormal[c(11:165)], 
                                                     BetaMatrix = BetaMatrix_T1[, c(11:165)], MvalueMatrix = MvalueMatrix_T1[, c(11:165)], 
                                                     ListofProbes = rownames(Output_FLonly_SDeviation_0.18_18$BetaMatrix_SD_Filtered), 
                                                     ClinicalFile = ClinicalFile_T1[c(11:165), ], FigureGenerate = "No", 
                                                     PNGorPDF = "png", ImageName = "None")
set.seed(1234)
Output_FLonly_SDeviation_0.16_18_clust <- Clustering(TypeofClustering = "InfiniumClust", Purity = purity_OICR_withNormal[c(11:165)], 
                                                     BetaMatrix = BetaMatrix_T1[, c(11:165)], MvalueMatrix = MvalueMatrix_T1[, c(11:165)], 
                                                     ListofProbes = rownames(Output_FLonly_SDeviation_0.16_18$BetaMatrix_SD_Filtered), 
                                                     ClinicalFile = ClinicalFile_T1[c(11:165), ], FigureGenerate = "No", 
                                                     PNGorPDF = "png", ImageName = "None")


# 11 March 2020
ViolinPlotsMehtylation(BetaMatrix = BetaMatrix_T1, 
                       MvalueMatrix = MvalueMatrix_T1, 
                       ClinicalFile = ClinicalFile_T1, 
                       CategoryToVisualize = "TYPE")


# 8 April 2020
cat("\n Running SD of 0.25 with no sex xsomes, no SNPs, of FL only samples =  214  probes")
Output_point2SD_All <- SDeviation(BetaMatrix = BetaMatrix_T1, 
                                  MvalueMatrix = MvalueMatrix_T1, 
                                  CutOff = 0.25)
Output_Clustering_point2SD_All <- Clustering(BetaMatrix = BetaMatrix_T1, 
                                             MvalueMatrix = MvalueMatrix_T1, 
                                             AnnotationFile = AnnotationFile,
                                             ListofProbes = rownames(Output_point2SD_All$BetaMatrix_SD_Filtered), 
                                             ClinicalFile = ClinicalFile_T1, 
                                             TumorPurity = TumorPurity,
                                             FigureGenerate = "No", 
                                             PNGorPDF = "png")

# Pick top 1,000 most variable probes
Output_SDeviation_1000 <- SDeviation(BetaMatrix = BetaMatrix_T1, 
                                     MvalueMatrix = MvalueMatrix_T1, 
                                     NumberofProbes = 1000)
Output_Clustering_1000 <- Clustering(BetaMatrix = BetaMatrix_T1, 
                                     MvalueMatrix = MvalueMatrix_T1, 
                                     AnnotationFile = AnnotationFile,
                                     ListofProbes = rownames(Output_SDeviation_1000$BetaMatrix_SD_Filtered), 
                                     ClinicalFile = ClinicalFile_T1, 
                                     TumorPurity = TumorPurity,
                                     FigureGenerate = "No", 
                                     PNGorPDF = "png")
Output_Clustering_1000_lables #<- rpmmClass_numbers

# Pick top 5,000 most variable probes
Output_SDeviation_5000_All <- SDeviation(BetaMatrix = BetaMatrix_T1, 
                                         MvalueMatrix = MvalueMatrix_T1, 
                                         NumberofProbes = 5000)
Output_Clustering_5000_All <- Clustering(BetaMatrix = BetaMatrix_T1, 
                                         MvalueMatrix = MvalueMatrix_T1, 
                                         AnnotationFile = AnnotationFile,
                                         ListofProbes = rownames(Output_SDeviation_5000_All$BetaMatrix_SD_Filtered), 
                                         ClinicalFile = ClinicalFile_T1, 
                                         TumorPurity = TumorPurity,
                                         FigureGenerate = "No", 
                                         PNGorPDF = "png")


# Pick top 10,000 most variable probes
Output_SDeviation_10000_All <- SDeviation(BetaMatrix = BetaMatrix_T1, 
                                          MvalueMatrix = MvalueMatrix_T1, 
                                          NumberofProbes = 10000)
Output_Clustering_10000_All <- Clustering(BetaMatrix = BetaMatrix_T1, 
                                          MvalueMatrix = MvalueMatrix_T1, 
                                          AnnotationFile = AnnotationFile,
                                          ListofProbes = rownames(Output_SDeviation_10000_All$BetaMatrix_SD_Filtered), 
                                          ClinicalFile = ClinicalFile_T1, 
                                          TumorPurity = TumorPurity,
                                          FigureGenerate = "No", 
                                          PNGorPDF = "png")

# Pick top 20,000 most variable probes
Output_SDeviation_20000_All <- SDeviation(BetaMatrix = BetaMatrix_T1, 
                                          MvalueMatrix = MvalueMatrix_T1, 
                                          NumberofProbes = 20000)
Output_Clustering_20000_All <- Clustering(BetaMatrix = BetaMatrix_T1, 
                                          MvalueMatrix = MvalueMatrix_T1, 
                                          AnnotationFile = AnnotationFile,
                                          ListofProbes = rownames(Output_SDeviation_20000_All$BetaMatrix_SD_Filtered), 
                                          ClinicalFile = ClinicalFile_T1, 
                                          TumorPurity = TumorPurity,
                                          FigureGenerate = "No", 
                                          PNGorPDF = "png")




# FL only samples
# Pick top 1,000 most variable probes
Output_SDeviation_1000_FLonly <- SDeviation(BetaMatrix = BetaMatrix_T1[, 11:165], 
                                            MvalueMatrix = MvalueMatrix_T1[, 11:165], 
                                            NumberofProbes = 1000)
Output_Clustering_1000_FLonly <- Clustering(BetaMatrix = BetaMatrix_T1[, 11:165], 
                                            MvalueMatrix = MvalueMatrix_T1[, 11:165], 
                                            AnnotationFile = AnnotationFile,
                                            ListofProbes = rownames(Output_SDeviation_1000_FLonly$BetaMatrix_SD_Filtered), 
                                            ClinicalFile = ClinicalFile_T1[11:165, ], 
                                            TumorPurity = TumorPurity[11:165, ],
                                            FigureGenerate = "No", 
                                            PNGorPDF = "png")
Output_SDeviation_1000_FLonly_labels #<- rpmmClass_numbers

# FL only samples
# Pick top 5,000 most variable probes
Output_SDeviation_5000_FLonly <- SDeviation(BetaMatrix = BetaMatrix_T1[, 11:165], 
                                            MvalueMatrix = MvalueMatrix_T1[, 11:165], 
                                            NumberofProbes = 5000)
Output_Clustering_5000_FLonly <- Clustering(BetaMatrix = BetaMatrix_T1[, 11:165], 
                                            MvalueMatrix = MvalueMatrix_T1[, 11:165], 
                                            AnnotationFile = AnnotationFile,
                                            ListofProbes = rownames(Output_SDeviation_5000_FLonly$BetaMatrix_SD_Filtered), 
                                            ClinicalFile = ClinicalFile_T1[11:165, ], 
                                            TumorPurity = TumorPurity[11:165, ],
                                            FigureGenerate = "No", 
                                            PNGorPDF = "png")
# FL only samples
# Pick top 10,000 most variable probes
Output_SDeviation_10000_FLonly <- SDeviation(BetaMatrix = BetaMatrix_T1[, 11:165], 
                                             MvalueMatrix = MvalueMatrix_T1[, 11:165], 
                                             NumberofProbes = 10000)
Output_Clustering_10000_FLonly <- Clustering(BetaMatrix = BetaMatrix_T1[, 11:165], 
                                             MvalueMatrix = MvalueMatrix_T1[, 11:165], 
                                             AnnotationFile = AnnotationFile,
                                             ListofProbes = rownames(Output_SDeviation_10000_FLonly$BetaMatrix_SD_Filtered), 
                                             ClinicalFile = ClinicalFile_T1[11:165, ], 
                                             TumorPurity = TumorPurity[11:165, ],
                                             FigureGenerate = "No", 
                                             PNGorPDF = "png")
# FL only samples
# Pick top 20,000 most variable probes
Output_SDeviation_20000_FLonly <- SDeviation(BetaMatrix = BetaMatrix_T1[, 11:165], 
                                             MvalueMatrix = MvalueMatrix_T1[, 11:165], 
                                             NumberofProbes = 20000)
Output_Clustering_20000_FLonly <- Clustering(BetaMatrix = BetaMatrix_T1[, 11:165], 
                                             MvalueMatrix = MvalueMatrix_T1[, 11:165], 
                                             AnnotationFile = AnnotationFile,
                                             ListofProbes = rownames(Output_SDeviation_20000_FLonly$BetaMatrix_SD_Filtered), 
                                             ClinicalFile = ClinicalFile_T1[11:165, ], 
                                             TumorPurity = TumorPurity[11:165, ],
                                             FigureGenerate = "No", 
                                             PNGorPDF = "png")





# Pick  most variable probes with >0.25 SD = 203 probes
Output_point25_FLonly <- SDeviation(BetaMatrix = BetaMatrix_T1[, 11:165], 
                                    MvalueMatrix = MvalueMatrix_T1[, 11:165], 
                                    CutOff = 0.2)
Output_Clustering_point25_FLonly  <- Clustering(BetaMatrix = BetaMatrix_T1[, 11:165], 
                                                MvalueMatrix = MvalueMatrix_T1[, 11:165], 
                                                AnnotationFile = AnnotationFile,
                                                ListofProbes = rownames(Output_point25_FLonly$BetaMatrix_SD_Filtered), 
                                                ClinicalFile = ClinicalFile_T1[11:165, ], 
                                                TumorPurity = TumorPurity[11:165, ],
                                                FigureGenerate = "No", 
                                                PNGorPDF = "png")


# 22 April 2020 - Robert recommended to go back to InfiniumClust

# Pick top 10,000 most variable probes
Output_SDeviation_10000_All <- SDeviation(BetaMatrix = BetaMatrix_T1, 
                                          MvalueMatrix = MvalueMatrix_T1, 
                                          NumberofProbes = 10000)

# calculating purity
library(InfiniumPurify)
data(iDMC) 
probes <- iDMC[["DLBC"]]
probes.true <- names(probes[probes == T])
beta.sel <- BetaMatrix_T1[row.names(BetaMatrix_T1) %in% probes.true,]
purity_OICR_withNormal <- getPurity(tumor.data = beta.sel, tumor.type = "DLBC")

TumorPurity10Jan2020 <- readRDS(file = "Purity_281probes_10Jan2020.rds")
dim(TumorPurity10Jan2020) # 170   1

# Pick top 5,000 most variable probes
Output_SDeviation_5000_All <- SDeviation(BetaMatrix = BetaMatrix_T1, 
                                         MvalueMatrix = MvalueMatrix_T1, 
                                         NumberofProbes = 5000)


Output_InfiniumClustering_5000_All <- Clustering(TypeofClustering = "InfiniumClust",
                                                 BetaMatrix = BetaMatrix_T1, 
                                                 MvalueMatrix = MvalueMatrix_T1, 
                                                 AnnotationFile = AnnotationFile,
                                                 ListofProbes = rownames(Output_SDeviation_5000_All$BetaMatrix_SD_Filtered), 
                                                 ClinicalFile = ClinicalFile_T1, 
                                                 TumorPurity = purity_OICR_withNormal,
                                                 FigureGenerate = "No", 
                                                 PNGorPDF = "png")

# Trying to using RNAseq data for clusters  
# These were generated in SNF_ASilva_10Feb2020_extended
# RNASeqCountMatrixMatched <- readRDS(file = paste0("RNASeqCountMatrix_Integrative.rds"))
# dim(RNASeqCountMatrixMatched) # 57820   132
# MvalueMatrixMatched <- readRDS(file = paste0("MvalueMatrix_Integrative.rds"))
# dim(MvalueMatrixMatched) # 595564    132
# ClinicalFileMatched <- readRDS(file = paste0("ClinicalFile_Integrative.rds"))
# dim(ClinicalFileMatched) # 132  25
# BetaMatrixMatched <- readRDS(file = paste0("BetaMatrix_Integrative.rds"))
# dim(BetaMatrixMatched) # 595564  132
# QCcombinedMatched <- readRDS(file = paste0("QCcombined_Integrative.rds"))
# dim(QCcombinedMatched) # 132 12
TumorPurity <- readRDS(file = paste0("Purity_281probes_10Jan2020.rds"))
PurityMatched <- as.matrix(TumorPurity[match(colnames(MvalueMatrixMatched), rownames(TumorPurity)), ])
rownames(PurityMatched) <- rownames(TumorPurity)[match(colnames(MvalueMatrixMatched), rownames(TumorPurity))]
dim(PurityMatched) # 132   1

set.seed(1234)
ClusteringSNF <- SNFClustering(RNAseqCountMatrix = RNASeqCountMatrixMatched, 
                               MvalueMatrix = MvalueMatrixMatched,
                               BetaMatrix = BetaMatrixMatched,
                               AnnotationFile = AnnotationFileEPIC,
                               TumorPurity = PurityMatched, 
                               ClinicalFile = ClinicalFileMatched, 
                               QCMatrix = QCcombinedMatched,
                               SurvivalFile = SurvivalFile,
                               CufoffUniqMapReads = 30000000, 
                               CutoffRRNAcontam = 40,
                               MethylationFeatureSelectionMethod = "MAD",
                               MethylationFeatureSelectionCutOff = NA,
                               MethylationFeatureSelectionNumberofProbes = 595564,
                               MethylationNormalizationMethod = NA, 
                               RNAseqFeatureSelectionMethod = "edgeR",
                               RNAseqFeatureSelectionCutOff = NA,
                               RNAseqFeatureSelectionNumberofProbes = NA,
                               RNAseqNormalizationMethod = NA,
                               NumberOfClusters = 1,
                               ImageName = "1ClusterModel",
                               ProduceImages = "No")

RNAseqCountsAllSamplesFilteredProbes <- readRDS(file = "33_RNASeqCountMatrixMatchedAllSamplesFilteredProbes.rds")
dim(RNAseqCountsAllSamplesFilteredProbes) # 41220   132


Match132with170 <- match(colnames(RNAseqCountsAllSamplesFilteredProbes), colnames(BetaMatrix_T1))
ClusterLabels2to4 <- readRDS(file = "InfiniumClustering2to4.rds")
ClusterLabels2to4[[2]][[2]][Match132with170]

RNASeqCountMatrixMatchedMethylation <- readRDS(file = "RNASeqCountMatrixMatchedMethylationProbes.rds")
# several zero column sums detected, remove those
which(colSums(RNASeqCountMatrixMatchedMethylation) == 0)

DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = as.matrix(RNASeqCountMatrixMatchedMethylation[, -which(colSums(RNASeqCountMatrixMatchedMethylation) == 0)]), 
                                                                   ContrastColumnName = "CLUSTER", 
                                                                   ClinicalFile = ClinicalFileMatched, 
                                                                   ClusterLabels = (ClusterLabels2to4[[2]][[2]][Match132with170])[-which(colSums(RNASeqCountMatrixMatchedMethylation) == 0)], 
                                                                   ProduceImages = "Yes", 
                                                                   PNGorPDF = "png") 
length(DifferentialExpressionRNAseqOutput$AllDEgenes) #  
#   logFCcutoff <- 2 # logFC; FDRcutoff <- 0.05 # FDR
# Contrast: one-two has 87 ensembl IDs (+ve)logFC and significant.
# Contrast: one-two has 40 ensembl IDs (-ve)logFC and significant

ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")
matchRNAseq <- match(substr(DifferentialExpressionRNAseqOutput$AllDEIDs$ensemblGeneId, 1, 15), ENSMBLid$gene)
DEgeneTable <- cbind(DifferentialExpressionRNAseqOutput$AllDEIDs, ENSMBLid$name[matchRNAseq])
length(unique(ENSMBLid$name[matchRNAseq])) # 55 unique genes
colnames(DEgeneTable)[7] <- "ENSMBLgene"
write.csv(DEgeneTable, file = "DEgenes_InfiniumClustContrast_14March2020.csv")


RNAseqDEgenes <- unique(ENSMBLid$name[matchRNAseq])



# Summarizing clinical characteristics by cluster
# FLIPI
median(SurvivalFile_OrderedbyBetaMatrixPatients$AGE_AT_DIAGNOSIS[which(ClusterLabels2to4[[2]][[2]][which(ClinicalFile_T1$TYPE == "FL")] == 1)])
median(SurvivalFile_OrderedbyBetaMatrixPatients$AGE_AT_DIAGNOSIS[which(ClusterLabels2to4[[2]][[2]][which(ClinicalFile_T1$TYPE == "FL")] == 2)])

table(SurvivalFile_OrderedbyBetaMatrixPatients$B_SYMPTOMS[which(ClusterLabels2to4[[2]][[2]][which(ClinicalFile_T1$TYPE == "FL")] == 1)])
table(SurvivalFile_OrderedbyBetaMatrixPatients$B_SYMPTOMS[which(ClusterLabels2to4[[2]][[2]][which(ClinicalFile_T1$TYPE == "FL")] == 2)])


table(SurvivalFile_OrderedbyBetaMatrixPatients$FLIPI_BINARY[which(ClusterLabels2to4[[2]][[2]][which(ClinicalFile_T1$TYPE == "FL")] == 1)])
table(SurvivalFile_OrderedbyBetaMatrixPatients$FLIPI_BINARY[which(ClusterLabels2to4[[2]][[2]][which(ClinicalFile_T1$TYPE == "FL")] == 2)])


# 3 May 2020
Output_Differential_probes_LimAdv <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1[which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL"), ], 
                                                             MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")], 
                                                             BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4, 5) == "FL")],
                                                             ProbeOrGene = "Probe", 
                                                             ContrastColumnName = "STAGE", 
                                                             RGChannelSet = Output_Data_1$RGChannelSet, 
                                                             SampleSheet = Output_Data_1$SampleSheet, 
                                                             ProduceImages = "No", PNGorPDF = "png")

dim(top_P_Coef11) #  179659      probes show difference between stages
matchDiffProbesSTAGEBeta <- match(rownames(top_P_Coef11), rownames(BetaMatrix_T1))
dim(BetaMatrix_T1[matchDiffProbesSTAGEBeta, ]) # 179659    170


Output_InfiniumClustering_179659_All <- Clustering(TypeofClustering = "InfiniumClust",
                                                   BetaMatrix = BetaMatrix_T1, 
                                                   MvalueMatrix = MvalueMatrix_T1, 
                                                   AnnotationFile = AnnotationFile,
                                                   ListofProbes = rownames(top_P_Coef11), 
                                                   ClinicalFile = ClinicalFile_T1, 
                                                   TumorPurity = purity_OICR_withNormal,
                                                   FigureGenerate = "No", 
                                                   PNGorPDF = "png")

# check if 5000 (SD) probes are found in 179659 STAGE probes
match5000with179659 <- match(rownames(Output_SDeviation_5000_All$BetaMatrix_SD_Filtered), rownames(top_P_Coef11))
length(which(is.na(rownames(top_P_Coef11)[match5000with179659]) == TRUE))
#5000-1899 = 3101 entries match between the two



Output_Visuals_Relation_to_Island_STAGE_2 <- ProportionVisualization(CategoryToVisualize = "UCSC_RefGene_Group", 
                                                                     BetaMatrix = BetaMatrix_T1, 
                                                                     AnnotationFile = AnnotationFile, 
                                                                     ClinicalFile = ClinicalFile_T1, 
                                                                     ClinicalCategory = "STAGE", 
                                                                     PlotWithinCategories="Yes", 
                                                                     FigureGenerate="Yes", 
                                                                     PNGorPDF="png", 
                                                                     ImageName = "UCSC_RefGene_Group_TYPE_9May2020")

BoxPlot(BetaMatrix = BetaMatrix_T1, 
        ClinicalFile  = ClinicalFile_T1, 
        CategoryToVisualize  = "TYPE",
        PNGorPDF = "png")

ViolinPlotsMethylation(BetaMatrix = BetaMatrix_T1, 
                       MvalueMatrix = MvalueMatrix_T1, 
                       ClinicalFile = ClinicalFile_T1, 
                       CategoryToVisualize = "STAGE", 
                       ClusterLabels = NA,
                       ProduceImages = "Yes",
                       PNGorPDF = "png") 

# ------------------------------------------ -
# 11 May 2020

DiffMethylatedRegionsInfiniumClusters <- DiffMethylatedRegions(Method = "DMRcate", 
                                                               BetaMatrix = BetaMatrix_T1, 
                                                               MvalueMatrix = MvalueMatrix_T1, 
                                                               ContrastColumnName = "CLUSTER", 
                                                               ClinicalFile = ClinicalFile_T1, 
                                                               ClusterLabels = ClusterLabels2to4[[2]][[2]],
                                                               AnnotationFile = AnnotationFile, 
                                                               ProduceImages = "No", 
                                                               DMR = 1, ExpressionFile = NA)


# 12 May 2020
# 17646 unique genes from 31740 DMRs were identified based on output of DMRCate  
# See this with RNAseq data
length(DiffMethylatedRegionsInfiniumClusters$GenesFallingInRegions$`one-two`)
DMRcateGenes <- DiffMethylatedRegionsInfiniumClusters$GenesFallingInRegions$`one-two`


ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")
matchRNAseq <- match(substr(rownames(RNASeqCountMatrixMatched), 1, 15), ENSMBLid$gene)

# match DMRcate genes with RNAseq genes
matchRNAseqDMRcate <- match(DMRcateGenes, ENSMBLid$name[matchRNAseq])
RNAseqGenesFromDMRcate <- ENSMBLid$name[matchRNAseq][matchRNAseqDMRcate[!is.na(matchRNAseqDMRcate)]]
# length(RNAseqGenesFromDMRcate) # 16134

# get RNAseq counts for methylation genes 
match16134 <- match(RNAseqGenesFromDMRcate, ENSMBLid$name)
matchRNAseqMethylationGenes <- match(ENSMBLid$gene[match16134], 
                                     substr(rownames(RNASeqCountMatrixMatched), 1, 15))
# matchRNAseqMethylationGenes[!is.na(matchRNAseqMethylationGenes)] # remove NA values
dim(RNASeqCountMatrixMatched[matchRNAseqMethylationGenes[!is.na(matchRNAseqMethylationGenes)], ]) # 16133  132
RNASeqCountMatrixMatchedMethylation <- RNASeqCountMatrixMatched[matchRNAseqMethylationGenes[!is.na(matchRNAseqMethylationGenes)], ]
# saveRDS(RNASeqCountMatrixMatchedMethylation, file = "RNASeqCountMatrixMatchedMethylationProbes.rds")
RNASeqCountMatrixMatchedMethylationNorm <- SNFtool::standardNormalization(log2(1 + RNASeqCountMatrixMatchedMethylation))
# [,orderColumnsRNAseq]

# order patients by cluster 

Match132with170 <- match(colnames(RNASeqCountMatrixMatchedMethylationNorm), colnames(BetaMatrix_T1))

orderColumnsRNAseq <- c(which(ClusterLabels2to4[[2]][[2]][Match132with170] == 1),
                        which(ClusterLabels2to4[[2]][[2]][Match132with170] == 2))
Cluster = factor(ClusterLabels2to4[[2]][[2]][Match132with170])
Disease = factor(ClinicalFile$TYPE[Match132with170])
Stage = factor(ClinicalFile$STAGE[Match132with170])
Sex = factor(ClinicalFile$SEX[Match132with170])
Translocation = factor(ClinicalFile$TRANSLOC_14_18[Match132with170])
TypeBiopsy =  factor(ClinicalFile$TYPE_BIOPSY[Match132with170])
SiteBiopsy =  factor(ClinicalFile$SITE_BIOPSY[Match132with170])

annotation_colRNAseq = data.frame(
  Cluster = factor(Cluster[orderColumnsRNAseq]),
  Disease = factor(Disease[orderColumnsRNAseq]),
  Stage = factor(Stage[orderColumnsRNAseq]),
  Sex = factor(Sex[orderColumnsRNAseq]),
  Translocation = factor(Translocation[orderColumnsRNAseq]),
  TypeBiopsy = factor(TypeBiopsy[orderColumnsRNAseq]),
  SiteBiopsy = factor(SiteBiopsy[orderColumnsRNAseq]))
rownames(annotation_colRNAseq) = colnames(RNASeqCountMatrixMatchedMethylationNorm[, orderColumnsRNAseq])

ann_colorsRNAseq = list(
  Cluster = c("1" = "#4363d8", "2" = "#f58231"), 
  Disease = c(DLBCL = "#d6604d", FL = "#66bd63", RLN = "#4575b4"), #option 5
  Stage = c(ADVANCED = "#762a83", LIMITED = "#c2a5cf"),
  Sex = c("F"="#b35806", "M"="#fdb863"),
  Translocation = c("0"="#e0e0e0","1"="#878787"),
  TypeBiopsy = c("TISSUE" = "#a6dba0", "CORE"="#878787"),
  SiteBiopsy = c("LN" = "#f1b6da" , "EN" = "#c51b7d"))



## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# Plot the heatmap
# log transformed counts
pheatmap::pheatmap(log(RNASeqCountMatrixMatchedMethylation+0.01)[, orderColumnsRNAseq],
                   show_colnames = T, 
                   show_rownames = F,
                   fontface = "italic", 
                   legend = T, #scale ="row", 
                   annotation_colors = ann_colorsRNAseq, 
                   border_color = "black", 
                   cluster_cols = F,
                   cluster_rows = T,
                   annotation_col = annotation_colRNAseq, 
                   color =  rev(morecols(50)))

# normalized counts by SNFnormlized function
pheatmap::pheatmap(RNASeqCountMatrixMatchedMethylationNorm[c(1:5), orderColumnsRNAseq],
                   show_colnames = T, 
                   show_rownames = T,
                   fontface = "italic", 
                   legend = T, scale ="row", 
                   annotation_colors = ann_colorsRNAseq, 
                   border_color = "black", 
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   annotation_col = annotation_colRNAseq, 
                   color =  rev(morecols(50)))


# Compare RNAseq DE genes and DMCRcate DE genes
matchRNAseqDE_DMRcateDE <- match(RNAseqGenesFromDMRcate, RNAseqDEgenes)
length(which(!is.na(matchRNAseqDE_DMRcateDE))) #54

RNAseqDEgenes[matchRNAseqDE_DMRcateDE[which(!is.na(matchRNAseqDE_DMRcateDE))]]


# ------------------------------------------ -
# 15 May 2020
# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) 


# add clinical file with cluster labels
ClinicalFile_T1CLUSTERs <- cbind(ClinicalFile_T1, InfiniumClustLabels$Cluster)
colnames(ClinicalFile_T1CLUSTERs)[26] <-"CLUSTER"
dim(ClinicalFile_T1CLUSTERs) # 170  26
Output_Differential_probes_InfiniumClust <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1CLUSTERs, 
                                                                    MvalueMatrix = MvalueMatrix_T1, 
                                                                    BetaMatrix = BetaMatrix_T1,
                                                                    ProbeOrGene = "Probe", 
                                                                    ContrastColumnName = "CLUSTER", 
                                                                    RGChannelSet = Output_Data_1$RGChannelSet, 
                                                                    SampleSheet = Output_Data_1$SampleSheet, 
                                                                    ProduceImages = "No", PNGorPDF = "png")

head(top_P) # from within DifferentialMethylation function
#               logFC    AveExpr        t      P.Value    adj.P.Val        B
# cg11090352 3.408936 -1.2430269 16.01931 4.382190e-36 2.609874e-30 70.73779
# cg09568691 2.922059 -0.7210139 14.22248 5.797192e-31 1.726300e-25 59.32127
dim(top_P) # 595564      6
matchMethProbesGenes <- match(rownames(top_P), AnnotationFile$V1)
methylationResults <- data.table(top_P, 
                                 Probe = rownames(top_P),
                                 GeneName = sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Name[matchMethProbesGenes]))




# no zero colSums: -which(colSums(RNAseqCountMatrixMatched) == 0)
Match132with170 <- match(colnames(RNAseqCountMatrixMatched), InfiniumClustLabels$ID)
DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = as.matrix(RNAseqCountMatrixMatched), 
                                                                   ContrastColumnName = "CLUSTER", 
                                                                   ClinicalFile = ClinicalFile_T1[Match132with170, ], 
                                                                   ClusterLabels = InfiniumClustLabels$Cluster[Match132with170], 
                                                                   ProduceImages = "Yes", 
                                                                   PNGorPDF = "png") 

head(topTagsPvalue[[1]])
#       logFC      logCPM         F       PValue          FDR      ensemblGeneId
# 1  6.802819  0.60750734 101.43487 4.827170e-18 2.791069e-13  ENSG00000229314.4
# 2 -8.138583  3.75965306  88.40850 1.048977e-16 3.032593e-12 ENSG00000163092.15
dim(topTagsPvalue[[1]]) # 57820     6
# topTagsPvalue <- topTagsPvalue[[1]][topTagsPvalue[[1]]$PValue < 0.01, ]

ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")
RNAseqtoGene <- match(substr(topTagsPvalue[[1]]$ensemblGeneId, 1, 15), ENSMBLid$gene)
# ENSMBLid$name[RNAseqtoGene]
RNAseqresults <- data.table(topTagsPvalue[[1]], GeneName = ENSMBLid$name[RNAseqtoGene])
compareRNAseqMethylationGenes <- match(RNAseqresults$GeneName, methylationResults$GeneName)

matchedMethylationResults <- methylationResults[compareRNAseqMethylationGenes[! is.na(compareRNAseqMethylationGenes)], ]
dim(matchedMethylationResults) #  19851     8
matchedRNAseqresults <- RNAseqresults[match(matchedMethylationResults$GeneName, RNAseqresults$GeneName), ]
dim(matchedRNAseqresults) # 19851     7

# colour vector to determine if both pvalues are significant or not
colVector <- rep("NoSig", nrow(matchedMethylationResults))
for(i in 1:nrow(matchedMethylationResults)) {
  if((matchedMethylationResults$adj.P.Val[i] < 0.05) && (matchedRNAseqresults$PValue[i] < 0.05)) {
    colVector[i] <- "Sig"
  }
}
table(colVector)
# NoSig   Sig 
# 17553  2298


combinedResults <- data.table(matchedMethylationResults, matchedRNAseqresults, colVector)
colnames(combinedResults)[1:8] <- paste0("Methylation", colnames(combinedResults)[1:8])
colnames(combinedResults)[9:15] <- paste0("RNAseq", colnames(combinedResults)[9:15])


p2 <- ggplot(combinedResults, aes(x = MethylationlogFC, y = RNAseqlogFC, color=factor(colVector))) +
  geom_point(size = 2) +
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "Methylation FC")+
  theme_bw() + 
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold") )+
  scale_color_manual(values = c("#bababa", "#b2182b"),
                     name = "Significance") +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))

plot(matchedMethylationResults$logFC, matchedRNAseqresults$logFC, col = colVector+1)



# Determine genes with inverse correlation 
inverseCorVector <- rep("NoSig", nrow(matchedMethylationResults))
for(i in 1:nrow(matchedMethylationResults)) {
  if(colVector[i] == "Sig") {
    if((matchedMethylationResults$logFC[i] > 0) && (matchedRNAseqresults$logFC[i] < 0)) {
      inverseCorVector[i] <- "PosMeth,NegExp"
    }
    
    if((matchedMethylationResults$logFC[i] < 0) && (matchedRNAseqresults$logFC[i] > 0)) {
      inverseCorVector[i] <- "NegMeth,PosExp"
    }
  }
}
table(inverseCorVector)
# NegMeth,PosExp          NoSig        PosMeth,NegExp 
# 749                     18855                   247

inverseCorVectorGenes <- data.table(combinedResults$MethylationGeneName, inverseCorVector)
inverseCorVectorGenes[which(inverseCorVector == "NegMeth,PosExp"), ]
inverseCorVectorGenes[which(inverseCorVector == "PosMeth,NegExp"), ]$V1


# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) 
DiffMethylatedRegionsInfiniumClusters <- DiffMethylatedRegions(Method = "DMRcate", 
                                                               BetaMatrix = BetaMatrix_T1, 
                                                               MvalueMatrix = MvalueMatrix_T1, 
                                                               ContrastColumnName = "CLUSTER", 
                                                               ClinicalFile = ClinicalFile_T1, 
                                                               ClusterLabels = ClusterLabels2to4[[2]][[2]],
                                                               AnnotationFile = AnnotationFile, 
                                                               ProduceImages = "No", 
                                                               DMR = 1, ExpressionFile = NA)

# ------------------------------------------ -
# 11 May 2020
# running Robert' script for generating cell type
require("devtools")
install_github("brentp/celltypes450")
library(celltypes450)
# Get CpGs for cell deconvolution (takes a long time to run)
# data(IDOLOptimizedCpGs)
RGChannelSet <- readRDS("2_RGChannelSet_25June2020.rds")
countsEPIC <- FlowSorted.Blood.EPIC::estimateCellCounts2(
                                    rgSet = RGChannelSet,
                                    compositeCellType = "Blood",
                                    processMethod = "preprocessNoob",
                                    probeSelect = "IDOL",
                                    cellTypes = c("CD8T", "CD4T", "NK", "Bcell","Mono", "Neu"),
                                    referencePlatform = "IlluminaHumanMethylationEPIC",
                                    referenceset = NULL,
                                    IDOLOptimizedCpGs = IDOLOptimizedCpGs,
                                    returnAll = TRUE,
                                    meanPlot = TRUE,
                                    verbose = FALSE)
names(countsEPIC) # "counts"         "compTable"      "normalizedData"
range(countsEPIC$counts) # -1.387779e-17  6.893340e-01
dim(countsEPIC$counts) # 170   6
range(countsEPIC$compTable) # 2.751438e-56 3.071209e+04
dim(countsEPIC$compTable) # 865859     11


  EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020 <- readRDS(file = paste0(MethylationDirPath, "/EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020.rds"))
dim(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020) # 170 6
colnames(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020)
# "CD8T"  "CD4T"  "NK"    "Bcell" "Mono"  "Neu"
match(rownames(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020), InfiniumClustLabels$ID) # match


# Add clusters to the table
Table <- EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020 %>%
  data.frame() %>% mutate(SAMPLE_ID = row.names(.)) %>%
  left_join(ClinicalFile_T1[, c("SAMPLE_ID", "STAGE", "TYPE", "SEX", "SITE_BIOPSY", 
                                "TYPE_BIOPSY", "INSTITUTION", "TRANSLOC_14_18", "EPIC_QC")])
Table <- cbind(Table, factor(InfiniumClustLabels$Cluster))
colnames(Table)[16] <- "CLUSTER"
Table[, which(colnames(Table) == "TRANSLOC_14_18")] <- factor(ClinicalFile_T1$TRANSLOC_14_18)

EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020 <- readRDS(file = paste0(MethylationDirPath, "/EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020.rds"))

CellTypesFLonly <- EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020[which(substr(rownames(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020), 4, 5) == "FL"), ]
estimateCellCountsMethylation(CellTypesFile = CellTypesFLonly,
                              RGChannelSet = NA,
                              ClinicalFile = ClinicalFile_T1[which(ClinicalFile_T1$TYPE == "FL"), ], 
                              ClusterLabels = InfiniumClustLabels$Cluster[which(ClinicalFile_T1$TYPE == "FL")],
                              FigureGenerate = "Yes",
                              PNGorPDF = "png")



# plotting by type
ClusterComparisonOptions <- list(c("ADVANCED", "LIMITED"))
Table %>%
  reshape2::melt(id.vars = c("SAMPLE_ID", "STAGE", "TYPE", "SEX", "SITE_BIOPSY", 
                             "TYPE_BIOPSY", "INSTITUTION", "TRANSLOC_14_18", "EPIC_QC",
                             "CLUSTER")) %>%
  ggpubr::ggboxplot(x = "variable", y = "value", fill = "STAGE",
                    add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                    xlab = "Cell type",
                    font.label = list(size = 20, color = "black"), 
                    palette = c("#d6604d", "#66bd63", "#4575b4")) +  
  ggtitle("Cell type vs. fraction of tumor") +
  stat_compare_means(comparisons = ClusterComparisonOptions) + 
  # Add pairwise comparisons p-value
  stat_compare_means(aes(group = STAGE))

# plotting by cluster
ClusterComparisonOptions <- list(c("1","2"))
Table %>%
  reshape2::melt(id.vars = c("SAMPLE_ID", "STAGE", "TYPE", "SEX", "SITE_BIOPSY", 
                             "TYPE_BIOPSY", "INSTITUTION", "TRANSLOC_14_18", "EPIC_QC",
                             "CLUSTER")) %>%
  ggpubr::ggboxplot(x = "variable", y = "value", fill = "CLUSTER",
                    add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                    xlab = "Cell type",
                    font.label = list(size = 20, color = "black"), 
                    palette = c("#4363d8", "#f58231")) +
  ggtitle("Cell type vs. fraction of tumor - all samples") +
  stat_compare_means(comparisons = ClusterComparisonOptions) + 
  # Add pairwise comparisons p-value
  stat_compare_means(aes(group = CLUSTER))


ClusterComparisonOptions <- list(c("0", "1"))
Table %>%
  reshape2::melt(id.vars = c("SAMPLE_ID", "STAGE", "TYPE", "SEX", "SITE_BIOPSY", 
                             "TYPE_BIOPSY", "INSTITUTION", "TRANSLOC_14_18", "EPIC_QC",
                             "CLUSTER")) %>%
  filter(TYPE == "FL") %>%
  filter(! is.na(TRANSLOC_14_18)) %>%
  ggpubr::ggboxplot(x = "variable", y = "value", fill = factor("TRANSLOC_14_18"),
                    add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                    xlab = "Cell type",
                    font.label = list(size = 20, color = "black"), 
                    palette = c("#e0e0e0","#878787")) +
  ggtitle("Cell type vs. fraction of tumor - FL only samples") +
  stat_compare_means(comparisons = ClusterComparisonOptions) + 
  # Add pairwise comparisons p-value
  stat_compare_means(aes(group = TRANSLOC_14_18))


# plotting by SITE_BIOPSY
ClusterComparisonOptions <- list(c("0", "1"))
Table %>%
  reshape2::melt(id.vars = c("SAMPLE_ID", "STAGE", "TYPE", "SEX", "SITE_BIOPSY", 
                             "TYPE_BIOPSY", "INSTITUTION", "TRANSLOC_14_18", "EPIC_QC",
                             "CLUSTER")) %>%
  filter(TYPE == "FL") %>%
  filter(! is.na(TRANSLOC_14_18)) %>%
  ggpubr::ggboxplot(x = "variable", y = "value", fill = "TRANSLOC_14_18",
                    add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                    xlab = "Cell type",
                    font.label = list(size = 20, color = "black"), 
                    palette = c("#e0e0e0","#878787")) +  
  ggtitle("Cell type vs. fraction of tumor") +
  stat_compare_means(comparisons = ClusterComparisonOptions) + 
  # Add pairwise comparisons p-value
  stat_compare_means(aes(group = TRANSLOC_14_18))



# plotting by SITE_BIOPSY
ClusterComparisonOptions <- list(c("Fair", "Good", "Poor"))
Table %>%
  reshape2::melt(id.vars = c("SAMPLE_ID", "STAGE", "TYPE", "SEX", "SITE_BIOPSY", 
                             "TYPE_BIOPSY", "INSTITUTION", "TRANSLOC_14_18", "EPIC_QC",
                             "CLUSTER")) %>%
  ggpubr::ggboxplot(x = "variable", y = "value", fill = "EPIC_QC",
                    add = "boxplot", ylab = "Fraction of tumor represented by cell type", 
                    xlab = "Cell type",
                    font.label = list(size = 20, color = "black"), 
                    palette = c('#4363d8', '#f58231', '#911eb4')) +  
  ggtitle("Cell type vs. fraction of tumor") +
  stat_compare_means(comparisons = ClusterComparisonOptions) + 
  # Add pairwise comparisons p-value
  stat_compare_means(aes(group = EPIC_QC))


# I was interpreting probably that 0.689 would mean estimated 68.9% of population B-cells, 
# but could be wrong. You could correlate B-cell estimate with purity to explore this. 
TumorPurity <- readRDS(file = paste0("Purity_281probes_10Jan2020.rds"))
identical(rownames(TumorPurity), rownames(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020)) # TRUE
cor(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020[,4], TumorPurity$purity) # 0.5590701
cor.test(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020[,4], TumorPurity$purity)
# 	Pearson's product-moment correlation

#data:  EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020[, 4] and TumorPurity$purity
#t = 8.7398, df = 168, p-value = 2.309e-15
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.4460937 0.6545087
#sample estimates:
#  cor 
#0.5590701 

library(corrplot)
library(RColorBrewer)
M <- cor(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020, TumorPurity) 
corrplot(M, type = "full",
         col = rev(brewer.pal(n = 8, name = "RdYlBu")) )


library("ggpubr")
p1 <- ggscatter(data.frame(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020, TumorPurity), 
                x = "CD8T", y = "purity", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "CD8T normalized counts", ylab = "Purity")
p2 <- ggscatter(data.frame(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020, TumorPurity), 
                x = "CD4T", y = "purity", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "CD4T normalized counts", ylab = "Purity")
p3 <- ggscatter(data.frame(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020, TumorPurity), 
                x = "NK", y = "purity", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "NK normalized counts", ylab = "Purity")
p4 <- ggscatter(data.frame(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020, TumorPurity), 
                x = "Bcell", y = "purity", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Bcell normalized counts", ylab = "Purity")
p5 <- ggscatter(data.frame(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020, TumorPurity), 
                x = "Mono", y = "purity", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Mono normalized counts", ylab = "Purity")
p6 <- ggscatter(data.frame(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020, TumorPurity), 
                x = "Neu", y = "purity", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Neu normalized counts", ylab = "Purity")

gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2)

# 13 May 2020
# Karin sent the new QC file for RNAseq FL_TGL_STAR_logQC_2020-05-13_summary_KI.xlsx
QCKarin13May2020 <- read.csv(file = "/Users/as/Desktop/FLOMICS/RNAseq/ExtendedStudy2020/FL_TGL_STAR_logQC_2020-05-13_summary_KI.csv", header = T)
dim(QCKarin13May2020) # 136  58

# old QCfile
QCcombinedMatched <- readRDS(file = paste0("QCcombined_Integrative.rds"))
dim(QCcombinedMatched) # 132 12

matchQCcombinedQCKarin <- match(QCcombinedMatched$SAMPLE_ID, QCKarin13May2020$SAMPLE_ID)
QCKarin13May2020132 <- QCKarin13May2020[matchQCcombinedQCKarin, ]
dim(QCKarin13May2020132) # 132  58
QCKarin13May2020132Combined <- cbind(QCKarin13May2020132, QCcombinedMatched$rrna_contam, ClusterLabels2to4[[2]][[2]][Match132with170])
colnames(QCKarin13May2020132Combined)[59] <- "rrna_contam"
colnames(QCKarin13May2020132Combined)[60] <- "CLUSTER_InfiniumClust"
write.csv(QCKarin13May2020132Combined, file = "FL_TGL_STAR_logQC_2020-05-13_summary_KI_ClusterContamAdded.csv")


#  14 May 2020
# add clinical file with cluster labels
ClinicalFile_T1CLUSTERs <- cbind(ClinicalFile_T1, ClusterLabels2to4[[2]][[2]])
colnames(ClinicalFile_T1CLUSTERs)[26] <-"CLUSTER"
dim(ClinicalFile_T1CLUSTERs) # 170  26
Output_Visuals_Relation_to_Island_CLUSTER <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", 
                                                                     BetaMatrix = BetaMatrix_T1, 
                                                                     AnnotationFile = AnnotationFile, 
                                                                     ClinicalFile = ClinicalFile_T1CLUSTERs, 
                                                                     ClinicalCategory = "STAGE", 
                                                                     PlotWithinCategories="Yes", 
                                                                     FigureGenerate="Yes", 
                                                                     PNGorPDF="png", 
                                                                     ImageName = "UCSC_RefGene_Group_TYPE_9May2020")

# 18 May 2020
QCMatrix <- read.csv("FL_TGL_STAR_logQC_2020-05-13_summary_KI_ClusterContamAdded.csv")
dim(QCMatrix) # 132  61
QCMatrix$Uniquely.mapped.reads..
QCMatrix$Uniquely.mapped
QCMatrix$X..of.reads.mapped.to.multiple.loci
QCMatrix$rrna_contam

TumorPurity <- readRDS(file = paste0("Purity_281probes_10Jan2020.rds"))

RNAseqCountMatrix <- read.csv("2020-05-19_FL_136_samples_featureCounts_hg19_raw.csv", 
                              row.names = 1)
# fixing missing "_T1" in FL samples 
colnames(RNAseqCountMatrix)[which(substr(colnames(RNAseqCountMatrix), 10, 10) == "")] <-
  paste0(colnames(RNAseqCountMatrix)[which(substr(colnames(RNAseqCountMatrix), 10, 10) == "")], "_T1")

# match RNAseq count matrix to QCMatrix
RNAseqCountMatrixMatched <- RNAseqCountMatrix[, match(QCMatrix$SAMPLE_ID, colnames(RNAseqCountMatrix))]

QCanalysisRNAseq<- QCRNAseq(RNAseqCountMatrix = RNAseqCountMatrixMatched, 
                            QCMatrix = QCMatrix, 
                            BetaMatrix = BetaMatrix_T1,  
                            MvalueMatrix = MvalueMatrix_T1, 
                            TumorPurity = TumorPurity,
                            ClinicalFile = ClinicalFile_T1,
                            SurvivalFile = SurvivalFile,
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
                            ProduceImages = "Yes")



# 21 May 2020 - using new matrix
# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) 
# add clinical file with cluster labels
ClinicalFile_T1CLUSTERs <- cbind(ClinicalFile_T1, ClusterLabels2to4[[2]][[2]])
colnames(ClinicalFile_T1CLUSTERs)[26] <-"CLUSTER"
dim(ClinicalFile_T1CLUSTERs) # 170  26
Output_Differential_probes_InfiniumClust <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1CLUSTERs, 
                                                                    MvalueMatrix = MvalueMatrix_T1, 
                                                                    BetaMatrix = BetaMatrix_T1,
                                                                    ProbeOrGene = "Probe", 
                                                                    ContrastColumnName = "CLUSTER", 
                                                                    RGChannelSet = Output_Data_1$RGChannelSet, 
                                                                    SampleSheet = Output_Data_1$SampleSheet, 
                                                                    ProduceImages = "No", PNGorPDF = "png")

head(top_P) # from within DifferentialMethylation function
#               logFC    AveExpr        t      P.Value    adj.P.Val        B
# cg11090352 3.408936 -1.2430269 16.01931 4.382190e-36 2.609874e-30 70.73779
# cg09568691 2.922059 -0.7210139 14.22248 5.797192e-31 1.726300e-25 59.32127
dim(top_P) # 595564      6
matchMethProbesGenes <- match(rownames(top_P), AnnotationFile$V1)
methylationResults <- data.table(top_P, 
                                 Probe = rownames(top_P),
                                 GeneName = sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Name[matchMethProbesGenes]))




# no zero colSums: which(colSums(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered) == 0)
MatchRNAseqwithBeta <- match(colnames(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), colnames(BetaMatrix_T1))
DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = as.matrix(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), 
                                                                   ContrastColumnName = "CLUSTER", 
                                                                   ClinicalFile = QCanalysisRNAseq$ClinicalFileMatchedSampleFiltered, 
                                                                   ClusterLabels = ClusterLabels2to4[[2]][[2]][MatchRNAseqwithBeta], 
                                                                   ProduceImages = "Yes", 
                                                                   PNGorPDF = "png") 

# Contrast: one-two has 147 ensembl IDs (+ve)logFC and significant.
# Contrast: one-two has 55 ensembl IDs (-ve)logFC and significant.
head(topTagsPvalue[[1]])
#       logFC      logCPM         F       PValue          FDR      ensemblGeneId
# 1  6.802819  0.60750734 101.43487 4.827170e-18 2.791069e-13  ENSG00000229314.4
# 2 -8.138583  3.75965306  88.40850 1.048977e-16 3.032593e-12 ENSG00000163092.15
dim(topTagsPvalue[[1]]) # 34749     6
# topTagsPvalue <- topTagsPvalue[[1]][topTagsPvalue[[1]]$PValue < 0.01, ]

ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")
RNAseqtoGene <- match(substr(topTagsPvalue[[1]]$ensemblGeneId, 1, 15), ENSMBLid$gene)
# ENSMBLid$name[RNAseqtoGene]
RNAseqresults <- data.table(topTagsPvalue[[1]], GeneName = ENSMBLid$name[RNAseqtoGene])
compareRNAseqMethylationGenes <- match(RNAseqresults$GeneName, methylationResults$GeneName)

matchedMethylationResults <- methylationResults[compareRNAseqMethylationGenes[! is.na(compareRNAseqMethylationGenes)], ]
dim(matchedMethylationResults) #   16941     8
matchedRNAseqresults <- RNAseqresults[match(matchedMethylationResults$GeneName, RNAseqresults$GeneName), ]
dim(matchedRNAseqresults) # 16941     7

# colour vector to determine if both pvalues are significant or not
colVector <- rep("NoSig", nrow(matchedMethylationResults))
for(i in 1:nrow(matchedMethylationResults)) {
  if((matchedMethylationResults$adj.P.Val[i] < 0.01) && (matchedRNAseqresults$PValue[i] < 0.01)) {
    colVector[i] <- "Sig"
  }
}
colVector <- relevel(factor(colVector), ref= "Sig")
table(colVector)
# Sig NoSig 
# 1691 15250 


combinedResults <- data.table(matchedMethylationResults, matchedRNAseqresults, colVector)
colnames(combinedResults)[1:8] <- paste0("Methylation", colnames(combinedResults)[1:8])
colnames(combinedResults)[9:15] <- paste0("RNAseq", colnames(combinedResults)[9:15])


p2 <- ggplot(combinedResults, aes(x = MethylationlogFC, y = RNAseqlogFC, color=factor(colVector))) +
  geom_point(size = 2) +
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "Methylation FC")+
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

plot(matchedMethylationResults$logFC, matchedRNAseqresults$logFC, col = c(colVector)+1)



# Determine genes with inverse correlation 
inverseCorVector <- as.vector(colVector)
for(i in 1:nrow(matchedMethylationResults)) {
  if(colVector[i] == "Sig") {
    if((matchedMethylationResults$logFC[i] > 0) && (matchedRNAseqresults$logFC[i] < 0)) {
      inverseCorVector[i] <- "PosMeth,NegExp"
    }
    
    if((matchedMethylationResults$logFC[i] < 0) && (matchedRNAseqresults$logFC[i] > 0)) {
      inverseCorVector[i] <- "NegMeth,PosExp"
    }
  }
}
table(inverseCorVector)
# NegMeth,PosExp          NoSig     PosMeth,NegExp            Sig 
# 564                     15250                 81           1046

combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName
combinedResultsNegMethPosExp <- data.table(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ])

inverseCorVectorGenes[which(inverseCorVector == "PosMeth,NegExp"), ]
combinedResultsPosMethNegExp <- data.table(combinedResults[which(inverseCorVector == "PosMeth,NegExp"), ])

combinedResultsCorrelation <-data.table(combinedResults, "correlation" = inverseCorVector)


p3 <- ggplot(combinedResultsCorrelation[- which(inverseCorVector == "NoSig"), ], aes(x = MethylationlogFC, y = RNAseqlogFC, color=factor(correlation))) +
  geom_point(size = 2) +
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "Methylation FC")+
  theme_bw() + 
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold") )+
  scale_color_manual(values = c("#de2d26", "#41ab5d", "#4292c6"),
                     name = "Significance") +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))



# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) 
DiffMethylatedRegionsInfiniumClusters <- DiffMethylatedRegions(Method = "DMRcate", 
                                                               BetaMatrix = BetaMatrix_T1, 
                                                               MvalueMatrix = MvalueMatrix_T1, 
                                                               ContrastColumnName = "CLUSTER", 
                                                               ClinicalFile = ClinicalFile_T1, 
                                                               ClusterLabels = ClusterLabels2to4[[2]][[2]],
                                                               AnnotationFile = AnnotationFile, 
                                                               ProduceImages = "No", 
                                                               DMR = 1, ExpressionFile = NA)

# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) 


# add clinical file with cluster labels
ClinicalFile_T1CLUSTERs <- cbind(ClinicalFile_T1, ClusterLabels2to4[[2]][[2]])
colnames(ClinicalFile_T1CLUSTERs)[26] <-"CLUSTER"
dim(ClinicalFile_T1CLUSTERs) # 170  26
Output_Differential_probes_InfiniumClust <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1CLUSTERs, 
                                                                    MvalueMatrix = MvalueMatrix_T1, 
                                                                    BetaMatrix = BetaMatrix_T1,
                                                                    ProbeOrGene = "Probe", 
                                                                    ContrastColumnName = "CLUSTER", 
                                                                    RGChannelSet = Output_Data_1$RGChannelSet, 
                                                                    SampleSheet = Output_Data_1$SampleSheet, 
                                                                    ProduceImages = "No", PNGorPDF = "png")

head(top_P) # from within DifferentialMethylation function
#               logFC    AveExpr        t      P.Value    adj.P.Val        B
# cg11090352 3.408936 -1.2430269 16.01931 4.382190e-36 2.609874e-30 70.73779
# cg09568691 2.922059 -0.7210139 14.22248 5.797192e-31 1.726300e-25 59.32127
dim(top_P) # 595564      6
matchMethProbesGenes <- match(rownames(top_P), AnnotationFile$V1)
methylationResults <- data.table(top_P, 
                                 Probe = rownames(top_P),
                                 GeneName = sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Name[matchMethProbesGenes]))




# no zero colSums: -which(colSums(RNASeqCountMatrixMatched) == 0)
DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = as.matrix(RNASeqCountMatrixMatched), 
                                                                   ContrastColumnName = "CLUSTER", 
                                                                   ClinicalFile = ClinicalFileMatched, 
                                                                   ClusterLabels = ClusterLabels2to4[[2]][[2]][Match132with170], 
                                                                   ProduceImages = "Yes", 
                                                                   PNGorPDF = "png") 

head(topTagsPvalue[[1]])
#       logFC      logCPM         F       PValue          FDR      ensemblGeneId
# 1  6.802819  0.60750734 101.43487 4.827170e-18 2.791069e-13  ENSG00000229314.4
# 2 -8.138583  3.75965306  88.40850 1.048977e-16 3.032593e-12 ENSG00000163092.15
dim(topTagsPvalue[[1]]) # 57820     6
# topTagsPvalue <- topTagsPvalue[[1]][topTagsPvalue[[1]]$PValue < 0.01, ]

ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")
RNAseqtoGene <- match(substr(topTagsPvalue[[1]]$ensemblGeneId, 1, 15), ENSMBLid$gene)
# ENSMBLid$name[RNAseqtoGene]
RNAseqresults <- data.table(topTagsPvalue[[1]], GeneName = ENSMBLid$name[RNAseqtoGene])
compareRNAseqMethylationGenes <- match(RNAseqresults$GeneName, methylationResults$GeneName)

matchedMethylationResults <- methylationResults[compareRNAseqMethylationGenes[! is.na(compareRNAseqMethylationGenes)], ]
dim(matchedMethylationResults) #  19851     8
matchedRNAseqresults <- RNAseqresults[match(matchedMethylationResults$GeneName, RNAseqresults$GeneName), ]
dim(matchedRNAseqresults) # 19851     7

# colour vector to determine if both pvalues are significant or not
colVector <- rep("NoSig", nrow(matchedMethylationResults))
for(i in 1:nrow(matchedMethylationResults)) {
  if((matchedMethylationResults$adj.P.Val[i] < 0.01) && (matchedRNAseqresults$PValue[i] < 0.01)) {
    colVector[i] <- "Sig"
  }
}
table(colVector)
# NoSig   Sig 
# 17553  2298


combinedResults <- data.table(matchedMethylationResults, matchedRNAseqresults, colVector)
colnames(combinedResults)[1:8] <- paste0("Methylation", colnames(combinedResults)[1:8])
colnames(combinedResults)[9:15] <- paste0("RNAseq", colnames(combinedResults)[9:15])


p2 <- ggplot(combinedResults, aes(x = MethylationlogFC, y = RNAseqlogFC, color=factor(colVector))) +
  geom_point(size = 2) +
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "Methylation FC")+
  theme_bw() + 
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold") )+
  scale_color_manual(values = c("#bababa", "#b2182b"),
                     name = "Significance") +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))

plot(matchedMethylationResults$logFC, matchedRNAseqresults$logFC, col = colVector+1)



# Determine genes with inverse correlation 
inverseCorVector <- rep("NoSig", nrow(matchedMethylationResults))
for(i in 1:nrow(matchedMethylationResults)) {
  if(colVector[i] == "Sig") {
    if((matchedMethylationResults$logFC[i] > 0) && (matchedRNAseqresults$logFC[i] < 0)) {
      inverseCorVector[i] <- "PosMeth,NegExp"
    }
    
    if((matchedMethylationResults$logFC[i] < 0) && (matchedRNAseqresults$logFC[i] > 0)) {
      inverseCorVector[i] <- "NegMeth,PosExp"
    }
  }
}
table(inverseCorVector)
# NegMeth,PosExp          NoSig        PosMeth,NegExp 
# 749                     18855                   247

inverseCorVectorGenes <- data.table(combinedResults$MethylationGeneName, inverseCorVector)
inverseCorVectorGenes[which(inverseCorVector == "NegMeth,PosExp"), ]
inverseCorVectorGenes[which(inverseCorVector == "PosMeth,NegExp"), ]$V1


# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) 
DiffMethylatedRegionsInfiniumClusters <- DiffMethylatedRegions(Method = "DMRcate", 
                                                               BetaMatrix = BetaMatrix_T1, 
                                                               MvalueMatrix = MvalueMatrix_T1, 
                                                               ContrastColumnName = "CLUSTER", 
                                                               ClinicalFile = ClinicalFile_T1, 
                                                               ClusterLabels = ClusterLabels2to4[[2]][[2]],
                                                               AnnotationFile = AnnotationFile, 
                                                               ProduceImages = "No", 
                                                               DMR = 1, ExpressionFile = NA)


# 18 May 2020
QCMatrix <- read.csv("FL_TGL_STAR_logQC_2020-05-13_summary_KI_ClusterContamAdded.csv")
dim(QCMatrix) # 132  61
QCMatrix$Uniquely.mapped.reads..
QCMatrix$Uniquely.mapped
QCMatrix$X..of.reads.mapped.to.multiple.loci
QCMatrix$rrna_contam

TumorPurity <- readRDS(file = paste0("Purity_281probes_10Jan2020.rds"))

RNAseqCountMatrix <- read.csv("2020-05-19_FL_136_samples_featureCounts_hg19_raw.csv", 
                              row.names = 1)
# fixing missing "_T1" in FL samples 
colnames(RNAseqCountMatrix)[which(substr(colnames(RNAseqCountMatrix), 10, 10) == "")] <-
  paste0(colnames(RNAseqCountMatrix)[which(substr(colnames(RNAseqCountMatrix), 10, 10) == "")], "_T1")

# match RNAseq count matrix to QCMatrix
RNAseqCountMatrixMatched <- RNAseqCountMatrix[, match(QCMatrix$SAMPLE_ID, colnames(RNAseqCountMatrix))]

QCanalysisRNAseq<- QCRNAseq(RNAseqCountMatrix = RNAseqCountMatrixMatched, 
                            QCMatrix = QCMatrix, 
                            BetaMatrix = BetaMatrix_T1,  
                            MvalueMatrix = MvalueMatrix_T1, 
                            TumorPurity = TumorPurity,
                            ClinicalFile = ClinicalFile_T1,
                            SurvivalFile = SurvivalFile,
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
                            ProduceImages = "Yes")



# 21 May 2020 - using new matrix
# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change)  probe level
# add clinical file with cluster labels
ClinicalFile_T1CLUSTERs <- cbind(ClinicalFile_T1, ClusterLabels2to4[[2]][[2]])
colnames(ClinicalFile_T1CLUSTERs)[26] <-"CLUSTER"
dim(ClinicalFile_T1CLUSTERs) # 170  26
Output_Differential_probes_InfiniumClust <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1CLUSTERs, 
                                                                    MvalueMatrix = MvalueMatrix_T1, 
                                                                    BetaMatrix = BetaMatrix_T1,
                                                                    ProbeOrGene = "Probe", 
                                                                    ContrastColumnName = "CLUSTER", 
                                                                    RGChannelSet = Output_Data_1$RGChannelSet, 
                                                                    SampleSheet = Output_Data_1$SampleSheet, 
                                                                    ProduceImages = "No", PNGorPDF = "png")

head(top_P) # from within DifferentialMethylation function
#               logFC    AveExpr        t      P.Value    adj.P.Val        B
# cg11090352 3.408936 -1.2430269 16.01931 4.382190e-36 2.609874e-30 70.73779
# cg09568691 2.922059 -0.7210139 14.22248 5.797192e-31 1.726300e-25 59.32127
dim(top_P) # 595564      6
matchMethProbesGenes <- match(rownames(top_P), AnnotationFile$V1)
methylationResults <- data.table(top_P, 
                                 Probe = rownames(top_P),
                                 GeneName = sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Name[matchMethProbesGenes]))




# no zero colSums: which(colSums(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered) == 0)
MatchRNAseqwithBeta <- match(colnames(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), colnames(BetaMatrix_T1))
DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = as.matrix(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), 
                                                                   ContrastColumnName = "CLUSTER", 
                                                                   ClinicalFile = QCanalysisRNAseq$ClinicalFileMatchedSampleFiltered, 
                                                                   ClusterLabels = ClusterLabels2to4[[2]][[2]][MatchRNAseqwithBeta], 
                                                                   ProduceImages = "Yes", 
                                                                   PNGorPDF = "png") 

# Contrast: one-two has 147 ensembl IDs (+ve)logFC and significant.
# Contrast: one-two has 55 ensembl IDs (-ve)logFC and significant.
head(topTagsPvalue[[1]])
#       logFC      logCPM         F       PValue          FDR      ensemblGeneId
# 1  6.802819  0.60750734 101.43487 4.827170e-18 2.791069e-13  ENSG00000229314.4
# 2 -8.138583  3.75965306  88.40850 1.048977e-16 3.032593e-12 ENSG00000163092.15
dim(topTagsPvalue[[1]]) # 34749     6
# topTagsPvalue <- topTagsPvalue[[1]][topTagsPvalue[[1]]$PValue < 0.01, ]

ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")
RNAseqtoGene <- match(substr(topTagsPvalue[[1]]$ensemblGeneId, 1, 15), ENSMBLid$gene)
# ENSMBLid$name[RNAseqtoGene]
RNAseqresults <- data.table(topTagsPvalue[[1]], GeneName = ENSMBLid$name[RNAseqtoGene])
compareRNAseqMethylationGenes <- match(RNAseqresults$GeneName, methylationResults$GeneName)

matchedMethylationResults <- methylationResults[compareRNAseqMethylationGenes[! is.na(compareRNAseqMethylationGenes)], ]
dim(matchedMethylationResults) #   16941     8
matchedRNAseqresults <- RNAseqresults[match(matchedMethylationResults$GeneName, RNAseqresults$GeneName), ]
dim(matchedRNAseqresults) # 16941     7

# colour vector to determine if both pvalues are significant or not
colVector <- rep("NoSig", nrow(matchedMethylationResults))
for(i in 1:nrow(matchedMethylationResults)) {
  if((matchedMethylationResults$adj.P.Val[i] < 0.01) && (matchedRNAseqresults$PValue[i] < 0.01)) {
    colVector[i] <- "Sig"
  }
}
colVector <- relevel(factor(colVector), ref= "Sig")
table(colVector)
# Sig NoSig 
# 1691 15250 


combinedResults <- data.table(matchedMethylationResults, matchedRNAseqresults, colVector)
colnames(combinedResults)[1:8] <- paste0("Methylation", colnames(combinedResults)[1:8])
colnames(combinedResults)[9:15] <- paste0("RNAseq", colnames(combinedResults)[9:15])


p2 <- ggplot(combinedResults, aes(x = MethylationlogFC, y = RNAseqlogFC, color=factor(colVector))) +
  geom_point(size = 2) +
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "Methylation FC")+
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

plot(matchedMethylationResults$logFC, matchedRNAseqresults$logFC, col = c(colVector)+1)



# Determine genes with inverse correlation 
inverseCorVector <- as.vector(colVector)
for(i in 1:nrow(matchedMethylationResults)) {
  if(colVector[i] == "Sig") {
    if((matchedMethylationResults$logFC[i] > 0) && (matchedRNAseqresults$logFC[i] < 0)) {
      inverseCorVector[i] <- "PosMeth,NegExp"
    }
    
    if((matchedMethylationResults$logFC[i] < 0) && (matchedRNAseqresults$logFC[i] > 0)) {
      inverseCorVector[i] <- "NegMeth,PosExp"
    }
  }
}
table(inverseCorVector)
# NegMeth,PosExp          NoSig     PosMeth,NegExp            Sig 
# 564                     15250                 81           1046

combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName
combinedResultsNegMethPosExp <- data.table(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ])

inverseCorVectorGenes[which(inverseCorVector == "PosMeth,NegExp"), ]
combinedResultsPosMethNegExp <- data.table(combinedResults[which(inverseCorVector == "PosMeth,NegExp"), ])

combinedResultsCorrelation <-data.table(combinedResults, "correlation" = inverseCorVector)


p3 <- ggplot(combinedResultsCorrelation[- which(inverseCorVector == "NoSig"), ], aes(x = MethylationlogFC, y = RNAseqlogFC, color=factor(correlation))) +
  geom_point(size = 2) +
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "Methylation FC")+
  theme_bw() + 
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold") )+
  scale_color_manual(values = c("#de2d26", "#41ab5d", "#4292c6"),
                     name = "Significance") +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))



# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) 
DiffMethylatedRegionsInfiniumClusters <- DiffMethylatedRegions(Method = "DMRcate", 
                                                               BetaMatrix = BetaMatrix_T1, 
                                                               MvalueMatrix = MvalueMatrix_T1, 
                                                               ContrastColumnName = "CLUSTER", 
                                                               ClinicalFile = ClinicalFile_T1, 
                                                               ClusterLabels = ClusterLabels2to4[[2]][[2]],
                                                               AnnotationFile = AnnotationFile, 
                                                               ProduceImages = "No", 
                                                               DMR = 1, ExpressionFile = NA)

# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) 


# add clinical file with cluster labels
ClinicalFile_T1CLUSTERs <- cbind(ClinicalFile_T1, ClusterLabels2to4[[2]][[2]])
colnames(ClinicalFile_T1CLUSTERs)[26] <-"CLUSTER"
dim(ClinicalFile_T1CLUSTERs) # 170  26
Output_Differential_probes_InfiniumClust <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1CLUSTERs, 
                                                                    MvalueMatrix = MvalueMatrix_T1, 
                                                                    BetaMatrix = BetaMatrix_T1,
                                                                    ProbeOrGene = "Probe", 
                                                                    ContrastColumnName = "CLUSTER", 
                                                                    RGChannelSet = Output_Data_1$RGChannelSet, 
                                                                    SampleSheet = Output_Data_1$SampleSheet, 
                                                                    ProduceImages = "No", PNGorPDF = "png")

head(top_P) # from within DifferentialMethylation function
#               logFC    AveExpr        t      P.Value    adj.P.Val        B
# cg11090352 3.408936 -1.2430269 16.01931 4.382190e-36 2.609874e-30 70.73779
# cg09568691 2.922059 -0.7210139 14.22248 5.797192e-31 1.726300e-25 59.32127
dim(top_P) # 595564      6
matchMethProbesGenes <- match(rownames(top_P), AnnotationFile$V1)
methylationResults <- data.table(top_P, 
                                 Probe = rownames(top_P),
                                 GeneName = sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Name[matchMethProbesGenes]))




# no zero colSums: -which(colSums(RNASeqCountMatrixMatched) == 0)
DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = as.matrix(RNASeqCountMatrixMatched), 
                                                                   ContrastColumnName = "CLUSTER", 
                                                                   ClinicalFile = ClinicalFileMatched, 
                                                                   ClusterLabels = ClusterLabels2to4[[2]][[2]][Match132with170], 
                                                                   ProduceImages = "Yes", 
                                                                   PNGorPDF = "png") 

head(topTagsPvalue[[1]])
#       logFC      logCPM         F       PValue          FDR      ensemblGeneId
# 1  6.802819  0.60750734 101.43487 4.827170e-18 2.791069e-13  ENSG00000229314.4
# 2 -8.138583  3.75965306  88.40850 1.048977e-16 3.032593e-12 ENSG00000163092.15
dim(topTagsPvalue[[1]]) # 57820     6
# topTagsPvalue <- topTagsPvalue[[1]][topTagsPvalue[[1]]$PValue < 0.01, ]

ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")
RNAseqtoGene <- match(substr(topTagsPvalue[[1]]$ensemblGeneId, 1, 15), ENSMBLid$gene)
# ENSMBLid$name[RNAseqtoGene]
RNAseqresults <- data.table(topTagsPvalue[[1]], GeneName = ENSMBLid$name[RNAseqtoGene])
compareRNAseqMethylationGenes <- match(RNAseqresults$GeneName, methylationResults$GeneName)

matchedMethylationResults <- methylationResults[compareRNAseqMethylationGenes[! is.na(compareRNAseqMethylationGenes)], ]
dim(matchedMethylationResults) #  19851     8
matchedRNAseqresults <- RNAseqresults[match(matchedMethylationResults$GeneName, RNAseqresults$GeneName), ]
dim(matchedRNAseqresults) # 19851     7

# colour vector to determine if both pvalues are significant or not
colVector <- rep("NoSig", nrow(matchedMethylationResults))
for(i in 1:nrow(matchedMethylationResults)) {
  if((matchedMethylationResults$adj.P.Val[i] < 0.01) && (matchedRNAseqresults$PValue[i] < 0.01)) {
    colVector[i] <- "Sig"
  }
}
table(colVector)
# NoSig   Sig 
# 17553  2298


combinedResults <- data.table(matchedMethylationResults, matchedRNAseqresults, colVector)
colnames(combinedResults)[1:8] <- paste0("Methylation", colnames(combinedResults)[1:8])
colnames(combinedResults)[9:15] <- paste0("RNAseq", colnames(combinedResults)[9:15])


p2 <- ggplot(combinedResults, aes(x = MethylationlogFC, y = RNAseqlogFC, color=factor(colVector))) +
  geom_point(size = 2) +
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "Methylation FC")+
  theme_bw() + 
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold") )+
  scale_color_manual(values = c("#bababa", "#b2182b"),
                     name = "Significance") +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))

plot(matchedMethylationResults$logFC, matchedRNAseqresults$logFC, col = colVector+1)



# Determine genes with inverse correlation 
inverseCorVector <- rep("NoSig", nrow(matchedMethylationResults))
for(i in 1:nrow(matchedMethylationResults)) {
  if(colVector[i] == "Sig") {
    if((matchedMethylationResults$logFC[i] > 0) && (matchedRNAseqresults$logFC[i] < 0)) {
      inverseCorVector[i] <- "PosMeth,NegExp"
    }
    
    if((matchedMethylationResults$logFC[i] < 0) && (matchedRNAseqresults$logFC[i] > 0)) {
      inverseCorVector[i] <- "NegMeth,PosExp"
    }
  }
}
table(inverseCorVector)
# NegMeth,PosExp          NoSig        PosMeth,NegExp 
# 749                     18855                   247

inverseCorVectorGenes <- data.table(combinedResults$MethylationGeneName, inverseCorVector)
inverseCorVectorGenes[which(inverseCorVector == "NegMeth,PosExp"), ]
inverseCorVectorGenes[which(inverseCorVector == "PosMeth,NegExp"), ]$V1









# 25 May 2020
# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change)  probe level
# Use new count matrix STAR
# QCMatrix <- read.csv("FL_TGL_STAR_logQC_2020-05-13_summary_KI_ClusterContamAdded.csv")
# dim(QCMatrix) # 132  61
# TumorPurity <- readRDS(file = paste0("Purity_281probes_10Jan2020.rds"))
#  RNAseqCountMatrix <- read.delim("STAR_quantmode_counts_matrix_FL_136_patients_rcd25May2020.txt", row.names = 1)
# dim(RNAseqCountMatrix) # 57820   136
# fixing missing "_T1" in FL samples 
# colnames(RNAseqCountMatrix)[which(substr(colnames(RNAseqCountMatrix), 10, 10) == "")] <-
#   paste0(colnames(RNAseqCountMatrix)[which(substr(colnames(RNAseqCountMatrix), 10, 10) == "")], "_T1")

# match RNAseq count matrix to QCMatrix
# RNAseqCountMatrixMatched <- RNAseqCountMatrix[, match(QCMatrix$SAMPLE_ID, colnames(RNAseqCountMatrix))]
# dim(RNAseqCountMatrixMatched) # 57820   132
# write.csv(RNAseqCountMatrixMatched, file = "STAR_quantmode_counts_matrix_MATCHED_FL_132_patients_rcd25May2020.csv")
# RNAseqCountMatrixMatched <- read.csv("STAR_quantmode_counts_matrix_MATCHED_FL_132_patients_rcd25May2020.csv", 
#                                      row.names = 1)
dim(RNAseqCountMatrixMatched) # 57820   132

QCanalysisRNAseq<- QCRNAseq(RNAseqCountMatrix = RNAseqCountMatrixMatched, 
                            QCMatrix = RNAseqQCFile, 
                            BetaMatrix = BetaMatrix_T1,  
                            MvalueMatrix = MvalueMatrix_T1, 
                            TumorPurity = TumorPurity,
                            ClinicalFile = ClinicalFile_T1,
                            SurvivalFile = SurvivalFile,
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
                            ProduceImages = "No")



# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change)  probe level
# add clinical file with cluster labels
ClinicalFile_T1CLUSTERs <- cbind(ClinicalFile_T1, InfiniumClustLabels$Cluster)
colnames(ClinicalFile_T1CLUSTERs)[26] <-"CLUSTER"
dim(ClinicalFile_T1CLUSTERs) # 170  26
Output_Differential_probes_InfiniumClust <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1CLUSTERs, 
                                                                    MvalueMatrix = MvalueMatrix_T1, 
                                                                    BetaMatrix = BetaMatrix_T1,
                                                                    ProbeOrGene = "Probe", 
                                                                    ContrastColumnName = "CLUSTER", 
                                                                    RGChannelSet = Output_Data_1$RGChannelSet, 
                                                                    SampleSheet = Output_Data_1$SampleSheet, 
                                                                    ProduceImages = "No", PNGorPDF = "png")

head(top_P) # from within DifferentialMethylation function
#               logFC    AveExpr        t      P.Value    adj.P.Val        B
# cg11090352 3.408936 -1.2430269 16.01931 4.382190e-36 2.609874e-30 70.73779
# cg09568691 2.922059 -0.7210139 14.22248 5.797192e-31 1.726300e-25 59.32127
dim(top_P) # 595564      6
matchMethProbesGenes <- match(rownames(top_P), AnnotationFile$V1)
methylationResults <- data.table(top_P, 
                                 Probe = rownames(top_P),
                                 GeneName = sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Name[matchMethProbesGenes]))




# no zero colSums: which(colSums(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered) == 0)
MatchRNAseqwithBeta <- match(colnames(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), colnames(BetaMatrix_T1))
DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = as.matrix(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), 
                                                                   ContrastColumnName = "CLUSTER", 
                                                                   ClinicalFile = QCanalysisRNAseq$ClinicalFileMatchedSampleFiltered, 
                                                                   ClusterLabels = InfiniumClustLabels$Cluster[MatchRNAseqwithBeta], 
                                                                   ProduceImages = "Yes", 
                                                                   PNGorPDF = "png") 

# Contrast: one-two has 116 ensembl IDs (+ve)logFC and significant. 
# Contrast: one-two has 50 ensembl IDs (-ve)logFC and significant. 
head(topTagsPvalue[[1]])
#       logFC      logCPM         F       PValue          FDR      ensemblGeneId
# 1  6.111903  0.8395439 46.25519 4.308221e-09 0.0001160263 ENSG00000125144
# 2  5.600054  4.2377716 44.45528 7.360204e-09 0.0001160263 ENSG00000196611
dim(topTagsPvalue[[1]]) # 31528     6
# topTagsPvalue <- topTagsPvalue[[1]][topTagsPvalue[[1]]$PValue < 0.01, ]

ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")
RNAseqtoGene <- match(substr(topTagsPvalue[[1]]$ensemblGeneId, 1, 15), ENSMBLid$gene)
# ENSMBLid$name[RNAseqtoGene]
RNAseqresults <- data.table(topTagsPvalue[[1]], GeneName = ENSMBLid$name[RNAseqtoGene])
compareRNAseqMethylationGenes <- match(RNAseqresults$GeneName, methylationResults$GeneName)

matchedMethylationResults <- methylationResults[compareRNAseqMethylationGenes[! is.na(compareRNAseqMethylationGenes)], ]
dim(matchedMethylationResults) #   16292     8
matchedRNAseqresults <- RNAseqresults[match(matchedMethylationResults$GeneName, RNAseqresults$GeneName), ]
dim(matchedRNAseqresults) # 16292     7

# colour vector to determine if both pvalues are significant or not
colVector <- rep("NoSig", nrow(matchedMethylationResults))
for(i in 1:nrow(matchedMethylationResults)) {
  if((matchedMethylationResults$adj.P.Val[i] < 0.05) && (matchedRNAseqresults$PValue[i] < 0.05)) {
    colVector[i] <- "Sig"
  }
}
colVector <- relevel(factor(colVector), ref= "Sig")
table(colVector)
# Sig NoSig 
# 3702 12590 


combinedResults <- data.table(matchedMethylationResults, matchedRNAseqresults, colVector)
colnames(combinedResults)[1:8] <- paste0("Methylation", colnames(combinedResults)[1:8])
colnames(combinedResults)[9:15] <- paste0("RNAseq", colnames(combinedResults)[9:15])


p2 <- ggplot(combinedResults, aes(x = MethylationlogFC, y = RNAseqlogFC, color=factor(colVector))) +
  geom_point(size = 2) +
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "Methylation FC")+
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

plot(matchedMethylationResults$logFC, matchedRNAseqresults$logFC, col = c(colVector)+1)



# Determine genes with inverse correlation 
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
table(inverseCorVector)
# NegMeth,NegExp  NegMeth,PosExp          NoSig     PosMeth,NegExp   PosMeth,PosExp 
#            192           1328          12590                300           1882 


A <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,PosExp"), ]$RNAseqGeneName)
B <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,NegExp"), ]$RNAseqGeneName)
C <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,NegExp"), ]$RNAseqGeneName)
D <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName)
E <-  unfactor(combinedResults[which(inverseCorVector == "NoSig"), ]$RNAseqGeneName)
length(A) <- length(B) <- length(C) <- length(D) <- length(E)

geneList <- cbind(A,B,C,D,E)
colnames(geneList) <- c("PosMethPosExp", "PosMethNegExp", "NegMethNegExp", "NegMethPosExp", "NoSig")
dim(geneList) # 12590     5
write.csv(geneList, file = "FCgeneList_allMethProbes_allGenes_1June2020.csv")

# combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ] %>% arrange(desc(abs(MethylationlogFC)))
# combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ] %>% arrange(desc(abs(RNAseqlogFC)))
# write.csv(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName,
#           file = "NegMethPosExp_25May2020.csv")
# combinedResultsNegMethPosExp <- data.table(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ])
A <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,PosExp"), ]$RNAseqGeneName)
B <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,NegExp"), ]$RNAseqGeneName)
C <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,NegExp"), ]$RNAseqGeneName)
D <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName)
E <-  unfactor(combinedResults[which(inverseCorVector == "NoSig"), ]$RNAseqGeneName)

gprofilerPosMethPosExp <- gprofiler2::gost(query = list(A), user_threshold = 0.01, 
                                           as_short_link = TRUE)
gprofilerPosMethNegExp <- gprofiler2::gost(query = list(B), user_threshold = 0.01,
                                           as_short_link = TRUE)
gprofilerNegMethNegExp <- gprofiler2::gost(query = list(C), user_threshold = 0.01,
                                           as_short_link = TRUE)
gprofilerNegMethPosExp <- gprofiler2::gost(query = list(D), user_threshold = 0.01,
                                           as_short_link = TRUE)
gprofilerNoSig <- gprofiler2::gost(query = list(E), user_threshold = 0.01,
                                   as_short_link = TRUE)


combinedResultsCorrelation <- data.table(combinedResults, "correlation" = inverseCorVector)


p3 <- ggplot2::ggplot(combinedResultsCorrelation[- which(inverseCorVector == "NoSig"), ], 
                      aes(x = MethylationlogFC, y = RNAseqlogFC, color=factor(correlation))) +
  geom_point(size = 2) +
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "Methylation FC")+
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
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))

# order by methylationFC
combinedResults %>% arrange(desc(abs(MethylationlogFC)))  %>% filter(colVector == "Sig")
# order by RNAseqFC
combinedResults %>% arrange(desc(abs(RNAseqlogFC)))  %>% filter(colVector == "Sig")





# 27 May 2020
# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change)  probe level
# Smae analysis as above, but only using "island" probes
# QCMatrix <- read.csv("FL_TGL_STAR_logQC_2020-05-13_summary_KI_ClusterContamAdded.csv")
# dim(QCMatrix) # 132  61
# TumorPurity <- readRDS(file = paste0("Purity_281probes_10Jan2020.rds"))
#  RNAseqCountMatrix <- read.delim("STAR_quantmode_counts_matrix_FL_136_patients_rcd25May2020.txt", row.names = 1)
# dim(RNAseqCountMatrix) # 57820   136
# fixing missing "_T1" in FL samples 
# colnames(RNAseqCountMatrix)[which(substr(colnames(RNAseqCountMatrix), 10, 10) == "")] <-
#   paste0(colnames(RNAseqCountMatrix)[which(substr(colnames(RNAseqCountMatrix), 10, 10) == "")], "_T1")

# match RNAseq count matrix to QCMatrix
# RNAseqCountMatrixMatched <- RNAseqCountMatrix[, match(QCMatrix$SAMPLE_ID, colnames(RNAseqCountMatrix))]
# dim(RNAseqCountMatrixMatched) # 57820   132
# write.csv(RNAseqCountMatrixMatched, file = "STAR_quantmode_counts_matrix_MATCHED_FL_132_patients_rcd25May2020.csv")
# RNAseqCountMatrixMatched <- read.csv("STAR_quantmode_counts_matrix_MATCHED_FL_132_patients_rcd25May2020.csv", 
#                                      row.names = 1)
dim(RNAseqCountMatrixMatched) # 57820   132

QCanalysisRNAseq<- QCRNAseq(RNAseqCountMatrix = RNAseqCountMatrixMatched, 
                            QCMatrix = RNAseqQCFile, 
                            BetaMatrix = BetaMatrix_T1,  
                            MvalueMatrix = MvalueMatrix_T1, 
                            TumorPurity = TumorPurity,
                            ClinicalFile = ClinicalFile_T1,
                            SurvivalFile = SurvivalFile,
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
                            ProduceImages = "No")



# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change)  probe level
# add clinical file with cluster labels
ClinicalFile_T1CLUSTERs <- cbind(ClinicalFile_T1, InfiniumClustLabels$Cluster)
colnames(ClinicalFile_T1CLUSTERs)[26] <-"CLUSTER"
dim(ClinicalFile_T1CLUSTERs) # 170  26
Output_Differential_probes_InfiniumClust <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1CLUSTERs, 
                                                                    MvalueMatrix = MvalueMatrix_T1, 
                                                                    BetaMatrix = BetaMatrix_T1,
                                                                    ProbeOrGene = "Probe", 
                                                                    ContrastColumnName = "CLUSTER", 
                                                                    RGChannelSet = Output_Data_1$RGChannelSet, 
                                                                    SampleSheet = Output_Data_1$SampleSheet, 
                                                                    ProduceImages = "No", PNGorPDF = "png")

head(top_P) # from within DifferentialMethylation function
#               logFC    AveExpr        t      P.Value    adj.P.Val        B
# cg11090352 3.408936 -1.2430269 16.01931 4.382190e-36 2.609874e-30 70.73779
# cg09568691 2.922059 -0.7210139 14.22248 5.797192e-31 1.726300e-25 59.32127
dim(top_P) # 595564      6



matchMethProbesGenes <- match(AnnotationFile %>% filter(Relation_to_Island == "Island") %>% pull("V1"), rownames(top_P))
length(match(rownames(top_P)[which((! is.na(matchMethProbesGenes)) == TRUE)], AnnotationFile$V1)) # 114772

methylationResults <- data.table(top_P[which((! is.na(matchMethProbesGenes)) == TRUE), ], 
                                 Probe = rownames(top_P[which((! is.na(matchMethProbesGenes)) == TRUE), ]),
                                 GeneName = sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Name[match(rownames(top_P)[which((! is.na(matchMethProbesGenes)) == TRUE)], AnnotationFile$V1)]))
# double check
# AnnotationFile$UCSC_RefGene_Name[grep("cg20689228", AnnotationFile$V1)]


# no zero colSums: which(colSums(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered) == 0)
MatchRNAseqwithBeta <- match(colnames(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), colnames(BetaMatrix_T1))
DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = as.matrix(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), 
                                                                   ContrastColumnName = "CLUSTER", 
                                                                   ClinicalFile = QCanalysisRNAseq$ClinicalFileMatchedSampleFiltered, 
                                                                   ClusterLabels = InfiniumClustLabels$Cluster[MatchRNAseqwithBeta], 
                                                                   ProduceImages = "Yes", 
                                                                   PNGorPDF = "png") 

# Contrast: one-two has 116 ensembl IDs (+ve)logFC and significant. 
# Contrast: one-two has 50 ensembl IDs (-ve)logFC and significant. 
head(topTagsPvalue[[1]])
#       logFC      logCPM         F       PValue          FDR      ensemblGeneId
# 1  6.111903  0.8395439 46.25519 4.308221e-09 0.0001160263 ENSG00000125144
# 2  5.600054  4.2377716 44.45528 7.360204e-09 0.0001160263 ENSG00000196611
dim(topTagsPvalue[[1]]) # 31528     6
# topTagsPvalue <- topTagsPvalue[[1]][topTagsPvalue[[1]]$PValue < 0.01, ]

ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")
RNAseqtoGene <- match(substr(topTagsPvalue[[1]]$ensemblGeneId, 1, 15), ENSMBLid$gene)
# ENSMBLid$name[RNAseqtoGene]
RNAseqresults <- data.table(topTagsPvalue[[1]], GeneName = ENSMBLid$name[RNAseqtoGene])
compareRNAseqMethylationGenes <- match(RNAseqresults$GeneName, methylationResults$GeneName)

matchedMethylationResults <- methylationResults[compareRNAseqMethylationGenes[! is.na(compareRNAseqMethylationGenes)], ]
dim(matchedMethylationResults) #   13189     8
matchedRNAseqresults <- RNAseqresults[match(matchedMethylationResults$GeneName, RNAseqresults$GeneName), ]
dim(matchedRNAseqresults) # 13189     7

# colour vector to determine if both pvalues are significant or not
colVector <- rep("NoSig", nrow(matchedMethylationResults))
for(i in 1:nrow(matchedMethylationResults)) {
  if((matchedMethylationResults$adj.P.Val[i] < 0.05) && (matchedRNAseqresults$PValue[i] < 0.05)) {
    colVector[i] <- "Sig"
  }
}
colVector <- relevel(factor(colVector), ref= "Sig")
table(colVector)
# Sig NoSig 
# 3192  9997  (previous results # 1606 14686 )



combinedResults <- data.table(matchedMethylationResults, matchedRNAseqresults, colVector)
dim(combinedResults) # 13189    16
colnames(combinedResults)[1:8] <- paste0("Methylation", colnames(combinedResults)[1:8])
colnames(combinedResults)[9:15] <- paste0("RNAseq", colnames(combinedResults)[9:15])


p2 <- ggplot(combinedResults, aes(x = MethylationlogFC, y = RNAseqlogFC, color=factor(colVector))) +
  geom_point(size = 2) +
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "Methylation FC (Island only)")+
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

plot(matchedMethylationResults$logFC, matchedRNAseqresults$logFC, col = c(colVector)+1)



# Determine genes with inverse correlation 
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
table(inverseCorVector)
# NegMeth,NegExp NegMeth,PosExp          NoSig PosMeth,NegExp PosMeth,PosExp 
# 179           1087           9997            284           1642


AA <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,PosExp"), ]$RNAseqGeneName)
BB <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,NegExp"), ]$RNAseqGeneName)
CC <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,NegExp"), ]$RNAseqGeneName)
DD <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName)
EE <-  unfactor(combinedResults[which(inverseCorVector == "NoSig"), ]$RNAseqGeneName)
length(AA) <- length(BB) <- length(CC) <- length(DD) <- length(EE)

geneList <- cbind(AA,BB,CC,DD,EE)
colnames(geneList) <- c("PosMethPosExp", "PosMethNegExp", "NegMethNegExp", "NegMethPosExp", "NoSig")
dim(geneList) # 9997    5
# write.csv(geneList, file = "FCgeneList_IslandMethProbes_allGenes_1June2020.csv")


# combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName
# combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ] %>% arrange(desc(abs(MethylationlogFC)))
# combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ] %>% arrange(desc(abs(RNAseqlogFC)))
# write.csv(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName,
#          file = "NegMethPosExp_25May2020.csv")
# combinedResultsNegMethPosExp <- data.table(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ])
# gprofilerNegMethPosExp <- gprofiler2::gost(query = list(unfactor(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName)), 
#                                            user_threshold = 0.01)


AA <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,PosExp"), ]$RNAseqGeneName)
BB <- unfactor(combinedResults[which(inverseCorVector == "PosMeth,NegExp"), ]$RNAseqGeneName)
CC <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,NegExp"), ]$RNAseqGeneName)
DD <- unfactor(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName)
EE <-  unfactor(combinedResults[which(inverseCorVector == "NoSig"), ]$RNAseqGeneName)

gprofilerPosMethPosExpAA <- gprofiler2::gost(query = list(AA), user_threshold = 0.01, 
                                           as_short_link = TRUE)
gprofilerPosMethNegExpBB <- gprofiler2::gost(query = list(BB), user_threshold = 0.01,
                                           as_short_link = TRUE)
gprofilerNegMethNegExpCC <- gprofiler2::gost(query = list(CC), user_threshold = 0.01,
                                           as_short_link = TRUE)
gprofilerNegMethPosExpDD <- gprofiler2::gost(query = list(DD), user_threshold = 0.01,
                                           as_short_link = TRUE)
gprofilerNoSigEE <- gprofiler2::gost(query = list(EE), user_threshold = 0.01,
                                   as_short_link = TRUE)


combinedResultsCorrelation <- data.table(combinedResults, "correlation" = inverseCorVector)
dim(combinedResultsCorrelation) # 13189    17


p3 <- ggplot2::ggplot(combinedResultsCorrelation[- which(inverseCorVector == "NoSig"), ], 
                      aes(x = MethylationlogFC, y = RNAseqlogFC, color=factor(correlation))) +
  geom_point(size = 2) +
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "Methylation FC (Island only)")+
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
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))

# order by methylationFC
combinedResults %>% arrange(desc(abs(MethylationlogFC)))  %>% filter(colVector == "Sig")
# order by RNAseqFC
combinedResults %>% arrange(desc(abs(RNAseqlogFC)))  %>% filter(colVector == "Sig")









# 26 May 2020
# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) via DMRcate
set.seed(1234)
DiffMethylatedRegionsInfiniumClusters <- DiffMethylatedRegions(Method = "DMRcate", 
                                                               BetaMatrix = BetaMatrix_T1, 
                                                               MvalueMatrix = MvalueMatrix_T1, 
                                                               ContrastColumnName = "CLUSTER", 
                                                               ClinicalFile = ClinicalFile_T1, 
                                                               ClusterLabels = InfiniumClustLabels$Cluster,
                                                               AnnotationFile = AnnotationFile, 
                                                               ProduceImages = "No", 
                                                               DMR = 1, ExpressionFile = NA)

# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) 

saveDetails <- matrix(NA, nrow = 1, ncol = 2)
for(i in 1:nrow(DiffMethylatedRegionsInfiniumClusters$GRangesObject[[1]]@elementMetadata)) {
  if (DiffMethylatedRegionsInfiniumClusters$GRangesObject[[1]]@elementMetadata[i,]$Stouffer < 0.05) {
    geneNames <- unlist(strsplit(DiffMethylatedRegionsInfiniumClusters$GRangesObject[[1]]@elementMetadata[i,]$overlapping.genes, ","))
    Info <- cbind(geneNames, rep(DiffMethylatedRegionsInfiniumClusters$GRangesObject[[1]]@elementMetadata[i, ]$meandiff, length(geneNames)))
    saveDetails <- rbind(saveDetails, Info)
  }
}

saveDetails <- saveDetails[- c(1),]
saveDetails <- as.data.frame(saveDetails)
colnames(saveDetails) <- c("Gene", "MeanDiff")
head(saveDetails)
dim(saveDetails) # 155955      2
saveDetails <- saveDetails[- which(is.na(saveDetails$Gene) == TRUE), ]
dim(saveDetails) # 153287      2

DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = as.matrix(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), 
                                                                   ContrastColumnName = "CLUSTER", 
                                                                   ClinicalFile = QCanalysisRNAseq$ClinicalFileMatchedSampleFiltered, 
                                                                   ClusterLabels = ClusterLabels2to4[[2]][[2]][MatchRNAseqwithBeta], 
                                                                   ProduceImages = "Yes", 
                                                                   PNGorPDF = "png") 

# Contrast: one-two has 116 ensembl IDs (+ve)logFC and significant. 
# Contrast: one-two has 50 ensembl IDs (-ve)logFC and significant. 
head(topTagsPvalue[[1]])
#       logFC      logCPM         F       PValue          FDR      ensemblGeneId
# 1  6.111903  0.8395439 46.25519 4.308221e-09 0.0001160263 ENSG00000125144
# 2  5.600054  4.2377716 44.45528 7.360204e-09 0.0001160263 ENSG00000196611
dim(topTagsPvalue[[1]]) # 31528     6
# topTagsPvalue <- topTagsPvalue[[1]][topTagsPvalue[[1]]$PValue < 0.01, ]

ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")
RNAseqtoGene <- match(substr(topTagsPvalue[[1]]$ensemblGeneId, 1, 15), ENSMBLid$gene)
# ENSMBLid$name[RNAseqtoGene]
RNAseqresults <- data.table(topTagsPvalue[[1]], GeneName = ENSMBLid$name[RNAseqtoGene])
compareRNAseqDMRcateGenes <- match(unfactor(RNAseqresults$GeneName), saveDetails$Gene)



matchedDMRcateResults <- saveDetails[compareRNAseqDMRcateGenes[! is.na(compareRNAseqDMRcateGenes)], ]
dim(matchedDMRcateResults) #   6623    2
matchedRNAseqresults <- RNAseqresults[match(matchedDMRcateResults$Gene, RNAseqresults$GeneName), ]
dim(matchedRNAseqresults) # 6623    7


# colour vector to determine if both pvalues are significant or not
colVectorDMRcate <- rep("NoSig", nrow(matchedDMRcateResults))
for(i in 1:nrow(matchedDMRcateResults)) {
  if(matchedRNAseqresults$PValue[i] < 0.05) {
    colVectorDMRcate[i] <- "Sig"
  }
}
colVectorDMRcate <- relevel(factor(colVectorDMRcate), ref= "Sig")
table(colVectorDMRcate)
# Sig NoSig 
#  1696  4927 

# Determine genes with inverse correlation 
inverseCorVectorDMRcate <- as.vector(colVectorDMRcate)
for(i in 1:nrow(matchedDMRcateResults)) {
  if(colVectorDMRcate[i] == "Sig") {
    if((unfactor(matchedDMRcateResults$MeanDiff[i]) > 0) && (matchedRNAseqresults$logFC[i] < 0)) {
      inverseCorVectorDMRcate[i] <- "PosMeth,NegExp"
    }
    
    if((unfactor(matchedDMRcateResults$MeanDiff[i]) < 0) && (matchedRNAseqresults$logFC[i] > 0)) {
      inverseCorVectorDMRcate[i] <- "NegMeth,PosExp"
    }
    
    if((unfactor(matchedDMRcateResults$MeanDiff[i]) < 0) && (matchedRNAseqresults$logFC[i] < 0)) {
      inverseCorVectorDMRcate[i] <- "NegMeth,NegExp"
    }
    if((unfactor(matchedDMRcateResults$MeanDiff[i]) > 0) && (matchedRNAseqresults$logFC[i] > 0)) {
      inverseCorVectorDMRcate[i] <- "PosMeth,PosExp"
    }
  }
}
table(inverseCorVectorDMRcate)
# NegMeth,NegExp     NegMeth,PosExp    NoSig    PosMeth,NegExp   PosMeth,PosExp 
# 84                 622               4927     151              839 


combinedResults <- data.table(matchedDMRcateResults, matchedRNAseqresults, colVectorDMRcate, inverseCorVectorDMRcate)
colnames(combinedResults)[1:2] <- paste0("DMRcate", colnames(combinedResults)[1:2])
colnames(combinedResults)[3:9] <- paste0("RNAseq", colnames(combinedResults)[3:9])
combinedResults$DMRcateMeanDiff <- round(as.numeric(unfactor(combinedResults$DMRcateMeanDiff)), 4)
typeof(combinedResults)
class(combinedResults)
dim(combinedResults) # 6623   11


DMRcateAA <- unfactor(combinedResults[which(inverseCorVectorDMRcate == "PosMeth,PosExp"), ]$RNAseqGeneName)
DMRcateBB <- unfactor(combinedResults[which(inverseCorVectorDMRcate == "PosMeth,NegExp"), ]$RNAseqGeneName)
DMRcateCC <- unfactor(combinedResults[which(inverseCorVectorDMRcate == "NegMeth,NegExp"), ]$RNAseqGeneName)
DMRcateDD <- unfactor(combinedResults[which(inverseCorVectorDMRcate == "NegMeth,PosExp"), ]$RNAseqGeneName)
DMRcateEE <-  unfactor(combinedResults[which(inverseCorVectorDMRcate == "NoSig"), ]$RNAseqGeneName)
length(DMRcateAA) <- length(DMRcateBB) <- length(DMRcateCC) <- length(DMRcateDD) <- length(DMRcateEE)

geneList <- cbind(DMRcateAA, DMRcateBB, DMRcateCC, DMRcateDD, DMRcateEE)
colnames(geneList) <- c("PosMethPosExp", "PosMethNegExp", "NegMethNegExp", "NegMethPosExp", "NoSig")
dim(geneList) # 4927    5
 write.csv(geneList, file = "FCgeneList_DMRcateAllMethProbes_allGenes_1June2020.csv")


# combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName
# combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ] %>% arrange(desc(abs(MethylationlogFC)))
# combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ] %>% arrange(desc(abs(RNAseqlogFC)))
# write.csv(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName,
#          file = "NegMethPosExp_25May2020.csv")
# combinedResultsNegMethPosExp <- data.table(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ])
# gprofilerNegMethPosExp <- gprofiler2::gost(query = list(unfactor(combinedResults[which(inverseCorVector == "NegMeth,PosExp"), ]$RNAseqGeneName)), 
#                                            user_threshold = 0.01)


DMRcateAA <- unfactor(combinedResults[which(inverseCorVectorDMRcate == "PosMeth,PosExp"), ]$RNAseqGeneName)
DMRcateBB <- unfactor(combinedResults[which(inverseCorVectorDMRcate == "PosMeth,NegExp"), ]$RNAseqGeneName)
DMRcateCC <- unfactor(combinedResults[which(inverseCorVectorDMRcate == "NegMeth,NegExp"), ]$RNAseqGeneName)
DMRcateDD <- unfactor(combinedResults[which(inverseCorVectorDMRcate == "NegMeth,PosExp"), ]$RNAseqGeneName)
DMRcateEE <-  unfactor(combinedResults[which(inverseCorVectorDMRcate == "NoSig"), ]$RNAseqGeneName)

gprofilerPosMethPosExpDMRcateAA <- gprofiler2::gost(query = list(DMRcateAA), user_threshold = 0.01, 
                                             as_short_link = TRUE)
gprofilerPosMethNegExpDMRcateBB <- gprofiler2::gost(query = list(DMRcateBB), user_threshold = 0.01,
                                             as_short_link = TRUE)
gprofilerNegMethNegExpDMRcateCC <- gprofiler2::gost(query = list(DMRcateCC), user_threshold = 0.01,
                                             as_short_link = TRUE)
gprofilerNegMethPosExpDMRcateDD <- gprofiler2::gost(query = list(DMRcateDD), user_threshold = 0.01,
                                             as_short_link = TRUE)
gprofilerNoSigDMRcateEE <- gprofiler2::gost(query = list(DMRcateEE), user_threshold = 0.01,
                                     as_short_link = TRUE)


p2 <- ggplot(combinedResults, aes(x = DMRcateMeanDiff, y = RNAseqlogFC, color=factor(colVectorDMRcate))) +
  geom_point(size = 2) +
  xlim(-0.12, 0.24) + 
  scale_y_continuous(name = "RNAseq FC") +
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
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))

p3 <- ggplot(combinedResults[which(combinedResults$colVectorDMRcate == "Sig"), ], aes(x = DMRcateMeanDiff, y = RNAseqlogFC, color=factor(inverseCorVectorDMRcate))) +
  geom_point(size = 2) +
  xlim(-0.12, 0.24) + 
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "DMRcate Mean Difference")+
  theme_bw() + 
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold") ) +
  scale_color_manual(values = c("#4292c6", "#de2d26", "#41ab5d", "#9970ab"),
                     name = "Significance")  +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))









# 1 June 2020
# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) via DMRcate
# But set a cut off for DMRcate because mean differences resulting are very close to 0
# Stouffer < 0.01 is used (instead of Stouffer < 0.05 like earlier). 
set.seed(1234)
DiffMethylatedRegionsInfiniumClusters <- DiffMethylatedRegions(Method = "DMRcate", 
                                                               BetaMatrix = BetaMatrix_T1, 
                                                               MvalueMatrix = MvalueMatrix_T1, 
                                                               ContrastColumnName = "CLUSTER", 
                                                               ClinicalFile = ClinicalFile_T1, 
                                                               ClusterLabels = InfiniumClustLabels$Cluster,
                                                               AnnotationFile = AnnotationFile, 
                                                               ProduceImages = "No", 
                                                               DMR = 1, ExpressionFile = NA)

# Creating a plot of logFC RNAseq (probe level) vs methylation logFC (log fold change) 

saveDetails <- matrix(NA, nrow = 1, ncol = 2)
for(i in 1:nrow(DiffMethylatedRegionsInfiniumClusters$GRangesObject[[1]]@elementMetadata)) {
  if (DiffMethylatedRegionsInfiniumClusters$GRangesObject[[1]]@elementMetadata[i,]$Stouffer < 0.01) {
    geneNames <- unlist(strsplit(DiffMethylatedRegionsInfiniumClusters$GRangesObject[[1]]@elementMetadata[i,]$overlapping.genes, ","))
    Info <- cbind(geneNames, rep(DiffMethylatedRegionsInfiniumClusters$GRangesObject[[1]]@elementMetadata[i, ]$meandiff, length(geneNames)))
    saveDetails <- rbind(saveDetails, Info)
  }
}

saveDetails <- saveDetails[-c(1), ]
saveDetails <- as.data.frame(saveDetails)
colnames(saveDetails) <- c("Gene", "MeanDiff")
head(saveDetails)
dim(saveDetails) # 145192       2
saveDetails <- saveDetails[- which(is.na(saveDetails$Gene) == TRUE), ]
dim(saveDetails) # 142696       2

DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = as.matrix(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), 
                                                                   ContrastColumnName = "CLUSTER", 
                                                                   ClinicalFile = QCanalysisRNAseq$ClinicalFileMatchedSampleFiltered, 
                                                                   ClusterLabels = InfiniumClustLabels$Cluster[MatchRNAseqwithBeta], 
                                                                   ProduceImages = "No", 
                                                                   PNGorPDF = "png") 

# Contrast: one-two has 116 ensembl IDs (+ve)logFC and significant. 
# Contrast: one-two has 50 ensembl IDs (-ve)logFC and significant. 
head(topTagsPvalue[[1]])
#       logFC      logCPM         F       PValue          FDR      ensemblGeneId
# 1  6.111903  0.8395439 46.25519 4.308221e-09 0.0001160263 ENSG00000125144
# 2  5.600054  4.2377716 44.45528 7.360204e-09 0.0001160263 ENSG00000196611
dim(topTagsPvalue[[1]]) # 31528     6
# topTagsPvalue <- topTagsPvalue[[1]][topTagsPvalue[[1]]$PValue < 0.01, ]

ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")
RNAseqtoGene <- match(substr(topTagsPvalue[[1]]$ensemblGeneId, 1, 15), ENSMBLid$gene) # 31528 ENSEMBLIDs from RNAseq included
# ENSMBLid$name[RNAseqtoGene]
RNAseqresults <- data.table(topTagsPvalue[[1]], GeneName = ENSMBLid$name[RNAseqtoGene])
compareRNAseqDMRcateGenes <- match(unfactor(RNAseqresults$GeneName), saveDetails$Gene)



matchedDMRcateResults <- saveDetails[compareRNAseqDMRcateGenes[! is.na(compareRNAseqDMRcateGenes)], ]
dim(matchedDMRcateResults) #   6327    2
matchedRNAseqresults <- RNAseqresults[match(matchedDMRcateResults$Gene, RNAseqresults$GeneName), ]
dim(matchedRNAseqresults) # 6327    2; # 5981 unique entries

# colour vector to determine if both pvalues are significant or not
colVectorDMRcate <- rep("NoSig", nrow(matchedDMRcateResults))
for(i in 1:nrow(matchedDMRcateResults)) {
  if(matchedRNAseqresults$PValue[i] < 0.05) {
    colVectorDMRcate[i] <- "Sig"
  }
}
colVectorDMRcate <- relevel(factor(colVectorDMRcate), ref= "Sig")
table(colVectorDMRcate)
# Sig NoSig 
# 1626  4701 

# Determine genes with inverse correlation 
inverseCorVectorDMRcate <- as.vector(colVectorDMRcate)
for(i in 1:nrow(matchedDMRcateResults)) {
  if(colVectorDMRcate[i] == "Sig") {
    if((unfactor(matchedDMRcateResults$MeanDiff[i]) > 0) && (matchedRNAseqresults$logFC[i] < 0)) {
      inverseCorVectorDMRcate[i] <- "PosMeth,NegExp"
    }
    
    if((unfactor(matchedDMRcateResults$MeanDiff[i]) < 0) && (matchedRNAseqresults$logFC[i] > 0)) {
      inverseCorVectorDMRcate[i] <- "NegMeth,PosExp"
    }
    
    if((unfactor(matchedDMRcateResults$MeanDiff[i]) < 0) && (matchedRNAseqresults$logFC[i] < 0)) {
      inverseCorVectorDMRcate[i] <- "NegMeth,NegExp"
    }
    if((unfactor(matchedDMRcateResults$MeanDiff[i]) > 0) && (matchedRNAseqresults$logFC[i] > 0)) {
      inverseCorVectorDMRcate[i] <- "PosMeth,PosExp"
    }
  }
}
table(inverseCorVectorDMRcate)
# NegMeth,NegExp NegMeth,PosExp          NoSig PosMeth,NegExp PosMeth,PosExp 
# 84            580           4701            149            813 


combinedResults <- data.table(matchedDMRcateResults, matchedRNAseqresults, colVectorDMRcate, inverseCorVectorDMRcate)
colnames(combinedResults)[1:2] <- paste0("DMRcate", colnames(combinedResults)[1:2])
colnames(combinedResults)[3:9] <- paste0("RNAseq", colnames(combinedResults)[3:9])
combinedResults$DMRcateMeanDiff <- round(as.numeric(unfactor(combinedResults$DMRcateMeanDiff)), 4)
typeof(combinedResults)
class(combinedResults)
dim(combinedResults) # 6327    11



p2 <- ggplot(combinedResults[abs(combinedResults$DMRcateMeanDiff) > 0.01, ], aes(x = DMRcateMeanDiff, y = RNAseqlogFC, color=factor(colVectorDMRcate))) +
  geom_point(size = 2) +
  xlim(-0.12, 0.24) + 
  scale_y_continuous(name = "RNAseq FC") +
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
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))


combinedResultsSig <- combinedResults[which(combinedResults$colVectorDMRcate == "Sig"), ]
dim(combinedResultsSig) # 1626   11
p3 <- ggplot(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ], aes(x = DMRcateMeanDiff, y = RNAseqlogFC, color=factor(inverseCorVectorDMRcate))) +
  geom_point(size = 2) +
  xlim(-0.12, 0.24) + 
  scale_y_continuous(name = "RNAseq FC") +
  labs(x = "DMRcate Mean Difference")+
  theme_bw() + 
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold") ) +
  scale_color_manual(values = c("#4292c6", "#de2d26", "#41ab5d", "#9970ab"),
                     name = "Significance")  +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))

DMRcateAAmeanDiff <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "PosMeth,PosExp"), ]$RNAseqGeneName)
DMRcateBBmeanDiff <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "PosMeth,NegExp"), ]$RNAseqGeneName)
DMRcateCCmeanDiff <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "NegMeth,NegExp"), ]$RNAseqGeneName)
DMRcateDDmeanDiff <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "NegMeth,PosExp"), ]$RNAseqGeneName)
DMRcateEEmeanDiff <-  unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "NoSig"), ]$RNAseqGeneName)
length(DMRcateBBmeanDiff) <- length(DMRcateCCmeanDiff) <- length(DMRcateDDmeanDiff) <- length(DMRcateAAmeanDiff)

geneList <- cbind(DMRcateAAmeanDiff, DMRcateBBmeanDiff, DMRcateCCmeanDiff, DMRcateDDmeanDiff)
colnames(geneList) <- c("PosMethPosExp", "PosMethNegExp", "NegMethNegExp", "NegMethPosExp")
dim(geneList) # 745   4
write.csv(geneList, file = "FCgeneList_DMRcate_MeanDiffGreat0.01_AllMethProbes_allGenes_15June2020.csv")



DMRcateAAmeanDiff <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "PosMeth,PosExp"), ]$RNAseqGeneName)
DMRcateBBmeanDiff <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "PosMeth,NegExp"), ]$RNAseqGeneName)
DMRcateCCmeanDiff <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "NegMeth,NegExp"), ]$RNAseqGeneName)
DMRcateDDmeanDiff <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "NegMeth,PosExp"), ]$RNAseqGeneName)

gprofilerPosMethPosExpDMRcateAAmeanDiff <- gprofiler2::gost(query = list(DMRcateAAmeanDiff), 
                                                            user_threshold = 0.01, 
                                                            as_short_link = TRUE)
gprofilerPosMethNegExpDMRcateBBmeanDiff <- gprofiler2::gost(query = list(DMRcateBBmeanDiff), 
                                                            user_threshold = 0.01,
                                                            as_short_link = TRUE)
gprofilerNegMethNegExpDMRcateCCmeanDiff <- gprofiler2::gost(query = list(DMRcateCCmeanDiff), 
                                                            user_threshold = 0.01,
                                                            as_short_link = TRUE)
gprofilerNegMethPosExpDMRcateDDmeanDiff <- gprofiler2::gost(query = list(DMRcateDDmeanDiff), 
                                                            user_threshold = 0.01,
                                                            as_short_link = TRUE)


# 15 June 2020
# Trying Robert's script
# Plot logFC chance of DMRcate vs RNAseq using P values
# When presenting pathways across different contrasts, a plot such as the one we used in synergy paper is useful.
# See https://github.com/kridel-lab/TAZ_LEN_SYNERGY/blob/master/LOCAL_SCRIPTS/FINAL_R_ANALYSIS/figure_scripts/Figure3_associated_script line 299 onwards.
PosMethPosExp <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "PosMeth,PosExp"), ]$RNAseqGeneName)
PosMethNegExp <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "PosMeth,NegExp"), ]$RNAseqGeneName)
NegMethNegExp <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "NegMeth,NegExp"), ]$RNAseqGeneName)
NegMethPosExp <- unfactor(combinedResultsSig[abs(combinedResultsSig$DMRcateMeanDiff) > 0.01, ][which(inverseCorVectorDMRcate == "NegMeth,PosExp"), ]$RNAseqGeneName)

ContrastsPosExp <- list(PosMethPosExp, NegMethPosExp)
names(ContrastsPosExp) <- c("PosMethPosExp", "NegMethPosExp")

ContrastsNegExp <- list(PosMethNegExp, NegMethNegExp)
names(ContrastsNegExp) <- c("PosMethNegExp", "NegMethNegExp")



gProfContrastsPosExp <- list()
gProfContrastsNegExp <- list()
for (i in 1:length(ContrastsPosExp)) {
gProfContrastsPosExp[[i]] <- gprofiler2::gost(query = ContrastsPosExp[i], 
                            user_threshold = 0.01)  %>%
                            .$result %>%
                            arrange(p_value) %>%
                            mutate(log.p = -log10(p_value)) %>%
                            mutate(term.name = paste(term_id, term_name, sep = "\n")) %>%
                            mutate(term.name = fct_reorder(term.name, p_value, .desc = TRUE)) %>%
                            dplyr::select(term_name, term_size, query_size, intersection_size, log.p) %>% 
                            mutate(condition = names(ContrastsPosExp[i]))
}

for (i in 1:length(ContrastsNegExp)) {
gProfContrastsNegExp[[i]] <- gprofiler2::gost(query = ContrastsNegExp[i], 
                                              user_threshold = 0.01)  %>%
                                              .$result %>%
                                              arrange(p_value) %>%
                                              mutate(log.p = -log10(p_value)) %>%
                                              mutate(term.name = paste(term_id, term_name, sep = "\n")) %>%
                                              mutate(term.name = fct_reorder(term.name, p_value, .desc = TRUE)) %>%
                                              dplyr::select(term_name, term_size, query_size, intersection_size, log.p) %>% 
                                              mutate(condition = names(ContrastsNegExp[i]))
}

posExp = data.frame(do.call(rbind, gProfContrastsPosExp))
NegExp = data.frame(do.call(rbind, gProfContrastsNegExp))

# got top 50 terms based on p value for each contrast
termNamesPosExp <- posExp %>%
                   group_by(term_name) %>% top_n(1, log.p) %>%
                   ungroup() %>% top_n(50, log.p) %>% .$term_name

termNamesNegExp <- NegExp %>%
                   group_by(term_name) %>% top_n(1, log.p) %>%
                   ungroup() %>% top_n(50, log.p) %>% .$term_name

# expanded.df_up <- pathways_up %>% tidyr::expand(term_name, condition)


plotposExp <- posExp %>%
              filter(term_name %in% termNamesPosExp) %>%
              ggplot(aes(x = factor(condition), y = term_name, fill = log.p, colour = factor(condition) )) + 
              geom_tile(colour = "white") +
              labs(x = "Category", y = "gProfiler Term Name", fill = "-log10(P value)") +
              # geom_text(aes(label = intersection_size), size = 2) +
              theme_classic() +
              scale_fill_gradientn(values = c(1.0, 0.75, 0.5, 0.25, 0), 
                                   colours = c("#67000d", "#a50f15", "#cb181d", "#ef3b2c", "white"), 
                                   na.value = "white") +  #facet_wrap(~ cell_line, scales = "free_x") + theme_classic() +
              theme(axis.title.x = element_blank(), 
                    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.text.y = element_text(size = 6),
                    panel.background = element_rect(colour = "black", size = 0.5), 
                    strip.text.x = element_text(size = 8), strip.background = element_blank())


plotNegExp <- NegExp %>%
              filter(term_name %in% termNamesNegExp) %>%
              ggplot(aes(x = factor(condition), y = term_name, fill = log.p, colour = factor(condition) )) + 
              geom_tile(colour = "white") +
              labs(x = "Category", y = "gProfiler Term Name", fill = "-log10(P value)") +
              theme_classic() +
              # geom_text(aes(label = intersection_size), size = 2) +
              scale_fill_gradientn(values = c(1.0, 0.75, 0.5, 0.25, 0), 
                                   colours = c("#023858", "#045a8d", "#0570b0", "#74a9cf", "white"), 
                                   na.value = "white") +
              theme(axis.title.x = element_blank(), 
                    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.text.y = element_text(size = 6),
                    panel.background = element_rect(colour = "black", size = 0.5), 
                    strip.text.x = element_text(size = 8), strip.background = element_blank())


gridExtra::grid.arrange(plotposExp, plotNegExp, ncol = 2)





# 25 May 2020
# Recreate RNAseq heatmap 14_Clustering_RPMM_InfiniumClust

# Pick top 5,000 most variable probes
Output_SDeviation_5000_All <- SDeviation(BetaMatrix = BetaMatrix_T1, 
                                         MvalueMatrix = MvalueMatrix_T1, 
                                         NumberofProbes = 5000)

Output_InfiniumClustering_5000_All <- Clustering(TypeofClustering = "InfiniumClust",
                                                 BetaMatrix = BetaMatrix_T1, 
                                                 MvalueMatrix = MvalueMatrix_T1, 
                                                 AnnotationFile = AnnotationFile,
                                                 ListofProbes = rownames(Output_SDeviation_5000_All$BetaMatrix_SD_Filtered), 
                                                 ClinicalFile = ClinicalFile_T1, 
                                                 TumorPurity = TumorPurity,
                                                 FigureGenerate = "No", 
                                                 PNGorPDF = "png")


RNASeqCountMatrixMatched = as.matrix(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered)
dim(RNASeqCountMatrixMatched) # 31528    58

ClusterLabels <- InfiniumClustLabels$Cluster
RNASeqCountMatrixMatched = as.matrix(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNormalized)
dim(RNASeqCountMatrixMatched) # 31528    58


as.matrix(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered)
RNASeqCountMatrixOriginalFeatureFiltered <- RNAseqCountMatrix132[match(rownames(as.matrix(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered)), rownames(RNAseqCountMatrix132)), ]
dim(RNASeqCountMatrixOriginalFeatureFiltered) # 31528   132

# ftp://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/
# https://bioinformatics.stackexchange.com/questions/4942/finding-gene-length-using-ensembl-id
txdb <- GenomicFeatures::makeTxDbFromGFF(file = 'Homo_sapiens.GRCh37.87.gff3',
                                         format = "gff3",
                                         dataSource = "ensembl",
                                         organism = "Homo sapiens")

txdb <- GenomicFeatures::makeTxDbFromEnsembl(organism="Homo sapiens",
                                             server="ensembldb.ensembl.org",
                                             username="anonymous", password=NULL, port=0L,
                                             tx_attrib=NULL)

## get gene information
all.genes <- GenomicFeatures::genes(txdb)
show(all.genes)

## import your list of gene names
my.genes <- c('ENSG00000141510','ENSG00000184571','ENSG00000011007')

## get the length of each of those genes
my.genes.lengths <- width(all.genes[my.genes])
## put the names back on the lengths
names(my.genes.lengths) <- my.genes

## print lengths
print(my.genes.lengths)

#cat("\n Saving data file ")
save.image(paste0("MethylOutput_", paste(format(Sys.time(), "%d%b%Y")), ".RData"))
# saveRDS(Output_Differential_probes_STAGE_7, "Output_Differential_probes_STAGE_7.RDS")








# 2 June 2020

TargetSeqBC <- data.table::fread(file = paste0(TargetedDNAseqDirPath, "/BC_Cancer_capseq_data.csv"))
dim(TargetSeqBC) # 380  10
table(TargetSeqBC$Mutation_Type)
#  Frame_Shift_Del   Frame_Shift_Ins      In_Frame_Del Missense_Mutation Nonsense_Mutation       Splice_Site 
#  14                 3                 1               300                41                21 
table(TargetSeqBC$Protein_Change)
length(unique(TargetSeqBC$SAMPLE_ID)) # 31
length(unique(TargetSeqBC$Hugo_Symbol)) # 66
# saveRDS(TargetSeqBC, file = "BC_Cancer_capseq_data.RDS")


# Creat an empty matrix with 0 for saving muation data
TargetSeqBCmutations01 <- matrix(0, nrow = length(unique(TargetSeqBC$SAMPLE_ID)), 
                                 ncol = length(unique(TargetSeqBC$Hugo_Symbol)))
dim(TargetSeqBCmutations01) # 31 66
rownames(TargetSeqBCmutations01) <- unique(TargetSeqBC$SAMPLE_ID)
colnames(TargetSeqBCmutations01) <- unique(TargetSeqBC$Hugo_Symbol)


# In the 0 matrix, only write 1 if the muation is present for protein coding (i.e., ignore if NA)
# On 3 June 2020, Robert asked to include NA entries as well. As a result, altered code. This code i marked by **
for (i in 1:length(unique(TargetSeqBC$SAMPLE_ID))) {
  for (j in 1:nrow(filter(TargetSeqBC, SAMPLE_ID == unique(TargetSeqBC$SAMPLE_ID)[i]))) {
    # ** if(! is.na(filter(TargetSeqBC, SAMPLE_ID == unique(TargetSeqBC$SAMPLE_ID)[i])$Protein_Change[j])) { 
    # ** if a protein change is NOT present; e.g., if a Splice_Site is present, don't use that
      colName <- match(filter(TargetSeqBC, SAMPLE_ID == unique(TargetSeqBC$SAMPLE_ID)[i])$Hugo_Symbol[j], colnames(TargetSeqBCmutations01))
      TargetSeqBCmutations01[i, colName] <- 1
    # ** }
  }
}


TargetSeqBCmutations01ClusterLabs <- InfiniumClustLabels$Cluster[match(rownames(TargetSeqBCmutations01), colnames(BetaMatrix_T1))]
TargetSeqBCmutations01ClusterLabs <- data.frame(cbind(TargetSeqBCmutations01ClusterLabs, TargetSeqBCmutations01))
colnames(TargetSeqBCmutations01ClusterLabs) <- c("Cluster", colnames(TargetSeqBCmutations01))
dim(TargetSeqBCmutations01ClusterLabs) # 31 67
head(TargetSeqBCmutations01ClusterLabs)
# write.csv(TargetSeqBCmutations01ClusterLabs, file = "BC_Cancer_capseq_data_with01ClusterLabs_AS_2June2020.csv")
# TargetSeqBCmutations01ClusterLabs <- read.csv("BC_Cancer_capseq_data_with01ClusterLabs_AS_2June2020.csv", row.names = 1)
# saveRDS(TargetSeqBCmutations01ClusterLabs, file = "BC_Cancer_capseq_data_with01ClusterLabs_AS_2June2020.rds")

chisq.test(y = TargetSeqBCmutations01ClusterLabs$Cluster, x = TargetSeqBCmutations01ClusterLabs$CREBBP)
# X-squared = 7.7149e-32, df = 1, p-value = 1

chisq.test(y = TargetSeqBCmutations01ClusterLabs$Cluster, x = TargetSeqBCmutations01ClusterLabs$EZH2)
# X-squared = 0.58941, df = 1, p-value = 0.4426

chisq.test(y = TargetSeqBCmutations01ClusterLabs$Cluster, x = TargetSeqBCmutations01ClusterLabs$BCL2)
# X-squared = 5.0748, df = 1, p-value = 0.02428

chisq.test(y = TargetSeqBCmutations01ClusterLabs$Cluster, x = TargetSeqBCmutations01ClusterLabs$SOCS1)
# X-squared = 1.3875e-30, df = 1, p-value = 1


TargetSeqBCmutations01ClusterLabsSignficance <- matrix(NA, ncol = 9, nrow = ncol(TargetSeqBCmutations01ClusterLabs))
colnames(TargetSeqBCmutations01ClusterLabsSignficance) <- c("Xsquared", "XsquaredP.Value", 
                                                            "XsquaredRes0-C1", "XsquaredRes1-C1", 
                                                            "XsquaredRes0-C2", "XsquaredRes1-C2", 
                                                            "fisherTOddsRatio", "fisherTP.Value",
                                                            "Frequency")
rownames(TargetSeqBCmutations01ClusterLabsSignficance) <- colnames(TargetSeqBCmutations01ClusterLabs)
for (k in 1:ncol(TargetSeqBCmutations01ClusterLabs)) {
  if(! all(TargetSeqBCmutations01ClusterLabs[, k] == 0)) { # only if all entries are not zero
    chisqTest  <- chisq.test(x = TargetSeqBCmutations01ClusterLabs[, k], y = TargetSeqBCmutations01ClusterLabs$Cluster)
    TargetSeqBCmutations01ClusterLabsSignficance[k, 1] <- chisqTest$statistic
    TargetSeqBCmutations01ClusterLabsSignficance[k, 2] <- chisqTest$p.value
    TargetSeqBCmutations01ClusterLabsSignficance[k, 3] <- chisqTest$residuals[1, 1]
    TargetSeqBCmutations01ClusterLabsSignficance[k, 4] <- chisqTest$residuals[2, 1]
    TargetSeqBCmutations01ClusterLabsSignficance[k, 5] <- chisqTest$residuals[1, 2]
    TargetSeqBCmutations01ClusterLabsSignficance[k, 6] <- chisqTest$residuals[2, 2]
    fisherTest <- fisher.test(x = TargetSeqBCmutations01ClusterLabs[, k], y = TargetSeqBCmutations01ClusterLabs$Cluster)
    TargetSeqBCmutations01ClusterLabsSignficance[k, 7] <- fisherTest$estimate
    TargetSeqBCmutations01ClusterLabsSignficance[k, 8] <- fisherTest$p.value
    TargetSeqBCmutations01ClusterLabsSignficance[k, 9] <- (length(which(TargetSeqBCmutations01ClusterLabs[, k] == 1)) / 
                                                            length(TargetSeqBCmutations01ClusterLabs[, k]))
  }
}

# remove first row as it is for cluster against cluster results
TargetSeqBCmutations01ClusterLabsSignficance <- TargetSeqBCmutations01ClusterLabsSignficance [-1,]
dim(TargetSeqBCmutations01ClusterLabsSignficance) # 66  9
head(TargetSeqBCmutations01ClusterLabsSignficance)
which(TargetSeqBCmutations01ClusterLabsSignficance[, 2] < 0.05) # BCL2
which(TargetSeqBCmutations01ClusterLabsSignficance[, 8] < 0.05) # BCL2
# write.csv(round(TargetSeqBCmutations01ClusterLabsSignficance, 4), file = "MutationSignificance.csv")
# 4 June 2020
BCL2mutCluster <- table(x = TargetSeqBCmutations01ClusterLabs$BCL2, y = TargetSeqBCmutations01ClusterLabs$Cluster)
fisher.test(BCL2mutCluster[c(2, 1), ])


# Arrange clinical file by BC 31 cases - to look at BCL2 translocations vs BCL2 mutations
BCL2mutBCL2trans <- table(ClinicalFile_T1[match(rownames(TargetSeqBCmutations01ClusterLabs), ClinicalFile_T1$SAMPLE_ID), ]$TRANSLOC_14_18, TargetSeqBCmutations01ClusterLabs$BCL2)
fisher.test(BCL2mutBCL2trans[c(2, 1), c(2, 1)])

# Arrange clinical file by BC 31 cases - to look at BCL2 translocations vs cluster labels
TransloCluster <- table(factor(ClinicalFile_T1[match(rownames(TargetSeqBCmutations01ClusterLabs), ClinicalFile_T1$SAMPLE_ID), ]$TRANSLOC_14_18), TargetSeqBCmutations01ClusterLabs$Cluster)
fisher.test(TransloCluster[c(2,1),])

# Arrange clinical file by BC 31 cases - to look at stage  vs BCL2 mutations
BCL2mutStage <- table(TargetSeqBCmutations01ClusterLabs$BCL2, ClinicalFile_T1[match(rownames(TargetSeqBCmutations01ClusterLabs), ClinicalFile_T1$SAMPLE_ID), ]$STAGE)
fisher.test(BCL2mutStage) # error









# 3 June 2020 - mutation calling from RNA seq file from Karin
# 4 June 2020 - updated based on Robert's recommendation to filter genes (69) absed on panels we submitted to BC 
# 8 June 202 - Use the new matrix Karin sent. Slack message " I missed one patient before so now they are included otherwise shouldn't be a big difference"
# link to panel: https://universityhealthnetwork-my.sharepoint.com/:x:/g/personal/robert_kridel_uhn_ca/EZoqOwj0zEVEv09Na1h02ZwBQMhTMGMKkN_vWhncg-XHpQ?e=iKfvnV
# 15 June 2020 - Update based on Teams meet with Robert 10 June 2020; remove "downstream", "upstream" 
# "exonic unknowns"; "UTR3" and "UTR5" 
#RNAseqMutationCalls <- data.table::fread(file = paste0(RNAseqDirPath, "/2020-06-02_opossum_variant_FL_rna-seq_filtered.txt"))
# RNAseqMutationCalls <- data.table::fread(file = paste0(RNAseqDirPath, "/2020-06-08_opossum_variant_FL_rna-seq_filtered.txt"))
dim(RNAseqMutationCalls) # 150423     73 (previously 99948    73)
head(RNAseqMutationCalls)
table(RNAseqMutationCalls$Variant_Type)
# SNP 
# 150423  
table(RNAseqMutationCalls$Variant_Classification)
#                  downstream .     exonic frameshift_deletion    exonic frameshift_insertion 
#  15427                           3121                           1540  
# exonic nonframeshift_deletion    exonic nonframeshift_insertion       exonic nonsynonymous_SNV 
# 2215                             44                           24383  
# exonic stopgain                 exonic stoploss                 exonic unknown 
# 2607                            172                            430 
# splicing .                      upstream .                         UTR3 . 
# 1835                            2794                          86367  
# UTR5 . 
# 9488 

length(unique(RNAseqMutationCalls$Hugo_Symbol)) # 9829
length(unique(RNAseqMutationCalls$SAMPLE_ID)) # 128
RNAseqMutationCalls$Hugo_Symbol[which(RNAseqMutationCalls$Variant_Classification == "exonic unknown")]
# write.csv(RNAseqMutationCalls$Hugo_Symbol[which(RNAseqMutationCalls$Variant_Classification == "exonic unknown")]
#, file = "MutationSignificance.csv")

# gprofilerExonicUnknown <- gprofiler2::gost(query = list(RNAseqMutationCalls$Hugo_Symbol[which(RNAseqMutationCalls$Variant_Classification == "exonic unknown")]), 
#                                           user_threshold = 0.01,
#                                           as_short_link = TRUE)

# Remove UTRs based on Robert's recommendation from June 2nd on slack (group convo Karin, Robert)
RNAseqMutationCallsFiltered <- RNAseqMutationCalls %>% filter(Variant_Classification != "UTR5 .") %>% 
                                                       filter(Variant_Classification != "UTR3 .") %>% 
                                                       filter(Variant_Classification != "downstream .") %>% 
                                                       filter(Variant_Classification != "upstream .") %>% 
                                                       filter(Variant_Classification != "exonic unknown") 

dim(RNAseqMutationCallsFiltered) # 35917    73
length(unique(RNAseqMutationCallsFiltered$Hugo_Symbol)) # 7544
length(unique(RNAseqMutationCallsFiltered$SAMPLE_ID)) # 124 
table(RNAseqMutationCallsFiltered$Variant_Classification)
# exonic frameshift_deletion    exonic frameshift_insertion  exonic nonframeshift_deletion exonic nonframeshift_insertion 
# 3121                           1540                           2215                              44 
# exonic nonsynonymous_SNV      exonic stopgain              exonic stoploss                     splicing . 
# 24383                           2607                             172                           1835  


# Creat an empty matrix with 0 for saving muation data
RNAseqMutationCallsFiltered01 <- matrix(0, nrow = length(unique(RNAseqMutationCallsFiltered$SAMPLE_ID)), 
                                      ncol = length(unique(RNAseqMutationCallsFiltered$Hugo_Symbol)))
FrequencyofRNAseqMutations <- RNAseqMutationCallsFiltered01 # create a table to keep track of frequency of mutations
dim(RNAseqMutationCallsFiltered01) #  124 7544
dim(FrequencyofRNAseqMutations) # 124 7544
rownames(RNAseqMutationCallsFiltered01) <- rownames(FrequencyofRNAseqMutations) <- unique(RNAseqMutationCallsFiltered$SAMPLE_ID)
colnames(RNAseqMutationCallsFiltered01) <- colnames(FrequencyofRNAseqMutations) <- unique(RNAseqMutationCallsFiltered$Hugo_Symbol)


# In the 0 matrix, only write 1 if the muation is present 
for (i in 1:length(unique(RNAseqMutationCallsFiltered$SAMPLE_ID))) {
  cat("\n Running i:", i)
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
dim(RNAseqMutationCallsFiltered01) #  124 7544
RNAseqMutationCallsFiltered01ClusterLabs <- InfiniumClustLabels$Cluster[match(rownames(RNAseqMutationCallsFiltered01), colnames(BetaMatrix_T1))]
FrequencyofRNAseqMutationsClusterLabs <- data.frame(cbind(RNAseqMutationCallsFiltered01ClusterLabs, FrequencyofRNAseqMutations))
RNAseqMutationCallsFiltered01ClusterLabs <- data.frame(cbind(RNAseqMutationCallsFiltered01ClusterLabs, RNAseqMutationCallsFiltered01))
dim(RNAseqMutationCallsFiltered01ClusterLabs) #  124 7545
dim(FrequencyofRNAseqMutationsClusterLabs) #  124 7545
colnames(RNAseqMutationCallsFiltered01ClusterLabs) <- c("Cluster", colnames(RNAseqMutationCallsFiltered01ClusterLabs)[-1])
colnames(FrequencyofRNAseqMutationsClusterLabs) <- c("Cluster", colnames(FrequencyofRNAseqMutationsClusterLabs)[-1])
dim(FrequencyofRNAseqMutationsClusterLabs) # 124 7545
dim(RNAseqMutationCallsFiltered01ClusterLabs) # 124 7545
head(RNAseqMutationCallsFiltered01ClusterLabs)
# write.csv(RNAseqMutationCallsFiltered01ClusterLabs, file = "2020-06-22_opossum_variant_FL_rna-seq_filtered_with01ClusterLabs_AS_22June2020.csv")
# saveRDS(RNAseqMutationCallsFiltered01ClusterLabs, file = "2020-06-22_opossum_variant_FL_rna-seq_filtered_with01ClusterLabs_AS_22June2020.rds")

# RNAseqMutationCallsFiltered01ClusterLabs <- readRDS(file = paste0(RNAseqDirPath, "/2020-06-22_opossum_variant_FL_rna-seq_filtered_with01ClusterLabs_AS_22June2020.rds"))
dim(RNAseqMutationCallsFiltered01ClusterLabs) # 124 7545

# Check entries of FrequencyofRNAseqMutations and RNAseqMutationCallsFiltered01
head(colnames(FrequencyofRNAseqMutations)) # "CFLAR"    "CCDC109B" "CREBBP"   "LYPLA2"   "BTBD7"    "TTC19"
which(FrequencyofRNAseqMutationsClusterLabs[, 2] > 0)  # for first gene: CFLAR
which(FrequencyofRNAseqMutationsClusterLabs[, 100] > 0)

head(colnames(RNAseqMutationCallsFiltered01))
which(RNAseqMutationCallsFiltered01ClusterLabs[, 2] == 1) 
which(RNAseqMutationCallsFiltered01ClusterLabs[, 100] == 1) 


# calculating for all patients
RNAseqMutationCalls01ClusterLabsSignficance <- matrix(NA, ncol = 11, nrow = ncol(RNAseqMutationCallsFiltered01ClusterLabs))
colnames(RNAseqMutationCalls01ClusterLabsSignficance) <- c("Xsquared", "XsquaredP.Value", 
                                                           "XsquaredRes0-C1", "XsquaredRes1-C1", 
                                                           "XsquaredRes0-C2", "XsquaredRes1-C2", 
                                                           "fisherTOddsRatio", "fisherTP.Value",
                                                           "MinOneMutation",
                                                           "FrequencyMinOneMutation", 
                                                           "NumAllMutations")
rownames(RNAseqMutationCalls01ClusterLabsSignficance) <- colnames(RNAseqMutationCallsFiltered01ClusterLabs)
dim(RNAseqMutationCalls01ClusterLabsSignficance) #  7545    11
for (k in 1:ncol(RNAseqMutationCallsFiltered01ClusterLabs)) { # for each gene/mutation
  cat("\nk is:", k)
  if(! all(RNAseqMutationCallsFiltered01ClusterLabs[, k] == 0, na.rm = T)) { # only if all entries are not zero
    Table <- table(x = RNAseqMutationCallsFiltered01ClusterLabs[, k], 
                   y = RNAseqMutationCallsFiltered01ClusterLabs$Cluster)
    chisqTest  <- chisq.test(Table[c(2, 1), ])
    fisherTest <- fisher.test(Table[c(2, 1), ])
    
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
head(RNAseqMutationCalls01ClusterLabsSignficance)
# remove first row as it is for cluster against cluster results
RNAseqMutationCalls01ClusterLabsSignficance <- RNAseqMutationCalls01ClusterLabsSignficance [-1, ]
dim(RNAseqMutationCalls01ClusterLabsSignficance) # 7544    11
head(RNAseqMutationCalls01ClusterLabsSignficance)

names(which(RNAseqMutationCalls01ClusterLabsSignficance[, 2] < 0.05)) # XsquaredP.Value; 60; 
# "CFLAR"   "GDI2"    "CDK13"   "NFX1"    "EZH2"    "ZFAND5"  "LUC7L3"  "FRG1"    "FAM193A" "IRF5"    "SP110"   "CERS5"   "SRP9"    "UBXN4"  
#  "DLAT"    "COL6A3"  "IGFBP7"  "UBXN7"   "UGP2"    "BCL2"    "AFF1"    "RPL37A"  "NBPF16"  "CORO1C"  "SLC11A2" "SGK1"    "TPT1"    "SKAP1"  
#  "FRG1B"   "SPSB1"   "ZNF669"  "RAPGEF5" "NAE1"    "MRPL30"  "CAPZA2"  "IFI30"   "RPS15A"  "DIP2A"   "RBBP8"   "CENPF"   "VEZT"    "DLG1"   
#  "RBL1"    "IRF2BP2" "MARCH1"  "COL3A1"  "METTL2A" "TMC6"    "EBF1"    "ROCK2"   "POLG"    "RNF19A"  "ROCK1"   "MYH9"    "USP6NL"  "CBLB"   
#  "NPIPA1"  "SRPK1"   "RAD50"   "RASA2"
length(unique(which(RNAseqMutationCalls01ClusterLabsSignficance[, 2] < 0.05))) # 60; based on XsquaredP.Value
which(RNAseqMutationCalls01ClusterLabsSignficance[, 8] < 0.05) # 
length(which(RNAseqMutationCalls01ClusterLabsSignficance[, 8] < 0.05)) # 110; fisherTP.Value
names(which(RNAseqMutationCalls01ClusterLabsSignficance[, 8] < 0.05)) 
RNAseqMutationCalls01ClusterLabsSignficance[which(RNAseqMutationCalls01ClusterLabsSignficance[, 8] < 0.05), ] # fisherTP.Value
# write.csv(round(RNAseqMutationCalls01ClusterLabsSignficance[which(RNAseqMutationCalls01ClusterLabsSignficance[, 8] < 0.05), ], 4), file = "MutationSignificance.csv")
RNAseqMutationCalls01ClusterLabsSignficance[which(RNAseqMutationCalls01ClusterLabsSignficance[, 8] < 0.05), ]
# look at overlap in elements
length(intersect(names(which(RNAseqMutationCalls01ClusterLabsSignficance[, 8] < 0.05)), 
                 names(which(RNAseqMutationCalls01ClusterLabsSignficance[, 2] < 0.05))))


EZH2Cluster <- table(RNAseqMutationCallsFiltered01ClusterLabs$EZH2, RNAseqMutationCallsFiltered01ClusterLabs$Cluster)
fisher.test(EZH2Cluster[c(2,1), ])


# Filter based on the panels sent to BC
# https://universityhealthnetwork-my.sharepoint.com/:x:/g/personal/robert_kridel_uhn_ca/EZoqOwj0zEVEv09Na1h02ZwBQMhTMGMKkN_vWhncg-XHpQ?e=iKfvnV
PM_FL_PANEL <- data.frame(read.csv(file = "PM_FL_PANEL_4June2020.csv", header = 1))
matchPMwithRNAseqMutCalls <- match(c(unfactor(PM_FL_PANEL$X..Refseq[1:70]), "BCL6", "BCL2"), rownames(RNAseqMutationCalls01ClusterLabsSignficance))
RNAseqMutationCalls01ClusterLabsSignficancePMPanelFiltered <- RNAseqMutationCalls01ClusterLabsSignficance[matchPMwithRNAseqMutCalls[! is.na(matchPMwithRNAseqMutCalls)], ]
dim(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFiltered) # 63  9
which(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFiltered[, 2] < 0.05) # EBF1 EZH2 IRF5 SGK1 BCL2
which(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFiltered[, 8] < 0.05) # EBF1 EZH2 IRF5 SGK1 BCL2 
# write.csv(round(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFiltered, 4), file = "MutationSignificance.csv")
RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlotting <- data.frame(
           RNAseqMutationCalls01ClusterLabsSignficancePMPanelFiltered, 
           NegLogOdds = -log2(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFiltered[, 7]),
           NegLogPval = -log10(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFiltered[, 8]))

library(ggplot2)
# Basic dot plot

# plot P value against odds ratio for 60 matching genes of PM_FL_PANEL
ggplot(data = RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlotting, 
       aes(x = NegLogOdds, y = NegLogPval, size = FrequencyMinOneMutation)) + 
       geom_point() + scale_y_continuous(name = "-log10(P value)") +
       # geom_text(aes(label = rownames(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlotting)), 
       #               hjust = 0, vjust = 0) + 
       ggrepel::geom_label_repel(aes(label = rownames(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlotting)),
                                 box.padding   = 0.35, 
                                 point.padding = 0.5,
                                 segment.color = 'grey50') +
       labs(x = "-log2(odds ratio)") +
       theme_bw() + 
       theme(text = element_text(size=20), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())





# subset mutation frequencies based on patients from T1 (81 patients)
# get the 81 patients
RNASeqCountMatrixMatched <- data.frame(readRDS(file = paste0(RNAseqDirPath, "2020-06-18STAR_quantmode_counts_matrix_FL_132_patients.RDS")))
dim(RNASeqCountMatrixMatched) # 57820   132
dim(RNAseqQCFile) # 132  60
RNAseqQC18June2020T1Samples81 <- QCRNAseq(RNAseqCountMatrix = RNASeqCountMatrixMatched, 
                                          QCMatrix = RNAseqQCFile, 
                                          BetaMatrix = BetaMatrix_T1,  
                                          MvalueMatrix = MvalueMatrix_T1, 
                                          TumorPurity = TumorPurity,
                                          ClinicalFile = ClinicalFile_T1,
                                          SurvivalFile = SurvivalFile,
                                          RNAseqSampleCufoffUniqMapReadCount = 10000000, 
                                          RNAseqSampleCufoffUniqMapReadPercentage = 50,
                                          RNAseqSampleCufoffReadsMappedMultipleLoci = 20,
                                          RNAseqSampleCutoffRRNAcontam = 40, 
                                          RNAseqFeatureSelectionMethod = "edgeR",
                                          RNAseqFeatureSelectionCutOff = NA,
                                          RNAseqFeatureSelectionNumberofProbes = NA,
                                          RNAseqNormalizationMethod = "edgeR",
                                          ImageName = "35_QCRNAseq_T2",
                                          PNGorPDF = "png",
                                          ProduceImages = "Yes")
dim(RNAseqQC18June2020T1Samples81$QCMatrixMatchedSampleFiltered) # 81 60
RNAseqQC18June2020T1Samples81$QCMatrixMatchedSampleFiltered$SAMPLE_ID

# Filter to tier from initial matrix 
# dim(RNAseqMutationCallsFiltered01ClusterLabs); initial matrix
RNAseqMutationCallsFiltered01ClusterLabsT1 <- RNAseqMutationCallsFiltered01ClusterLabs[match(RNAseqQC18June2020T1Samples81$QCMatrixMatchedSampleFiltered$SAMPLE_ID, 
                                                                                           rownames(RNAseqMutationCallsFiltered01ClusterLabs)), ]
dim(RNAseqMutationCallsFiltered01ClusterLabsT1) # 81 7545
FrequencyofRNAseqMutationsClusterLabsT1 <- FrequencyofRNAseqMutationsClusterLabs[match(RNAseqQC18June2020T1Samples81$QCMatrixMatchedSampleFiltered$SAMPLE_ID, 
                                                                                          rownames(FrequencyofRNAseqMutationsClusterLabs)), ]
dim(FrequencyofRNAseqMutationsClusterLabsT1) # 81 7545
# write.csv(RNAseqMutationCallsFiltered01ClusterLabsT1, file = "2020-06-08_opossum_variant_FL_rna-seq_filtered_with01ClusterLabs_T1OnlySAMPLES_AS_16June2020_.csv")

RNAseqMutationCalls01ClusterLabsSignficanceT1 <- matrix(NA, ncol = 11, nrow = ncol(RNAseqMutationCallsFiltered01ClusterLabsT1))
colnames(RNAseqMutationCalls01ClusterLabsSignficanceT1) <- c("Xsquared", "XsquaredP.Value", 
                                                           "XsquaredRes0-C1", "XsquaredRes1-C1", 
                                                           "XsquaredRes0-C2", "XsquaredRes1-C2", 
                                                           "fisherTOddsRatio", "fisherTP.Value",
                                                           "MinOneMutation",
                                                           "FrequencyMinOneMutation", 
                                                           "NumAllMutations")
rownames(RNAseqMutationCalls01ClusterLabsSignficanceT1) <- colnames(RNAseqMutationCallsFiltered01ClusterLabsT1)
dim(RNAseqMutationCalls01ClusterLabsSignficanceT1) #  7545    11
for (k in 1:ncol(RNAseqMutationCallsFiltered01ClusterLabsT1)) { # for each gene/mutation
  cat("\nk is:", k)
  if(! all(RNAseqMutationCallsFiltered01ClusterLabsT1[, k] == 0, na.rm = T)) { # only if all entries are not zero
    Table <- table(x = RNAseqMutationCallsFiltered01ClusterLabsT1[, k], 
                   y = RNAseqMutationCallsFiltered01ClusterLabsT1$Cluster)
    chisqTest  <- chisq.test(Table[c(2, 1), ])
    fisherTest <- fisher.test(Table[c(2, 1), ])
    
    RNAseqMutationCalls01ClusterLabsSignficanceT1[k, 1] <- chisqTest$statistic
    RNAseqMutationCalls01ClusterLabsSignficanceT1[k, 2] <- chisqTest$p.value
    RNAseqMutationCalls01ClusterLabsSignficanceT1[k, 3] <- chisqTest$residuals[1, 1]
    RNAseqMutationCalls01ClusterLabsSignficanceT1[k, 4] <- chisqTest$residuals[2, 1]
    RNAseqMutationCalls01ClusterLabsSignficanceT1[k, 5] <- chisqTest$residuals[1, 2]
    RNAseqMutationCalls01ClusterLabsSignficanceT1[k, 6] <- chisqTest$residuals[2, 2]
    RNAseqMutationCalls01ClusterLabsSignficanceT1[k, 7] <- fisherTest$estimate
    RNAseqMutationCalls01ClusterLabsSignficanceT1[k, 8] <- fisherTest$p.value
    RNAseqMutationCalls01ClusterLabsSignficanceT1[k, 9] <- length(which(RNAseqMutationCallsFiltered01ClusterLabsT1[, k] == 1))
    RNAseqMutationCalls01ClusterLabsSignficanceT1[k, 10] <- (length(which(RNAseqMutationCallsFiltered01ClusterLabsT1[, k] == 1)) / 
                                                             length(RNAseqMutationCallsFiltered01ClusterLabsT1[, k]))
    RNAseqMutationCalls01ClusterLabsSignficanceT1[k, 11] <- (sum(FrequencyofRNAseqMutationsClusterLabsT1[, k], na.rm = TRUE))
  }
}
head(RNAseqMutationCalls01ClusterLabsSignficanceT1)
# remove first row as it is for cluster against cluster results
RNAseqMutationCalls01ClusterLabsSignficanceT1 <- RNAseqMutationCalls01ClusterLabsSignficanceT1 [-1, ]
dim(RNAseqMutationCalls01ClusterLabsSignficanceT1) # 7544    11
head(RNAseqMutationCalls01ClusterLabsSignficanceT1)

names(which(RNAseqMutationCalls01ClusterLabsSignficanceT1[, 2] < 0.05)) # XsquaredP.Value; ; 
length(unique(which(RNAseqMutationCalls01ClusterLabsSignficanceT1[, 2] < 0.05))) # 34; based on XsquaredP.Value
which(RNAseqMutationCalls01ClusterLabsSignficanceT1[, 8] < 0.05) # 
length(which(RNAseqMutationCalls01ClusterLabsSignficanceT1[, 8] < 0.05)) # 62; fisherTP.Value
names(which(RNAseqMutationCalls01ClusterLabsSignficanceT1[, 8] < 0.05)) 
RNAseqMutationCalls01ClusterLabsSignficanceT1[which(RNAseqMutationCalls01ClusterLabsSignficanceT1[, 8] < 0.05), ] # fisherTP.Value
# write.csv(round(RNAseqMutationCalls01ClusterLabsSignficanceT1[which(RNAseqMutationCalls01ClusterLabsSignficanceT1[, 8] < 0.05), ], 4), file = "MutationSignificance.csv")
RNAseqMutationCalls01ClusterLabsSignficanceT1[which(RNAseqMutationCalls01ClusterLabsSignficanceT1[, 8] < 0.05), ]
# look at overlap in elements
length(intersect(names(which(RNAseqMutationCalls01ClusterLabsSignficanceT1[, 8] < 0.05)), 
                 names(which(RNAseqMutationCalls01ClusterLabsSignficanceT1[, 2] < 0.05))))


EZH2Cluster <- table(RNAseqMutationCallsFiltered01ClusterLabsT1$EZH2, RNAseqMutationCallsFiltered01ClusterLabsT1$Cluster)
fisher.test(EZH2Cluster[c(2,1), ])


# Filter based on the panels sent to BC
# https://universityhealthnetwork-my.sharepoint.com/:x:/g/personal/robert_kridel_uhn_ca/EZoqOwj0zEVEv09Na1h02ZwBQMhTMGMKkN_vWhncg-XHpQ?e=iKfvnV
PM_FL_PANEL <- data.frame(read.csv(file = "PM_FL_PANEL_4June2020.csv", header = 1))
matchPMwithRNAseqMutCalls <- match(c(unfactor(PM_FL_PANEL$X..Refseq[1:70]), "BCL6", "BCL2"), rownames(RNAseqMutationCalls01ClusterLabsSignficanceT1))
RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredT1 <- RNAseqMutationCalls01ClusterLabsSignficanceT1[matchPMwithRNAseqMutCalls[! is.na(matchPMwithRNAseqMutCalls)], ]
dim(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredT1) # 63  11
which(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredT1[, 2] < 0.05) # CARD11   EBF1   IRF5   BCL2
which(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredT1[, 8] < 0.05) #  BTG2 CARD11   EBF1   IRF5   SGK1   BCL2
# write.csv(round(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredT1, 4), file = "MutationSignificance.csv")
RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlottingT1 <- data.frame(
  RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredT1, 
  NegLogOdds = -log2(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredT1[, 7]),
  NegLogPval = -log10(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredT1[, 8]))

library(ggplot2)
# Basic dot plot

# plot P value against odds ratio for 60 matching genes of PM_FL_PANEL
ggplot(data = RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlottingT1, 
       aes(x = NegLogOdds, y = NegLogPval, size = FrequencyMinOneMutation)) + 
  geom_point() + scale_y_continuous(name = "-log10(P value)") +
  # geom_text(aes(label = rownames(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlotting)), 
  #               hjust = 0, vjust = 0) + 
  ggrepel::geom_label_repel(aes(label = rownames(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlottingT1)),
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            segment.color = 'grey50') +
  labs(x = "-log2(odds ratio)") +
  theme_bw() + 
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())








# compare RNA expression of FOXP1 between C1 and C1 in FOXP1; EZH2 
# with Tier 1 patients: 58
RNAseqQC18June2020T1Samples81 <- QCRNAseq(RNAseqCountMatrix = RNASeqCountMatrixMatched, 
                                          QCMatrix = RNAseqQCFile, 
                                          BetaMatrix = BetaMatrix_T1,  
                                          MvalueMatrix = MvalueMatrix_T1, 
                                          TumorPurity = TumorPurity,
                                          ClinicalFile = ClinicalFile_T1,
                                          SurvivalFile = SurvivalFile,
                                          RNAseqSampleCufoffUniqMapReadCount = 10000000, 
                                          RNAseqSampleCufoffUniqMapReadPercentage = 70,
                                          RNAseqSampleCufoffReadsMappedMultipleLoci = 20,
                                          RNAseqSampleCutoffRRNAcontam = 40, 
                                          RNAseqFeatureSelectionMethod = "edgeR",
                                          RNAseqFeatureSelectionCutOff = NA,
                                          RNAseqFeatureSelectionNumberofProbes = NA,
                                          RNAseqNormalizationMethod = "edgeR",
                                          ImageName = "35_QCRNAseq_T1",
                                          PNGorPDF = "png",
                                          ProduceImages = "Yes")


RNAseqQC18June2020T2Samples104 <- QCRNAseq(RNAseqCountMatrix = RNASeqCountMatrixMatched, 
                                           QCMatrix = RNAseqQCFile, 
                                           BetaMatrix = BetaMatrix_T1,  
                                           MvalueMatrix = MvalueMatrix_T1, 
                                           TumorPurity = TumorPurity,
                                           ClinicalFile = ClinicalFile_T1,
                                           SurvivalFile = SurvivalFile,
                                           RNAseqSampleCufoffUniqMapReadCount = 10000000, 
                                           RNAseqSampleCufoffUniqMapReadPercentage = 50,
                                           RNAseqSampleCufoffReadsMappedMultipleLoci = 20,
                                           RNAseqSampleCutoffRRNAcontam = 40, 
                                           RNAseqFeatureSelectionMethod = "edgeR",
                                           RNAseqFeatureSelectionCutOff = NA,
                                           RNAseqFeatureSelectionNumberofProbes = NA,
                                           RNAseqNormalizationMethod = "edgeR",
                                           ImageName = "35_QCRNAseq_T2",
                                           PNGorPDF = "png",
                                           ProduceImages = "Yes")

ClinicalFile_T1Cluster <- data.frame(ClinicalFile_T1, "CLUSTER" = InfiniumClustLabels$Cluster)
RNAseqQC18June2020T3Samples132 <- QCRNAseq(RNAseqCountMatrix = RNASeqCountMatrixMatched, 
                                           QCMatrix = RNAseqQCFile, 
                                           BetaMatrix = BetaMatrix_T1,  
                                           MvalueMatrix = MvalueMatrix_T1, 
                                           TumorPurity = TumorPurity,
                                           ClinicalFile = ClinicalFile_T1Cluster,
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
                                           ProduceImages = "Yes")

RNASeqCountMatrixNorm <- RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered
RNAseqtoGene <- match(rownames(RNASeqCountMatrixNorm), ENSMBLid$gene)

RNASeqCountMatrixNorm <- data.frame("GeneName" = ENSMBLid$name[RNAseqtoGene], RNASeqCountMatrixNorm)
dim(RNASeqCountMatrixNorm) # 29829    82

FOXP1exp <- data.frame(Cluster = factor(InfiniumClustLabels$Cluster[match(colnames(RNASeqCountMatrixNorm), InfiniumClustLabels$ID)])[-1], 
                       Expression = as.numeric(t(RNASeqCountMatrixNorm[which(RNASeqCountMatrixNorm$GeneName == "FOXP1"), ])[-1, ]),
                       Mutation = factor(RNAseqMutationCallsFiltered01ClusterLabs$FOXP1[match(colnames(RNASeqCountMatrixNorm)[-1], 
                                                                                              rownames(RNAseqMutationCallsFiltered01ClusterLabs))]))
rownames(FOXP1exp) <- names(t(RNASeqCountMatrixNorm[which(RNASeqCountMatrixNorm$GeneName == "FOXP1"), ])[-1, ])
FOXP1exp <- FOXP1exp[- which(is.na(FOXP1exp$Mutation) == TRUE), ]


ComparisonOptionsFOXP1 <- list(names(table(FOXP1exp$Mutation)))
FOXP1Plot <- ggpubr::ggboxplot(FOXP1exp, x = "Mutation", y = "Expression", fill = "Cluster",
                                add = "boxplot", ylab = "Normalized BCL2 RNAseq expression", 
                                font.label = list(size = 20, color = "black"), 
                                palette = c('#4363d8', '#f58231')) +
                                ggtitle("Cluster vs. BCL2 expression") +
                                stat_compare_means(comparisons = ComparisonOptionsFOXP1) + 
                                # Add pairwise comparisons p-value
                                stat_compare_means(aes(group = Cluster))

ComparisonOptionsFOXP1 <- list(names(table(FOXP1exp$Cluster)))
FOXP1Plot <- ggpubr::ggboxplot(FOXP1exp, x = "Cluster", y = "Expression", fill = "Cluster",
                               add = "boxplot", ylab = "Normalized FOXP1 RNAseq expression", 
                               font.label = list(size = 20, color = "black"), 
                               palette = c('#4363d8', '#f58231')) +
                               ggtitle("Cluster vs. FOXP1 expression") +
                               stat_compare_means(comparisons = ComparisonOptionsFOXP1) + 
                               # Add pairwise comparisons p-value
                               stat_compare_means(paired = FALSE)
                                


# 10 June 2020
DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = as.matrix(QCanalysisRNAseq$RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered), 
                                                                   ContrastColumnName = "CLUSTER", 
                                                                   ClinicalFile = QCanalysisRNAseq$ClinicalFileMatchedSampleFiltered, 
                                                                   ClusterLabels = InfiniumClustLabels$Cluster[MatchRNAseqwithBeta], 
                                                                   ProduceImages = "Yes", 
                                                                   PNGorPDF = "png") 

# Contrast: one-two has 116 ensembl IDs (+ve)logFC and significant. 
# Contrast: one-two has 50 ensembl IDs (-ve)logFC and significant. 

DifferentialExpressionRNAseqOutput$AllDEgenes
ENSMBLid <- read.csv(file = "hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv")

# # looking at all 166 DE gene IDs
RNAseqtoGene <- match(substr(DifferentialExpressionRNAseqOutput$AllDEgenes, 1, 15), ENSMBLid$gene)
unique(ENSMBLid$name[RNAseqtoGene[! is.na(RNAseqtoGene)]]) 
length(unique(ENSMBLid$name[RNAseqtoGene[! is.na(RNAseqtoGene)]])) # 123
table(ENSMBLid$type[RNAseqtoGene[! is.na(RNAseqtoGene)]])

# looking at 116 gene IDs
RNAseqtoGeneUP <- match(DifferentialExpressionRNAseqOutput$TopTagsEachContrastUP[[1]], ENSMBLid$gene)
length(unique(ENSMBLid$name[RNAseqtoGeneUP[! is.na(RNAseqtoGeneUP)]])) # 90
table(ENSMBLid$type[RNAseqtoGeneUP[! is.na(RNAseqtoGeneUP)]])
write.csv(data.frame("gene" = ENSMBLid$name[RNAseqtoGeneUP[! is.na(RNAseqtoGeneUP)]],
                     "type" = ENSMBLid$type[RNAseqtoGeneUP[! is.na(RNAseqtoGeneUP)]]), 
          "T1_GeneSignificance_PosFC_June2020.csv")
gprofilerRNAseqtoGeneUP <- gprofiler2::gost(query = list(ENSMBLid$name[RNAseqtoGeneUP[! is.na(RNAseqtoGeneUP)]]), 
                                                    user_threshold = 0.01,
                                                    as_short_link = TRUE)


# looking at 50 gene IDs 
RNAseqtoGeneDown <- match(DifferentialExpressionRNAseqOutput$TopTagsEachContrastDOWN[[1]], ENSMBLid$gene)
length(unique(ENSMBLid$name[RNAseqtoGeneDown[! is.na(RNAseqtoGeneDown)]])) # 33
table(ENSMBLid$type[RNAseqtoGeneDown[! is.na(RNAseqtoGeneDown)]])
write.csv(data.frame("gene" = ENSMBLid$name[RNAseqtoGeneDown[! is.na(RNAseqtoGeneDown)]],
                     "type" = ENSMBLid$type[RNAseqtoGeneDown[! is.na(RNAseqtoGeneDown)]]), 
          "T1_GeneSignificance_NegFC_June2020.csv")
gprofilerRNAseqtoGeneDown <- gprofiler2::gost(query = list(ENSMBLid$name[RNAseqtoGeneDown[! is.na(RNAseqtoGeneDown)]]), 
                                            user_threshold = 0.01,
                                            as_short_link = TRUE)











### 17 June 2020
PieChartCreation <- PieCharts(BetaMatrix = Output_SDeviation_5000_All$BetaMatrix_SD_Filtered, 
                              AnnotationFile = AnnotationFile, 
                              ProduceImages = "Yes", 
                              PNGorPDF = "png", 
                              ImageName = "AllProbes")

Output_SDeviation_5000_All <- SDeviation(BetaMatrix = BetaMatrix_T1, 
                                         MvalueMatrix = MvalueMatrix_T1, 
                                         NumberofProbes = 5000)


BoxPlotsClusters <- BoxPlotsMethylation(BetaMatrix = Output_SDeviation_5000_All$BetaMatrix_SD_Filtered, 
                                        MvalueMatrix = Output_SDeviation_5000_All$MvalueMatrix_SD_Filtered, 
                                        ClinicalFile = ClinicalFile_T1, 
                                        CategoryToVisualize = "EPIC_QC", 
                                        ClusterLabels = InfiniumClustLabels$Cluster, 
                                        TumorPurity = TumorPurity,
                                        PlotClustersWithinCategories = "Yes", 
                                        ProduceImages = "Yes",
                                        PNGorPDF = "png")


BoxPlotsClusters <- BoxPlotsMethylation(BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix), 4, 5) == "FL")], 
                                        MvalueMatrix = MvalueMatrix_T1[, which(substr(colnames(MvalueMatrix_T1), 4, 5) == "FL")], 
                                        ClinicalFile = ClinicalFile_T1[which(substr(ClinicalFile_T1$SAMPLE_ID, 4, 5) == "FL"), ],
                                        CategoryToVisualize = "CLUSTER", 
                                        ClusterLabels = InfiniumClustLabels$Cluster[which(substr(colnames(BetaMatrix), 4, 5) == "FL")],
                                        TumorPurity = TumorPurity,
                                        PlotClustersWithinCategories = NA, 
                                        ProduceImages = "Yes",
                                        PNGorPDF = "png")







# Pick top 5000 and all probes
Output_SDeviation_5000_All <- SDeviation(BetaMatrix = BetaMatrix_T1, 
                                         MvalueMatrix = MvalueMatrix_T1, 
                                         NumberofProbes = 5000)

Output_InfiniumClustering_5000_All <- Clustering(TypeofClustering = "InfiniumClust",
                                                 BetaMatrix = BetaMatrix_T1, 
                                                 MvalueMatrix = MvalueMatrix_T1, 
                                                 AnnotationFile = AnnotationFile,
                                                 ListofProbes = rownames(Output_SDeviation_5000_All$BetaMatrix_SD_Filtered), 
                                                 ClinicalFile = ClinicalFile_T1, 
                                                 TumorPurity = TumorPurity,
                                                 FigureGenerate = "No", 
                                                 PNGorPDF = "png")


MethylationDensityPlotCluster <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                       ClusterLabels = NA, 
                                                       PlotWithinAnnotationCategories = "Yes", 
                                                       ClinicalCategoryToVisualize = "TYPE", 
                                                       BetaMatrix = BetaMatrix_T1, 
                                                       AnnotationFile = AnnotationFile, 
                                                       ClinicalFile = ClinicalFile_T1, 
                                                       SampleSheet = NA, 
                                                       FigureGenerate = "Yes", 
                                                       ImageName = "Relation_to_Island", 
                                                       PNGorPDF = "png")



ClinicalFile_T1Cluster <- data.frame(ClinicalFile_T1, "CLUSTER" = InfiniumClustLabels$Cluster)
HeatPlotAll <- HeatPlot(BetaMatrix = BetaMatrix_T1, 
                       ClinicalFile = ClinicalFile_T1, 
                       CategoryToVisualize = "TYPE", 
                       PNGorPDF = "png")


Output_Visuals_Relation_to_Island_CLUSTER <- ProportionVisualization(CategoryToVisualize = "Relation_to_Island", 
                                                                     ClusterLabels = NA, 
                                                                     BetaMatrix = BetaMatrix_T1, 
                                                                     AnnotationFile = AnnotationFile, 
                                                                     ClinicalFile = ClinicalFile_T1Cluster, 
                                                                     ClinicalCategory = "CLUSTER", 
                                                                     PlotWithinCategories = "Yes", 
                                                                     FigureGenerate = "Yes", 
                                                                     PNGorPDF = "png", 
                                                                     ImageName = "Relation_to_Island_June2020")

Output_Visuals_Relation_to_Island_TYPE <- ProportionVisualization(CategoryToVisualize = NA, 
                                                                     ClusterLabels = InfiniumClustLabels$Cluster, 
                                                                     BetaMatrix = BetaMatrix_T1, 
                                                                     AnnotationFile = AnnotationFile, 
                                                                     ClinicalFile = ClinicalFile_T1Cluster, 
                                                                     ClinicalCategory = "CLUSTER", 
                                                                     PlotWithinCategories = "No", 
                                                                     FigureGenerate = "Yes", 
                                                                     PNGorPDF = "png", 
                                                                     ImageName = "CLUSTER_June2020")





# Reading new RNAseq QC file from Karin - 18 June 2020 and adding rRNA contamination data 
# RNAseqQC <- data.table::fread("~/Users/anjalisilva/Desktop/UHN/FLOMICS/RNAseq/ExtendedStudy2020/FL_TGL_STAR_logQC_2020-06-18_summary_KI.csv")
# dim(RNAseqQC) # 136  59
# saveRDS(RNAseqQC, file = "FL_TGL_STAR_logQC_2020-06-18_summary_KI.RDS")
# RNAseqQC <- readRDS("FL_TGL_STAR_logQC_2020-06-18_summary_KI.RDS")
# matchColNamesRNAseqBeta <- match(RNAseqQC$SAMPLE_ID, colnames(BetaMatrix_T1))
# length(matchColNamesRNAseqBeta[! is.na(matchColNamesRNAseqBeta)]) # 132
# matchColNamesRNAseqBeta <- matchColNamesRNAseqBeta[! is.na(matchColNamesRNAseqBeta)] # remove NA

# QCcombinedMatched <- readRDS(file = paste0("QCcombined_Integrative.rds"))
# dim(QCcombinedMatched) # 132 12
# matchColNamesBetarRNAc <- match(colnames(BetaMatrix_T1)[matchColNamesRNAseqBeta], QCcombinedMatched$SAMPLE_ID)
# 
# matchColNamesRNAseq <- match(colnames(BetaMatrix_T1)[matchColNamesRNAseqBeta], RNAseqQC$SAMPLE_ID)
# RNAseqQC132patients <- data.frame(RNAseqQC[matchColNamesRNAseq, c(2:59)], 
#                                  rRNAcontam = QCcombinedMatched$rrna_contam[matchColNamesBetarRNAc], 
#                                   Cluster = InfiniumClustLabels$Cluster[matchColNamesRNAseqBeta])
# write.csv(RNAseqQC132patients, "FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")
# saveRDS(RNAseqQC132patients, "FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.RDS")


# RNAseqCountMatrix <- data.table::fread("/Users/anjalisilva/Desktop/UHN/FLOMICS/RNAseq/ExtendedStudy2020/2020-06-18STAR_quantmode_counts_matrix_FL_136_patients.txt")
# RNAseqCountMatrixToSave <- RNAseqCountMatrix[, -1]
# rownames(RNAseqCountMatrixToSave) <- RNAseqCountMatrix[, 1]$gene
# saveRDS(RNAseqCountMatrixToSave, file = "2020-06-18STAR_quantmode_counts_matrix_FL_136_patients.RDS")

# QC_data_TGL13_11Feb2020 <- data.table::fread("/Users/anjalisilva/Desktop/UHN/FLOMICS/RNAseq/ExtendedStudy2020/QC_data_TGL13_11Feb2020.csv")
# dim(QC_data_TGL13_11Feb2020) # 205  25
# saveRDS(QC_data_TGL13_11Feb2020, file = "QC_data_TGL13_11Feb2020.RDS")

# InfiniumClustLabels2 <- read.csv("InfiniumClustLabels2.csv")
# saveRDS(InfiniumClustLabels2[, c(2:3)], file = "InfiniumClustLabels2.RDS")








# Working with new RNAseq data from Karin - 18 June 2020

# RNAseqCountMatrix18June2020 <- data.frame(readRDS(file = paste0(RNAseqDirPath, "/2020-06-18STAR_quantmode_counts_matrix_FL_136_patients.RDS")))
# dim(RNAseqCountMatrix18June2020) # 57820   136
# range(RNAseqCountMatrix18June2020) # 0 9693321
# correct column names in count matrix 
 colnames(RNAseqCountMatrix18June2020) <- 
   c(colnames(RNAseqCountMatrix18June2020)[1:10], 
   paste0(substr(colnames(RNAseqCountMatrix18June2020[, which(substr(colnames(RNAseqCountMatrix18June2020), 4, 5) == "FL")]), 1, 9), "_T1"),
   colnames(RNAseqCountMatrix18June2020)[136])
# match names with 132 samples in RNAseqQC (i.e., in beta matrix )
# matchNames <- match(RNAseqQCFile$SAMPLE_ID, colnames(RNAseqCountMatrix18June2020))
# length(matchNames[! is.na(matchNames)]) #132
# RNAseqCountMatrix18June2020Filtered <- RNAseqCountMatrix18June2020[, matchNames[! is.na(matchNames)]]
# dim(RNAseqCountMatrix18June2020Filtered) # 57820   132
# rownames(RNAseqCountMatrix18June2020Filtered) <- RNAseqCountMatrix[, 1]$gene
# saveRDS(RNAseqCountMatrix18June2020Filtered, file = "2020-06-18STAR_quantmode_counts_matrix_FL_132_patients.RDS")

RNASeqCountMatrixMatched <- data.frame(readRDS(file = paste0(RNAseqDirPath, "2020-06-18STAR_quantmode_counts_matrix_FL_132_patients.RDS")))
dim(RNASeqCountMatrixMatched) # 57820   132
dim(RNAseqQCFile) # 132  60

ClinicalFile_T1Cluster <- data.frame(ClinicalFile_T1, "CLUSTER" = InfiniumClustLabels$Cluster)
ClinicalFile_T1ClusterRNAseq <- ClinicalFile_T1Cluster[match(colnames(RNASeqCountMatrixMatched),
                                      ClinicalFile_T1Cluster$SAMPLE_ID), ]
dim(ClinicalFile_T1ClusterRNAseq) # 132  26

exprsRNAseq <- edgeR::DGEList(counts = RNASeqCountMatrixMatched, 
                              group = factor(ClinicalFile_T1ClusterRNAseq$CLUSTER))
exprsRNAseq$counts <- exprsRNAseq$counts
RNASeqCountMatrixMatchedSampleFilteredFeatureFiltered <- exprsRNAseq$counts
# perform the TMM normalization and display the normalization factors
Data2RNAseqFilterededgeRNormFactors <- edgeR::calcNormFactors(exprsRNAseq)
# plotMDS(Data2RNAseqFilterededgeRNorm)
# Generate matrix
Data2RNAseqFilterededgeRNorm <- edgeR::cpm(Data2RNAseqFilterededgeRNormFactors, 
                                           normalized.lib.sizes = TRUE, log = TRUE)
# dim(Data2RNAseqFilterededgeRNorm) # 34749    58



RNAseqQC18June2020T1Samples81 <- QCRNAseq(RNAseqCountMatrix = RNASeqCountMatrixMatched, 
                                          QCMatrix = RNAseqQCFile, 
                                          BetaMatrix = BetaMatrix_T1,  
                                          MvalueMatrix = MvalueMatrix_T1, 
                                          TumorPurity = TumorPurity,
                                          ClinicalFile = ClinicalFile_T1,
                                          SurvivalFile = SurvivalFile,
                                          RNAseqSampleCufoffUniqMapReadCount = 10000000, 
                                          RNAseqSampleCufoffUniqMapReadPercentage = 70,
                                          RNAseqSampleCufoffReadsMappedMultipleLoci = 20,
                                          RNAseqSampleCutoffRRNAcontam = 40, 
                                          RNAseqFeatureSelectionMethod = "edgeR",
                                          RNAseqFeatureSelectionCutOff = NA,
                                          RNAseqFeatureSelectionNumberofProbes = NA,
                                          RNAseqNormalizationMethod = "edgeR",
                                           ImageName = "35_QCRNAseq_T1",
                                          PNGorPDF = "png",
                                          ProduceImages = "Yes")

RNAseqQC18June2020T2Samples104 <- QCRNAseq(RNAseqCountMatrix = RNASeqCountMatrixMatched, 
                                          QCMatrix = RNAseqQCFile, 
                                          BetaMatrix = BetaMatrix_T1,  
                                          MvalueMatrix = MvalueMatrix_T1, 
                                          TumorPurity = TumorPurity,
                                          ClinicalFile = ClinicalFile_T1,
                                          SurvivalFile = SurvivalFile,
                                          RNAseqSampleCufoffUniqMapReadCount = 10000000, 
                                          RNAseqSampleCufoffUniqMapReadPercentage = 50,
                                          RNAseqSampleCufoffReadsMappedMultipleLoci = 20,
                                          RNAseqSampleCutoffRRNAcontam = 40, 
                                          RNAseqFeatureSelectionMethod = "edgeR",
                                          RNAseqFeatureSelectionCutOff = NA,
                                          RNAseqFeatureSelectionNumberofProbes = NA,
                                          RNAseqNormalizationMethod = "edgeR",
                                          ImageName = "35_QCRNAseq_T2",
                                          PNGorPDF = "png",
                                          ProduceImages = "Yes")

ClinicalFile_T1Cluster <- data.frame(ClinicalFile_T1, "CLUSTER" = InfiniumClustLabels$Cluster)
RNAseqQC18June2020T3Samples132 <- QCRNAseq(RNAseqCountMatrix = RNASeqCountMatrixMatched, 
                                           QCMatrix = RNAseqQCFile, 
                                           BetaMatrix = BetaMatrix_T1,  
                                           MvalueMatrix = MvalueMatrix_T1, 
                                           TumorPurity = TumorPurity,
                                           ClinicalFile = ClinicalFile_T1Cluster,
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
                                           ProduceImages = "Yes")
