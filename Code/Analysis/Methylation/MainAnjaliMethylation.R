# Updated on 25 June 2020
# Author: Anjali Silva

# Source files
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


# Set pathways:
MethylationDirPath <- "~/Desktop/UHN/FLOMICS/Methylation/Pipeline"
RNAseqDirPath <- "~/Desktop/UHN/FLOMICS/RNAseq/ExtendedStudy2020/"
TargetedDNAseqDirPath <- "~/Desktop/UHN/FLOMICS/TargetedSequencing"
CodeDirPath <- "/Volumes/GoogleDrive/My Drive/UHN/FLOMICS/FLOMICS-Anjali/Methylation/Pipeline"


# Uploading needed datasets
# Annotation file sent by Alberto on 9 Nov 2018
# library(readr)
# AnnotationFile <- data.table::fread(file = "Ann_IlluminaHumanMethylationEPICanno.ilm10b2.hg19.csv")
AnnotationFile <-  readRDS(file = paste0(MethylationDirPath, "/AnnotationFile_EPIC.rds"))
dim(AnnotationFile) # 866836     47

# SurvivalFile <- data.table::fread(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/Data_Survival/clinical_data_rcd17Jan2020.txt")
# cat("\n Dimension of SurvivalFile is", dim(SurvivalFile)) # 496 23
SurvivalFile <- readRDS(file = paste0(MethylationDirPath, "/SurvivalFile_rcd17Jan2020.rds"))
dim(SurvivalFile) # 496  23


# From PipelineStep1
sheet <- readRDS(file = paste0(MethylationDirPath, "/1_MethylArraySheet_updSamples_Ordered_T1.rds"))
dim(sheet) #  170   8

# Reading clinical file, should be in current directory
# Altered location of clinical file on 10 Oct 2019 
# ClinicalFile <- data.table::fread(file = "/Users/anjalisilva/Desktop/UHN/FLOMICS/Data_Clinical/sample_annotations_rcd16Oct2019.txt")
# The following file was generated from PipelineStep1
ClinicalFile_T1 <- readRDS(file = paste0(MethylationDirPath, "/1_ClinicalFile_rcd16Oct2019_updSamples_Ordered_T1.rds"))
dim(ClinicalFile_T1) # 170  25

BetaMatrix_T1 <- readRDS(file = paste0(MethylationDirPath, "/2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.rds"))
dim(BetaMatrix_T1) # 595564    170
range(BetaMatrix_T1) # 3.957731e-05 9.999547e-01

MvalueMatrix_T1 <- readRDS(file = paste0(MethylationDirPath, "/2_MvalueMatrix_updSamples_Ordered_T1_FilteredProbes.rds"))
dim(MvalueMatrix_T1) #   595564    170
range(MvalueMatrix_T1) # -7.93023  7.72350  

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


InfiniumClustLabels <- readRDS(file = paste0(MethylationDirPath, "/InfiniumClustLabels2.RDS"))
dim(InfiniumClustLabels) # 170   2
range(InfiniumClustLabels$Cluster) # 1 2

EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020 <- readRDS(file = paste0(MethylationDirPath, "/EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020.rds"))
dim(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020) # 170 6
colnames(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020)
# "CD8T"  "CD4T"  "NK"    "Bcell" "Mono"  "Neu"
range(EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020)
# -1.387779e-17  6.893340e-01



ENSMBLid <- read.csv(file = paste0(RNAseqDirPath, "/hg19_genes_using_hg38_identifiers_no_dups_Karin_04March2020.csv"))
head(ENSMBLid$gene) #  ENSG00000000003 ENSG00000000005


RNAseqCountMatrixMatched <- data.frame(readRDS(file = paste0(RNAseqDirPath, "2020-06-18STAR_quantmode_counts_matrix_FL_132_patients.RDS")))
dim(RNAseqCountMatrixMatched) # 57820   132
range(RNAseqCountMatrixMatched) # 0 9693321

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
