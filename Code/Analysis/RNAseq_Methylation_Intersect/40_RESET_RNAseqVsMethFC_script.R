# 13 August 2020
# This is a script, not a funciton
# Using files from reset tools on http://ciriellolab.org/reset/reset.html, carry out analysis RNAseq vs Methylation 
# data falling in TSS (transcription start site).
# Author: Anjali Silva



# Resource to detect Epigenetically Silenced and Enhanced Targets
DirectoryResetScripts <- "/OtherScripts/reset/R"
source(paste0(getwd(), paste0(DirectoryResetScripts,"/eventScore.R")))
source(paste0(getwd(), paste0(DirectoryResetScripts,"/methStatus.R")))
source(paste0(getwd(), paste0(DirectoryResetScripts,"/reset.R")))
source(paste0(getwd(), paste0(DirectoryResetScripts,"/shared_functions.R")))
source(paste0(getwd(), paste0(DirectoryResetScripts,"/methNorSel.R")))

load(paste0(getwd(), "/OtherScripts/reset/R/promoter-probes-list.rdata"))
View(promoter.probes.list)
head(promoter.probes.list)
names(promoter.probes.list)
# [1] "ProbID"             "index"              "Gene_Name"          "Coordinate_Pos"    
# [5] "Gene_Fantom"        "Chromosome"         "Promoter_Pos_start" "Promoter_Pos_End"  
# [9] "Promoter_No"        "Entrezgene_ID"      "Hgnc_ID"            "Uniprot_ID"      


# select TSS probes
AnnotationFileEdited <- EPICAnnotationFile
matchBetaProbes <- match(rownames(BetaMatrix_T1), AnnotationFileEdited$V1)
AnnotationFileEditedBetaMatrix <- AnnotationFileEdited[matchBetaProbes, ]
dim(AnnotationFileEditedBetaMatrix) # 595564     47
AnnotationFileEditedBetaMatrix$UCSC_RefGene_Group <- sub("\\;.*", "", AnnotationFileEditedBetaMatrix$UCSC_RefGene_Group)
promoterProbes <- AnnotationFileEditedBetaMatrix %>% filter(substr(UCSC_RefGene_Group, 1, 3) == "TSS") %>% pull("V1")
matchTSSProbes <- match(promoterProbes, rownames(BetaMatrix_T1))
length(matchTSSProbes) # 114897



ClinicalFile_T1Cluster <- data.frame(ClinicalFile_T1, "CLUSTER" = InfiniumClustLabels$Cluster)
RNAseqQC18June2020T3Samples132 <- QCRNAseq(RNAseqCountMatrix = RNAseqCountMatrixMatched, 
                                           QCMatrix = RNAseqQCFile, 
                                           RNAseqAnnotationFile = ENSMBLid, 
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
                                           ProduceImages = "No")

RNAseqT3ProteinCode <- RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations[which(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations$type == "protein_coding"), ]
dim(RNAseqT3ProteinCode) # 17190   135
rownames(RNAseqT3ProteinCode) <- RNAseqT3ProteinCode$name
RNAseqT3ProteinCode <- RNAseqT3ProteinCode[, c(4:134)]
dim(RNAseqT3ProteinCode) # 17190   131
# RNAseq counts for T3 (131 patients )


# Run first function from RESET - 
methNorSelOutput <- methNorSel(normal.mtx = BetaMatrix_T1[matchTSSProbes, c(166:170)],
                               probe.list = promoter.probes.list)
# this function first extract the data regarding the selected probes, and later based on the beta-values prepare two seperate normal-sample data-sets:
## 1. hypermethylated probes in normal condition: Enh.probes
## 2. hypomethylated probes in normal condition: Sil.probes
names(methNorSelOutput) # "normal.sil.probes" "normal.enh.probes"
typeof(methNorSelOutput)
dim(methNorSelOutput$normal.sil.probes) # 45229     5
class(methNorSelOutput$normal.sil.probes)
dim(methNorSelOutput$normal.enh.probes) # 2798    5


silProbesFL <- gsub(".*@", "", rownames(methNorSelOutput$normal.sil.probes))
length(silProbesFL) # 45229
write.csv(silProbesFL, file = "RESET_HypomethylatedProbesNormalCondition_Enhancers.csv")
enhProbesFL <- gsub(".*@", "", rownames(methNorSelOutput$normal.enh.probes))
length(enhProbesFL) # 2798
write.csv(enhProbesFL, file = "RESET_HypermethylatedProbesNormalCondition_Silencers.csv")


# Match samples between RNAseq and methylation - 131 tumor samples
matchRNAseqBeta <- match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized),
                         colnames(BetaMatrix_T1))

# Run reset function on all tumor samples - enhancers
resetResultsENH <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                         meth.tumor = BetaMatrix_T1[, matchRNAseqBeta[1:131]],
                         transcriptome = RNAseqT3ProteinCode,
                         methylation.event = c('enh'),
                         FDR.permutation.no = 100, 
                         seed = 100)

names(resetResultsENH)
dim(resetResultsENH$normal.meth)  # 2798    5
dim(resetResultsENH$meth.tumor.all) # 2798  131
dim(resetResultsENH$meth.tumor.status.all) # 2798  131
dim(resetResultsENH$beta.dist.normal.specific) #  2798    4
dim(resetResultsENH$beta.dist.normal.universal) # NULL
dim(resetResultsENH$Score.report) #  2798    5
dim(resetResultsENH$transcriptome) # 2798  131
dim(resetResultsENH$transcriptome.norm.dis) # 2798  131
dim(resetResultsENH$meth.tumor.matched) # 2798  131
dim(resetResultsENH$meth.tumor.status.matched) #  2798  131
dim(resetResultsENH$FDR.res) #  585   2
dim(resetResultsENH$permutation.res) #  2798  102
dim(resetResultsENH$score.cutoff)

# Sadegh Saghafinia <Sadegh.saghafinia@epfl.ch>
# In general score above 1.5 is what you can trust. So if the FDR is significant enough, 
# then maybe consider 1.5 as a good cut.off.


resetScoreENH <- data.frame(resetResultsENH$Score.report[which(! is.na(resetResultsENH$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENH)){
  resetScoreENH[i, 6] <- resetResultsENH$FDR.res$FDR[which(resetScoreENH$Score[i] == resetResultsENH$FDR.res$observed.scores)]
}
head(resetScoreENH)
dim(resetScoreENH)
range(resetScoreENH$Score) # 0.02586767 1.29846829
View(resetScoreENH)
write.csv(resetScoreENH, file = "RESET_HypomethylatedProbes_Enhancers.csv")


# Run reset function on all 131 tumor samples - silencers
resetResultsSIL <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                         meth.tumor = BetaMatrix_T1[, matchRNAseqBeta[1:131]],
                         transcriptome = RNAseqT3ProteinCode,
                         methylation.event = c('sil'),
                         FDR.permutation.no = 100, 
                         seed = 100)


names(resetResultsSIL)
dim(resetResultsSIL$normal.meth)  # 2798    5
dim(resetResultsSIL$meth.tumor.all) # 2798  131
dim(resetResultsSIL$meth.tumor.status.all) # 2798  131
dim(resetResultsSIL$beta.dist.normal.specific) #  2798    4
dim(resetResultsSIL$beta.dist.normal.universal) # NULL
dim(resetResultsSIL$Score.report) #  2798    5
dim(resetResultsSIL$transcriptome) # 2798  131
dim(resetResultsSIL$transcriptome.norm.dis) # 2798  131
dim(resetResultsSIL$meth.tumor.matched) # 2798  131
dim(resetResultsSIL$meth.tumor.status.matched) #  2798  131
dim(resetResultsSIL$FDR.res) #  585   2
dim(resetResultsSIL$permutation.res) #  2798  102
dim(resetResultsSIL$score.cutoff)


resetScoreSIL <- data.frame(resetResultsSIL$Score.report[which(! is.na(resetResultsSIL$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSIL)) {
  resetScoreSIL[i, 6] <- resetResultsSIL$FDR.res$FDR[which(resetScoreSIL$Score[i] == resetResultsSIL$FDR.res$observed.scores)]
}
head(resetScoreSIL)
range(resetScoreSIL$Score) # 0.002308231 1.269213329
View(resetScoreSIL)





# Run reset function on all tumor samples - enhancers - FL only
resetResultsENHFLonly <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                               meth.tumor = BetaMatrix_T1[, matchRNAseqBeta[11:131]],
                               transcriptome = RNAseqT3ProteinCode[, c(11:131)],
                               methylation.event = c('enh'),
                               FDR.permutation.no = 100, 
                               seed = 100)


resetScoreENHFLonly <- data.frame(resetResultsENHFLonly$Score.report[which(! is.na(resetResultsENHFLonly$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFLonly)) {
  resetScoreENHFLonly[i, 6] <- resetResultsENHFLonly$FDR.res$FDR[which(resetScoreENHFLonly$Score[i] == resetResultsENHFLonly$FDR.res$observed.scores)]
}
head(resetScoreENHFLonly)
range(resetScoreENHFLonly$Score) # 0.004389081 1.150200828
View(resetScoreENHFLonly)



# Run reset function on all tumor samples - silencers - FL only
resetResultsSILFLonly <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                               meth.tumor = BetaMatrix_T1[, matchRNAseqBeta[11:131]],
                               transcriptome = RNAseqT3ProteinCode[, c(11:131)],
                               methylation.event = c('sil'),
                               FDR.permutation.no = 100, 
                               seed = 100)


resetScoreSILFLonly <- data.frame(resetResultsSILFLonly$Score.report[which(! is.na(resetResultsSILFLonly$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFLonly)) {
  resetScoreSILFLonly[i, 6] <- resetResultsSILFLonly$FDR.res$FDR[which(resetScoreSILFLonly$Score[i] == resetResultsSILFLonly$FDR.res$observed.scores)]
}
head(resetScoreSILFLonly)
range(resetScoreSILFLonly$Score) # 0.004389081 1.150200828
View(resetScoreSILFLonly)
