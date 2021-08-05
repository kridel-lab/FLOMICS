# Updated 3 Aug 2021
# 1 October 2020
# This is a script, not a funciton
# Author: Anjali Silva

# For images, the MethylMix.R from original package, function MethylMix_PlotModel()
# was used for enhancing. 



#### Download needed packages ####

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("MethylMix")
library("MethylMix")


#### Example #####

# Optional register cluster to run in parallel
library(doParallel)
cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)

# Methylation data for ovarian cancer
# cancerSite <- "OV"
# targetDirectory <- paste0(getwd(), "/")

# Downloading methylation data
# METdirectories <- Download_DNAmethylation(cancerSite, targetDirectory, TRUE)

# Processing methylation data
# METProcessedData <- MethylMix::Preprocess_DNAmethylation(cancerSite, METdirectories)

# Saving methylation processed data
# saveRDS(METProcessedData, file = paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))

# Clustering methylation data
# res <- ClusterProbes(METProcessedData[[1]], METProcessedData[[2]])

data(METcancer)
dim(METcancer) # 14 251
range(METcancer) # 0.00000 0.95195

data(METnormal)
dim(METnormal) # 14  4
range(METnormal) # 0.055919 0.796010

data(GEcancer)
dim(GEcancer) # 14 251
range(GEcancer) # -6.2051  5.7490

identical(rownames(METcancer), rownames(GEcancer))


MethylMixResults <- MethylMix::MethylMix(METcancer, GEcancer, METnormal,
                                         listOfGenes = rownames(METcancer),
                                         NoNormalMode = TRUE)
# MethylMixResults <- readRDS("MethylMixStep2.rds")
names(MethylMixResults)
MethylMixResults$MethylationDrivers

# Saving methylation clustered data
# toSave <- list(METcancer = res[[1]], METnormal = res[[2]], ProbeMapping = res$ProbeMapping)
# saveRDS(toSave, file = paste0(targetDirectory, "MET_", cancerSite, "_Clustered.rds"))

stopCluster(cl)
# saveRDS(testing, file = "MethylMixStep2.rds")


#### Applying to our data - Read needed data ####

# MethylationDirPath <- "~/Desktop/UHN/FLOMICS/Methylation/Pipeline"
# BetaMatrix_T1 <- readRDS(file = paste0(MethylationDirPath, "/2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.rds"))


#### Select TSS probes for downstream analysis ####

AnnotationFileEdited <- EPICAnnotationFile
matchBetaProbes <- match(rownames(BetaMatrix_T1), AnnotationFileEdited$V1)
AnnotationFileEditedBetaMatrix <- AnnotationFileEdited[matchBetaProbes, ]
dim(AnnotationFileEditedBetaMatrix) # 595564     47
AnnotationFileEditedBetaMatrix$UCSC_RefGene_Group <- sub("\\;.*", "", AnnotationFileEditedBetaMatrix$UCSC_RefGene_Group)
promoterProbes <- AnnotationFileEditedBetaMatrix %>% 
  dplyr::filter(substr(UCSC_RefGene_Group, 1, 3) == "TSS") %>% 
  dplyr::pull("V1")
matchTSSProbes <- match(promoterProbes, rownames(BetaMatrix_T1))
length(matchTSSProbes) # 114897


#### select RNAseq samples (T3) - edge R normalized ####
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







#### Applying to our data - TSS; using Karin's Kallisto RNAseq data ####

#testing <- MethylMix::ClusterProbes(MET_Cancer = BetaMatrix_T1[, c(165:170)],
#              MET_Normal = BetaMatrix_T1[, c(1:164)])
#testing <- readRDS("MethylMixr.rds")
#names(testing)
#testing$MET_Cancer_Clustered
#testing$MET_Normal_Clustered
#testing$ProbeMapping


#  Methylation data probes to genes
EPICAnnotationFile <- data.frame(EPICAnnotationFile,
                                 UCSC_RefGene_NameOnly = sub("\\;.*", "", EPICAnnotationFile$UCSC_RefGene_Name))
dim(EPICAnnotationFile) # 866836     48
BetaTSSMatrixT1 <- BetaMatrix_T1[matchTSSProbes, ]
matchTSSProbesMatch <- match(rownames(BetaTSSMatrixT1), EPICAnnotationFile$V1)
rownames(BetaTSSMatrixT1) <- EPICAnnotationFile$UCSC_RefGene_NameOnly[matchTSSProbesMatch]
if(length(which(is.na(rownames(BetaTSSMatrixT1)) == TRUE)) > 0) {
  BetaTSSMatrixT1 <- BetaTSSMatrixT1[- which(is.na(rownames(BetaTSSMatrixT1)) == TRUE), ]
}
dim(BetaTSSMatrixT1) # 114897    170


# match methylation samples with RNAseq samples
matchBetaRNAseq <- match(colnames(RNAseqNormalizedCountDataTab[, c(11:131)]), colnames(BetaTSSMatrixT1))

# RNAseq data ENEMBLidS to genes
matchRNAseqENSMBLid <- match(rownames(RNAseqNormalizedCountDataTab), ENSMBLid$gene)

RNAseqNormalizedCountDataTabExtended <- data.frame(RNAseqNormalizedCountDataTab, 
                                                   chr = ENSMBLid$chr[matchRNAseqENSMBLid], 
                                                   name = ENSMBLid$name[matchRNAseqENSMBLid],
                                                   type = ENSMBLid$type[matchRNAseqENSMBLid])
# select protein coding only genes
RNAseqProtCodT3Mat <- RNAseqNormalizedCountDataTabExtended[which((RNAseqNormalizedCountDataTabExtended$type == "protein_coding") == TRUE), ]
dim(RNAseqProtCodT3Mat) # 17405   135
# select those not in sex chromosomes
RNAseqProtCodT3Mat <- RNAseqProtCodT3Mat[! (substr(RNAseqProtCodT3Mat$chr, 4, 5) == "X" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "Y" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "M"), ]
dim(RNAseqProtCodT3Mat) # 16739   135
table(RNAseqProtCodT3Mat$chr)
rownames(RNAseqProtCodT3Mat) <- RNAseqProtCodT3Mat$name
RNAseqProtCodT3Mat <- RNAseqProtCodT3Mat[, c(1:132)]
dim(RNAseqProtCodT3Mat) # 16726   132

# check colnames
identical(colnames(RNAseqProtCodT3Mat[, c(11:131)]), colnames(BetaTSSMatrixT1[, matchBetaRNAseq]))


set.seed(1234)
MethylMixResultsTSS <- MethylMix::MethylMix(METcancer = BetaTSSMatrixT1[, matchBetaRNAseq], 
                                            GEcancer = as.matrix(RNAseqProtCodT3Mat[, c(11:131)]), 
                                            METnormal = BetaTSSMatrixT1[, c(166:170)],
                                            NoNormalMode = TRUE)
# saveRDS(MethylMixResultsTSS, file = "41_MethylMixResultsTSS.rds")
MethylMixResultsTSS <- readRDS(file = "41_MethylMixResultsTSS.rds")

names(MethylMixResultsTSS) # 6 names
MethylMixResultsTSS$MethylationDrivers  # 1025; genes identified as transcriptionally predictive and differentially methylated 
length(MethylMixResultsTSS$MethylationDrivers) # 1025
MethylMixResultsTSS$NrComponents # number of methylation states found for each driver gene.
MethylMixResultsTSS$MixtureStates # DM-values for each driver gene (is this  mean??)
length(MethylMixResultsTSS$MixtureStates) # 1025

unlist(MethylMixResultsTSS$MixtureStates)[order(abs(unlist(MethylMixResultsTSS$MixtureStates)), decreasing = TRUE)]
length(unlist(MethylMixResultsTSS$MixtureStates)[order(abs(unlist(MethylMixResultsTSS$MixtureStates)), decreasing = TRUE)]) # 1085
# look at reverse order
tail(unlist(MethylMixResultsTSS$MixtureStates)[order(unlist(MethylMixResultsTSS$MixtureStates), decreasing = TRUE)])



MethylMixResultsTSS$MethylationStates # DM-values for all driver genes (rows) and all samples (columns)
MethylMixResultsTSS$Classifications # Matrix with integers indicating to which mixture component
# each cancer sample was assigned to, for each gene.
MethylMixResultsTSS$Models # Beta mixture model parameters for each driver gene.
names(MethylMixResultsTSS$Models)
MethylMixResultsTSS$Models$ZCCHC24$llike

# Plot the most famous methylated gene for glioblastoma
MethylMixResultsTSS$MixtureStates$MTSS1
plots <- MethylMix::MethylMix_PlotModel(GeneName = "AFF3", 
                                        MixtureModelResults = MethylMixResultsTSS, 
                                        METcancer = BetaTSSMatrixT1[, matchBetaRNAseq],
                                        GEcancer = RNAseqNormalizedCountDataTabMatrix[, c(11:131)],
                                        METnormal = BetaTSSMatrixT1[, c(166:170)])
plots$MixtureModelPlot
plots$CorrelationPlot


# my own plotting
GeneName <- "RNF180"


GeneNameMyPlot <- data.frame(normalBeta = c(BetaTSSMatrixT1[, c(166:170)][which(rownames(BetaTSSMatrixT1[, c(166:170)]) == GeneName), ][1, ], rep(NA, 116)),
                             tumorBeta = BetaTSSMatrixT1[, matchBetaRNAseq][which(rownames(BetaTSSMatrixT1[, matchBetaRNAseq]) == GeneName), ][1, ],
                             normalRNAseq = c(RNAseqProtCodT3Mat[which(rownames(RNAseqProtCodT3Mat) == GeneName), 132], rep(NA, 120)),
                             tumorRNAseq = RNAseqProtCodT3Mat[which(rownames(RNAseqProtCodT3Mat[, c(11:131)]) == GeneName), c(11:131)])

ggplot2::ggplot(data = GeneNameMyPlot) + 
  geom_point(mapping = aes(x = normalBeta, y = normalRNAseq)) +
  labs(title = "RNF180",
       x = "normalized methylation beta", 
       y = "normalized RNAseq count",
       colour = "Event") +
  theme_bw() +
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))


GeneNameMyPlot %>%
  reshape2::melt() %>%
  ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                    ylab = "Value", 
                    font.label = list(size = 20, color = "black")) +
  ggtitle("Type vs. Value") +
  stat_compare_means() + stat_n_text()



all(which((MethylMixResultsTSS$MethylationDrivers %in% rownames(RNAseqNormalizedCountDataTabMatrix)) == TRUE)) # 1025




# selecting entries with no spurious clusters
numberLength <- 1:1025
sizeSpuriousClusterExclude <- 2
for (i in seq(along = numberLength)) {
  # cat("\n Table ", table(MethylMixResultsTSS$Classifications[i, ]))
  if(any((table(MethylMixResultsTSS$Classifications[i, ]) < sizeSpuriousClusterExclude) == TRUE)) {
    print(i)
    cat(rownames(MethylMixResultsTSS$Classifications)[i])
  }
}




# Try with edgeR normalized RNAseq count matrix - Error
identical(colnames(BetaTSSMatrixT1[, matchBetaRNAseq]), 
          colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations[, c(15:135)]))
RNAseqCountsNormEdgeR <- as.matrix(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations[, c(15:135)])
dim(RNAseqCountsNormEdgeR) # 44219   121
rownames(RNAseqCountsNormEdgeR) <- RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations$name
identical(colnames(BetaTSSMatrixT1[, matchBetaRNAseq]), 
          colnames(RNAseqCountsNormEdgeR)) # TRUE

set.seed(1234)
MethylMixResultsTSSedgeR <- MethylMix::MethylMix(METcancer = BetaTSSMatrixT1[, matchBetaRNAseq], 
                                                 GEcancer = RNAseqCountsNormEdgeR, 
                                                 METnormal = BetaTSSMatrixT1[, c(166:170)],
                                                 NoNormalMode = TRUE)

names(MethylMixResultsTSSedgeR) # 6 names
MethylMixResultsTSSedgeR$MethylationDrivers  # 1025; genes identified as transcriptionally predictive and differentially methylated 
length(MethylMixResultsTSSedgeR$MethylationDrivers) # 1025
MethylMixResultsTSSedgeR$NrComponents # number of methylation states found for each driver gene.
MethylMixResultsTSSedgeR$MixtureStates # DM-values for each driver gene (is this  mean??)
MethylMixResultsTSSedgeR$MethylationStates # DM-values for all driver genes (rows) and all samples (columns)
MethylMixResultsTSSedgeR$Classifications # Matrix with integers indicating to which mixture component
# each cancer sample was assigned to, for each gene.
MethylMixResultsTSSedgeR$Models # Beta mixture model parameters for each driver gene.
names(MethylMixResultsTSSedgeR$Models)
MethylMixResultsTSSedgeR$Models$ZCCHC24$llike
# Error in { : 
#    task 1041 failed - "variable lengths differ (found for 'METcancer[Genes[i], ]')"
# save.image("41_MethylMixl_RNAseqVsMethFC_script.RData")


#### Compare with RESET and MethylMix TSS data ####
# MethylMix has 1025 genes

# compare with RESET tumor hypomethylated cases (normal enhanced - hypermethylated)
# RESET had 540 genes

MethylMixResultsTSS$MethylationDrivers[which(MethylMixResultsTSS$MethylationDrivers %in% resetScoreENHFLonlyKallisto$Gene == TRUE)]
length(MethylMixResultsTSS$MethylationDrivers[which(MethylMixResultsTSS$MethylationDrivers %in% resetScoreENHFLonlyKallisto$Gene == TRUE)])
hypomethylatedTumorRESET <- MethylMixResultsTSS$MethylationDrivers[which(MethylMixResultsTSS$MethylationDrivers %in% resetScoreENHFLonlyKallisto$Gene == TRUE)]
# 16 genes found overalp
# "ENO1"    "PRC1"    "PLCL2"   "FAM8A1"  "IGSF11"  "PRKACA"  "MFGE8"   "TNIP1"   "NSUN4"   "CD81"   
# "MSANTD3" "NUP62"   "MTSS1"   "GRIP1"   "HEPHL1"  "IL17B" 

RESEThypoTable <- resetScoreENHFLonlyKallisto[match(hypomethylatedTumorRESET, resetScoreENHFLonlyKallisto$Gene), ]
RESEThypoTable[order(RESEThypoTable$Score, decreasing = TRUE), ]


MethylMixResultsTSS$MixtureStates$PLCL2    


testingGene <- "IGSF11"
testingDataFrame <- data.frame(Values = c(resetResultsENHFLonlyKallisto$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonlyKallisto$normal.meth)) == testingGene)[1], ],
                                          resetResultsENHFLonlyKallisto$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonlyKallisto$meth.tumor.all)) == testingGene)[1], ],
                                          resetResultsENHFLonlyKallisto$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonlyKallisto$transcriptome)) == testingGene)[1], ]),
                               Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
dim(testingDataFrame) # 247 2
head(testingDataFrame)

testingDataFrame %>%  # methylation only
  filter(Type != "TumorRNAseq") %>%
  ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                    ylab = "Values", 
                    font.label = list(size = 20, color = "black")) +
  ggtitle(paste0("Type vs. Value: ", testingGene)) +
  ggpubr::stat_compare_means() + EnvStats::stat_n_text()


# plotting RNAseq vs methylation value - SPO11, FCRLB
GeneToPlot <- data.frame(Beta = t(BetaTSSMatrixT1[which(rownames(BetaTSSMatrixT1) == testingGene), matchBetaRNAseq]),
                         Event = t(ENHFLmeth.tumor.status.allKallisto[which(ENHFLmeth.tumor.status.allKallisto$Gene == testingGene), c(1:121)]),
                         RNAseqCounts = t(RNAseqProtCodT3Mat[which(rownames(RNAseqProtCodT3Mat) == testingGene), c(11:131)]))

ggplot2::ggplot(GeneToPlot, aes(x = Beta.NSUN4, y = NSUN4, 
                                color = factor(cg13488501.p4.NSUN4))) +
  geom_point(size = 2)+
  labs(title = "cg13488501.p4.NSUN4",
       x = "normalized methylation beta", 
       y = "normalized RNAseq count",
       colour = "Event") +
  theme_bw() +
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))




# compare with RESET tumor hypermethylated cases (normal enhanced - hypomethylated)
# RESET had 5621 genes
length(MethylMixResultsTSS$MethylationDrivers[which(MethylMixResultsTSS$MethylationDrivers %in% resetScoreSILFLonlyKallisto$Gene == TRUE)])
head(MethylMixResultsTSS$MethylationDrivers[which(MethylMixResultsTSS$MethylationDrivers %in% resetScoreSILFLonlyKallisto$Gene == TRUE)])
tail(MethylMixResultsTSS$MethylationDrivers[which(MethylMixResultsTSS$MethylationDrivers %in% resetScoreSILFLonlyKallisto$Gene == TRUE)])
hypermethylatedTumorRESET <- MethylMixResultsTSS$MethylationDrivers[which(MethylMixResultsTSS$MethylationDrivers %in% resetScoreSILFLonlyKallisto$Gene == TRUE)]
# 312 genes found
# "RBM38"  "NES"    "GTF3C6" "SERF2"  "YWHAH"  "GCLC"  


RESEThyperTable <- resetScoreSILFLonlyKallisto[match(hypermethylatedTumorRESET, resetScoreSILFLonlyKallisto$Gene), ]
RESEThyperTable <- RESEThyperTable[order(RESEThyperTable$Score, decreasing = TRUE), ]
RESEThyperTable[which(RESEThyperTable$Gene == "FAM134B"), ]


testingGene <- "BMP7"
MethylMixResultsTSS$MixtureStates$HOXB4                                

plots <- MethylMix::MethylMix_PlotModel(GeneName = testingGene, 
                                        MixtureModelResults = MethylMixResultsTSS, 
                                        METcancer = BetaTSSMatrixT1[, matchBetaRNAseq],
                                        GEcancer = as.matrix(RNAseqProtCodT3Mat[, c(11:131)]),
                                        METnormal = BetaTSSMatrixT1[, c(166:170)])


testingDataFrame <- data.frame(Values = c(resetResultsSILFLonlyKallisto$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonlyKallisto$normal.meth)) == testingGene)[1], ],
                                          resetResultsSILFLonlyKallisto$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonlyKallisto$meth.tumor.all)) == testingGene)[1], ],
                                          resetResultsSILFLonlyKallisto$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonlyKallisto$transcriptome)) == testingGene)[1], ]),
                               Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
dim(testingDataFrame) # 247 2
head(testingDataFrame)

testingDataFrame %>%  # methylation only
  filter(Type != "TumorRNAseq") %>%
  ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                    ylab = "Values", 
                    font.label = list(size = 20, color = "black")) +
  ggtitle(paste0("Type vs. Value: ", testingGene)) +
  ggpubr::stat_compare_means() + EnvStats::stat_n_text()


# plotting RNAseq vs methylation value - SPO11, FCRLB
GeneToPlot <- data.frame(Beta = t(BetaTSSMatrixT1[which(rownames(BetaTSSMatrixT1) == testingGene), matchBetaRNAseq]),
                         Event = t(SILFLmeth.tumor.status.allKallisto[which(SILFLmeth.tumor.status.allKallisto$Gene == testingGene), c(1:121)]),
                         RNAseqCounts = t(RNAseqProtCodT3Mat[which(rownames(RNAseqProtCodT3Mat) == testingGene), c(11:131)]))

ggplot2::ggplot(GeneToPlot, aes(x =  Beta.BMP7, y = BMP7, 
                                color = factor(Event.cg14839404.p1.BMP7))) +
  geom_point(size = 2)+
  labs(title = "cg14839404.p1.BMP7",
       x = "normalized methylation beta", 
       y = "normalized RNAseq count",
       colour = "Event") +
  theme_bw() +
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))




# Saving all genes from MethylMix, RESET overlaps
length(hypomethylatedTumorRESET) <- 
  length(hypermethylatedTumorRESET) <-
  length(MethylMixResultsTSS$MethylationDrivers)

dataFrameAll <- data.frame(MethylMixGenes = MethylMixResultsTSS$MethylationDrivers,
                           hypomethylatedTumorRESET = hypomethylatedTumorRESET,
                           hypermethylatedTumorRESET = hypermethylatedTumorRESET)
# write.csv(dataFrameAll, "MethylMixRESEToverlap_October2020.csv")





# Plot the most famous methylated gene for glioblastoma

# selecting entries with no spurious clusters for hypermethylatedTumorRESET
numberLength <- 1:length(hypermethylatedTumorRESET)
sizeSpuriousClusterExclude <- 120
for (i in seq(along = numberLength)) {
  entryNumber <- which(rownames(MethylMixResultsTSS$Classifications) == hypermethylatedTumorRESET[i])
  if(any((table(MethylMixResultsTSS$Classifications[entryNumber, ]) < sizeSpuriousClusterExclude) == TRUE)) {
    print(i)
    cat(rownames(MethylMixResultsTSS$Classifications)[entryNumber])
  }
}


table(MethylMixResultsTSS$Classifications[which(rownames(MethylMixResultsTSS$Classifications) == "RGMB"), ])




MethylMixResultsTSS$MixtureStates$SLC35D2
plots <- MethylMix::MethylMix_PlotModel(GeneName = "SLC35D2", 
                                        MixtureModelResults = MethylMixResultsTSS, 
                                        METcancer = BetaTSSMatrixT1[, matchBetaRNAseq],
                                        GEcancer = RNAseqNormalizedCountDataTabMatrix[, c(11:131)],
                                        METnormal = BetaTSSMatrixT1[, c(166:170)])
plots$MixtureModelPlot
plots$CorrelationPlot










#### gProfiler - Pathway analysis MethylMix TSS data ####

# gprofiler analysis of  all MethylMix 1025 genes
MethylMix1025 <- GProfilerAnalysis(GeneIDs = list(MethylMixResultsTSS$MethylationDrivers),
                                       Organism = "hsapiens",
                                       OrderedQuery = TRUE,
                                       PvalAlphaLevel = 0.01,
                                       PositiveorNegFC = NA, # with respect to expression
                                       ConditionName = "All",
                                       ProduceImages = "Yes", 
                                       PNGorPDF = "png")
MethylMix1025$shortLink




# gprofiler analysis of genes from methylmix overlapping with hypomethylatedTumorRESET 
hypoTumRESET <- GProfilerAnalysis(GeneIDs = list(hypomethylatedTumorRESET),
                             Organism = "hsapiens",
                             OrderedQuery = TRUE,
                             PvalAlphaLevel = 0.01,
                             PositiveorNegFC = NA, # with respect to expression
                             ConditionName = "All",
                             ProduceImages = "Yes", 
                             PNGorPDF = "png")
hypoTumRESET$shortLink



# gprofiler analysis of genes from methylmix overlapping with hypermethylatedTumorRESET 
hyperTumRESET <- GProfilerAnalysis(GeneIDs = list(hypermethylatedTumorRESET),
                                  Organism = "hsapiens",
                                  OrderedQuery = TRUE,
                                  PvalAlphaLevel = 0.01,
                                  PositiveorNegFC = NA, # with respect to expression
                                  ConditionName = "All",
                                  ProduceImages = "Yes", 
                                  PNGorPDF = "png")
hyperTumRESET$shortLink




# gprofiler on GMT file for all MethylMix 1025 genes
RNAseqDirPath <- "~/Desktop/UHN/FLOMICS/RNAseq/ExtendedStudy2020/"
customid <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, "GSEA/Staudt_signatureDB_030920.gmt"))
# Your custom annotations ID is gp__v8Ll_32EB_QGE
# You can use this ID as an 'organism' name in all the related enrichment tests against this custom source.
# Just use: gost(my_genes, organism = 'gp__v8Ll_32EB_QGE')

testingGMT <- gprofiler2::gost(query = list(MethylMixResultsTSS$MethylationDrivers), 
                               organism = customid,
                               ordered_query = FALSE,
                               user_threshold = 0.01,
                               as_short_link = TRUE)

names(testingGMT)
# https://biit.cs.ut.ee/gplink/l/_mkGUD-mT5



# gprofiler on GMT file for methylmix overlapping with hypomethylatedTumorRESET 
GMThypoTumRESET <- gprofiler2::gost(query = list(hypomethylatedTumorRESET),
                                     organism = customid,
                                     ordered_query = FALSE,
                                     user_threshold = 0.01,
                                     as_short_link = TRUE)
GMThypoTumRESET


# gprofiler on GMT file for methylmix overlapping with hypermethylatedTumorRESET 
GMThyperTumRESET <- gprofiler2::gost(query = list(hypermethylatedTumorRESET),
                                    organism = customid,
                                    ordered_query = FALSE,
                                    user_threshold = 0.01,
                                    as_short_link = TRUE)
GMThyperTumRESET




#### rGREAT - Pathway analysis MethylMix TSS data ####  

BiocManager::install("rGREAT")
library(rGREAT)

BiocManager::install("circlize")
library(circlize)

matchingEntries <- match(MethylMixResultsTSS$MethylationDrivers, annGRCh37hg19Oct2020Gene$Gene)
annGRCh37hg19Oct2020Gene$Gene[matchingEntries]
length(annGRCh37hg19Oct2020Gene$Gene[matchingEntries]) # 1025

rGREATMethDrivers1025 <- data.frame(chr = paste0("chr", annGRCh37hg19Oct2020Gene$chr[matchingEntries]), 
                                    start = as.numeric(annGRCh37hg19Oct2020Gene$start[matchingEntries]),
                                    end =  as.numeric(annGRCh37hg19Oct2020Gene$end[matchingEntries]))
                                    
rGREATMethDrivers1025 <- rGREAT::submitGreatJob(rGREATMethDrivers1025)

tbrGREATMethDrivers1025 <- rGREAT::getEnrichmentTables(rGREATMethDrivers1025)
# GO molecular function
head(tbrGREATMethDrivers1025$`GO Molecular Function`)
GOmolecularFunction <- data.frame(tbrGREATMethDrivers1025$`GO Molecular Function`)
GOmolecularFunction <- GOmolecularFunction[order(GOmolecularFunction$Binom_Adjp_BH), ]
GOmolecularFunction[GOmolecularFunction$Binom_Adjp_BH < 0.01, c(1, 2, 9)] 


GOBiologicalProcess <- data.frame(tbrGREATMethDrivers1025$`GO Biological Process`)
GOBiologicalProcess <- GOBiologicalProcess[order(GOBiologicalProcess$Binom_Adjp_BH), ]
GOBiologicalProcess[GOBiologicalProcess$Binom_Adjp_BH < 0.01, c(1, 2, 9)] 



GOCellularComponent <- data.frame(tbrGREATMethDrivers1025$`GO Cellular Component`)
GOCellularComponent <- GOCellularComponent[order(GOCellularComponent$Binom_Adjp_BH), ]
GOCellularComponent[GOCellularComponent$Binom_Adjp_BH < 0.01, c(1, 2, 9)] 



availableOntologies(rGREATMethDrivers1025)
# [1] "GO Molecular Function"     "GO Biological Process"    
# [3] "GO Cellular Component"     "Mouse Phenotype"          
# [5] "Mouse Phenotype Single KO" "Human Phenotype"          
# [7] "Ensembl Genes" 

resrGREATMethDrivers1025 <- rGREAT::plotRegionGeneAssociationGraphs(rGREATMethDrivers1025)
plotRegionGeneAssociationGraphs(rGREATMethDrivers1025, type = 1)
plotRegionGeneAssociationGraphs(rGREATMethDrivers1025, 
                                ontology = "GO Molecular Function", 
                                termID = "GO:0005515")







#### Applying to our data - All probes; using Karin's Kallisto RNAseq data ####

#testing <- MethylMix::ClusterProbes(MET_Cancer = BetaMatrix_T1[, c(165:170)],
#              MET_Normal = BetaMatrix_T1[, c(1:164)])
#testing <- readRDS("MethylMixr.rds")
#names(testing)
#testing$MET_Cancer_Clustered
#testing$MET_Normal_Clustered
#testing$ProbeMapping

#  Methylation data probes to genes
EPICAnnotationFile <- data.frame(EPICAnnotationFile,
                                 UCSC_RefGene_NameOnly = sub("\\;.*", "", EPICAnnotationFile$UCSC_RefGene_Name))
dim(EPICAnnotationFile) # 866836     48
matchAllProbesMatch <- match(rownames(BetaMatrix_T1), EPICAnnotationFile$V1)
BetaMatrix_T1Names <- BetaMatrix_T1
rownames(BetaMatrix_T1Names) <- EPICAnnotationFile$UCSC_RefGene_NameOnly[matchAllProbesMatch]
if(length(which(is.na(rownames(BetaMatrix_T1Names)) == TRUE)) > 0) {
  BetaMatrix_T1Names <- BetaMatrix_T1Names[- which(is.na(rownames(BetaMatrix_T1Names)) == TRUE), ]
}
dim(BetaMatrix_T1Names) # 595564    170

if(length(which(rownames(BetaMatrix_T1Names) == "")) > 0) {
  BetaMatrix_T1Names <- BetaMatrix_T1Names[- which(rownames(BetaMatrix_T1Names) == ""), ]
}
dim(BetaMatrix_T1Names) # 428129    170


# match methylation samples with RNAseq samples
matchBetaRNAseq <- match(colnames(RNAseqNormalizedCountDataTab[, c(11:131)]), colnames(BetaTSSMatrixT1))

# RNAseq data ENEMBLidS to genes
matchRNAseqENSMBLid <- match(rownames(RNAseqNormalizedCountDataTab), ENSMBLid$gene)

RNAseqNormalizedCountDataTabExtended <- data.frame(RNAseqNormalizedCountDataTab, 
                                                   name = ENSMBLid$name[matchRNAseqENSMBLid],
                                                   type = ENSMBLid$type[matchRNAseqENSMBLid])
# select protein coding only genes
RNAseqProtCodT3Mat <- RNAseqNormalizedCountDataTabExtended[which((RNAseqNormalizedCountDataTabExtended$type == "protein_coding") == TRUE), c(1:132)]
dim(RNAseqProtCodT3Mat) # 17405   132
rownames(RNAseqProtCodT3Mat) <- RNAseqNormalizedCountDataTabExtended[which((RNAseqNormalizedCountDataTabExtended$type == "protein_coding") == TRUE), 133]


# check colnames
identical(colnames(RNAseqProtCodT3Mat[, c(11:131)]), colnames(BetaMatrix_T1Names[, matchBetaRNAseq]))
# TRUE

set.seed(1234)
MethylMixResultsAll <- MethylMix::MethylMix(METcancer = BetaMatrix_T1Names[, matchBetaRNAseq], 
                                            GEcancer = as.matrix(RNAseqProtCodT3Mat[, c(11:131)]), 
                                            METnormal = BetaMatrix_T1Names[, c(166:170)],
                                            NoNormalMode = TRUE)
# saveRDS(MethylMixResultsAll, file = "MethylMixResultsAll.rds")
MethylMixResultsAll <- readRDS(file = "41_MethylMixResultsAll.rds")

names(MethylMixResultsAll) # 6 names
MethylMixResultsAll$MethylationDrivers  # 1034; genes identified as transcriptionally predictive and differentially methylated 
length(MethylMixResultsAll$MethylationDrivers) # 1034
MethylMixResultsAll$NrComponents # number of methylation states found for each driver gene.
MethylMixResultsAll$MixtureStates # DM-values for each driver gene (is this  mean??)
MethylMixResultsAll$MethylationStates # DM-values for all driver genes (rows) and all samples (columns)
MethylMixResultsAll$Classifications # Matrix with integers indicating to which mixture component
# each cancer sample was assigned to, for each gene.
MethylMixResultsAll$Models # Beta mixture model parameters for each driver gene.
names(MethylMixResultsAll$Models)
MethylMixResultsAll$Models$ZCCHC24$llike


unlist(MethylMixResultsAll$MixtureStates)[order(abs(unlist(MethylMixResultsAll$MixtureStates)), decreasing = TRUE)]
length(unlist(MethylMixResultsAll$MixtureStates)[order(abs(unlist(MethylMixResultsAll$MixtureStates)), decreasing = TRUE)]) #  1093
head(unlist(MethylMixResultsAll$MixtureStates)[order(unlist(MethylMixResultsAll$MixtureStates), decreasing = TRUE)])
# SPRED1     ZIK1    ACSL6    EYA2   RNF180    NTSR1             ASCL4
# 0.4449486 0.4416508 0.4332334 0.4195578 0.4153790 0.4121466 

# look at reverse order
tail(unlist(MethylMixResultsAll$MixtureStates)[order(unlist(MethylMixResultsAll$MixtureStates), decreasing = TRUE)])
# CARD11      HTR3A    PRDM151      PAX8    PRR15L1    PPFIA11 
# -0.2942115 -0.3209366 -0.3270504 -0.3299306 -0.3694646 -0.7949000

MethylMixResultsAll$MixtureStates$ASCL4


# Plot the most famous methylated gene for glioblastoma
plots <- MethylMix::MethylMix_PlotModel(GeneName = "CARD11", 
                                        MixtureModelResults = MethylMixResultsAll, 
                                        METcancer = BetaMatrix_T1Names[, matchBetaRNAseq],
                                        GEcancer = RNAseqNormalizedCountDataTabMatrix[, c(11:131)],
                                        METnormal = BetaMatrix_T1Names[, c(166:170)])
plots$MixtureModelPlot
plots$CorrelationPlot


all(which((MethylMixResultsTSS$MethylationDrivers %in% rownames(RNAseqNormalizedCountDataTabMatrix)) == TRUE)) # 1025


# selecting entries with no spurious clusters
numberLength <- 1:1034
sizeSpuriousClusterExclude <- 6
for (i in seq(along = numberLength)) {
  # cat("\n Table ", table(MethylMixResultsTSS$Classifications[i, ]))
  if(any((table(MethylMixResultsAll$Classifications[i, ]) < sizeSpuriousClusterExclude) == TRUE)) {
    print(i)
    cat(rownames(MethylMixResultsAll$Classifications)[i])
  }
}

table(MethylMixResultsAll$Classifications[which(rownames(MethylMixResultsAll$Classifications) == "SPRED1"), ])


#### gProfiler - Pathway analysis MethylMix TSS data ####

# gprofiler analysis of  all MethylMix 1025 genes
MethylMix1034 <- GProfilerAnalysis(GeneIDs = list(MethylMixResultsAll$MethylationDrivers),
                                   Organism = "hsapiens",
                                   OrderedQuery = TRUE,
                                   PvalAlphaLevel = 0.01,
                                   PositiveorNegFC = NA, # with respect to expression
                                   ConditionName = "All",
                                   ProduceImages = "Yes", 
                                   PNGorPDF = "png")
MethylMix1034$shortLink
# "https://biit.cs.ut.ee/gplink/l/3dLtcb1nS9"


# gprofiler on GMT file for all MethylMix 1025 genes
RNAseqDirPath <- "~/Desktop/UHN/FLOMICS/RNAseq/ExtendedStudy2020/"
customid <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, "GSEA/Staudt_signatureDB_030920.gmt"))
# Your custom annotations ID is gp__v8Ll_32EB_QGE
# You can use this ID as an 'organism' name in all the related enrichment tests against this custom source.
# Just use: gost(my_genes, organism = 'gp__v8Ll_32EB_QGE')

testingGMTAll <- gprofiler2::gost(query = list(MethylMixResultsAll$MethylationDrivers), 
                               organism = customid,
                               ordered_query = FALSE,
                               user_threshold = 0.01,
                               as_short_link = TRUE)

testingGMTAll
# https://biit.cs.ut.ee/gplink/l/IBxQ1MdKT8


# gprofiler on GMT file for methylmix overlapping with hypomethylatedTumorRESET 
GMThypoTumRESET <- gprofiler2::gost(query = list(hypomethylatedTumorRESET),
                                    organism = customid,
                                    ordered_query = FALSE,
                                    user_threshold = 0.01,
                                    as_short_link = TRUE)
GMThypoTumRESET


# gprofiler on GMT file for methylmix overlapping with hypermethylatedTumorRESET 
GMThyperTumRESET <- gprofiler2::gost(query = list(hypermethylatedTumorRESET),
                                     organism = customid,
                                     ordered_query = FALSE,
                                     user_threshold = 0.01,
                                     as_short_link = TRUE)
GMThyperTumRESET





#### Run Enrichr - All probes ####

#  Enrichr - tool for analysing gene sets and returns any enrichment of common annotated biological features
# https://maayanlab.cloud/Enrichr/

  
dbs <- enrichR::listEnrichrDbs()
dbs$libraryName

# based on slack discussion with Robert 30 June 2020
dbsSelected <- c("ChEA_2016", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
                 "ENCODE_Histone_Modifications_2015", "ENCODE_TF_ChIP-seq_2015",
                 "Epigenomics_Roadmap_HM_ChIP-seq")
enrichrResults <- enrichR::enrichr(as.vector(MethylMixResultsAll$MethylationDrivers), 
                                   databases = dbsSelected)
enrichr1 <- data.frame(enrichrResults[1]) %>%
  dplyr::select(term = ChEA_2016.Term, 
                adj.P.val = ChEA_2016.Adjusted.P.value) %>%
  dplyr::mutate(minus.log10.P.val = -log10(adj.P.val)) %>%
  dplyr::mutate(database = "ChEA_2016") %>%
  dplyr::filter(grepl("Human", term)) %>%
  head(10)

enrichr2 <- data.frame(enrichrResults[2]) %>%
  dplyr::select(term = ENCODE_and_ChEA_Consensus_TFs_from_ChIP.X.Term, 
                adj.P.val = ENCODE_and_ChEA_Consensus_TFs_from_ChIP.X.Adjusted.P.value) %>%
  dplyr::mutate(minus.log10.P.val = -log10(adj.P.val)) %>%
  dplyr::mutate(database = paste0(dbsSelected[2])) %>%
  head(10)
enrichr3 <- data.frame(enrichrResults[3]) %>%
  dplyr::select(term = ENCODE_Histone_Modifications_2015.Term, 
                adj.P.val = ENCODE_Histone_Modifications_2015.Adjusted.P.value) %>%
  dplyr::mutate(minus.log10.P.val = -log10(adj.P.val)) %>%
  dplyr::mutate(database = paste0(dbsSelected[3])) %>%
  dplyr::filter(grepl("hg19", term)) %>%
  head(10)

enrichr4 <- data.frame(enrichrResults[4]) %>%
  dplyr::select(term = ENCODE_TF_ChIP.seq_2015.Term, 
                adj.P.val = ENCODE_TF_ChIP.seq_2015.Adjusted.P.value) %>%
  dplyr::mutate(minus.log10.P.val = -log10(adj.P.val)) %>%
  dplyr::mutate(database = paste0(dbsSelected[4])) %>%
  dplyr::filter(grepl("hg19", term)) %>%
  head(10)

enrichr5
  dplyr::select(term = Epigenomics_Roadmap_HM_ChIP.seq.Term, 
                adj.P.val = Epigenomics_Roadmap_HM_ChIP.seq.Adjusted.P.value) %>%
  dplyr::mutate(minus.log10.P.val = -log10(adj.P.val)) %>%
  dplyr::mutate(database = paste0(dbsSelected[5])) %>%
  head(10)

  
#### Rerun 24 November 2020 using edgeR normalized (non-log) data #####################
#### select RNAseq samples (T3) - edge R normalized ####

AllClustersLab <- read.csv(file = "InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv")
SNFClusterLabs <- AllClustersLab$SNFClust[- which(is.na(AllClustersLab$SNFClust) == TRUE)]
  
  
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
  
  
#### Applying to our data - TSS; using edgeR norm (non-log) RNAseq data - 101 Samples ####
  
  #testing <- MethylMix::ClusterProbes(MET_Cancer = BetaMatrix_T1[, c(165:170)],
  #              MET_Normal = BetaMatrix_T1[, c(1:164)])
  #testing <- readRDS("MethylMixr.rds")
  #names(testing)
  #testing$MET_Cancer_Clustered
  #testing$MET_Normal_Clustered
  #testing$ProbeMapping
  
  
  #  Methylation data probes to genes
  EPICAnnotationFile <- data.frame(EPICAnnotationFile,
                                   UCSC_RefGene_NameOnly = sub("\\;.*", "", EPICAnnotationFile$UCSC_RefGene_Name))
  dim(EPICAnnotationFile) # 866836     48
  BetaTSSMatrixT1 <- BetaMatrix_T1[matchTSSProbes, ]
  matchTSSProbesMatch <- match(rownames(BetaTSSMatrixT1), EPICAnnotationFile$V1)
  rownames(BetaTSSMatrixT1) <- EPICAnnotationFile$UCSC_RefGene_NameOnly[matchTSSProbesMatch]
  if(length(which(is.na(rownames(BetaTSSMatrixT1)) == TRUE)) > 0) {
    BetaTSSMatrixT1 <- BetaTSSMatrixT1[- which(is.na(rownames(BetaTSSMatrixT1)) == TRUE), ]
  }
  dim(BetaTSSMatrixT1) # 114897    170
  
  
  # match methylation samples with RNAseq samples
  matchBetaRNAseq <- match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized[, c(11:131)]), 
                           colnames(BetaTSSMatrixT1))
  
  
  # RNA seq
  dim(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized)
  # 44219   132
  range(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized)
  # 0.00 92597.52
  range(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedLOG)
  # -2.72230 16.49869
  
  # RNAseq data ENEMBLidS to genes
  matchRNAseqENSMBLid <- match(rownames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized), ENSMBLid$gene)
  RNAseqNormalizedCountDataTabExtended <- data.frame(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized, 
                                                     ENSEMBLID = rownames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized),
                                                     chr = ENSMBLid$chr[matchRNAseqENSMBLid], 
                                                     name = ENSMBLid$name[matchRNAseqENSMBLid],
                                                     type = ENSMBLid$type[matchRNAseqENSMBLid])
  # remove NA values
  RNAseqNormalizedCountDataTabExtended <- RNAseqNormalizedCountDataTabExtended[- which(is.na(RNAseqNormalizedCountDataTabExtended$type) == TRUE), ]
  dim(RNAseqNormalizedCountDataTabExtended) # 31467   136
  
  # select protein coding only genes
  RNAseqProtCodT3Mat <- RNAseqNormalizedCountDataTabExtended[which((RNAseqNormalizedCountDataTabExtended$type == "protein_coding") == TRUE), ]
  dim(RNAseqProtCodT3Mat) # 17190   136
  # select those not in sex chromosomes
  RNAseqProtCodT3Mat <- RNAseqProtCodT3Mat[! (substr(RNAseqProtCodT3Mat$chr, 4, 5) == "X" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "Y" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "M"), ]
  dim(RNAseqProtCodT3Mat) # 17190   136
  table(RNAseqProtCodT3Mat$chr)
  
  matchSNF101 <- match(AllClustersLab$ID[- which(is.na(AllClustersLab$SNFClust) == TRUE)], colnames(RNAseqProtCodT3Mat))
  RNAseqProtCodT3MatName <- RNAseqProtCodT3Mat[, matchSNF101]
  rownames(RNAseqProtCodT3MatName) <- RNAseqProtCodT3Mat$name
  range(RNAseqProtCodT3MatName) # 0.00 17766.33
  dim(RNAseqProtCodT3MatName) # 17190   101

  matchSNF101 <- match(AllClustersLab$ID[- which(is.na(AllClustersLab$SNFClust) == TRUE)], colnames(BetaTSSMatrixT1))
  
  
  # check colnames
  identical(colnames(RNAseqProtCodT3MatName), colnames(BetaTSSMatrixT1[, matchSNF101])) # TRUE
  
  
  set.seed(1234)
  MethylMixResultsTSSSNF <- MethylMix::MethylMix(METcancer = BetaTSSMatrixT1[, matchSNF101], 
                                              GEcancer = as.matrix(RNAseqProtCodT3MatName), 
                                              METnormal = BetaTSSMatrixT1[, c(166:170)],
                                              NoNormalMode = TRUE)
  # saveRDS(MethylMixResultsTSSSNF, file = "41_MethylMixResultsTSS_SNF_Nov2020.rds")
  MethylMixResultsTSSSNF <- readRDS(file = "41_MethylMixResultsTSS_SNF_Nov2020.rds")
  
  names(MethylMixResultsTSSSNF) # 6 names
  MethylMixResultsTSSSNF$MethylationDrivers  # 343; genes identified as transcriptionally predictive and differentially methylated 
  length(MethylMixResultsTSSSNF$MethylationDrivers) # 343
  MethylMixResultsTSSSNF$NrComponents # number of methylation states found for each driver gene.
  MethylMixResultsTSSSNF$MixtureStates # DM-values for each driver gene (is this  mean??)
  length(MethylMixResultsTSSSNF$MixtureStates) # 343
  
  unlist(MethylMixResultsTSSSNF$MixtureStates)[order(abs(unlist(MethylMixResultsTSSSNF$MixtureStates)), decreasing = TRUE)]
  length(unlist(MethylMixResultsTSSSNF$MixtureStates)[order(abs(unlist(MethylMixResultsTSSSNF$MixtureStates)), decreasing = TRUE)]) # 462
  # look at reverse order
  tail(unlist(MethylMixResultsTSSSNF$MixtureStates)[order(unlist(MethylMixResultsTSSSNF$MixtureStates), decreasing = TRUE)])
  
  
  
  MethylMixResultsTSSSNF$MethylationStates # DM-values for all driver genes (rows) and all samples (columns)
  MethylMixResultsTSSSNF$Classifications # Matrix with integers indicating to which mixture component
  # each cancer sample was assigned to, for each gene.
  MethylMixResultsTSSSNF$Models # Beta mixture model parameters for each driver gene.
  names(MethylMixResultsTSSSNF$Models)
  MethylMixResultsTSSSNF$Models$ZCCHC24$llike
  
  
  
  table(MethylMixResultsTSSSNF$NrComponents)
  # 1   2   3 
  # 225 117   1
  
  # Most DM value (with 1 component) if - =  tumor hypomethylated compared to normal
  # Most DM values (with 1 component) if + = tumor hypermethylated compared to normal
  # Save entries with only one component and see if tumor hyper or hypo methylated
  genesToLoop <- which(MethylMixResultsTSSSNF$NrComponents == 1)
  emptySumVectorSNF <- vector()
  for(i in 1:length(genesToLoop)) { 
    cat("\n", i)
    if(MethylMixResultsTSSSNF$MixtureStates[[genesToLoop[i]]] < 0) {
      emptySumVectorSNF[i] <- -1
    } else if(MethylMixResultsTSSSNF$MixtureStates[[genesToLoop[i]]] > 0) {
      emptySumVectorSNF[i] <- 1
    } 
  }
  
  table(emptySumVectorSNF)
  # emptySumVector
  # -1  1 
  # 167  58
  which(emptySumVectorSNF == -1)
  
  
  

  # Plot the most famous methylated gene for glioblastoma
  GeneName <- "EPHA4"
  MethylMixResultsTSSSNF$MixtureStates$GPRC5C
  
  # look at stage
  matchClinical101 <- match(colnames(BetaTSSMatrixT1[, matchSNF101]), ClinicalFile_T1$SAMPLE_ID)
  table(MethylMixResultsTSSSNF$Classifications[which(rownames(MethylMixResultsTSSSNF$Classifications) == GeneName), ],
        ClinicalFile_T1$STAGE[matchClinical101])
  chisq.test(table(MethylMixResultsTSSSNF$Classifications[which(rownames(MethylMixResultsTSSSNF$Classifications) == GeneName), ],
                   ClinicalFile_T1$STAGE[matchClinical101]))
  plots <- MethylMix::MethylMix_PlotModel(GeneName = GeneName, 
                                          MixtureModelResults = MethylMixResultsTSSSNF, 
                                          METcancer = BetaTSSMatrixT1[, matchSNF101],
                                          GEcancer = as.matrix(RNAseqProtCodT3MatName),
                                          METnormal = BetaTSSMatrixT1[, c(166:170)])
  plots$MixtureModelPlot
  plots$CorrelationPlot
  
  
  # my own plotting
  normalBeta <- as.vector(c(BetaTSSMatrixT1[, c(166:170)][which(rownames(BetaTSSMatrixT1[, c(166:170)]) == GeneName), ][1, ], rep(NA, 96)))
  tumorBeta <- as.vector(BetaTSSMatrixT1[, matchSNF101][which(rownames(BetaTSSMatrixT1[, matchSNF101]) == GeneName), ][1, ])
  # tumorRNAseq <- RNAseqProtCodT3MatNameC1[which(rownames(RNAseqProtCodT3MatNameC1) == GeneName), ]
  GeneNameMyPlot <- data.frame(
    normalBeta = normalBeta,
    tumorBeta = tumorBeta
    # tumorRNAseq = RNAseqProtCodT3MatName
  )
  
  GeneNameMyPlot %>%
    reshape2::melt() %>%
    ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                      ylab = "Value", 
                      font.label = list(size = 20, color = "black")) +
    ggtitle(paste0(GeneName)) +
    ggpubr::stat_compare_means() + 
    EnvStats::stat_n_text()
  
  
  

  ggplot2::ggplot(data = GeneNameMyPlot) + 
    geom_point(mapping = aes(x = normalBeta, y = normalRNAseq)) +
    labs(title = "RNF180",
         x = "normalized methylation beta", 
         y = "normalized RNAseq count",
         colour = "Event") +
    theme_bw() +
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold")) +
    theme(aspect.ratio = 1, 
          legend.position = "right", 
          panel.background = element_rect(colour = "black", size=1.5),  
          axis.title =  element_text(face = "bold"))
  

  
  all(which((MethylMixResultsTSSSNF$MethylationDrivers %in% rownames(RNAseqNormalizedCountDataTabMatrix)) == TRUE)) # 1025
  
  
  
  
  # selecting entries with no spurious clusters
  
  # selecting entries with no spurious clusters
  sizeSpuriousClusterExclude <- 10
  saveGeneNamesAll <- NA
  cat("\n Methylation drivers that should be kept:")
  for (i in 1:nrow(MethylMixResultsTSSSNF$Classifications)) {
    #cat("\n i = ", i, "\n")
    if(length(unique(MethylMixResultsTSSSNF$Classifications[i, ])) > 1) { # if more than one component
      if(all((table(MethylMixResultsTSSSNF$Classifications[i, ]) > sizeSpuriousClusterExclude) == TRUE)) {
        saveGeneNamesAll <- c(saveGeneNamesAll, rownames(MethylMixResultsTSSSNF$Classifications)[i])
      }
    }
  }
  
  # Arranging such genes with no suprious clusters by DM values
  matchNoSpurious <- match(saveGeneNamesAll[-1], names(MethylMixResultsTSSSNF$MixtureStates))
  saveTopResults <- list()
  for (i in 1:length(matchNoSpurious)) {
    if(length(MethylMixResultsTSSSNF$MixtureStates[[matchNoSpurious[i]]]) > 1) {
      saveTopResults[[i]] <- max(abs(MethylMixResultsTSSSNF$MixtureStates[[matchNoSpurious[i]]]))
    } else {
      saveTopResults[[i]] <- MethylMixResultsTSSSNF$MixtureStates[[matchNoSpurious[i]]]
    }
  }
  names(saveTopResults) <- names(MethylMixResultsTSSSNF$MixtureStates)[matchNoSpurious]
  # sort values from highest to lowest
  sort(unlist(saveTopResults), decreasing = TRUE)
  
  
  gprofiler101 <- GProfilerAnalysis(GeneIDs = list(MethylMixResultsTSSSNF$MethylationDrivers),
                                         Organism = "hsapiens",
                                         OrderedQuery = TRUE,
                                         PvalAlphaLevel = 0.01,
                                         PositiveorNegFC = NA, # with respect to expression
                                         ConditionName = "101SNFsamples",
                                         ProduceImages = "Yes", 
                                         PNGorPDF = "png")  
  gprofiler101$shortLink  
  # https://biit.cs.ut.ee/gplink/l/vcH1_SFNQ2
  dim(gprofiler101$TableAllValues) # 54
  
  
  gprofilersameC1C2 <- GProfilerAnalysis(GeneIDs = list(sameC1C2),
                                         Organism = "hsapiens",
                                         OrderedQuery = TRUE,
                                         PvalAlphaLevel = 0.01,
                                         PositiveorNegFC = NA, # with respect to expression
                                         ConditionName = "sameC1C2",
                                         ProduceImages = "Yes", 
                                         PNGorPDF = "png")  
  
  
  
  
#### Applying to our data - TSS; using edgeR norm (non-log) RNAseq data - C1 Samples (18Nov2020) ####
  
  
  #  Methylation data probes to genes
  EPICAnnotationFile <- data.frame(EPICAnnotationFile,
                                   UCSC_RefGene_NameOnly = sub("\\;.*", "", EPICAnnotationFile$UCSC_RefGene_Name))
  dim(EPICAnnotationFile) # 866836     48
  BetaTSSMatrixT1 <- BetaMatrix_T1[matchTSSProbes, ]
  matchTSSProbesMatch <- match(rownames(BetaTSSMatrixT1), EPICAnnotationFile$V1)
  rownames(BetaTSSMatrixT1) <- EPICAnnotationFile$UCSC_RefGene_NameOnly[matchTSSProbesMatch]
  if(length(which(is.na(rownames(BetaTSSMatrixT1)) == TRUE)) > 0) {
    BetaTSSMatrixT1 <- BetaTSSMatrixT1[- which(is.na(rownames(BetaTSSMatrixT1)) == TRUE), ]
  }
  dim(BetaTSSMatrixT1) # 114897    170
  
  
  # match methylation samples with RNAseq samples
  matchBetaRNAseq <- match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized[, c(11:131)]), 
                           colnames(BetaTSSMatrixT1))
  
  
  # RNA seq
  dim(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized)
  # 44219   132
  range(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized)
  # 0.00 92597.52
  range(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedLOG)
  # -2.72230 16.49869
  
  # RNAseq data ENEMBLidS to genes
  matchRNAseqENSMBLid <- match(rownames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized), ENSMBLid$gene)
  RNAseqNormalizedCountDataTabExtended <- data.frame(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized, 
                                                     ENSEMBLID = rownames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized),
                                                     chr = ENSMBLid$chr[matchRNAseqENSMBLid], 
                                                     name = ENSMBLid$name[matchRNAseqENSMBLid],
                                                     type = ENSMBLid$type[matchRNAseqENSMBLid])
  # remove NA values
  RNAseqNormalizedCountDataTabExtended <- RNAseqNormalizedCountDataTabExtended[- which(is.na(RNAseqNormalizedCountDataTabExtended$type) == TRUE), ]
  dim(RNAseqNormalizedCountDataTabExtended) # 31467   136
  
  # select protein coding only genes
  RNAseqProtCodT3Mat <- RNAseqNormalizedCountDataTabExtended[which((RNAseqNormalizedCountDataTabExtended$type == "protein_coding") == TRUE), ]
  dim(RNAseqProtCodT3Mat) # 17190   136
  # select those not in sex chromosomes
  RNAseqProtCodT3Mat <- RNAseqProtCodT3Mat[! (substr(RNAseqProtCodT3Mat$chr, 4, 5) == "X" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "Y" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "M"), ]
  dim(RNAseqProtCodT3Mat) # 17190   136
  table(RNAseqProtCodT3Mat$chr)
  
  matchSNF101 <- match(AllClustersLab$ID[- which(is.na(AllClustersLab$SNFClust) == TRUE)], colnames(RNAseqProtCodT3Mat))
  RNAseqProtCodT3MatName <- RNAseqProtCodT3Mat[, matchSNF101]
  rownames(RNAseqProtCodT3MatName) <- RNAseqProtCodT3Mat$name
  range(RNAseqProtCodT3MatName) # 0.00 17766.33
  dim(RNAseqProtCodT3MatName) # 17190   101
  
  
  AllClustersLab <- read.csv(file = "InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv")
  
  matchSNF101C1 <- match(AllClustersLab$ID[which(AllClustersLab$SNFClust == 1)], colnames(BetaTSSMatrixT1))
  matchSNF101C1RNaseq <- match(AllClustersLab$ID[which(AllClustersLab$SNFClust == 1)], colnames(RNAseqProtCodT3MatName))
  
  
  # check colnames
  identical(colnames(RNAseqProtCodT3MatName), colnames(BetaTSSMatrixT1[, matchSNF101])) # TRUE
  
  
  set.seed(1234)
  MethylMixResultsTSSSNFC1 <- MethylMix::MethylMix(METcancer = BetaTSSMatrixT1[, matchSNF101C1], 
                                                 GEcancer = as.matrix(RNAseqProtCodT3MatName[, matchSNF101C1RNaseq]), 
                                                 METnormal = BetaTSSMatrixT1[, c(166:170)],
                                                 NoNormalMode = TRUE)
  # saveRDS(MethylMixResultsTSSSNFC1, file = "41_MethylMixResultsTSSC1_SNF_Nov2020.rds")
  MethylMixResultsTSSSNFC1 <- readRDS(file = "41_MethylMixResultsTSSC1_SNF_Nov2020.rds")
  
  names(MethylMixResultsTSSSNFC1) # 6 names
  MethylMixResultsTSSSNFC1$MethylationDrivers  # 82; genes identified as transcriptionally predictive and differentially methylated 
  length(MethylMixResultsTSSSNFC1$MethylationDrivers) # 82
  MethylMixResultsTSSSNFC1$NrComponents # number of methylation states found for each driver gene.
  MethylMixResultsTSSSNFC1$MixtureStates # DM-values for each driver gene (is this  mean??)
  length(MethylMixResultsTSSSNFC1$MixtureStates) # 82
  
  unlist(MethylMixResultsTSSSNFC1$MixtureStates)[order(abs(unlist(MethylMixResultsTSSSNFC1$MixtureStates)), decreasing = TRUE)]
  length(unlist(MethylMixResultsTSSSNFC1$MixtureStates)[order(abs(unlist(MethylMixResultsTSSSNFC1$MixtureStates)), decreasing = TRUE)]) # 108
  # look at reverse order
  tail(unlist(MethylMixResultsTSSSNFC1$MixtureStates)[order(unlist(MethylMixResultsTSSSNFC1$MixtureStates), decreasing = TRUE)])

  
  MethylMixResultsTSSSNFC1$MethylationStates # DM-values for all driver genes (rows) and all samples (columns)
  MethylMixResultsTSSSNFC1$Classifications # Matrix with integers indicating to which mixture component
  # each cancer sample was assigned to, for each gene.
  MethylMixResultsTSSSNFC1$Models # Beta mixture model parameters for each driver gene.
  names(MethylMixResultsTSSSNFC1$Models)
  MethylMixResultsTSSSNFC1$Models$ZCCHC24$llike
  
  
  table(MethylMixResultsTSSSNFC1$NrComponents)
  # 1  2 
  # 56 26
  
  # Most DM value (with 1 component) if - =  tumor hypomethylated compared to normal
  # Most DM values (with 1 component) if + = tumor hypermethylated compared to normal
  # Save entries with only one component and see if tumor hyper or hypo methylated
  genesToLoop <- which(MethylMixResultsTSSSNFC1$NrComponents == 1)
  emptySumVectorC1 <- vector()
  for(i in 1:length(genesToLoop)) { 
    cat("\n", i)
    if(MethylMixResultsTSSSNFC1$MixtureStates[[genesToLoop[i]]] < 0) {
      emptySumVectorC1[i] <- -1
    } else if(MethylMixResultsTSSSNFC1$MixtureStates[[genesToLoop[i]]] > 0) {
      emptySumVectorC1[i] <- 1
    } 
  }
  
  table(emptySumVectorC1)
  # emptySumVector
  # -1  1 
  # 38 18

  
  
  # gProfiler analysis of all genes
  gprofilersameC1 <- GProfilerAnalysis(GeneIDs = list(MethylMixResultsTSSSNFC1$MethylationDrivers),
                                           Organism = "hsapiens",
                                           OrderedQuery = TRUE,
                                           PvalAlphaLevel = 0.05,
                                           PositiveorNegFC = NA, # with respect to expression
                                           ConditionName = "C1SNF",
                                           ProduceImages = "Yes", 
                                           PNGorPDF = "png")
  gprofilersameC1$shortLink
  # "https://biit.cs.ut.ee/gplink/l/Jy_lFejXRI"
  dim(gprofilersameC1$TableAllValues) # 6
  
  
  
  # Plot the most famous methylated gene for glioblastoma
  GeneName <- "TUBA3C"
  MethylMixResultsTSSSNFC1$MixtureStates$TUBA3C
  
  # look at stage
  matchClinicalC1 <- match(colnames(BetaTSSMatrixT1[, matchSNF101C1]), ClinicalFile_T1$SAMPLE_ID)
  table(MethylMixResultsTSSSNFC1$Classifications[which(rownames(MethylMixResultsTSSSNFC1$Classifications) == GeneName), ],
        ClinicalFile_T1$STAGE[matchClinicalC1])
  chisq.test(table(MethylMixResultsTSSSNFC1$Classifications[which(rownames(MethylMixResultsTSSSNFC1$Classifications) == GeneName), ],
                   ClinicalFile_T1$STAGE[matchClinicalC1]))
  plots <- MethylMix::MethylMix_PlotModel(GeneName = GeneName, 
                                          MixtureModelResults = MethylMixResultsTSSSNFC1, 
                                          METcancer = BetaTSSMatrixT1[, matchSNF101C1],
                                          GEcancer = as.matrix(RNAseqProtCodT3MatName[, matchSNF101C1RNaseq]),
                                          METnormal = BetaTSSMatrixT1[, c(166:170)])
  plots$MixtureModelPlot
  plots$CorrelationPlot
  
  
  # my own plotting
  RNAseqProtCodT3MatNameC1 <- as.matrix(RNAseqProtCodT3MatName[, matchSNF101C1RNaseq])
  normalBeta <- as.vector(c(BetaTSSMatrixT1[, c(166:170)][which(rownames(BetaTSSMatrixT1[, c(166:170)]) == GeneName), ][1, ], rep(NA, 48)))
  tumorBeta <- as.vector(BetaTSSMatrixT1[, matchSNF101C1][which(rownames(BetaTSSMatrixT1[, matchSNF101C1]) == GeneName), ][1, ])
  # tumorRNAseq <- RNAseqProtCodT3MatNameC1[which(rownames(RNAseqProtCodT3MatNameC1) == GeneName), ]
  GeneNameMyPlot <- data.frame(
    normalBeta = normalBeta,
    tumorBeta = tumorBeta
    # tumorRNAseq = tumorRNAseq
    )
  
  GeneNameMyPlot %>%
    reshape2::melt() %>%
    ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                      ylab = "Value", 
                      font.label = list(size = 20, color = "black")) +
    ggtitle(paste0(GeneName)) +
    stat_compare_means() + stat_n_text()
  
  ggplot2::ggplot(data = GeneNameMyPlot) + 
    geom_point(mapping = aes(x = normalBeta, y = normalRNAseq)) +
    labs(title = "RNF180",
         x = "normalized methylation beta", 
         y = "normalized RNAseq count",
         colour = "Event") +
    theme_bw() +
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold")) +
    theme(aspect.ratio = 1, 
          legend.position = "right", 
          panel.background = element_rect(colour = "black", size=1.5),  
          axis.title =  element_text(face = "bold"))
  
  
  all(which((MethylMixResultsTSSSNFC1$MethylationDrivers %in% rownames(RNAseqProtCodT3MatName)) == TRUE)) # 1025
  
  

  # selecting entries with no spurious clusters
  sizeSpuriousClusterExclude <- 10
  saveGeneNamesC1 <- NA
  cat("\n Methylation drivers that should be kept:")
  for (i in 1:nrow(MethylMixResultsTSSSNFC1$Classifications)) {
    #cat("\n i = ", i, "\n")
    if(length(unique(MethylMixResultsTSSSNFC1$Classifications[i, ])) > 1) { # if more than one component
      if(all((table(MethylMixResultsTSSSNFC1$Classifications[i, ]) > sizeSpuriousClusterExclude) == TRUE)) {
        saveGeneNamesC1 <- c(saveGeneNamesC1, rownames(MethylMixResultsTSSSNFC1$Classifications)[i])
        }
    }
  }
  
  # Arranging such genes with no suprious clusters by DM values
  matchNoSpurious <- match(saveGeneNamesC1[-1], names(MethylMixResultsTSSSNFC1$MixtureStates))
  saveTopResults <- list()
  for (i in 1:length(matchNoSpurious)) {
    if(length(MethylMixResultsTSSSNFC1$MixtureStates[[matchNoSpurious[i]]]) > 1) {
      saveTopResults[[i]] <- max(abs(MethylMixResultsTSSSNFC1$MixtureStates[[matchNoSpurious[i]]]))
    } else {
      saveTopResults[[i]] <- MethylMixResultsTSSSNFC1$MixtureStates[[matchNoSpurious[i]]]
    }
  }
  names(saveTopResults) <- names(MethylMixResultsTSSSNFC1$MixtureStates)[matchNoSpurious]
  # sort values from highest to lowest
  sort(unlist(saveTopResults), decreasing = TRUE)
  
  
#### Applying to our data - TSS; using edgeR norm (non-log) RNAseq data - C2 Samples (18Nov2020)  ####
  
  #testing <- MethylMix::ClusterProbes(MET_Cancer = BetaMatrix_T1[, c(165:170)],
  #              MET_Normal = BetaMatrix_T1[, c(1:164)])
  #testing <- readRDS("MethylMixr.rds")
  #names(testing)
  #testing$MET_Cancer_Clustered
  #testing$MET_Normal_Clustered
  #testing$ProbeMapping
  
  
  #  Methylation data probes to genes
  EPICAnnotationFile <- data.frame(EPICAnnotationFile,
                                   UCSC_RefGene_NameOnly = sub("\\;.*", "", EPICAnnotationFile$UCSC_RefGene_Name))
  dim(EPICAnnotationFile) # 866836     48
  BetaTSSMatrixT1 <- BetaMatrix_T1[matchTSSProbes, ]
  matchTSSProbesMatch <- match(rownames(BetaTSSMatrixT1), EPICAnnotationFile$V1)
  rownames(BetaTSSMatrixT1) <- EPICAnnotationFile$UCSC_RefGene_NameOnly[matchTSSProbesMatch]
  if(length(which(is.na(rownames(BetaTSSMatrixT1)) == TRUE)) > 0) {
    BetaTSSMatrixT1 <- BetaTSSMatrixT1[- which(is.na(rownames(BetaTSSMatrixT1)) == TRUE), ]
  }
  dim(BetaTSSMatrixT1) # 114897    170
  
  
  # match methylation samples with RNAseq samples
  matchBetaRNAseq <- match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized[, c(11:131)]), 
                           colnames(BetaTSSMatrixT1))
  
  
  # RNA seq
  dim(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized)
  # 44219   132
  range(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized)
  # 0.00 92597.52
  range(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedLOG)
  # -2.72230 16.49869
  
  # RNAseq data ENEMBLidS to genes
  matchRNAseqENSMBLid <- match(rownames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized), ENSMBLid$gene)
  RNAseqNormalizedCountDataTabExtended <- data.frame(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized, 
                                                     ENSEMBLID = rownames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized),
                                                     chr = ENSMBLid$chr[matchRNAseqENSMBLid], 
                                                     name = ENSMBLid$name[matchRNAseqENSMBLid],
                                                     type = ENSMBLid$type[matchRNAseqENSMBLid])
  # remove NA values
  RNAseqNormalizedCountDataTabExtended <- RNAseqNormalizedCountDataTabExtended[- which(is.na(RNAseqNormalizedCountDataTabExtended$type) == TRUE), ]
  dim(RNAseqNormalizedCountDataTabExtended) # 31467   136
  
  # select protein coding only genes
  RNAseqProtCodT3Mat <- RNAseqNormalizedCountDataTabExtended[which((RNAseqNormalizedCountDataTabExtended$type == "protein_coding") == TRUE), ]
  dim(RNAseqProtCodT3Mat) # 17190   136
  # select those not in sex chromosomes
  RNAseqProtCodT3Mat <- RNAseqProtCodT3Mat[! (substr(RNAseqProtCodT3Mat$chr, 4, 5) == "X" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "Y" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "M"), ]
  dim(RNAseqProtCodT3Mat) # 17190   136
  table(RNAseqProtCodT3Mat$chr)
  
  matchSNF101 <- match(AllClustersLab$ID[- which(is.na(AllClustersLab$SNFClust) == TRUE)], colnames(RNAseqProtCodT3Mat))
  RNAseqProtCodT3MatName <- RNAseqProtCodT3Mat[, matchSNF101]
  rownames(RNAseqProtCodT3MatName) <- RNAseqProtCodT3Mat$name
  range(RNAseqProtCodT3MatName) # 0.00 17766.33
  dim(RNAseqProtCodT3MatName) # 17190   101
  
  
  AllClustersLab <- read.csv(file = "InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv")
  
  matchSNF101C2 <- match(AllClustersLab$ID[which(AllClustersLab$SNFClust == 2)], colnames(BetaTSSMatrixT1))
  matchSNF101C2RNaseq <- match(AllClustersLab$ID[which(AllClustersLab$SNFClust == 2)], colnames(RNAseqProtCodT3MatName))
  
  
  # check colnames
  identical(colnames(RNAseqProtCodT3MatName[, matchSNF101C2RNaseq]), 
            colnames(BetaTSSMatrixT1[, matchSNF101C2])) # TRUE
  
  
  set.seed(1234)
  MethylMixResultsTSSSNFC2 <- MethylMix::MethylMix(METcancer = as.matrix(BetaTSSMatrixT1[, matchSNF101C2]), 
                                                   GEcancer = as.matrix(RNAseqProtCodT3MatName[, matchSNF101C2RNaseq]), 
                                                   METnormal = as.matrix(BetaTSSMatrixT1[, c(166:170)]),
                                                   NoNormalMode = TRUE)
  # saveRDS(MethylMixResultsTSSSNFC2, file = "41_MethylMixResultsTSSC2_SNF_Nov2020.rds")
  MethylMixResultsTSSSNFC2 <- readRDS(file = "41_MethylMixResultsTSSC2_SNF_Nov2020.rds")
  
  names(MethylMixResultsTSSSNFC2) # 6 names
  MethylMixResultsTSSSNFC2$MethylationDrivers  # 82; genes identified as transcriptionally predictive and differentially methylated 
  length(MethylMixResultsTSSSNFC2$MethylationDrivers) # 241
  MethylMixResultsTSSSNFC2$NrComponents # number of methylation states found for each driver gene.
  MethylMixResultsTSSSNFC2$MixtureStates # DM-values for each driver gene (is this  mean??)
  length(MethylMixResultsTSSSNFC2$MixtureStates) # 241
  
  
  # Most DM value (with 1 component) if - =  tumor hypomethylated compared to normal
  # Most DM values (with 1 component) if + = tumor hypermethylated compared to normal
  # Save entries with only one component and see if tumor hyper or hypo methylated
  genesToLoop <- which(MethylMixResultsTSSSNFC2$NrComponents == 1)
  emptySumVector <- vector()
  for(i in 1:length(genesToLoop)) { 
    cat("\n", i)
    if(MethylMixResultsTSSSNFC2$MixtureStates[[genesToLoop[i]]] < 0) {
      emptySumVector[i] <- -1
    } else if(MethylMixResultsTSSSNFC2$MixtureStates[[genesToLoop[i]]] > 0) {
      emptySumVector[i] <- 1
    } 
  }
  
  table(emptySumVector)
  # emptySumVector
  # -1  1 
  # 85 46
  
  unlist(MethylMixResultsTSSSNFC2$MixtureStates)[order(abs(unlist(MethylMixResultsTSSSNFC2$MixtureStates)), decreasing = TRUE)]
  length(unlist(MethylMixResultsTSSSNFC2$MixtureStates)[order(abs(unlist(MethylMixResultsTSSSNFC2$MixtureStates)), decreasing = TRUE)]) # 108
  # look at reverse order
  tail(unlist(MethylMixResultsTSSSNFC2$MixtureStates)[order(unlist(MethylMixResultsTSSSNFC2$MixtureStates), decreasing = TRUE)])
  
  
  
  MethylMixResultsTSSSNFC2$MethylationStates # DM-values for all driver genes (rows) and all samples (columns)
  MethylMixResultsTSSSNFC2$Classifications # Matrix with integers indicating to which mixture component
  # each cancer sample was assigned to, for each gene.
  MethylMixResultsTSSSNFC2$Models # Beta mixture model parameters for each driver gene.
  names(MethylMixResultsTSSSNFC2$Models)
  MethylMixResultsTSSSNFC2$Models$ZCCHC24$llike
  
  
  # gProfiler analysis of C1 genes
  gprofilersameC2MethylMix <- GProfilerAnalysis(GeneIDs = list(MethylMixResultsTSSSNFC2$MethylationDrivers),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = NA, # with respect to expression
                                                ConditionName = "allinC2",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
  gprofilersameC2MethylMix$shortLink
  # "https://biit.cs.ut.ee/gplink/l/-FEzEpRNRU"
  
  dim(gprofilersameC2MethylMix$TableAllValues) # 57  6
  
  
  # Plot the most famous methylated gene for glioblastoma
  GeneName <- "GPR25"
  MethylMixResultsTSSSNFC2$MixtureStates$GPR25
  sum(MethylMixResultsTSSSNFC2$MixtureStates$GPR25)
  
  # look at stage
  matchClinicalC2 <- match(colnames(BetaTSSMatrixT1[, matchSNF101C2]), ClinicalFile_T1$SAMPLE_ID)
  table(MethylMixResultsTSSSNFC2$Classifications[which(rownames(MethylMixResultsTSSSNFC2$Classifications) == GeneName), ],
        ClinicalFile_T1$STAGE[matchClinicalC2])
  chisq.test(table(MethylMixResultsTSSSNFC2$Classifications[which(rownames(MethylMixResultsTSSSNFC2$Classifications) == GeneName), ],
                   ClinicalFile_T1$STAGE[matchClinicalC2]))
  plots <- MethylMix::MethylMix_PlotModel(GeneName = GeneName, 
                                          MixtureModelResults = MethylMixResultsTSSSNFC2, 
                                          METcancer = BetaTSSMatrixT1[, matchSNF101C2],
                                          GEcancer = as.matrix(RNAseqProtCodT3MatName[, matchSNF101C2RNaseq]),
                                          METnormal = BetaTSSMatrixT1[, c(166:170)])
  plots$MixtureModelPlot
  plots$CorrelationPlot
  
  
  # my own plotting
  RNAseqProtCodT3MatNameC2 <- as.matrix(RNAseqProtCodT3MatName[, matchSNF101C2RNaseq])
  normalBeta <- as.vector(c(BetaTSSMatrixT1[, c(166:170)][which(rownames(BetaTSSMatrixT1[, c(166:170)]) == GeneName), ][1, ], rep(NA, 43)))
  tumorBeta <- as.vector(BetaTSSMatrixT1[, matchSNF101C2][which(rownames(BetaTSSMatrixT1[, matchSNF101C2]) == GeneName), ][1, ])
  # tumorRNAseq <- RNAseqProtCodT3MatNameC1[which(rownames(RNAseqProtCodT3MatNameC1) == GeneName), ]
  GeneNameMyPlot <- data.frame(
    normalBeta = normalBeta,
    tumorBeta = tumorBeta
    # tumorRNAseq = tumorRNAseq
  )
  
  GeneNameMyPlot %>%
    reshape2::melt() %>%
    ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                      ylab = "Value", 
                      font.label = list(size = 20, color = "black")) +
    ggtitle(paste0(GeneName)) +
    stat_compare_means() + stat_n_text()
  
  ggplot2::ggplot(data = GeneNameMyPlot) + 
    geom_point(mapping = aes(x = normalBeta, y = normalRNAseq)) +
    labs(title = "RNF180",
         x = "normalized methylation beta", 
         y = "normalized RNAseq count",
         colour = "Event") +
    theme_bw() +
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold")) +
    theme(aspect.ratio = 1, 
          legend.position = "right", 
          panel.background = element_rect(colour = "black", size=1.5),  
          axis.title =  element_text(face = "bold"))
  
  
  
  all(which((MethylMixResultsTSSSNFC2$MethylationDrivers %in% rownames(RNAseqProtCodT3MatName)) == TRUE)) # 1025
  
  
  
  # selecting entries with no spurious clusters
  sizeSpuriousClusterExclude <- 10
  saveGeneNamesC2 <- NA
  cat("\n Methylation drivers that should be kept:")
  for (i in 1:nrow(MethylMixResultsTSSSNFC2$Classifications)) {
    #cat("\n i = ", i, "\n")
    if(length(unique(MethylMixResultsTSSSNFC2$Classifications[i, ])) > 1) { # if more than one component
      if(all((table(MethylMixResultsTSSSNFC2$Classifications[i, ]) > sizeSpuriousClusterExclude) == TRUE)) {
        saveGeneNamesC2 <- c(saveGeneNamesC2, rownames(MethylMixResultsTSSSNFC2$Classifications)[i])
      }
    }
  }
  
  # Arranging such genes with no suprious clusters by DM values
  matchNoSpurious <- match(saveGeneNamesC2[-1], names(MethylMixResultsTSSSNFC2$MixtureStates))
  saveTopResults <- list()
  for (i in 1:length(matchNoSpurious)) {
    if(length(MethylMixResultsTSSSNFC2$MixtureStates[[matchNoSpurious[i]]]) > 1) {
      saveTopResults[[i]] <- max(abs(MethylMixResultsTSSSNFC2$MixtureStates[[matchNoSpurious[i]]]))
    } else {
      saveTopResults[[i]] <- MethylMixResultsTSSSNFC2$MixtureStates[[matchNoSpurious[i]]]
    }
  }
  names(saveTopResults) <- names(MethylMixResultsTSSSNFC2$MixtureStates)[matchNoSpurious]
  # sort values from highest to lowest
  sort(unlist(saveTopResults), decreasing = TRUE)
  
#### Looking at overalp #####
  
  mainlist <- list(All = MethylMixResultsTSSSNF$MethylationDrivers,
                   C1 = MethylMixResultsTSSSNFC1$MethylationDrivers,
                   C2 = MethylMixResultsTSSSNFC2$MethylationDrivers)
  
  OutputSNFClusters <- VennDiagramAnalysis(MainList = mainlist, 
                      Labels = c("All","C1","C2"),
                      FigureGenerate = "Yes", 
                      ImageName = "MethylMixResultsSNFclusters", 
                      PNGorPDF = "png")
  
  NumbComponents <- data.frame("All" = table(MethylMixResultsTSSSNF$NrComponents), 
                               "C1" = table(MethylMixResultsTSSSNFC2$NrComponents),
                               "C2" = table(MethylMixResultsTSSSNFC2$NrComponents))
  #manually adjust
  NumbComponents$C1.Freq <- c(int(table(MethylMixResultsTSSSNFC1$NrComponents)), 0)
  NumbComponents <- melt(NumbComponents)
  
  #1. summary barplot - adapted based on code from Karin I. 
  NumbComponentsPlot <- ggpubr::ggbarplot(NumbComponents, 
                                          x = "variable", y ="value", fill = "All.Var1",
                                          label = TRUE,   ylab = "Number of driver genes", 
                                          lab.pos = "in",
                                          xlab = "Contrast", 
                                          font.label = list(size = 20, color = "black")) +
                                          ggtitle("Comparison of All, SNF C1 and SNF C2") +
                                          rotate_x_text(90)
                            

  # Differences between C1 and C2
  diffC1C2 <- setdiff(MethylMixResultsTSSSNFC2$MethylationDrivers, 
                      MethylMixResultsTSSSNFC1$MethylationDrivers)
  length(diffC1C2) # 226
  
  # Same between C1 and C2
  sameC1C2 <- intersect(MethylMixResultsTSSSNFC2$MethylationDrivers, 
                      MethylMixResultsTSSSNFC1$MethylationDrivers)
  length(sameC1C2) # 15
  #"SOD1"   "ZNF677" "MKRN3"  "SPEF2"  "ZNF528" "ZNF667" "SKAP1"  "H3F3C"  "TMED3"  "NLRP2" 
  # "AFF3"   "SEC62"  "PSMA8"  "RASSF6" "ZNF608"
  
  gprofilersameC1C2 <- GProfilerAnalysis(GeneIDs = list(sameC1C2),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = NA, # with respect to expression
                                                ConditionName = "sameC1C2",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
  
#### Applying to our data - TSS; using edgeR norm (non-log) RNAseq data - 121 Samples ####
  
  #testing <- MethylMix::ClusterProbes(MET_Cancer = BetaMatrix_T1[, c(165:170)],
  #              MET_Normal = BetaMatrix_T1[, c(1:164)])
  #testing <- readRDS("MethylMixr.rds")
  #names(testing)
  #testing$MET_Cancer_Clustered
  #testing$MET_Normal_Clustered
  #testing$ProbeMapping
  
  
  #  Methylation data probes to genes
  EPICAnnotationFile <- data.frame(EPICAnnotationFile,
                                   UCSC_RefGene_NameOnly = sub("\\;.*", "", EPICAnnotationFile$UCSC_RefGene_Name))
  dim(EPICAnnotationFile) # 866836     48
  BetaTSSMatrixT1 <- BetaMatrix_T1[matchTSSProbes, ]
  matchTSSProbesMatch <- match(rownames(BetaTSSMatrixT1), EPICAnnotationFile$V1)
  rownames(BetaTSSMatrixT1) <- EPICAnnotationFile$UCSC_RefGene_NameOnly[matchTSSProbesMatch]
  if(length(which(is.na(rownames(BetaTSSMatrixT1)) == TRUE)) > 0) {
    BetaTSSMatrixT1 <- BetaTSSMatrixT1[- which(is.na(rownames(BetaTSSMatrixT1)) == TRUE), ]
  }
  dim(BetaTSSMatrixT1) # 114897    170
  
  
  # match methylation samples with RNAseq samples
  matchBetaRNAseq <- match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized[, c(11:131)]), 
                           colnames(BetaTSSMatrixT1))
  
  
  # RNA seq
  dim(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized)
  # 44219   132
  range(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized)
  # 0.00 92597.52
  range(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedLOG)
  # -2.72230 16.49869
  
  # RNAseq data ENEMBLidS to genes
  matchRNAseqENSMBLid <- match(rownames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized), ENSMBLid$gene)
  RNAseqNormalizedCountDataTabExtended <- data.frame(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized, 
                                                     ENSEMBLID = rownames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized),
                                                     chr = ENSMBLid$chr[matchRNAseqENSMBLid], 
                                                     name = ENSMBLid$name[matchRNAseqENSMBLid],
                                                     type = ENSMBLid$type[matchRNAseqENSMBLid])
  # remove NA values
  RNAseqNormalizedCountDataTabExtended <- RNAseqNormalizedCountDataTabExtended[- which(is.na(RNAseqNormalizedCountDataTabExtended$type) == TRUE), ]
  dim(RNAseqNormalizedCountDataTabExtended) # 31467   136
  
  # select protein coding only genes
  RNAseqProtCodT3Mat <- RNAseqNormalizedCountDataTabExtended[which((RNAseqNormalizedCountDataTabExtended$type == "protein_coding") == TRUE), ]
  dim(RNAseqProtCodT3Mat) # 17190   136
  # select those not in sex chromosomes
  RNAseqProtCodT3Mat <- RNAseqProtCodT3Mat[! (substr(RNAseqProtCodT3Mat$chr, 4, 5) == "X" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "Y" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "M"), ]
  dim(RNAseqProtCodT3Mat) # 17190   136
  table(RNAseqProtCodT3Mat$chr)
  
  # match samples between 
  BetaTSSMatrix_T1MatMatched121 <- match(colnames(RNAseqProtCodT3Mat), colnames(BetaTSSMatrixT1[, c(11:165)]))
  BetaTSSMatrix_T1MatMatched121 <- BetaTSSMatrix_T1MatMatched121[! is.na(BetaTSSMatrix_T1MatMatched121)]
  length(BetaTSSMatrix_T1MatMatched121) # 121
  dim(BetaTSSMatrixT1[, c(11:165)][, BetaTSSMatrix_T1MatMatched121]) # 114897    121
  BetaTSSMatrix121 <- BetaTSSMatrixT1[, c(11:165)][, BetaTSSMatrix_T1MatMatched121]
  dim(BetaTSSMatrix121) # 114897    121
  
  
  RNAseqProtCodT3MatMatMatched121 <- match(colnames(BetaTSSMatrix121), 
                                           colnames(RNAseqProtCodT3Mat))
  RNAseqProtCod121 <- RNAseqProtCodT3Mat[, RNAseqProtCodT3MatMatMatched121]
  dim(RNAseqProtCod121) #  17190   121
  rownames(RNAseqProtCod121) <- RNAseqProtCodT3Mat$name
  range(RNAseqProtCod121) # 0.00 17766.33
  dim(RNAseqProtCod121) # 17190   121
  
   
  # check colnames
  identical(colnames(RNAseqProtCod121), colnames(BetaTSSMatrix121)) # TRUE
  
  
  set.seed(1234)
  MethylMixResultsTSS121 <- MethylMix::MethylMix(METcancer = as.matrix(BetaTSSMatrix121), 
                                                 GEcancer = as.matrix(RNAseqProtCod121), 
                                                 METnormal = BetaTSSMatrixT1[, c(166:170)],
                                                 NoNormalMode = TRUE)
  # saveRDS(MethylMixResultsTSS121, file = "41_MethylMixResultsTSS_121_Dec2020.rds")
  # saveRDS(BetaTSSMatrix121, file = "BetaTSSMatrix121.rds")
  # saveRDS(BetaTSSMatrixT1[, c(166:170)], file = "BetaTSSMatrixRLN.rds")
  # saveRDS(RNAseqProtCod121, file = "RNAseqProtCod121.rds")
  # BetaTSSMatrixT1Match <- match(colnames(RNAseqProtCod121), colnames(BetaTSSMatrixT1))
  # dim(BetaTSSMatrixT1[, BetaTSSMatrixT1Match]) #  114897    121
  # BetaTSSMatrixT1[, c(166:170)]
  MethylMixResultsTSS121 <- readRDS(file = "41_MethylMixResultsTSS_121_Dec2020.rds")
  
  names(MethylMixResultsTSS121) # 6 names
  head(MethylMixResultsTSS121$NrComponents) # number of components per each driver gene
  head(MethylMixResultsTSS121$MethylationStates) # DM-values for each driver gene (is this  mean??)
  head(MethylMixResultsTSS121$MixtureStates) # Mean methylation of each component
  head(MethylMixResultsTSS121$Models) # likelihood, and other values of mixture beta model
  head(MethylMixResultsTSS121$Classifications) # Classification of each patient
  head(MethylMixResultsTSS121$MethylationDrivers)  # 370 ; genes identified as transcriptionally predictive and differentially methylated 
  
  length(MethylMixResultsTSS121$MethylationDrivers) # 370
  MethylMixResultsTSS121$NrComponents # number of methylation states found for each driver gene.
  length(MethylMixResultsTSS121$MixtureStates) # 370
  
  unlist(MethylMixResultsTSS121$MixtureStates)[order(abs(unlist(MethylMixResultsTSS121$MixtureStates)), decreasing = TRUE)]
  length(unlist(MethylMixResultsTSS121$MixtureStates)[order(abs(unlist(MethylMixResultsTSS121$MixtureStates)), decreasing = TRUE)]) # 462
  # look at reverse order
  tail(unlist(MethylMixResultsTSS121$MixtureStates)[order(unlist(MethylMixResultsTSS121$MixtureStates), decreasing = TRUE)])
  
  
  
  MethylMixResultsTSS121$MethylationStates # DM-values for all driver genes (rows) and all samples (columns)
  MethylMixResultsTSS121$Classifications # Matrix with integers indicating to which mixture component
  # each cancer sample was assigned to, for each gene.
  MethylMixResultsTSS121$Models # Beta mixture model parameters for each driver gene.
  names(MethylMixResultsTSS121$Models)
  MethylMixResultsTSS121$Models$ZCCHC24$llike
  
  
  table(MethylMixResultsTSS121$NrComponents)
  # 1   2   3 
  # 271  98   1
  
  # Most DM value (with 1 component) if - =  tumor hypomethylated compared to normal
  # Most DM values (with 1 component) if + = tumor hypermethylated compared to normal
  # Save entries with only one component and see if tumor hyper or hypo methylated
  genesToLoop <- which(MethylMixResultsTSS121$NrComponents == 1)
  emptySumVector121 <- vector()
  for(i in 1:length(genesToLoop)) { 
    cat("\n", i)
    if(MethylMixResultsTSS121$MixtureStates[[genesToLoop[i]]] < 0) {
      emptySumVector121[i] <- -1
    } else if(MethylMixResultsTSS121$MixtureStates[[genesToLoop[i]]] > 0) {
      emptySumVector121[i] <- 1
    } 
  }
  
  table(emptySumVector121)
  # emptySumVector
  # -1  1 
  # 200  71 
  
  MethylMixResultsTSS121HypoGenes <- 
    names(MethylMixResultsTSS121$NrComponents[genesToLoop][which(emptySumVector121 == -1)])
  sort(noquote(paste0(MethylMixResultsTSS121HypoGenes, sep = ",")))
  
  MethylMixResultsTSS121HyperGenes <- 
    names(MethylMixResultsTSS121$NrComponents[genesToLoop][which(emptySumVector121 == 1)])
  sort(noquote(paste0(MethylMixResultsTSS121HyperGenes, sep = ",")))
  
  
  
# from those with 2 component, get the those with hypomethylation
# genesToLoop2 <- which(MethylMixResultsTSS121$NrComponents == 1)
genesToLoop2 <- which(MethylMixResultsTSS121$NrComponents > 1)
pathNow <- getwd()  
plotList <- list()
for (i in 1:length(genesToLoop2)) { 
  cat("\n", i )
  GeneName <- names(genesToLoop2[i])
  plots <- MethylMix::MethylMix_PlotModel(GeneName = GeneName, 
                                          MixtureModelResults = MethylMixResultsTSS121, 
                                          METcancer = BetaTSSMatrix121,
                                          GEcancer = as.matrix(RNAseqProtCod121),
                                          METnormal = BetaTSSMatrixT1[, c(166:170)])
  plotList[[i]] <- plots[[1]]
}


for (i in 1:length(genesToLoop2)) {
  cat("\n", i )
  GeneName <- names(genesToLoop2[i])
  file_name = paste("iris_plot_", i, ".tiff", sep = "")
  grDevices::png(paste0(pathNow,"/img/MethylMix_n=121_more2components_12Feb2021/1Component/Density/41_MethylMix_", i,"_", GeneName, ".png"))
  print(plotList[[i]])
  grDevices::dev.off()
}
  

  
  
  
  
  
  
  # Plot the most famous methylated gene for glioblastoma
  GeneName <- "POU4F1"
  MethylMixResultsTSS121$MixtureStates$MKRN3
  
  # look at stage
  matchClinical121 <- match(colnames(BetaTSSMatrix121), ClinicalFile_T1$SAMPLE_ID)
  table(MethylMixResultsTSS121$Classifications[which(rownames(MethylMixResultsTSS121$Classifications) == GeneName), ],
        ClinicalFile_T1$STAGE[matchClinical121])
  chisq.test(table(MethylMixResultsTSS121$Classifications[which(rownames(MethylMixResultsTSS121$Classifications) == GeneName), ],
                   ClinicalFile_T1$STAGE[matchClinical121]))
  plots <- MethylMix::MethylMix_PlotModel(GeneName = GeneName, 
                                          MixtureModelResults = MethylMixResultsTSS121, 
                                          METcancer = BetaTSSMatrix121,
                                          GEcancer = as.matrix(RNAseqProtCod121),
                                          METnormal = BetaTSSMatrixT1[, c(166:170)])
  plots$MixtureModelPlot
  plots$CorrelationPlot

  
  
  
  # my own plotting - 121
  if(length(which(rownames(BetaTSSMatrixT1[, c(166:170)]) == GeneName)) == 1) {
    normalBeta <- as.vector(c(BetaTSSMatrixT1[, c(166:170)][which(rownames(BetaTSSMatrixT1[, c(166:170)]) == GeneName), ], rep(NA, 116)))
  } else {
    normalBeta <- as.vector(c(BetaTSSMatrixT1[, c(166:170)][which(rownames(BetaTSSMatrixT1[, c(166:170)]) == GeneName), ][1, ], rep(NA, 116)))
  }
  
  if(length(which(rownames(BetaTSSMatrix121) == GeneName)) == 1) {
    tumorBeta <- as.vector(BetaTSSMatrix121[which(rownames(BetaTSSMatrix121) == GeneName), ])
  } else {
    tumorBeta <- as.vector(BetaTSSMatrix121[which(rownames(BetaTSSMatrix121) == GeneName), ][1, ])
  }
  
  # tumorRNAseq <- RNAseqProtCodT3MatNameC1[which(rownames(RNAseqProtCodT3MatNameC1) == GeneName), ]
  GeneNameMyPlot <- data.frame(
    normalBeta = normalBeta,
    tumorBeta = tumorBeta)
  
  
  GeneNameMyPlot %>%
    reshape2::melt() %>%
    ggpubr::ggboxplot(x = "variable", 
                    y = "value", 
                    width = 0.8,
                    size = 1,
                    ylab = " Mean Beta Value", 
                    xlab = "",
                    font.label = list(size = 20, color = "black"),
                    color = "variable",
                    palette = c("#2166ac", "#5aae61"),
                    add = "jitter",
                    legend = "none") +
    ggtitle(paste0(GeneName)) +
    ggpubr::stat_compare_means() +
    EnvStats::stat_n_text()
  
  
  
  ggplot2::ggplot(data = GeneNameMyPlot) + 
    geom_point(mapping = aes(x = normalBeta, y = normalRNAseq)) +
    labs(title = "RNF180",
         x = "normalized methylation beta", 
         y = "normalized RNAseq count",
         colour = "Event") +
    theme_bw() +
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold")) +
    theme(aspect.ratio = 1, 
          legend.position = "right", 
          panel.background = element_rect(colour = "black", size=1.5),  
          axis.title =  element_text(face = "bold"))
  
  
  
  all(which((MethylMixResultsTSS121$MethylationDrivers %in% rownames(RNAseqNormalizedCountDataTabMatrix)) == TRUE)) # 1025
  
  
  
  
  # selecting entries with no spurious clusters
  
  # selecting entries with no spurious clusters
  sizeSpuriousClusterExclude <- 10
  saveGeneNamesAll121 <- NA
  cat("\n Methylation drivers that should be kept:")
  for (i in 1:nrow(MethylMixResultsTSS121$Classifications)) {
    #cat("\n i = ", i, "\n")
    if(length(unique(MethylMixResultsTSS121$Classifications[i, ])) > 1) { # if more than one component
      if(all((table(MethylMixResultsTSS121$Classifications[i, ]) > sizeSpuriousClusterExclude) == TRUE)) {
        saveGeneNamesAll121 <- c(saveGeneNamesAll121, rownames(MethylMixResultsTSS121$Classifications)[i])
      }
    }
  }
  
  # Arranging such genes with no suprious clusters by DM values
  matchNoSpurious <- match(saveGeneNamesAll121[-1], names(MethylMixResultsTSS121$MixtureStates))
  saveTopResults <- list()
  for (i in 1:length(matchNoSpurious)) {
    if(length(MethylMixResultsTSS121$MixtureStates[[matchNoSpurious[i]]]) > 1) {
      saveTopResults[[i]] <- max(abs(MethylMixResultsTSS121$MixtureStates[[matchNoSpurious[i]]]))
    } else {
      saveTopResults[[i]] <- MethylMixResultsTSS121$MixtureStates[[matchNoSpurious[i]]]
    }
  }
  names(saveTopResults) <- names(MethylMixResultsTSS121$MixtureStates)[matchNoSpurious]
  # sort values from highest to lowest
  sort(unlist(saveTopResults), decreasing = TRUE)
  
  
  
  
  gprofiler121 <- GProfilerAnalysis(GeneIDs = list( MethylMixResultsTSS121$MethylationDrivers ),
                                         Organism = "hsapiens",
                                         OrderedQuery = TRUE,
                                         PvalAlphaLevel = 0.01,
                                         PositiveorNegFC = NA, # with respect to expression
                                         ConditionName = "121FL",
                                         ProduceImages = "Yes", 
                                         PNGorPDF = "png")  
  
  gprofiler121$shortLink
  # https://biit.cs.ut.ee/gplink/l/WdH9jWhrSG
  
  dim(gprofiler121$TableAllValues) # 150   6
  
  
#### Applying to our data - TSS; using edgeR norm (non-log) RNAseq data - 121 Samples - Determine hyper or hypomethylated genes ####  
  
  # determine if genes are hyper or hypo methylated - wrote my own function using median
  # write a loop that iterate along genes
  
  geneMethylationMethylMix121 <- vector(mode = "list", length = length(MethylMixResultsTSS121$MethylationDrivers))
  names(geneMethylationMethylMix121) <- MethylMixResultsTSS121$MethylationDrivers
  
  for (gene in 1:length(MethylMixResultsTSS121$MethylationDrivers)) {
    print(gene)
    GeneName <- MethylMixResultsTSS121$MethylationDrivers[gene] # select gene
    
    matchBetaRNAseq <- match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized[, c(11:131)]), 
                             colnames(BetaMatrix_T1))
    if(length(which(rownames(BetaTSSMatrixT1[, c(166:170)]) == GeneName)) > 1) { # if gene match multiple probes
      normalBeta <- as.vector(c(BetaTSSMatrixT1[, c(166:170)][which(rownames(BetaTSSMatrixT1[, c(166:170)]) == GeneName), ][1, ], rep(NA, 116)))
      tumorBeta <- as.vector(BetaTSSMatrixT1[, matchBetaRNAseq][which(rownames(BetaTSSMatrixT1[, matchBetaRNAseq]) == GeneName), ][1, ])
    } else { # if gene match a single probe
      normalBeta <- as.vector(c(BetaTSSMatrixT1[, c(166:170)][which(rownames(BetaTSSMatrixT1[, c(166:170)]) == GeneName), ], rep(NA, 116)))
      tumorBeta <- as.vector(BetaTSSMatrixT1[, matchBetaRNAseq][which(rownames(BetaTSSMatrixT1[, matchBetaRNAseq]) == GeneName), ])
    }
    
    normalBetaMedian <- median(normalBeta, na.rm = TRUE)
    tumorBetaMedian <- median(tumorBeta, na.rm = TRUE)
    
    if (normalBetaMedian > tumorBetaMedian) {
      geneMethylationMethylMix121[[gene]] <- "TumorHypomethylated"
    } else if (normalBetaMedian < tumorBetaMedian) {
      geneMethylationMethylMix121[[gene]] <- "TumorHypermethylated"
    } else {
      geneMethylationMethylMix121[[gene]] <- "Unclear"
    }
  }
  
  UnlistGeneMethylationMethylMix121 <- unlist(geneMethylationMethylMix121)
  table(UnlistGeneMethylationMethylMix121)
  #TumorHypermethylated  TumorHypomethylated 
  #116                  254 
  (116 / (116 + 254)) * 100 # 31.35
  (254 / (116 + 254)) * 100 # 68.65
  
  
  
  
  
  
genesToLoop <- which(MethylMixResultsTSS121$NrComponents == 1)
  
# gProfiler analysis by direction   
# 121 - 1 Component results
MethylMixResultsTSS121HypoGenes <- 
    names(MethylMixResultsTSS121$NrComponents[genesToLoop][which(emptySumVector121 == -1)])

MethylMixResultsTSS121HyperGenes <- 
    names(MethylMixResultsTSS121$NrComponents[genesToLoop][which(emptySumVector121 == 1)])
  
# 121 - 2 Component results
HypoGenesComp2 <- c("AFF3", "ATP6V1G3", "CARD11", "CD1C", "CEACAM1", "CLNK", "DNAJB7", "FCRL2", 
                    "FBXO18", "FCRL3", "FRK", "GRIP1", "HEPHL1", "IGLL5", "LARP4B", "LEMD1", "LYPD6B", "MMP7", 
                    "MKRN3", "MYH15", "NLRP2", "P2RY14", "PHYHD1", "PSMA8", "RASSF6", "RGPD2", "SIT1", "SLC39A10", 
                    "SYT5", "SYT17", "TDRD1", "TAGAP", "TDGF1", "TVP23A", "VPREB3")
length(HypoGenesComp2) #35

HyperGenesComp2or3 <- c("ACADM", "AMOTL1", "ARHGAP21", "BEND4", "BMP7", "CARKD", "CAT", "CCSER1", "CHN2", "CNTLN", "EPHA4", "EYA2", "FAAH", "GABRB2", "HOXB4", 
                    "INO80", "INPP5A", "INADL", "MICU3", "NPTX2", "NR3C2", "OXR1", "PCCA", "PCDH7", "POU4F1", "PROX1", "RIC3", "RNF180", 
                    "ROBO1", "SDR42E1", "SKAP1", "SLC35D2", "SURF4", "UCHL1", "USP44", "VAV3", "ZNF471", "ZNF382", "ZCCHC14", "ZNF354C", 
                    "ZNF155", "TSPYL5", "ZNF470", "FAM169A", "ZFP28", "ZNF502", "ZNF677", "ZNF486", "ZNF772", "ZNF354B", "ZNF528", "ZNF813", 
                    "ZNF781", "ZNF667", "ZNF793", "ZNF234", "ZNF512B", "ZSCAN23", "ZNF135", "ZNF568", "ZFP37", "ZNF391", "ZIK1", "IGF2BP1")  
length(HyperGenesComp2or3) # 64
  
HypoGenes <- c(MethylMixResultsTSS121HypoGenes, HypoGenesComp2)
length(unique(HypoGenes)) # 235
write.csv(data.frame(HypoGenes), file = "41_MethylMix_n=121_HypoGenes235.csv")
HyperGenes <- c(MethylMixResultsTSS121HyperGenes, HyperGenesComp2or3)
length(unique(HyperGenes)) # 135
write.csv(data.frame(HyperGenes), file = "41_MethylMix_n=121_HyperGenes135.csv")




#### Applying to our data - TSS; using edgeR norm (non-log) RNAseq data - 121 Samples - Pathway Analysis ####


# gProfiler analysis of 1 component genes
gprofiler121OneCompGenes <- GProfilerAnalysis(GeneIDs = list(c(MethylMixResultsTSS121HypoGenes, 
                                                                     MethylMixResultsTSS121HyperGenes)),
                                                    Organism = "hsapiens",
                                                    OrderedQuery = FALSE,
                                                    PvalAlphaLevel = 0.01,
                                                    PositiveorNegFC = NA, # with respect to expression
                                                    ConditionName = "OneCompGenes",
                                                    ProduceImages = "Yes", 
                                                    PNGorPDF = "png")  
gprofiler121OneCompGenes$shortLink
# https://biit.cs.ut.ee/gplink/l/O1bHpzMpT2
dim(gprofiler121OneCompGenes$TableAllValues) # 171 6
customid56789 <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, "GSEA/GSEAgmtFiles/All.gmt")) # gp__nq4U_Dhqm_okk
gprofiler121OneCompGenes2 <- GProfilerAnalysis(GeneIDs = list(c(MethylMixResultsTSS121HypoGenes, 
                                                                   MethylMixResultsTSS121HyperGenes)),
                                                  Organism = "gp__nq4U_Dhqm_okk",
                                                  OrderedQuery = FALSE,
                                                  PvalAlphaLevel = 0.01,
                                                  PositiveorNegFC = NA, # with respect to expression
                                                  ConditionName = "OneCompGenes",
                                                  ProduceImages = "Yes", 
                                                  PNGorPDF = "png")  
gprofiler121OneCompGenes2$shortLink
# https://biit.cs.ut.ee/gplink/l/E2A8uUKxSv
dim(gprofiler121OneCompGenes2$TableAllValues) # 30 6


# gProfiler analysis of 2 component or more genes
gprofiler121MoreThan2CompGenes <- GProfilerAnalysis(GeneIDs = list(c(HypoGenesComp2, HyperGenesComp2or3)),
                                            Organism = "hsapiens",
                                            OrderedQuery = FALSE,
                                            PvalAlphaLevel = 0.01,
                                            PositiveorNegFC = NA, # with respect to expression
                                            ConditionName = "MoreThan2CompGenes",
                                            ProduceImages = "Yes", 
                                            PNGorPDF = "png")  
gprofiler121MoreThan2CompGenes$shortLink
# https://biit.cs.ut.ee/gplink/l/hlDV3HeRTB
dim(gprofiler121MoreThan2CompGenes$TableAllValues) # 21  6

customid56789 <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, "GSEA/GSEAgmtFiles/All.gmt")) # gp__nq4U_Dhqm_okk
gprofiler121MoreThan2CompGenes2 <- GProfilerAnalysis(GeneIDs = list(c(HypoGenesComp2, HyperGenesComp2or3)),
                                               Organism = "gp__nq4U_Dhqm_okk",
                                               OrderedQuery = FALSE,
                                               PvalAlphaLevel = 0.01,
                                               PositiveorNegFC = NA, # with respect to expression
                                               ConditionName = "MoreThan2CompGenes",
                                               ProduceImages = "Yes", 
                                               PNGorPDF = "png")  
gprofiler121MoreThan2CompGenes2$shortLink
# https://biit.cs.ut.ee/gplink/l/7m3ygWlnSF
dim(gprofiler121MoreThan2CompGenes2$TableAllValues) # 2 6









HypoGenes <- c(MethylMixResultsTSS121HypoGenes, HypoGenesComp2)
length(HypoGenes) # 235
HyperGenes <- c(MethylMixResultsTSS121HyperGenes, HyperGenesComp2or3)
length(HyperGenes) # 135



# gProfiler analysis of Hyper vs Hypo genes
gprofiler121HyperGenes <- GProfilerAnalysis(GeneIDs = list(HyperGenes),
                                  Organism = "hsapiens",
                                  OrderedQuery = FALSE,
                                  PvalAlphaLevel = 0.01,
                                  PositiveorNegFC = "positive", # with respect to expression
                                  ConditionName = "121FLHyperGenes",
                                  ProduceImages = "Yes", 
                                  PNGorPDF = "png")  
gprofiler121HyperGenes$shortLink
# https://biit.cs.ut.ee/gplink/l/Jp19RUGtSS
dim(gprofiler121HyperGenes$TableAllValues) # 59  6

gprofiler121HypoGenes <- GProfilerAnalysis(GeneIDs = list(HypoGenes),
                                           Organism = "hsapiens",
                                           OrderedQuery = FALSE,
                                           PvalAlphaLevel = 0.01,
                                           PositiveorNegFC = "negative", # with respect to expression
                                           ConditionName = "121FLHypoGenes",
                                           ProduceImages = "Yes", 
                                           PNGorPDF = "png")  
gprofiler121HypoGenes$shortLink
# https://biit.cs.ut.ee/gplink/l/9V2NJcpzQK
dim(gprofiler121HypoGenes$TableAllValues) # 24  6
  
  


# Trying Staudt signatures on Hyper and Hypo genes
customid56789 <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, "GSEA/Staudt_signatureDB_030920.gmt")) # 
gprofiler121HypoGenesStaudt <- GProfilerAnalysis(GeneIDs = list(HypoGenes),
                                           Organism = "gp__TBhL_fGa4_KUg",
                                           OrderedQuery = FALSE,
                                           PvalAlphaLevel = 0.01,
                                           PositiveorNegFC = "negative", # with respect to expression
                                           ConditionName = "121FLHypoGenesStaudt",
                                           ProduceImages = "Yes", 
                                           PNGorPDF = "png")  
gprofiler121HypoGenesStaudt$shortLink
# https://biit.cs.ut.ee/gplink/l/fXOFHwPKQT
dim(gprofiler121HypoGenesStaudt$TableAllValues) # 9 6
gprofiler121HyperGenesStaudt <- GProfilerAnalysis(GeneIDs = list(HyperGenes),
                                                 Organism = "gp__TBhL_fGa4_KUg",
                                                 OrderedQuery = FALSE,
                                                 PvalAlphaLevel = 0.01,
                                                 PositiveorNegFC = "positive", # with respect to expression
                                                 ConditionName = "121FLHyperGenesStaudt",
                                                 ProduceImages = "Yes", 
                                                 PNGorPDF = "png")  
gprofiler121HyperGenesStaudt$shortLink
# https://biit.cs.ut.ee/gplink/l/7nISUUWXQ4
dim(gprofiler121HyperGenesStaudt$TableAllValues) # No results




# Hallmarks gmt file
# customidHallmark <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, 
#                                                        "GSEA/GSEAgmtFiles/h.all.v7.2.symbols.gmt")) # gp__ufbn_z4Fu_r4M
# gprofiler121HypoGenesHallmark <- GProfilerAnalysis(GeneIDs = list(HypoGenes),
#                                                  Organism = "gp__ufbn_z4Fu_r4M",
#                                                  OrderedQuery = FALSE,
#                                                  PvalAlphaLevel = 0.05,
#                                                  PositiveorNegFC = "negative", # with respect to expression
#                                                  ConditionName = "121FLHypoGenesStaudt",
#                                                  ProduceImages = "Yes", 
#                                                  PNGorPDF = "png")  
# gprofiler121HypoGenesHallmark$shortLink
# "https://biit.cs.ut.ee/gplink/l/pD5XQVWCRv
# dim(gprofiler121HypoGenesHallmark$TableAllValues) # 3 6
# gprofiler121HyperGenesHallmark <- GProfilerAnalysis(GeneIDs = list(HyperGenes),
#                                                   Organism = "gp__ufbn_z4Fu_r4M",
#                                                   OrderedQuery = FALSE,
#                                                   PvalAlphaLevel = 0.05,
#                                                   PositiveorNegFC = "positive", # with respect to expression
#                                                   ConditionName = "121FLHyperGenesStaudt",
#                                                   ProduceImages = "Yes", 
#                                                   PNGorPDF = "png")  
# gprofiler121HyperGenesHallmark$shortLink
# https://biit.cs.ut.ee/gplink/l/F0VyGkqPSB
# dim(gprofiler121HyperGenesHallmark$TableAllValues) # 1 6




# Curated reactome v7.2 gmt file
# customidReactome <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, 
#                                                        "GSEA/GSEAgmtFiles/c2.cp.reactome.v7.2.symbols.gmt")) # 
# gprofiler121HypoGenesReactome <- GProfilerAnalysis(GeneIDs = list(HypoGenes),
#                                                    Organism = "gp__UZgc_0j9m_Xrc",
#                                                    OrderedQuery = FALSE,
#                                                    PvalAlphaLevel = 0.05,
#                                                    PositiveorNegFC = "negative", # with respect to expression
#                                                    ConditionName = "121FLHypoGenesStaudt",
#                                                    ProduceImages = "Yes", 
#                                                    PNGorPDF = "png")  
# gprofiler121HypoGenesReactome$shortLink
# https://biit.cs.ut.ee/gplink/l/CUXEANO8Q3
# dim(gprofiler121HypoGenesReactome$TableAllValues) # 2 6
# gprofiler121HyperGenesReactome <- GProfilerAnalysis(GeneIDs = list(HyperGenes),
#                                                     Organism = "gp__UZgc_0j9m_Xrc",
#                                                     OrderedQuery = FALSE,
#                                                     PvalAlphaLevel = 0.01,
#                                                     PositiveorNegFC = "positive", # with respect to expression
#                                                     ConditionName = "121FLHyperGenesStaudt",
#                                                     ProduceImages = "Yes", 
#                                                     PNGorPDF = "png")  
# gprofiler121HyperGenesReactome$shortLink
# https://biit.cs.ut.ee/gplink/l/YXuH3fCHR8
# dim(gprofiler121HyperGenesReactome$TableAllValues) # 1 6





# Curated wikipathways v7.2 gmt file
# customidWikipathways <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, 
#                                                        "GSEA/GSEAgmtFiles/c2.cp.wikipathways.v7.2.symbols.gmt")) # 
# gprofiler121HypoGenesWikipathways <- GProfilerAnalysis(GeneIDs = list(HypoGenes),
#                                                    Organism = "gp__oC6Y_RR6i_9Bk",
#                                                    OrderedQuery = FALSE,
#                                                    PvalAlphaLevel = 0.01,
#                                                    PositiveorNegFC = "negative", # with respect to expression
#                                                    ConditionName = "121FLHypoGenesStaudt",
#                                                    ProduceImages = "Yes", 
#                                                    PNGorPDF = "png")  
# gprofiler121HypoGenesWikipathways$shortLink
# https://biit.cs.ut.ee/gplink/l/RMPD55PYTu
# dim(gprofiler121HypoGenesWikipathways$TableAllValues) # 2 6
# gprofiler121HyperGenesWikipathways <- GProfilerAnalysis(GeneIDs = list(HyperGenes),
#                                                     Organism = "gp__oC6Y_RR6i_9Bk",
#                                                     OrderedQuery = FALSE,
#                                                     PvalAlphaLevel = 0.05,
#                                                     PositiveorNegFC = "positive", # with respect to expression
#                                                     ConditionName = "121FLHyperGenesStaudt",
#                                                     ProduceImages = "Yes", 
#                                                     PNGorPDF = "png")  
# gprofiler121HyperGenesWikipathways$shortLink
# https://biit.cs.ut.ee/gplink/l/3rKzXItLS9
# dim(gprofiler121HyperGenesWikipathways$TableAllValues) #1 6





# c5.go.bp.v7.2.symbols gmt file
# customidGObp <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, 
#                                                            "GSEA/GSEAgmtFiles/c5.go.bp.v7.2.symbols.gmt")) # 
# gprofiler121HypoGenesGObp <- GProfilerAnalysis(GeneIDs = list(HypoGenes),
#                                                        Organism = "gp__h3Eq_nNaX_k6E",
#                                                        OrderedQuery = FALSE,
#                                                        PvalAlphaLevel = 0.01,
#                                                        PositiveorNegFC = "negative", # with respect to expression
#                                                        ConditionName = "121FLHypoGenesStaudt",
#                                                        ProduceImages = "Yes", 
#                                                        PNGorPDF = "png")  
# gprofiler121HypoGenesGObp$shortLink
# https://biit.cs.ut.ee/gplink/l/BY1kscx7Tz
# dim(gprofiler121HypoGenesGObp$TableAllValues) # 1 6
# gprofiler121HyperGenesGObp <- GProfilerAnalysis(GeneIDs = list(HyperGenes),
#                                                         Organism = "gp__h3Eq_nNaX_k6E",
#                                                         OrderedQuery = FALSE,
#                                                         PvalAlphaLevel = 0.05,
#                                                         PositiveorNegFC = "positive", # with respect to expression
#                                                         ConditionName = "121FLHyperGenesStaudt",
#                                                         ProduceImages = "Yes", 
#                                                         PNGorPDF = "png")  
# gprofiler121HyperGenesGObp$shortLink
# No results to show
# dim(gprofiler121HyperGenesGObp$TableAllValues) # No results




# c6.all.v7.2.symbols
# customidC6 <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, 
#                                                    "GSEA/GSEAgmtFiles/c6.all.v7.2.symbols.gmt")) # 
# gprofiler121HypoGenesC6 <- GProfilerAnalysis(GeneIDs = list(HypoGenes),
#                                                Organism = "gp__sDkN_wGw3_bMU",
#                                                OrderedQuery = FALSE,
#                                                PvalAlphaLevel = 0.01,
#                                                PositiveorNegFC = "negative", # with respect to expression
#                                                ConditionName = "121FLHypoGenesStaudt",
#                                                ProduceImages = "Yes", 
#                                                PNGorPDF = "png")  
# gprofiler121HypoGenesC6$shortLink
# https://biit.cs.ut.ee/gplink/l/YPqgczSYQF
# dim(gprofiler121HypoGenesC6$TableAllValues) # 1 6
# gprofiler121HyperGenesC6 <- GProfilerAnalysis(GeneIDs = list(HyperGenes),
#                                                 Organism = "gp__sDkN_wGw3_bMU",
#                                                 OrderedQuery = FALSE,
#                                                 PvalAlphaLevel = 0.05,
#                                                 PositiveorNegFC = "positive", # with respect to expression
#                                                 ConditionName = "121FLHyperGenesStaudt",
#                                                 ProduceImages = "Yes", 
#                                                 PNGorPDF = "png")  
# gprofiler121HyperGenesC6$shortLink
# No results to show
# dim(gprofiler121HyperGenesC6$TableAllValues) # No results






# c7.all.v7.2.symbols
# customidC7 <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, 
#                                                  "GSEA/GSEAgmtFiles/c7.all.v7.2.symbols.gmt")) # 
# gprofiler121HypoGenesC7 <- GProfilerAnalysis(GeneIDs = list(HypoGenes),
#                                              Organism = "gp__CCGo_SCaE_BMA",
#                                              OrderedQuery = FALSE,
#                                              PvalAlphaLevel = 0.01,
#                                              PositiveorNegFC = "negative", # with respect to expression
#                                              ConditionName = "121FLHypoGenesStaudt",
#                                              ProduceImages = "Yes", 
#                                              PNGorPDF = "png")  
# gprofiler121HypoGenesC7$shortLink
# https://biit.cs.ut.ee/gplink/l/HFM83f0RS0
# dim(gprofiler121HypoGenesC7$TableAllValues) # 2 6
# gprofiler121HyperGenesC7 <- GProfilerAnalysis(GeneIDs = list(HyperGenes),
#                                               Organism = "gp__CCGo_SCaE_BMA",
#                                               OrderedQuery = FALSE,
#                                               PvalAlphaLevel = 0.01,
#                                               PositiveorNegFC = "positive", # with respect to expression
#                                               ConditionName = "121FLHyperGenesStaudt",
#                                               ProduceImages = "Yes", 
#                                               PNGorPDF = "png")  
# gprofiler121HyperGenesC7$shortLink
# https://biit.cs.ut.ee/gplink/l/-5JtglzQR8
# dim(gprofiler121HyperGenesC7$TableAllValues) # No results


# Trying Combined GMT file from GSEA pathways on Hyper and Hypo genes
customid56789 <- gprofiler2::upload_GMT_file(paste0(RNAseqDirPath, "GSEA/GSEAgmtFiles/All.gmt")) # gp__nq4U_Dhqm_okk
gprofiler121HypoGenesAllGSEA <- GProfilerAnalysis(GeneIDs = list(HypoGenes),
                                                 Organism = "gp__nq4U_Dhqm_okk",
                                                 OrderedQuery = FALSE,
                                                 PvalAlphaLevel = 0.01,
                                                 PositiveorNegFC = "negative", # with respect to expression
                                                 ConditionName = "121FLHypoGenesStaudt",
                                                 ProduceImages = "Yes", 
                                                 PNGorPDF = "png")  
gprofiler121HypoGenesAllGSEA$shortLink
# https://biit.cs.ut.ee/gplink/l/o91WBxJFSc
dim(gprofiler121HypoGenesAllGSEA$TableAllValues) # 7 6
gprofiler121HyperGenesAllGSEA <- GProfilerAnalysis(GeneIDs = list(HyperGenes),
                                                  Organism = "gp__nq4U_Dhqm_okk",
                                                  OrderedQuery = FALSE,
                                                  PvalAlphaLevel = 0.01,
                                                  PositiveorNegFC = "positive", # with respect to expression
                                                  ConditionName = "121FLHyperGenesStaudt",
                                                  ProduceImages = "Yes", 
                                                  PNGorPDF = "png")  
gprofiler121HyperGenesAllGSEA$shortLink
# https://biit.cs.ut.ee/gplink/l/1GNsI2LSTc
dim(gprofiler121HyperGenesAllGSEA$TableAllValues) # 2 6


#### Applying to our data - TSS; using edgeR norm (non-log) RNAseq data - 121 Samples - Bar Chart ####


# make a bar chart of hyper and hypo genes
HypoGenes <- c(MethylMixResultsTSS121HypoGenes, HypoGenesComp2)
length(HypoGenes) # 235
HyperGenes <- c(MethylMixResultsTSS121HyperGenes, HyperGenesComp2or3)
length(HyperGenes) # 135

MethylMix121AllDF <- data.frame(Gene = c(HypoGenes, HyperGenes),
                                MethylationState = c(rep(1, length(MethylMixResultsTSS121HypoGenes)), 
                                                     rep(2, length(HypoGenesComp2)), 
                                                     rep(1, length(MethylMixResultsTSS121HyperGenes)),  
                                                     rep(2, (length(HyperGenesComp2or3) - 1)), 
                                                     rep(3, 1)),
                                MethylationLevel = c(rep("Hypomethylated", length(HypoGenes)), 
                                                     rep("Hypermethylated", length(HyperGenes))))

MethylMix121AllDF <- data.frame(table(MethylMix121AllDF$MethylationState, 
                 MethylMix121AllDF$MethylationLevel))

ComponentColors <- c("1" = "#4d9221", "2" = "#e6f5d0", 
                "3" = "#c51b7d")


MethylMix121AllDFGroup <- MethylMix121AllDF %>%
  ggplot2::ggplot(aes(fill = Var1, y = Freq, x = reorder(Var2, - Freq))) +
  ggplot2::geom_bar(stat = "identity", width = 0.7) +
  #   ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.7) +
  labs(x = "",
       y = "Number of Genes",
       fill = "Number of \nMethylation \nStates") +
  #scale_y_continuous(labels = scales::percent) +
  scale_y_continuous(n.breaks = 10) +
  scale_fill_manual(values = ComponentColors) +
  theme_bw(base_rect_size = 0.6) + 
  theme(axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 0),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black")) +
  theme(aspect.ratio = 1, text = element_text(size = 15))

pathNow <- getwd()
ggplot2::ggsave(paste0(pathNow, "/img/41_MethylMixComponent121.", PNGorPDF),
                width = 5.5, height = 5.5)


  


#### Applying to our data - TSS; using edgeR norm (non-log) RNAseq data - 101 Samples C1 vs C2 ####
# 1 March 2021


#  Methylation data probes to genes
EPICAnnotationFile <- data.frame(EPICAnnotationFile,
                                 UCSC_RefGene_NameOnly = sub("\\;.*", "", EPICAnnotationFile$UCSC_RefGene_Name))
dim(EPICAnnotationFile) # 866836     48
BetaTSSMatrixT1 <- BetaMatrix_T1[matchTSSProbes, ]
matchTSSProbesMatch <- match(rownames(BetaTSSMatrixT1), EPICAnnotationFile$V1)
rownames(BetaTSSMatrixT1) <- EPICAnnotationFile$UCSC_RefGene_NameOnly[matchTSSProbesMatch]
if(length(which(is.na(rownames(BetaTSSMatrixT1)) == TRUE)) > 0) {
  BetaTSSMatrixT1 <- BetaTSSMatrixT1[- which(is.na(rownames(BetaTSSMatrixT1)) == TRUE), ]
}
dim(BetaTSSMatrixT1) # 114897    170


# match methylation samples with RNAseq samples
matchBetaRNAseq <- match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized[, c(11:131)]), 
                         colnames(BetaTSSMatrixT1))


# RNA seq
dim(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized)
# 44219   132
range(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized)
# 0.00 92597.52
range(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedLOG)
# -2.72230 16.49869

# RNAseq data ENEMBLidS to genes
matchRNAseqENSMBLid <- match(rownames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized), ENSMBLid$gene)
RNAseqNormalizedCountDataTabExtended <- data.frame(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized, 
                                                   ENSEMBLID = rownames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized),
                                                   chr = ENSMBLid$chr[matchRNAseqENSMBLid], 
                                                   name = ENSMBLid$name[matchRNAseqENSMBLid],
                                                   type = ENSMBLid$type[matchRNAseqENSMBLid])
# remove NA values
RNAseqNormalizedCountDataTabExtended <- RNAseqNormalizedCountDataTabExtended[- which(is.na(RNAseqNormalizedCountDataTabExtended$type) == TRUE), ]
dim(RNAseqNormalizedCountDataTabExtended) # 31467   136

# select protein coding only genes
RNAseqProtCodT3Mat <- RNAseqNormalizedCountDataTabExtended[which((RNAseqNormalizedCountDataTabExtended$type == "protein_coding") == TRUE), ]
dim(RNAseqProtCodT3Mat) # 17190   136
# select those not in sex chromosomes
RNAseqProtCodT3Mat <- RNAseqProtCodT3Mat[! (substr(RNAseqProtCodT3Mat$chr, 4, 5) == "X" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "Y" | substr(RNAseqProtCodT3Mat$chr, 4, 5) == "M"), ]
dim(RNAseqProtCodT3Mat) # 17190   136
table(RNAseqProtCodT3Mat$chr)


AllClustersLab <- read.csv(file = "InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv")
matchSNF101 <- match(AllClustersLab$ID[- which(is.na(AllClustersLab$SNFClust10Feb2021) == TRUE)], 
                     colnames(RNAseqProtCodT3Mat))
RNAseqProtCodT3MatName <- RNAseqProtCodT3Mat[, matchSNF101]
rownames(RNAseqProtCodT3MatName) <- RNAseqProtCodT3Mat$name
range(RNAseqProtCodT3MatName) #  0.00 17766.33
dim(RNAseqProtCodT3MatName) # 17190   101

matchSNF101 <- match(AllClustersLab$ID[- which(is.na(AllClustersLab$SNFClust10Feb2021) == TRUE)], 
                     colnames(BetaTSSMatrixT1))
BetaTSSMatrixT1Cancer <- BetaTSSMatrixT1[, matchSNF101]
dim(BetaTSSMatrixT1Cancer) # 114897    101
ClusterLabsSNF <- AllClustersLab$SNFClust10Feb2021[which(! is.na(AllClustersLab$SNFClust10Feb2021))]
length(ClusterLabsSNF) # 101

# check colnames
identical(colnames(RNAseqProtCodT3MatName), colnames(BetaTSSMatrixT1[, matchSNF101])) # TRUE


set.seed(1234)
MethylMixResultsTSSSNFC1vsC2 <- MethylMix::MethylMix(METcancer = BetaTSSMatrixT1Cancer[, which(ClusterLabsSNF == 1)], 
                                               GEcancer = as.matrix(RNAseqProtCodT3MatName)[, which(ClusterLabsSNF == 1)], 
                                               METnormal = BetaTSSMatrixT1Cancer[, which(ClusterLabsSNF == 2)],
                                               NoNormalMode = TRUE)
# saveRDS(MethylMixResultsTSSSNFC1vsC2, file = "41_MethylMixResultsTSSSNFC1vsC2_1March2021.rds")
MethylMixResultsTSSSNFC1vsC2 <- readRDS(file = "41_MethylMixResultsTSSSNFC1vsC2_1March2021.rds")

names(MethylMixResultsTSSSNFC1vsC2) # 6 names
MethylMixResultsTSSSNFC1vsC2$MethylationDrivers  # 249; genes identified as transcriptionally predictive and differentially methylated 
length(MethylMixResultsTSSSNFC1vsC2$MethylationDrivers) # 249
MethylMixResultsTSSSNFC1vsC2$NrComponents # number of methylation states found for each driver gene.
MethylMixResultsTSSSNFC1vsC2$MixtureStates # DM-values for each driver gene (is this  mean??)
length(MethylMixResultsTSSSNFC1vsC2$MixtureStates) # 249

unlist(MethylMixResultsTSSSNFC1vsC2$MixtureStates)[order(abs(unlist(MethylMixResultsTSSSNFC1vsC2$MixtureStates)), decreasing = TRUE)]
length(unlist(MethylMixResultsTSSSNFC1vsC2$MixtureStates)[order(abs(unlist(MethylMixResultsTSSSNFC1vsC2$MixtureStates)), decreasing = TRUE)]) # 377
# look at reverse order
tail(unlist(MethylMixResultsTSSSNFC1vsC2$MixtureStates)[order(unlist(MethylMixResultsTSSSNFC1vsC2$MixtureStates), decreasing = TRUE)])



MethylMixResultsTSSSNFC1vsC2$MethylationStates # DM-values for all driver genes (rows) and all samples (columns)
MethylMixResultsTSSSNFC1vsC2$Classifications # Matrix with integers indicating to which mixture component
# each cancer sample was assigned to, for each gene.
MethylMixResultsTSSSNFC1vsC2$Models # Beta mixture model parameters for each driver gene.
names(MethylMixResultsTSSSNFC1vsC2$Models)
MethylMixResultsTSSSNFC1vsC2$Models$ZCCHC24$llike



table(MethylMixResultsTSSSNFC1vsC2$NrComponents)
# 1   2   3 
# 123 124   2 

# Most DM value (with 1 component) if - =  tumor hypomethylated compared to normal
# Most DM values (with 1 component) if + = tumor hypermethylated compared to normal
# Save entries with only one component and see if tumor hyper or hypo methylated
genesToLoop <- which(MethylMixResultsTSSSNFC1vsC2$NrComponents == 1)
emptySumVectorSNF <- vector()
for(i in 1:length(genesToLoop)) { 
  cat("\n", i)
  if(MethylMixResultsTSSSNFC1vsC2$MixtureStates[[genesToLoop[i]]] < 0) {
    emptySumVectorSNF[i] <- -1
  } else if(MethylMixResultsTSSSNFC1vsC2$MixtureStates[[genesToLoop[i]]] > 0) {
    emptySumVectorSNF[i] <- 1
  } 
}

table(emptySumVectorSNF)
# emptySumVector
# -1  1 
# 81 42
which(emptySumVectorSNF == -1)




# Plot the most famous methylated gene for glioblastoma
GeneName <- "EPHA4"
MethylMixResultsTSSSNF$MixtureStates$GPRC5C

# look at stage
matchClinical101 <- match(colnames(BetaTSSMatrixT1[, matchSNF101]), ClinicalFile_T1$SAMPLE_ID)
table(MethylMixResultsTSSSNF$Classifications[which(rownames(MethylMixResultsTSSSNF$Classifications) == GeneName), ],
      ClinicalFile_T1$STAGE[matchClinical101])
chisq.test(table(MethylMixResultsTSSSNF$Classifications[which(rownames(MethylMixResultsTSSSNF$Classifications) == GeneName), ],
                 ClinicalFile_T1$STAGE[matchClinical101]))
plots <- MethylMix::MethylMix_PlotModel(GeneName = GeneName, 
                                        MixtureModelResults = MethylMixResultsTSSSNF, 
                                        METcancer = BetaTSSMatrixT1[, matchSNF101],
                                        GEcancer = as.matrix(RNAseqProtCodT3MatName),
                                        METnormal = BetaTSSMatrixT1[, c(166:170)])
plots$MixtureModelPlot
plots$CorrelationPlot


# my own plotting
normalBeta <- as.vector(c(BetaTSSMatrixT1[, c(166:170)][which(rownames(BetaTSSMatrixT1[, c(166:170)]) == GeneName), ][1, ], rep(NA, 96)))
tumorBeta <- as.vector(BetaTSSMatrixT1[, matchSNF101][which(rownames(BetaTSSMatrixT1[, matchSNF101]) == GeneName), ][1, ])
# tumorRNAseq <- RNAseqProtCodT3MatNameC1[which(rownames(RNAseqProtCodT3MatNameC1) == GeneName), ]
GeneNameMyPlot <- data.frame(
  normalBeta = normalBeta,
  tumorBeta = tumorBeta
  # tumorRNAseq = RNAseqProtCodT3MatName
)

GeneNameMyPlot %>%
  reshape2::melt() %>%
  ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                    ylab = "Value", 
                    font.label = list(size = 20, color = "black")) +
  ggtitle(paste0(GeneName)) +
  ggpubr::stat_compare_means() + 
  EnvStats::stat_n_text()




ggplot2::ggplot(data = GeneNameMyPlot) + 
  geom_point(mapping = aes(x = normalBeta, y = normalRNAseq)) +
  labs(title = "RNF180",
       x = "normalized methylation beta", 
       y = "normalized RNAseq count",
       colour = "Event") +
  theme_bw() +
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))



all(which((MethylMixResultsTSSSNF$MethylationDrivers %in% rownames(RNAseqNormalizedCountDataTabMatrix)) == TRUE)) # 1025





# selecting entries with no spurious clusters
sizeSpuriousClusterExclude <- 10
saveGeneNamesAll <- NA
cat("\n Methylation drivers that should be kept:")
for (i in 1:nrow(MethylMixResultsTSSSNF$Classifications)) {
  #cat("\n i = ", i, "\n")
  if(length(unique(MethylMixResultsTSSSNF$Classifications[i, ])) > 1) { # if more than one component
    if(all((table(MethylMixResultsTSSSNF$Classifications[i, ]) > sizeSpuriousClusterExclude) == TRUE)) {
      saveGeneNamesAll <- c(saveGeneNamesAll, rownames(MethylMixResultsTSSSNF$Classifications)[i])
    }
  }
}

# Arranging such genes with no suprious clusters by DM values
matchNoSpurious <- match(saveGeneNamesAll[-1], names(MethylMixResultsTSSSNF$MixtureStates))
saveTopResults <- list()
for (i in 1:length(matchNoSpurious)) {
  if(length(MethylMixResultsTSSSNF$MixtureStates[[matchNoSpurious[i]]]) > 1) {
    saveTopResults[[i]] <- max(abs(MethylMixResultsTSSSNF$MixtureStates[[matchNoSpurious[i]]]))
  } else {
    saveTopResults[[i]] <- MethylMixResultsTSSSNF$MixtureStates[[matchNoSpurious[i]]]
  }
}
names(saveTopResults) <- names(MethylMixResultsTSSSNF$MixtureStates)[matchNoSpurious]
# sort values from highest to lowest
sort(unlist(saveTopResults), decreasing = TRUE)


gprofiler101 <- GProfilerAnalysis(GeneIDs = list(MethylMixResultsTSSSNF$MethylationDrivers),
                                  Organism = "hsapiens",
                                  OrderedQuery = TRUE,
                                  PvalAlphaLevel = 0.01,
                                  PositiveorNegFC = NA, # with respect to expression
                                  ConditionName = "101SNFsamples",
                                  ProduceImages = "Yes", 
                                  PNGorPDF = "png")  
gprofiler101$shortLink  
# https://biit.cs.ut.ee/gplink/l/vcH1_SFNQ2
dim(gprofiler101$TableAllValues) # 54


gprofilersameC1C2 <- GProfilerAnalysis(GeneIDs = list(sameC1C2),
                                       Organism = "hsapiens",
                                       OrderedQuery = TRUE,
                                       PvalAlphaLevel = 0.01,
                                       PositiveorNegFC = NA, # with respect to expression
                                       ConditionName = "sameC1C2",
                                       ProduceImages = "Yes", 
                                       PNGorPDF = "png")  





#### Looking at overalp #####
  
  mainlist <- list(Al101 = MethylMixResultsTSSSNF$MethylationDrivers,
                   C1 = MethylMixResultsTSSSNFC1$MethylationDrivers,
                   C2 = MethylMixResultsTSSSNFC2$MethylationDrivers,
                   All121 = MethylMixResultsTSS121$MethylationDrivers)
  
  OutputSNFClusters121 <- VennDiagramAnalysis(MainList = mainlist, 
                                           Labels = c("All101", "C1", "C2", "All121"),
                                           FigureGenerate = "Yes", 
                                           ImageName = "MethylMixResultsSNFclusters", 
                                           PNGorPDF = "png")
  
  NumbComponents <- data.frame("All101" = table(MethylMixResultsTSSSNF$NrComponents), 
                               "C1" = table(MethylMixResultsTSSSNFC2$NrComponents),
                               "C2" = table(MethylMixResultsTSSSNFC2$NrComponents),
                               "All121" = table(MethylMixResultsTSS121$NrComponents))
  #manually adjust
  NumbComponents$C1.Freq <- c(int(table(MethylMixResultsTSSSNFC1$NrComponents)), 0)
  NumbComponents <- melt(NumbComponents)
  
  #1. summary barplot - adapted based on code from Karin I. 
  NumbComponentsPlot <- ggpubr::ggbarplot(NumbComponents, 
                                          x = "variable", y ="value", fill = "All121.Var1",
                                          label = TRUE,   ylab = "Number of driver genes", 
                                          lab.pos = "in",
                                          xlab = "Contrast", 
                                          font.label = list(size = 20, color = "black")) +
    ggtitle("Comparison of All (n = 101, 121), SNF C1 (n = 53) and SNF C2 (n = 48)") +
    rotate_x_text(90)
  
  
  # Differences between n = 101 and n = 121
  diffn101n121 <- setdiff(MethylMixResultsTSS121$MethylationDrivers, 
                          MethylMixResultsTSSSNF$MethylationDrivers)
  length(diffn101n121) # 123
  
  # Same between n = 101 and n = 121
  samen101n121 <- intersect(MethylMixResultsTSSSNF$MethylationDrivers, 
                        MethylMixResultsTSS121$MethylationDrivers)
  length(samen101n121) # 247

  
  gprofilersamen101n121 <- GProfilerAnalysis(GeneIDs = list(MethylMixResultsTSS121$MethylationDrivers),
                                         Organism = "hsapiens",
                                         OrderedQuery = TRUE,
                                         PvalAlphaLevel = 0.01,
                                         PositiveorNegFC = NA, # with respect to expression
                                         ConditionName = "n101n121",
                                         ProduceImages = "Yes", 
                                         PNGorPDF = "png")
  
  gprofilersamen101n121$shortLink
  # "https://biit.cs.ut.ee/gplink/l/PBUmBgolSF"
  
  
  
  gprofilerdiff123 <- GProfilerAnalysis(GeneIDs = list(diffn101n121),
                                             Organism = "hsapiens",
                                             OrderedQuery = TRUE,
                                             PvalAlphaLevel = 0.01,
                                             PositiveorNegFC = NA, # with respect to expression
                                             ConditionName = "diff123",
                                             ProduceImages = "Yes", 
                                             PNGorPDF = "png")
  gprofilerdiff123$shortLink
  # "https://biit.cs.ut.ee/gplink/l/QJeKeePKT7"
  dim(gprofilerdiff123$TableAllValues) # 65
  
#### Save ####
# save.image("41_MethylMixl_RNAseqVsMethFC_script.RData")
# [END]
