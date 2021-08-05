# Updated 3 Aug 2021
# 13 August 2020
# This is a script, not a funciton
# Using files from reset tools on http://ciriellolab.org/reset/reset.html, carry out analysis RNAseq vs Methylation 
# data falling in TSS (transcription start site).
# Author: Anjali Silva


#### Source files from RESET ####
# Resource to detect Epigenetically Silenced and Enhanced Targets
DirectoryResetScripts <- "/OtherScripts/reset/R"
source(paste0(getwd(), paste0(DirectoryResetScripts,"/eventScore.R")))
source(paste0(getwd(), paste0(DirectoryResetScripts,"/methStatus.R")))
source(paste0(getwd(), paste0(DirectoryResetScripts,"/FDRcal.R")))
source(paste0(getwd(), paste0(DirectoryResetScripts,"/reset.R")))
source(paste0(getwd(), paste0(DirectoryResetScripts,"/shared_functions.R")))
source(paste0(getwd(), paste0(DirectoryResetScripts,"/methNorSel.R")))

load(paste0(getwd(), "/OtherScripts/reset/R/promoter-probes-list.rdata"))
# View(promoter.probes.list)
head(promoter.probes.list)
names(promoter.probes.list)
# [1] "ProbID"             "index"              "Gene_Name"          "Coordinate_Pos"    
# [5] "Gene_Fantom"        "Chromosome"         "Promoter_Pos_start" "Promoter_Pos_End"  
# [9] "Promoter_No"        "Entrezgene_ID"      "Hgnc_ID"            "Uniprot_ID"  
dim(promoter.probes.list) # 221948     12

#### select TSS probes ####

AnnotationFileEdited <- EPICAnnotationFile
matchBetaProbes <- match(rownames(BetaMatrix_T1), AnnotationFileEdited$V1)
AnnotationFileEditedBetaMatrix <- AnnotationFileEdited[matchBetaProbes, ]
dim(AnnotationFileEditedBetaMatrix) # 595564     47
AnnotationFileEditedBetaMatrix$UCSC_RefGene_Group <- sub("\\;.*", "", 
                                                         AnnotationFileEditedBetaMatrix$UCSC_RefGene_Group)
promoterProbes <- AnnotationFileEditedBetaMatrix %>% 
  filter(substr(UCSC_RefGene_Group, 1, 3) == "TSS") %>% pull("V1")
matchTSSProbes <- match(promoterProbes, rownames(BetaMatrix_T1))
length(matchTSSProbes) # 114897


#### select RNAseq samples (T3) ####
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

RNAseqT3ProteinCode <- 
  RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations[which
  (RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations$type == "protein_coding"), ]
dim(RNAseqT3ProteinCode) # 17190   135
rownames(RNAseqT3ProteinCode) <- RNAseqT3ProteinCode$name
RNAseqT3ProteinCode <- RNAseqT3ProteinCode[, c(4:134)]
dim(RNAseqT3ProteinCode) # 17190   131
# RNAseq counts for T3 (131 patients )





#### select RNAseq samples (T2) ####
ClinicalFile_T1Cluster <- data.frame(ClinicalFile_T1, "CLUSTER" = InfiniumClustLabels$Cluster)
RNAseqQC18June2020T2Samples104 <- QCRNAseq(RNAseqCountMatrix = RNAseqCountMatrixMatched, 
                                           QCMatrix = RNAseqQCFile, 
                                           RNAseqAnnotationFile = ENSMBLid, 
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
                                           ProduceImages = "No")
RNAseqT2ProteinCode <- 
  RNAseqQC18June2020T2Samples104$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations[which
  (RNAseqQC18June2020T2Samples104$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations$type == "protein_coding"), ]
dim(RNAseqT2ProteinCode) # 15080   107
rownames(RNAseqT2ProteinCode) <- RNAseqT2ProteinCode$name
RNAseqT2ProteinCode <- RNAseqT2ProteinCode[, c(4:107)]
dim(RNAseqT2ProteinCode) # 15080   104
# RNAseq counts for T2 (104 patients )


RNAseqT2ProteinCodeNormalized <- 
  RNAseqQC18June2020T2Samples104$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized[which
  (RNAseqQC18June2020T2Samples104$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations$type == "protein_coding"), ]



#### select RNAseq samples (T1) ####
ClinicalFile_T1Cluster <- data.frame(ClinicalFile_T1, "CLUSTER" = InfiniumClustLabels$Cluster)
RNAseqQC18June2020T1Samples81 <- QCRNAseq(RNAseqCountMatrix = RNAseqCountMatrixMatched, 
                                          QCMatrix = RNAseqQCFile, 
                                          RNAseqAnnotationFile = ENSMBLid, 
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
RNAseqT1ProteinCode <- 
  RNAseqQC18June2020T1Samples81$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations[which
  (RNAseqQC18June2020T1Samples81$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations$type == "protein_coding"), ]
dim(RNAseqT1ProteinCode) # 15080   84
rownames(RNAseqT1ProteinCode) <- RNAseqT1ProteinCode$name
RNAseqT1ProteinCode <- RNAseqT1ProteinCode[, c(4:84)]
dim(RNAseqT1ProteinCode) #  14914    81





#### Run first function from RESET with RLN samples as normal ####
set.seed(1234)
methNorSelOutput <- methNorSel(normal.mtx = BetaMatrix_T1[matchTSSProbes, c(166:170)],
                               probe.list = promoter.probes.list)
# this function first extract the data regarding the selected probes, and 
# later based on the beta-values prepare two seperate normal-sample data-sets:
## 1. hypermethylated probes in normal condition: Enh.probes
## 2. hypomethylated probes in normal condition: Sil.probes
names(methNorSelOutput) # "normal.sil.probes" "normal.enh.probes"
typeof(methNorSelOutput)
dim(methNorSelOutput$normal.sil.probes) # 45229     5
dim(methNorSelOutput$normal.enh.probes) # 2798    5
class(methNorSelOutput$normal.sil.probes)
normal.sil.probes <- data.frame(methNorSelOutput$normal.sil.probes, 
                                rowMedians = rowMedians(methNorSelOutput$normal.sil.probes))

ggplot(normal.sil.probes, aes(x = c(1:45229), y = rowMedians)) + # plot rowMeans of sil.probes
  geom_point(size = 2) +
  labs(x = "probe", y = "Median Beta Value") +
  theme_bw() + 
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  scale_color_manual(values = c("#d6604d", "#66bd63", "#4575b4")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))
silProbesFL <- gsub(".*@", "", rownames(methNorSelOutput$normal.sil.probes))
length(silProbesFL) # 45229
# write.csv(silProbesFL, file = "RESET_HypermethylatedProbesNormalCondition_Enhancers.csv")
source("26_Gprofiler.R")
gprofilersilProbesFL <- GProfilerAnalysis(GeneIDs = list(silProbesFL),
                                          Organism = "hsapiens",
                                          OrderedQuery = TRUE,
                                          PvalAlphaLevel = 0.01,
                                          PositiveorNegFC = NA, # with respect to expression
                                          ConditionName = "All",
                                          ProduceImages = "Yes", 
                                          PNGorPDF = "png")
gprofilersilProbesFL$shortLink  # "https://biit.cs.ut.ee/gplink/l/ZLxq3tvKRC"



enhProbesFL <- gsub(".*@", "", rownames(methNorSelOutput$normal.enh.probes))
length(enhProbesFL) # 2798
# write.csv(enhProbesFL, file = "RESET_HypomethylatedProbesNormalCondition_Silencers.csv")
normal.enh.probes <- data.frame(methNorSelOutput$normal.enh.probes, 
                                rowMedians = rowMedians(methNorSelOutput$normal.enh.probes))
ggplot(normal.enh.probes, aes(x = c(1:2798), y = rowMedians)) + # plot rowMeans of enh.probes
  geom_point(size = 2) +
  labs(x = "probe", y = "Median Beta Value") +
  theme_bw() + 
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  scale_color_manual(values = c("#d6604d", "#66bd63", "#4575b4")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))
silProbesFL <- gsub(".*@", "", rownames(methNorSelOutput$normal.sil.probes))

gprofilerenhProbesFL <- GProfilerAnalysis(GeneIDs = list(enhProbesFL),
                                          Organism = "hsapiens",
                                          OrderedQuery = TRUE,
                                          PvalAlphaLevel = 0.01,
                                          PositiveorNegFC = NA, # with respect to expression
                                          ConditionName = "All",
                                          ProduceImages = "Yes", 
                                          PNGorPDF = "png")
gprofilerenhProbesFL$shortLink # "https://biit.cs.ut.ee/gplink/l/Aw1-UW-RS_"




#### RESET_HypermethylatedProbes_FL+DLBCL_(n = 131)_Enhancers ####

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
dim(resetResultsENH$score.cutoff)# NULL

# Sadegh Saghafinia <Sadegh.saghafinia@epfl.ch>
# In general score above 1.5 is what you can trust. So if the FDR is significant enough, 
# then maybe consider 1.5 as a good cut.off.


resetScoreENH <- data.frame(resetResultsENH$Score.report[which
                 (! is.na(resetResultsENH$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENH)){
  resetScoreENH[i, 6] <- resetResultsENH$FDR.res$FDR[which(resetScoreENH$Score[i] == 
                         resetResultsENH$FDR.res$observed.scores)]
}
head(resetScoreENH)
dim(resetScoreENH) # 585   6
length(unique(resetScoreENH$Gene)) # 476
range(resetScoreENH$Score) # 0.02586767 1.29846829
View(resetScoreENH)
write.csv(resetScoreENH, file = "RESET_HypermethylatedProbes_FLDLBCL_Enhancers.csv")

# Plot No.Methylation.Events vs Score
resetScoreENH %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()

# Plot No.Methylation.Events vs Score
resetScoreENH %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


# Plot top 100 Genes vs Score after selecting by FDR
resetScoreENH %>%  
  filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes (Enhancers) For FL+DLBCL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


#### RESET_HypomethylatedProbes_FL+DLBCL_(n = 131)_Silencers ####
# Run reset function on all 131 tumor samples - silencers

# Match samples between RNAseq and methylation - 131 tumor samples
matchRNAseqBeta <- 
  match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized),
  colnames(BetaMatrix_T1))

resetResultsSIL <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                         meth.tumor = BetaMatrix_T1[, matchRNAseqBeta[1:131]],
                         transcriptome = RNAseqT3ProteinCode,
                         methylation.event = c('sil'),
                         FDR.permutation.no = 100, 
                         seed = 100)


names(resetResultsSIL)
dim(resetResultsSIL$normal.meth)  # 45229    5
dim(resetResultsSIL$meth.tumor.all) # 45229  131
dim(resetResultsSIL$meth.tumor.status.all) # 45229  131
dim(resetResultsSIL$beta.dist.normal.specific) #  45229    4
dim(resetResultsSIL$beta.dist.normal.universal) # NULL
dim(resetResultsSIL$Score.report) #  45229    5
dim(resetResultsSIL$transcriptome) # 45229  131
dim(resetResultsSIL$transcriptome.norm.dis) # 45229  131
dim(resetResultsSIL$meth.tumor.matched) # 45229  131
dim(resetResultsSIL$meth.tumor.status.matched) #  45229  131
dim(resetResultsSIL$FDR.res) #  6611   2
dim(resetResultsSIL$permutation.res) #  45229  102
dim(resetResultsSIL$score.cutoff) # NULL


resetScoreSIL <- data.frame(resetResultsSIL$Score.report[which
                 (! is.na(resetResultsSIL$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSIL)) {
  resetScoreSIL[i, 6] <- resetResultsSIL$FDR.res$FDR[which(resetScoreSIL$Score[i] == 
                                           resetResultsSIL$FDR.res$observed.scores)]
}
head(resetScoreSIL)
dim(resetScoreSIL) # 6611    6
length(unique(resetScoreSIL$Gene)) # 3925
range(resetScoreSIL$Score) # 0.002308231 1.269213329
# View(resetScoreSIL)
# write.csv(resetScoreSIL, file = "RESET_HypomethylatedProbes_FLDLBCL_Silencers.csv")

# Plot No.Methylation.Events vs Score
resetScoreSIL %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()

# Plot No.Methylation.Events vs Score
resetScoreSIL %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score after selecting based on FDR
resetScoreSIL %>%  
  filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes (Silencers) For FL+DLBCL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


# Plot top 100 Genes vs Score after selecting based on FDR - coloured by FDR
resetScoreSIL %>%  
  filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes (Silencers) For FL+DLBCL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))

# Plot top 100 Genes vs Score after selecting based on FDR - coloured by No.Methylation.Events
resetScoreSIL %>%  
  filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = No.Methylation.Events)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes (Silencers) For FL+DLBCL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


# Plot top 100 Genes vs Score after selecting based on FDR - coloured by Promoter.no
resetScoreSIL %>%  
  filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = factor(Promoter.no))) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes (Silencers) For FL+DLBCL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


###########################################################################################
#### RESET_HypomethylatedProbes_FL_(T3 n = 131-10 = 121)_HypermethylatedNormal(Enhancers) ####

# Run reset function on all tumor samples - enhancers - FL only
resetResultsENHFLonly <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                               meth.tumor = BetaMatrix_T1[, matchRNAseqBeta[11:131]],
                               transcriptome = RNAseqT3ProteinCode[, c(11:131)],
                               methylation.event = c('enh'),
                               FDR.permutation.no = 100, 
                               seed = 100)


dim(resetResultsENHFLonly$FDR.res) #  588   2
head(resetResultsENHFLonly$meth.tumor.status.all)
ENHFLmeth.tumor.status.all <- data.frame(resetResultsENHFLonly$meth.tumor.status.all, 
                                         Gene = gsub(".*@", "", 
                                         rownames(resetResultsENHFLonly$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.all)
dim(ENHFLmeth.tumor.status.all) # 2798  122
dim(resetResultsENHFLonly$meth.tumor.all) # 2798  121
dim(resetResultsENHFLonly$transcriptome) # 2798  121

genelist <- c("SPO11", "FCRLB")
testingGene <- genelist[1]
ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.all$Gene == testingGene), ]
rowSums(ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.all$Gene == testingGene), c(1:121)])



plotting <- function() {
  # If one probe
  if(length(which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)) == 1) {
    
    
    testingDataFrame <- data.frame(Values = 
                                  c(resetResultsENHFLonly$normal.meth[which
                                  (gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ],
                                  resetResultsENHFLonly$meth.tumor.all[which
                                  (gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), ],
                                  resetResultsENHFLonly$transcriptome[which
                                  (gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), ]),
                                  Type = c(rep("RLNMethylation", 5), 
                                           rep("TumorMethylation", 121), 
                                           rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsENHFLonly$normal.meth[which
                                             (gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsENHFLonly$meth.tumor.all[which
                                             (gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsENHFLonly$transcriptome[which
                                             (gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    # only looking at samples with a status marked as 1
    hypoProbes <- which(ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.all$Gene == testingGene), ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLonly$normal.meth[which
                                                    (gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ],
                                                    resetResultsENHFLonly$meth.tumor.all[which
                                                    (gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), hypoProbes],
                                                    resetResultsENHFLonly$transcriptome[which
                                                    (gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), 
                                                  rep("TumorMethylation", length(hypoProbes)), 
                                                  rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) # 247 2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLonly$normal.meth[which
                                                   (gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ]),
                                      MeanMethTumor = mean(resetResultsENHFLonly$meth.tumor.all[which
                                                   (gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLonly$transcriptome[which
                                                   (gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsENHFLonly$normal.meth[which
                                          (gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ]),
                                         t(resetResultsENHFLonly$meth.tumor.all[which
                                          (gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsENHFLonly$transcriptome[which
                                          (gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 250   4
    
    
    testingDataFrame %>%  # highest score: cg20802515.p1.SPO11 methylation only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsENHFLonly$normal.meth[which
                                (gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsENHFLonly$meth.tumor.all[which
                                (gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsENHFLonly$transcriptome[which
                                (gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    # look at 2nd probe with highest score
    hypoProbes <- which(ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.all$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLonly$normal.meth[which
                                                  (gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)[2], ],
                                                  resetResultsENHFLonly$meth.tumor.all[which
                                                  (gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene)[2], hypoProbes],
                                                  resetResultsENHFLonly$transcriptome[which
                                                  (gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene)[2], hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), 
                                                  rep("TumorMethylation", length(hypoProbes)), 
                                                  rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) # 25  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLonly$normal.meth[which
                                      (gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsENHFLonly$meth.tumor.all[which
                                      (gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene)[2], hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLonly$transcriptome[which
                                      (gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene)[2], hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  }
}



resetScoreENHFLonly <- data.frame(resetResultsENHFLonly$Score.report[which
                                  (! is.na(resetResultsENHFLonly$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFLonly)) {
  resetScoreENHFLonly[i, 6] <- resetResultsENHFLonly$FDR.res$FDR[which
                               (resetScoreENHFLonly$Score[i] == resetResultsENHFLonly$FDR.res$observed.scores)]
}
head(resetScoreENHFLonly)
dim(resetScoreENHFLonly) # 588 6
length(unique(resetScoreENHFLonly$Gene)) # 489
range(resetScoreENHFLonly$Score) # 0.004389081 1.150200828
range(resetScoreENHFLonly$FDR) # 0.0540000 0.9513823
resetScoreENHFLonly[which(resetScoreENHFLonly$Gene == testingGene), ]

resetScoreENHFLonly[order(resetScoreENHFLonly$Score, decreasing = TRUE),]


# View(resetScoreENHFLonly)
# write.csv(resetScoreENHFLonly, file = "RESET_HypermethylatedProbes_FL_Enhancers.csv")

# Plot No.Methylation.Events vs Score
resetScoreENHFLonly %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()

# Plot No.Methylation.Events vs Score
resetScoreENHFLonly %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreENHFLonly %>%  
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreENHFLonly %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))





# plotting RNAseq vs methylation value - SPO11, FCRLB
# get the probe IDs
genelist <- c("SPO11", "FCRLB")
testingGene <- genelist[2]
probesToPlot <- resetScoreENHFLonly[which(resetScoreENHFLonly$Gene == testingGene), ]$ProbID
which(rownames(BetaMatrix_T1[, matchRNAseqBeta[11:131]]) == tolower(probesToPlot))
# 35143 54921 78854

# Get RNAseq counts
which(rownames(RNAseqT3ProteinCode) == testingGene) 
# 679

# Get normalized RNAseq values
RNAseqT3ProteinCodeNormalized <- 
  RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations[which
  (RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations$type == "protein_coding"), ]
dim(RNAseqT3ProteinCodeNormalized) # 17190   136
rownames(RNAseqT3ProteinCodeNormalized) <- RNAseqT3ProteinCodeNormalized$name
RNAseqT3ProteinCodeNormalized <- RNAseqT3ProteinCodeNormalized[, c(15:135)]
dim(RNAseqT3ProteinCodeNormalized) # 17190   121
range(RNAseqT3ProteinCodeNormalized) # -2.72230 14.11687

# getting the events
which(ENHFLmeth.tumor.status.all$Gene == testingGene)
# 937 939 941

GeneToPlot <- data.frame(Beta = (BetaMatrix_T1[which(rownames
                         (BetaMatrix_T1[, matchRNAseqBeta[11:131]]) == tolower(probesToPlot)), matchRNAseqBeta[11:131]]),
                         Event = t(ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.all$Gene == testingGene), c(1:121)]),
                         RNAseqCounts = t(RNAseqT3ProteinCode[which(rownames(RNAseqT3ProteinCode) == testingGene), c(11:131)]),
                         RNAseqNormalized = t(RNAseqT3ProteinCodeNormalized[which(rownames(RNAseqT3ProteinCode) == testingGene), ]))

ggplot2::ggplot(GeneToPlot, aes(x = Beta, y = FCRLB.1, 
                        color = factor(cg23641748.p17.FCRLB))) +
                geom_point(size = 2)+
                labs(title = "cg23641748 - FCRLB",
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







#### RESET_HypomethylatedProbes_FL_(T3 n = 131-10 = 121)_HypermethylatedNormal(Enhancers) with  FDR.permutation.no = 1000 ####
## Increase FDR.permutation.no to 1000 from 100

# Run reset function on all tumor samples - enhancers - FL only
resetResultsENHFLonly1000 <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                                   meth.tumor = BetaMatrix_T1[, matchRNAseqBeta[11:131]],
                                   transcriptome = RNAseqT3ProteinCode[, c(11:131)],
                                   methylation.event = c('enh'),
                                   FDR.permutation.no = 1000, 
                                   seed = 100)

dim(resetResultsENHFLonly1000$FDR.res) #  588   2
head(resetResultsENHFLonly1000$meth.tumor.status.all)
ENHFLmeth.tumor.status.all <- data.frame(resetResultsENHFLonly1000$meth.tumor.status.all, 
                                         Gene = gsub(".*@", "", 
                                         rownames(resetResultsENHFLonly1000$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.all)
dim(ENHFLmeth.tumor.status.all) # 2798  122
dim(resetResultsENHFLonly$meth.tumor.all) # 2798  121
dim(resetResultsENHFLonly$transcriptome) # 2798  121


resetScoreENHFLonly1000 <- data.frame(resetResultsENHFLonly1000$Score.report[which
                           (! is.na(resetResultsENHFLonly1000$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFLonly1000)) {
  resetScoreENHFLonly1000[i, 6] <- resetResultsENHFLonly1000$FDR.res$FDR[which
                                   (resetScoreENHFLonly1000$Score[i] ==
                                   resetResultsENHFLonly1000$FDR.res$observed.scores)]
}
head(resetScoreENHFLonly1000)
dim(resetScoreENHFLonly1000) # 588 6
length(unique(resetScoreENHFLonly1000$Gene)) # 489
range(resetScoreENHFLonly1000$Score) # 0.004389081 1.150200828
genelist <- c("SPO11", "FCRLB")
testingGene <- genelist[1]
resetScoreENHFLonly1000[which(resetScoreENHFLonly1000$Gene == testingGene), ]

#### RESET_HypermethylatedProbesFL_(T3 n = 131-10 = 121)_HypomethylatedNormal(Silencers) ####

# Match samples between RNAseq and methylation - 131 tumor samples
matchRNAseqBeta <-
  match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized),
  colnames(BetaMatrix_T1))

# Run reset function on all tumor samples - silencers - FL only
resetResultsSILFLonly <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                               meth.tumor = BetaMatrix_T1[, matchRNAseqBeta[11:131]],
                               transcriptome = RNAseqT3ProteinCode[, c(11:131)],
                               methylation.event = c('sil'),
                               FDR.permutation.no = 100, 
                               seed = 100)
dim(resetResultsSILFLonly$FDR.res) #  6452   2
SILFLmeth.tumor.status.all <- data.frame(resetResultsSILFLonly$meth.tumor.status.all, 
                                         Gene = gsub(".*@", "", 
                                         rownames(resetResultsSILFLonly$meth.tumor.status.all)))
head(SILFLmeth.tumor.status.all)
dim(SILFLmeth.tumor.status.all) # 45229   122
dim(resetResultsSILFLonly$meth.tumor.all) # 45229  121
dim(resetResultsSILFLonly$transcriptome) # 45229  121

genelist <- c("BMP3", "MME", "HLTF", "CAMK1")
testingGene <- genelist[1]
SILFLmeth.tumor.status.all[which(SILFLmeth.tumor.status.all$Gene == testingGene), ]
rowSums(SILFLmeth.tumor.status.all[which(SILFLmeth.tumor.status.all$Gene == testingGene), c(1:121)])




plotting <- function() {
  # If one probe
  if(length(which(gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene)) == 1) {
    testingDataFrame <- data.frame(Values = c(resetResultsSILFLonly$normal.meth[which
                                              (gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene), ],
                                              resetResultsSILFLonly$meth.tumor.all[which
                                              (gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene), ],
                                              resetResultsSILFLonly$transcriptome[which
                                              (gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene), ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsSILFLonly$normal.meth[which
                                (gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsSILFLonly$meth.tumor.all[which
                                (gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsSILFLonly$transcriptome[which
                                (gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    hyperProbes <- which(SILFLmeth.tumor.status.all[which(SILFLmeth.tumor.status.all$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsSILFLonly$normal.meth[which
                                                  (gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene)[2], ],
                                                  resetResultsSILFLonly$meth.tumor.all[which
                                                  (gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene)[2], hyperProbes],
                                                  resetResultsSILFLonly$transcriptome[which
                                                  (gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene)[2], hyperProbes]),
                                         Type = c(rep("RLNMethylation", 5), 
                                                  rep("TumorMethylation", length(hyperProbes)), 
                                                  rep("TumorRNAseq", length(hyperProbes))))
    dim(testingDataFrameStatus) # 33  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsSILFLonly$normal.meth[which
                                      (gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsSILFLonly$meth.tumor.all[which
                                      (gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene)[2], hyperProbes]),
                                      MeanTranscriptome = mean(resetResultsSILFLonly$transcriptome[which
                                      (gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene)[2], hyperProbes]))
    
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsSILFLonly$normal.meth[which
                                          (gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene), ]),
                                         t(resetResultsSILFLonly$meth.tumor.all[which
                                          (gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsSILFLonly$transcriptome[which
                                          (gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), 
                                            rep("TumorMethylation", 121), 
                                            rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247   4
    
    
    testingDataFrame %>%  # highest score: cg01941671@p1@BMP3 methylation only
      melt(id.vars = "Type")  %>%
      filter(variable == "cg01941671.p1.BMP3") %>% # highest score: cg01941671@p1@BMP3
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg01941671.p1.BMP3") %>% # highest score: cg01941671.p1.BMP3
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsSILFLonly$normal.meth[which
                                (gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsSILFLonly$meth.tumor.all[which
                                (gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsSILFLonly$transcriptome[which
                                (gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
  }
}

resetScoreSILFLonly <- data.frame(resetResultsSILFLonly$Score.report[which
                       (! is.na(resetResultsSILFLonly$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFLonly)) {
  resetScoreSILFLonly[i, 6] <- resetResultsSILFLonly$FDR.res$FDR[which
                               (resetScoreSILFLonly$Score[i] == resetResultsSILFLonly$FDR.res$observed.scores)]
}
head(resetScoreSILFLonly)
dim(resetScoreSILFLonly) #  6452  6
length(unique(resetScoreSILFLonly$Gene)) # 3848
range(resetScoreSILFLonly$Score) # 0.003460362 1.282929804
resetScoreSILFLonly[which(resetScoreSILFLonly$Gene == testingGene), ]

# View(resetScoreSILFLonly)
# write.csv(resetScoreSILFLonly, file = "RESET_HypomethylatedProbes_FL_Silencers.csv")

# Plot No.Methylation.Events vs Score
resetScoreSILFLonly %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()


# Plot No.Methylation.Events vs Score
resetScoreSILFLonly %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreSILFLonly %>% 
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreSILFLonly %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# plotting RNAseq vs methylation value - BMP3, FCRLB
# get the probe IDs
genelist <- c("BMP3", "MME")
testingGene <- genelist[2]
probesToPlot <- resetScoreSILFLonly[which(resetScoreSILFLonly$Gene == testingGene), ]$ProbID
which(rownames(BetaMatrix_T1[, matchRNAseqBeta[11:131]]) == tolower(probesToPlot)[8])
# 218098 125779 1060 97335 86432 63883 32900 99749

# Get RNAseq counts
which(rownames(RNAseqT3ProteinCode) == testingGene) 
# 679

# Get normalized RNAseq values
RNAseqT3ProteinCodeNormalized <- 
  RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations[which(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations$type == "protein_coding"), ]
dim(RNAseqT3ProteinCodeNormalized) # 17190   136
rownames(RNAseqT3ProteinCodeNormalized) <- RNAseqT3ProteinCodeNormalized$name
RNAseqT3ProteinCodeNormalized <- RNAseqT3ProteinCodeNormalized[, c(15:135)]
dim(RNAseqT3ProteinCodeNormalized) # 17190   121
range(RNAseqT3ProteinCodeNormalized) # -2.72230 14.11687

# getting the events
which(SILFLmeth.tumor.status.all$Gene == testingGene)
# 937 939 941

GeneToPlot <- data.frame(Beta = t(BetaMatrix_T1[c(218098, 125779, 1060, 97335, 86432, 63883, 32900, 99749), 
                                matchRNAseqBeta[11:131]]),
                         Event = t(SILFLmeth.tumor.status.all[which
                                 (SILFLmeth.tumor.status.all$Gene == testingGene), c(1:121)]),
                         RNAseqCounts = t(RNAseqT3ProteinCode[which
                                        (rownames(RNAseqT3ProteinCode) == testingGene), c(11:131)]),
                         RNAseqNormalized = t(RNAseqT3ProteinCodeNormalized[which
                                            (rownames(RNAseqT3ProteinCode) == testingGene), ]))

ggplot2::ggplot(GeneToPlot, aes(x = Beta.cg07862554, y = MME.1, 
                                color = factor(Event.cg07862554.p8.MME))) +
  geom_point(size = 2)+
  labs(title = "cg07862554 - MME",
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



#### SNF: RESET_HypomethylatedProbes_FL_(T3 n = 131-10 = 121)_HypermethylatedNormal(Enhancers) ####

# bind tumor + normal samples
identical(rownames(ENHFLmeth.tumor.status.all), 
          rownames(methNorSelOutput$normal.enh.probes)) # TRUE

dim(methNorSelOutput$normal.enh.probes) # 2798   5

# bind together tumor and normal for ENH (in normal)
EventMatrix <- ENHFLmeth.tumor.status.all[, c(1:121)]
dim(EventMatrix) # 2798  121
# rownames(EventMatrix) <- gsub("@.*", "", rownames(EventMatrix)) # error




TumorPurity <- readRDS(file = paste0(MethylationDirPath, "/Purity_281probes_10Jan2020.rds"))
dim(TumorPurity) # 170   1

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

# cluster only RNAseq protein coding genes
RNAseqQC18June2020T3Samples132Protein <- 
  RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations[which
  (RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations$type == "protein_coding"), ]
dim(RNAseqQC18June2020T3Samples132Protein) # 17190   135
RNAseqQC18June2020T3Samples132Protein <- as.matrix(RNAseqQC18June2020T3Samples132Protein[, c(4:135)])
dim(RNAseqQC18June2020T3Samples132Protein) # 17190   132


TargetedDNAseqDirPath <- "~/Desktop/UHN/FLOMICS/TargetedSequencing"
mutMergedT1Robert <- 
  read.csv(paste0(TargetedDNAseqDirPath, "/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl_25Aug2020.csv"), 
           row.names = 1)
dim(mutMergedT1Robert) # 55 138
colnames(mutMergedT1Robert)
class(mutMergedT1Robert)

# Fixed in mut.merged.df.T1.poor.cov.excl_25Aug2020 file 
# mutMergedT1Robert[, which(colnames(mutMergedT1Robert) == "LY_FL_179_T1")] # NA values
# mutMergedT1Robert[, which(colnames(mutMergedT1Robert) == "LY_FL_181_T1")] # NA values
# which(is.na(colSums(mutMergedT1Robert)) == TRUE) # LY_FL_179_T1 50; LY_FL_181_T1  51 
# mutMergedT1Robert[, c(50, 51)] <- 0
which(is.na(colSums(mutMergedT1Robert)) == TRUE) # named integer(0)




matchRNAseqMutation <- match(colnames(RNAseqQC18June2020T3Samples132Protein), colnames(mutMergedT1Robert))
dim(mutMergedT1Robert[, matchRNAseqMutation[! is.na(matchRNAseqMutation)]]) # 55 101

mutMergedT1Robert101 <- mutMergedT1Robert[, matchRNAseqMutation[! is.na(matchRNAseqMutation)]]
dim(mutMergedT1Robert101) # 55 101

matchDNAseqBeta <- match(colnames(mutMergedT1Robert101), colnames(EventMatrix))

matchDNAseqRNAseq <- match(colnames(mutMergedT1Robert101), colnames(RNAseqQC18June2020T3Samples132Protein))

matchDNAseqClinicalFile <- match(colnames(mutMergedT1Robert101), ClinicalFile_T1$SAMPLE_ID)

matchDNAseqTumorPurity <- match(colnames(mutMergedT1Robert101), rownames(TumorPurity))
identical(rownames(TumorPurity)[matchDNAseqTumorPurity], colnames(mutMergedT1Robert101))


# Combine mutations with BCL2 and BCL6 translocations
BCL2BCL6BreakapartPredictions <- data.table::fread(file = paste0(MethylationDirPath, 
                                                                 "/2020-07-05_BCL2_BCL6_rearrangements.txt"))
BCL2BCL6BreakapartPredictions <- data.frame(BCL2BCL6BreakapartPredictions)
dim(BCL2BCL6BreakapartPredictions) # 183   3
rownames(BCL2BCL6BreakapartPredictions) <- BCL2BCL6BreakapartPredictions$External_ID
BCL2BCL6BreakapartPredictions <- BCL2BCL6BreakapartPredictions[, -1]
matchMutationRearrangements <- match(colnames(mutMergedT1Robert101), rownames(BCL2BCL6BreakapartPredictions))
rearrangementsT1Robert108 <- BCL2BCL6BreakapartPredictions[matchMutationRearrangements, ]
identical(colnames(t(rearrangementsT1Robert108)), colnames(mutMergedT1Robert101)) #TRUE
mutRearrT1Robert101 <- rbind(mutMergedT1Robert101, t(rearrangementsT1Robert108))
dim(mutRearrT1Robert101) # 57 101

which(is.na(colSums(mutRearrT1Robert101)) == TRUE)
#LY_FL_244_T1 LY_FL_253_T1 LY_FL_264_T1 
#62           71           82
mutRearrT1Robert101[56:57, which(is.na(colSums(mutRearrT1Robert101)) == TRUE)]
# manual adjustment
mutRearrT1Robert101[56:57, 62] <- c(0, 0)
mutRearrT1Robert101[56:57, 71] <- c(1, 0)
mutRearrT1Robert101[56:57, 82] <- c(0, 0)

# check
mutRearrT1Robert101[56:57, c(62, 71, 82)] 



ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all <- 
                                  SNFClustering(RNAseqCountMatrix = as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq]), 
                                                MvalueMatrix = EventMatrix[, matchDNAseqBeta],
                                                BetaMatrix = EventMatrix[, matchDNAseqBeta], 
                                                TargetedSeqMatrix = mutRearrT1Robert101,
                                                MethylationAnnotationFile = EPICAnnotationFile,
                                                RNAseqAnnotationFile = ENSMBLid, 
                                                QCMatrix = RNAseqQC18June2020T3Samples132$QCMatrixMatchedSampleFiltered[matchDNAseqRNAseq, ],
                                                TumorPurity = as.matrix(TumorPurity[matchDNAseqTumorPurity, ], ncol = 1),
                                                ClinicalFile = ClinicalFile_T1[matchDNAseqClinicalFile, ], 
                                                SurvivalFile = SurvivalFile,
                                                NumberOfClusters = 2,
                                                ImageName = paste0("ENHFLmethTumorStatusAll", date()),
                                                PNGorPDF = "png",
                                                ProduceImages = "Yes") 



ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all <- ClusterLabels

table(ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all)

InfiniumClustLabels <- readRDS(file = paste0(MethylationDirPath, "/InfiniumClustLabels2.RDS"))
InfSNFMethEvents <- match(colnames(EventMatrix[, matchDNAseqBeta]), InfiniumClustLabels$ID)
identical(colnames(EventMatrix[, matchDNAseqBeta]), unfactor(InfiniumClustLabels$ID[InfSNFMethEvents]))
# TRUE
table(ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all, InfiniumClustLabels$Cluster[InfSNFMethEvents])
chisq.test(table(ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all, InfiniumClustLabels$Cluster[InfSNF5000SD]))


matchDNAseqBetaMat_T1 <- match(colnames(as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq])), colnames(BetaMatrix_T1))
ClinicalFile_T1SNFCluster2 <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all)
DifferentialMethylationProbesSNF <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1SNFCluster2, 
                                                            MvalueMatrix = MvalueMatrix_T1[, matchDNAseqBetaMat_T1], 
                                                            BetaMatrix = BetaMatrix_T1[, matchDNAseqBetaMat_T1],
                                                            ProbeOrGene = "Probe", 
                                                            ContrastColumnName = "CLUSTER", 
                                                            RGChannelSet = NA, 
                                                            SampleSheet = NA, 
                                                            ProduceImages = "Yes", 
                                                            PNGorPDF = "png")

ClinicalFile_T1SNFCluster2 <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all)
DifferentialMethylationProbesSNF <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1SNFCluster2, 
                                                            MvalueMatrix = OutputSDeviation5000All$MvalueMatrix_SD_Filtered[, matchDNAseqBetaMat_T1], 
                                                            BetaMatrix = OutputSDeviation5000All$BetaMatrix_SD_Filtered[, matchDNAseqBetaMat_T1],
                                                            ProbeOrGene = "Probe", 
                                                            ContrastColumnName = "CLUSTER", 
                                                            RGChannelSet = NA, 
                                                            SampleSheet = NA, 
                                                            ProduceImages = "Yes", 
                                                            PNGorPDF = "png")

# Only looking at TSS probes
ClinicalFile_T1SNFCluster2 <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all)
DifferentialMethylationProbesSNF <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1SNFCluster2, 
                                                            MvalueMatrix = MvalueMatrix_T1[matchTSSProbes, matchDNAseqBetaMat_T1], 
                                                            BetaMatrix = BetaMatrix_T1[matchTSSProbes, matchDNAseqBetaMat_T1],
                                                            ProbeOrGene = "Probe", 
                                                            ContrastColumnName = "CLUSTER", 
                                                            RGChannelSet = NA, 
                                                            SampleSheet = NA, 
                                                            ProduceImages = "Yes", 
                                                            PNGorPDF = "png")

ClinicalFile_T1SNFCluster2 <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all)
MethylationDensityOutput <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                   ClusterLabels = ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all, 
                                                   PlotWithinAnnotationCategories = "No", 
                                                   ClinicalCategoryToVisualize = "CLUSTER", 
                                                   BetaMatrix = BetaMatrix_T1[, matchDNAseqBetaMat_T1], 
                                                   AnnotationFile = EPICAnnotationFile, 
                                                   ClinicalFile = ClinicalFile_T1SNFCluster2, 
                                                   SampleSheet = NA, 
                                                   FigureGenerate = "Yes", 
                                                   ImageName = paste0("MethylationPlot_ENHFLmeth.tumor.status.all", ImageName), 
                                                   PNGorPDF = "png")

ClinicalFile_T1SNFCluster2 <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all)
# identical(colnames(RNAseqT3ProteinCodeNormalized[, matchSamplesNorm]), colnames(BetaMatrix_T1[, matchDNAseqBetaMat_T1]))
matchSamplesNorm <- match(colnames(as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq])), colnames(RNAseqT3ProteinCodeNormalized))
BoxPlotsClusters <- BoxPlotsMethylation(BetaMatrix = BetaMatrix_T1[, matchDNAseqBetaMat_T1], 
                                        MvalueMatrix = MvalueMatrix_T1[, matchDNAseqBetaMat_T1],
                                        RNAseqMatrix = as.matrix(RNAseqT3ProteinCodeNormalized[, matchSamplesNorm]), 
                                        ClinicalFile = ClinicalFile_T1SNFCluster2, 
                                        CategoryToVisualize = "CLUSTER", 
                                        ClusterLabels = ClinicalFile_T1SNFCluster2$CLUSTER, 
                                        TumorPurity = TumorPurity,
                                        PlotClustersWithinCategories = "Yes", 
                                        ProduceImages = "Yes",
                                        PNGorPDF = "png")


identical(colnames(BetaMatrix_T1[, matchDNAseqBeta]), colnames(as.matrix(RNAseqT3ProteinCodeNormalized)))

##### Connect SNF clusters (above) with Targeted Seq data: RESET_HypomethylatedProbes_FL_(T3 n = 131-10 = 121)_HypermethylatedNormal(Enhancers) ####
table(ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all)
#  1  2 
# 56 45 


# calculating for all patients
DNAseqMutSignficance <- matrix(NA, 
                               ncol = 10, 
                               nrow = nrow(mutRearrT1Robert101))
colnames(DNAseqMutSignficance) <- c("Xsquared", "XsquaredP.Value", 
                                    "XsquaredRes1-C1", "XsquaredRes0-C1", 
                                    "XsquaredRes1-C2", "XsquaredRes0-C2", 
                                    "fisherTOddsRatio", "fisherTP.Value",
                                    "NumbMutations",
                                    "FrequencyMinOneMutation")
rownames(DNAseqMutSignficance) <- rownames(mutRearrT1Robert101)
for (k in 1:ncol(mutRearrT1Robert101)) { # for each gene/mutation
  cat("\nk is:", k)
  if(! all(mutRearrT1Robert101[k, ] == 0, na.rm = T)) { # only if all entries are not zero
    Table <- table(MutationPresence = mutRearrT1Robert101[k, ], 
                   Cluster = ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all)
    Table <- Table[c(2, 1), ]
    # The reference group will be the first cell \n")
    chisqTest  <- stats::chisq.test(Table)
    fisherTest <- stats::fisher.test(Table)
    
    DNAseqMutSignficance[k, 1] <- chisqTest$statistic
    DNAseqMutSignficance[k, 2] <- chisqTest$p.value
    DNAseqMutSignficance[k, 3] <- chisqTest$residuals[1, 1]
    DNAseqMutSignficance[k, 4] <- chisqTest$residuals[2, 1]
    DNAseqMutSignficance[k, 5] <- chisqTest$residuals[1, 2]
    DNAseqMutSignficance[k, 6] <- chisqTest$residuals[2, 2]
    DNAseqMutSignficance[k, 7] <- fisherTest$estimate
    DNAseqMutSignficance[k, 8] <- fisherTest$p.value
    DNAseqMutSignficance[k, 9] <- length(which(mutRearrT1Robert101[, k] == 1))
    DNAseqMutSignficance[k, 10] <- (length(which(mutRearrT1Robert101[, k] == 1)) / 
                                      length(mutRearrT1Robert101[, k]))
  }
}



DNAseqMutSignficancePlottingT1 <- data.frame(
  DNAseqMutSignficance, 
  NegLogOdds = -log2(DNAseqMutSignficance[, 7]),
  NegLogPval = -log10(DNAseqMutSignficance[, 8]))






# Basic dot plot
# plot P value against odds ratio for 60 matching genes of PM_FL_PANEL
DNAseqMutSignficanceGGplotT1 <- ggplot2::ggplot(
  data = DNAseqMutSignficancePlottingT1, 
  aes(x = NegLogOdds, 
      y = NegLogPval, 
      size = FrequencyMinOneMutation)) + 
  ggplot2::geom_point() + 
  scale_y_continuous(name = "-log10(P value)") +
  # geom_text(aes(label = rownames(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlotting)), 
  # hjust = 0, vjust = 0) + 
  ggrepel::geom_label_repel(aes(label = rownames(DNAseqMutSignficancePlottingT1)),
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            segment.color = 'grey50') +
  labs(x = "-log2(odds ratio)") +
  theme_bw() + 
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "red")

# CREBBP mutation
CREBBPMutLabs <- table(ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all, mutRearrT1Robert101[, which(colnames(mutRearrT1Robert101) == "CREBBP")])[, c(2, 1)]
fisher.test(CREBBPMutLabs)


#### SNF: RESET_HypermethylatedProbesFL_(T3 n = 131-10 = 121)_HypomethylatedNormal(Silencers) ####

# bind tumor + normal samples
identical(rownames(SILFLmeth.tumor.status.all), 
          rownames(methNorSelOutput$normal.sil.probes)) # TRUE

dim(methNorSelOutput$normal.sil.probes) # 45229   5

# bind together tumor and normal for ENH (in normal)
EventMatrix <- SILFLmeth.tumor.status.all[, c(1:121)]
dim(EventMatrix) # 45229   121
# rownames(EventMatrix) <- gsub("@.*", "", rownames(EventMatrix)) # error


TumorPurity <- readRDS(file = paste0(MethylationDirPath, "/Purity_281probes_10Jan2020.rds"))
dim(TumorPurity) # 170   1

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

# cluster only RNAseq protein coding genes
RNAseqQC18June2020T3Samples132Protein <- RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations[which(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations$type == "protein_coding"), ]
dim(RNAseqQC18June2020T3Samples132Protein) # 17190   135
RNAseqQC18June2020T3Samples132Protein <- as.matrix(RNAseqQC18June2020T3Samples132Protein[, c(4:135)])
dim(RNAseqQC18June2020T3Samples132Protein) # 17190   132


TargetedDNAseqDirPath <- "~/Desktop/UHN/FLOMICS/TargetedSequencing"
mutMergedT1Robert <- read.csv(paste0(TargetedDNAseqDirPath, "/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl_25Aug2020.csv"), row.names = 1)
dim(mutMergedT1Robert) # 55 138
colnames(mutMergedT1Robert)
class(mutMergedT1Robert)

# Fixed in mut.merged.df.T1.poor.cov.excl_25Aug2020 file 
# mutMergedT1Robert[, which(colnames(mutMergedT1Robert) == "LY_FL_179_T1")] # NA values
# mutMergedT1Robert[, which(colnames(mutMergedT1Robert) == "LY_FL_181_T1")] # NA values
# which(is.na(colSums(mutMergedT1Robert)) == TRUE) # LY_FL_179_T1 50; LY_FL_181_T1  51 
# mutMergedT1Robert[, c(50, 51)] <- 0
which(is.na(colSums(mutMergedT1Robert)) == TRUE) # named integer(0)




matchRNAseqMutation <- match(colnames(RNAseqQC18June2020T3Samples132Protein), colnames(mutMergedT1Robert))
dim(mutMergedT1Robert[, matchRNAseqMutation[! is.na(matchRNAseqMutation)]]) # 55 101

mutMergedT1Robert101 <- mutMergedT1Robert[, matchRNAseqMutation[! is.na(matchRNAseqMutation)]]
dim(mutMergedT1Robert101) # 55 101

matchDNAseqBeta <- match(colnames(mutMergedT1Robert101), colnames(EventMatrix))

matchDNAseqRNAseq <- match(colnames(mutMergedT1Robert101), colnames(RNAseqQC18June2020T3Samples132Protein))

matchDNAseqClinicalFile <- match(colnames(mutMergedT1Robert101), ClinicalFile_T1$SAMPLE_ID)


# Combine mutations with BCL2 and BCL6 translocations
BCL2BCL6BreakapartPredictions <- data.table::fread(file = paste0(MethylationDirPath, "/2020-07-05_BCL2_BCL6_rearrangements.txt"))
BCL2BCL6BreakapartPredictions <- data.frame(BCL2BCL6BreakapartPredictions)
dim(BCL2BCL6BreakapartPredictions) # 183   3
rownames(BCL2BCL6BreakapartPredictions) <- BCL2BCL6BreakapartPredictions$External_ID
BCL2BCL6BreakapartPredictions <- BCL2BCL6BreakapartPredictions[, -1]
matchMutationRearrangements <- match(colnames(mutMergedT1Robert101), rownames(BCL2BCL6BreakapartPredictions))
rearrangementsT1Robert108 <- BCL2BCL6BreakapartPredictions[matchMutationRearrangements, ]
identical(colnames(t(rearrangementsT1Robert108)), colnames(mutMergedT1Robert101)) #TRUE
mutRearrT1Robert101 <- rbind(mutMergedT1Robert101, t(rearrangementsT1Robert108))
dim(mutRearrT1Robert101) # 57 101

which(is.na(colSums(mutRearrT1Robert101)) == TRUE)
#LY_FL_244_T1 LY_FL_253_T1 LY_FL_264_T1 
#62           71           82
mutRearrT1Robert101[56:57, which(is.na(colSums(mutRearrT1Robert101)) == TRUE)]
# manual adjustment
mutRearrT1Robert101[56:57, 62] <- c(0, 0)
mutRearrT1Robert101[56:57, 71] <- c(1, 0)
mutRearrT1Robert101[56:57, 82] <- c(0, 0)

# check
mutRearrT1Robert101[56:57, c(62, 71, 82)] 



ClusteringSNFMethEventRESET101_SILFLmeth.tumor.status.all <- 
  SNFClustering(RNAseqCountMatrix = as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq]), 
                MvalueMatrix = EventMatrix[, matchDNAseqBeta],
                BetaMatrix = EventMatrix[, matchDNAseqBeta], 
                TargetedSeqMatrix = mutRearrT1Robert101,
                MethylationAnnotationFile = EPICAnnotationFile,
                RNAseqAnnotationFile = ENSMBLid, 
                QCMatrix = RNAseqQC18June2020T3Samples132$QCMatrixMatchedSampleFiltered[matchDNAseqRNAseq, ],
                TumorPurity = as.matrix(TumorPurity[matchDNAseqTumorPurity, ], ncol = 1),
                ClinicalFile = ClinicalFile_T1[matchDNAseqClinicalFile, ], 
                SurvivalFile = SurvivalFile,
                NumberOfClusters = 2,
                ImageName = paste0("SILFLmethTumorStatusAll", date()),
                PNGorPDF = "png",
                ProduceImages = "Yes") 



ClusteringSNFMethEventRESET101_SILFLmeth.tumor.status.all <- ClusterLabels

table(ClusteringSNFMethEventRESET101_SILFLmeth.tumor.status.all)

InfiniumClustLabels <- readRDS(file = paste0(MethylationDirPath, "/InfiniumClustLabels2.RDS"))
InfSNFMethEvents <- match(colnames(EventMatrix[, matchDNAseqBeta]), InfiniumClustLabels$ID)
identical(colnames(EventMatrix[, matchDNAseqBeta]), unfactor(InfiniumClustLabels$ID[InfSNFMethEvents]))
# TRUE
table(ClusteringSNFMethEventRESET101_SILFLmeth.tumor.status.all, InfiniumClustLabels$Cluster[InfSNFMethEvents])
chisq.test(table(ClusteringSNFMethEventRESET101_SILFLmeth.tumor.status.all, InfiniumClustLabels$Cluster[InfSNF5000SD]))



ClinicalFile_T1SNFCluster <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventRESET101_SILFLmeth.tumor.status.all)
matchDNAseqBetaMat_T1 <- match(colnames(as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq])), colnames(BetaMatrix_T1))
DifferentialMethylationProbesSNF <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1SNFCluster, 
                                                            MvalueMatrix = MvalueMatrix_T1[, matchDNAseqBetaMat_T1], 
                                                            BetaMatrix = BetaMatrix_T1[, matchDNAseqBetaMat_T1],
                                                            ProbeOrGene = "Probe", 
                                                            ContrastColumnName = "CLUSTER", 
                                                            RGChannelSet = NA, 
                                                            SampleSheet = NA, 
                                                            ProduceImages = "Yes", 
                                                            PNGorPDF = "png")


ClinicalFile_T1SNFCluster <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventRESET101_SILFLmeth.tumor.status.all)
DifferentialMethylationProbesSNF <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1SNFCluster, 
                                                            MvalueMatrix = OutputSDeviation5000All$MvalueMatrix_SD_Filtered[, matchDNAseqBetaMat_T1], 
                                                            BetaMatrix = OutputSDeviation5000All$BetaMatrix_SD_Filtered[, matchDNAseqBetaMat_T1],
                                                            ProbeOrGene = "Probe", 
                                                            ContrastColumnName = "CLUSTER", 
                                                            RGChannelSet = NA, 
                                                            SampleSheet = NA, 
                                                            ProduceImages = "Yes", 
                                                            PNGorPDF = "png")

MethylationDensityOutput <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                   ClusterLabels = ClusteringSNFMethEventRESET101_SILFLmeth.tumor.status.all, 
                                                   PlotWithinAnnotationCategories = "No", 
                                                   ClinicalCategoryToVisualize = "CLUSTER", 
                                                   BetaMatrix = BetaMatrix_T1[, matchDNAseqBetaMat_T1], 
                                                   AnnotationFile = EPICAnnotationFile, 
                                                   ClinicalFile = ClinicalFile_T1SNFCluster, 
                                                   SampleSheet = NA, 
                                                   FigureGenerate = "Yes", 
                                                   ImageName = paste0("MethylationPlot_SILFLmeth.tumor.status.all", ImageName), 
                                                   PNGorPDF = "png")

rownames(RNAseqQC18June2020T3Samples132Protein)

DifferentialExpressionRNAseqOutput <- DifferentialExpressionRNAseq(RNAseqCountMatrix = RNAseqQC18June2020T3Samples132Protein[, c(11:132)], 
                                                                   ContrastColumnName = "TYPE", 
                                                                   AnnotationFileEnsemblToGenes = RNAseqAnnotationFile,
                                                                   ClinicalFile = ClinicalFile, 
                                                                   ClusterLabels = ClusterLabelsSNF[[Clusters]], 
                                                                   ProduceImages = "Yes", 
                                                                   PNGorPDF = "png")


ClinicalFile_T1SNFCluster <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventRESET101_SILFLmeth.tumor.status.all)
# identical(colnames(RNAseqT3ProteinCodeNormalized[, matchSamplesNorm]), colnames(BetaMatrix_T1[, matchDNAseqBetaMat_T1]))
matchSamplesNorm <- match(colnames(as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq])), colnames(RNAseqT3ProteinCodeNormalized))
BoxPlotsClusters <- BoxPlotsMethylation(BetaMatrix = BetaMatrix_T1[, matchDNAseqBetaMat_T1], 
                                        MvalueMatrix = MvalueMatrix_T1[, matchDNAseqBetaMat_T1],
                                        RNAseqMatrix = as.matrix(RNAseqT3ProteinCodeNormalized[, matchSamplesNorm]), 
                                        ClinicalFile = ClinicalFile_T1SNFCluster, 
                                        CategoryToVisualize = "CLUSTER", 
                                        ClusterLabels = ClinicalFile_T1SNFCluster2$CLUSTER, 
                                        TumorPurity = TumorPurity,
                                        PlotClustersWithinCategories = "Yes", 
                                        ProduceImages = "Yes",
                                        PNGorPDF = "png")



##### Connect SNF clusters (above) with Targeted Seq data: RESET_HypermethylatedProbesFL_(T3 n = 131-10 = 121)_HypomethylatedNormal(Silencers) ####
table(ClusteringSNFMethEventRESET101_SILFLmeth.tumor.status.all)
#  1  2 
# 61 40 


# calculating for all patients
DNAseqMutSignficance <- matrix(NA, 
                               ncol = 10, 
                               nrow = nrow(mutRearrT1Robert101))
colnames(DNAseqMutSignficance) <- c("Xsquared", "XsquaredP.Value", 
                                    "XsquaredRes1-C1", "XsquaredRes0-C1", 
                                    "XsquaredRes1-C2", "XsquaredRes0-C2", 
                                    "fisherTOddsRatio", "fisherTP.Value",
                                    "NumbMutations",
                                    "FrequencyMinOneMutation")
rownames(DNAseqMutSignficance) <- rownames(mutRearrT1Robert101)
for (k in 1:ncol(mutRearrT1Robert101)) { # for each gene/mutation
  cat("\nk is:", k)
  if(! all(mutRearrT1Robert101[k, ] == 0, na.rm = T)) { # only if all entries are not zero
    Table <- table(MutationPresence = mutRearrT1Robert101[k, ], 
                   Cluster = ClusteringSNFMethEventRESET101_SILFLmeth.tumor.status.all)
    Table <- Table[c(2, 1), ]
    # The reference group will be the first cell \n")
    chisqTest  <- stats::chisq.test(Table)
    fisherTest <- stats::fisher.test(Table)
    
    DNAseqMutSignficance[k, 1] <- chisqTest$statistic
    DNAseqMutSignficance[k, 2] <- chisqTest$p.value
    DNAseqMutSignficance[k, 3] <- chisqTest$residuals[1, 1]
    DNAseqMutSignficance[k, 4] <- chisqTest$residuals[2, 1]
    DNAseqMutSignficance[k, 5] <- chisqTest$residuals[1, 2]
    DNAseqMutSignficance[k, 6] <- chisqTest$residuals[2, 2]
    DNAseqMutSignficance[k, 7] <- fisherTest$estimate
    DNAseqMutSignficance[k, 8] <- fisherTest$p.value
    DNAseqMutSignficance[k, 9] <- length(which(mutRearrT1Robert101[, k] == 1))
    DNAseqMutSignficance[k, 10] <- (length(which(mutRearrT1Robert101[, k] == 1)) / 
                                      length(mutRearrT1Robert101[, k]))
  }
}



DNAseqMutSignficancePlottingT1 <- data.frame(
  DNAseqMutSignficance, 
  NegLogOdds = -log2(DNAseqMutSignficance[, 7]),
  NegLogPval = -log10(DNAseqMutSignficance[, 8]))






# Basic dot plot
# plot P value against odds ratio for 60 matching genes of PM_FL_PANEL
DNAseqMutSignficanceGGplotT1 <- ggplot2::ggplot(
  data = DNAseqMutSignficancePlottingT1, 
  aes(x = NegLogOdds, 
      y = NegLogPval, 
      size = FrequencyMinOneMutation)) + 
  ggplot2::geom_point() + 
  scale_y_continuous(name = "-log10(P value)") +
  # geom_text(aes(label = rownames(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlotting)), 
  # hjust = 0, vjust = 0) + 
  ggrepel::geom_label_repel(aes(label = rownames(DNAseqMutSignficancePlottingT1)),
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            segment.color = 'grey50') +
  labs(x = "-log2(odds ratio)") +
  theme_bw() + 
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "red")

# CREBBP mutation
CREBBPMutLabs <- table(ClusteringSNFMethEventRESET101_ENHFLmeth.tumor.status.all, mutRearrT1Robert101[, which(colnames(mutRearrT1Robert101) == "CREBBP")])[, c(2, 1)]
fisher.test(CREBBPMutLabs)



#### SNF: RESET_Hyper and Hypo methylated Probes FL_(T3 n = 131-10 = 121) Combined####


# bind together tumor and normal for ENH and SIL 
EventMatrix <- rbind(SILFLmeth.tumor.status.all[, c(1:121)], ENHFLmeth.tumor.status.all[, c(1:121)])
dim(EventMatrix) # 48027   121
length(unique(rownames(EventMatrix))) # 48027

length(which(rowSums(EventMatrix) == 0)) # 21256
EventMatrixNoZero <- EventMatrix[- which(rowSums(EventMatrix) == 0), ]


TumorPurity <- readRDS(file = paste0(MethylationDirPath, "/Purity_281probes_10Jan2020.rds"))
dim(TumorPurity) # 170   1

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

# cluster only RNAseq protein coding genes
RNAseqQC18June2020T3Samples132Protein <- RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations[which(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations$type == "protein_coding"), ]
dim(RNAseqQC18June2020T3Samples132Protein) # 17190   135
RNAseqQC18June2020T3Samples132Protein <- as.matrix(RNAseqQC18June2020T3Samples132Protein[, c(4:135)])
dim(RNAseqQC18June2020T3Samples132Protein) # 17190   132


TargetedDNAseqDirPath <- "~/Desktop/UHN/FLOMICS/TargetedSequencing"
mutMergedT1Robert <- read.csv(paste0(TargetedDNAseqDirPath, "/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl_25Aug2020.csv"), row.names = 1)
dim(mutMergedT1Robert) # 55 138
colnames(mutMergedT1Robert)
class(mutMergedT1Robert)

# Fixed in mut.merged.df.T1.poor.cov.excl_25Aug2020 file 
# mutMergedT1Robert[, which(colnames(mutMergedT1Robert) == "LY_FL_179_T1")] # NA values
# mutMergedT1Robert[, which(colnames(mutMergedT1Robert) == "LY_FL_181_T1")] # NA values
# which(is.na(colSums(mutMergedT1Robert)) == TRUE) # LY_FL_179_T1 50; LY_FL_181_T1  51 
# mutMergedT1Robert[, c(50, 51)] <- 0
which(is.na(colSums(mutMergedT1Robert)) == TRUE) # named integer(0)




matchRNAseqMutation <- match(colnames(RNAseqQC18June2020T3Samples132Protein), colnames(mutMergedT1Robert))
dim(mutMergedT1Robert[, matchRNAseqMutation[! is.na(matchRNAseqMutation)]]) # 55 101

mutMergedT1Robert101 <- mutMergedT1Robert[, matchRNAseqMutation[! is.na(matchRNAseqMutation)]]
dim(mutMergedT1Robert101) # 55 101

matchDNAseqBeta <- match(colnames(mutMergedT1Robert101), colnames(EventMatrix))

matchDNAseqRNAseq <- match(colnames(mutMergedT1Robert101), colnames(RNAseqQC18June2020T3Samples132Protein))

matchDNAseqClinicalFile <- match(colnames(mutMergedT1Robert101), ClinicalFile_T1$SAMPLE_ID)


# Combine mutations with BCL2 and BCL6 translocations
BCL2BCL6BreakapartPredictions <- data.table::fread(file = paste0(MethylationDirPath, "/2020-07-05_BCL2_BCL6_rearrangements.txt"))
BCL2BCL6BreakapartPredictions <- data.frame(BCL2BCL6BreakapartPredictions)
dim(BCL2BCL6BreakapartPredictions) # 183   3
rownames(BCL2BCL6BreakapartPredictions) <- BCL2BCL6BreakapartPredictions$External_ID
BCL2BCL6BreakapartPredictions <- BCL2BCL6BreakapartPredictions[, -1]
matchMutationRearrangements <- match(colnames(mutMergedT1Robert101), rownames(BCL2BCL6BreakapartPredictions))
rearrangementsT1Robert108 <- BCL2BCL6BreakapartPredictions[matchMutationRearrangements, ]
identical(colnames(t(rearrangementsT1Robert108)), colnames(mutMergedT1Robert101)) #TRUE
mutRearrT1Robert101 <- rbind(mutMergedT1Robert101, t(rearrangementsT1Robert108))
dim(mutRearrT1Robert101) # 57 101

which(is.na(colSums(mutRearrT1Robert101)) == TRUE)
#LY_FL_244_T1 LY_FL_253_T1 LY_FL_264_T1 
#62           71           82
mutRearrT1Robert101[56:57, which(is.na(colSums(mutRearrT1Robert101)) == TRUE)]
# manual adjustment
mutRearrT1Robert101[56:57, 62] <- c(0, 0)
mutRearrT1Robert101[56:57, 71] <- c(1, 0)
mutRearrT1Robert101[56:57, 82] <- c(0, 0)

# check
mutRearrT1Robert101[56:57, c(62, 71, 82)] 

# check if identical colnames
identical(colnames(as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq])), 
          colnames(EventMatrixNoZero[, matchDNAseqBeta]))
identical(colnames(EventMatrixNoZero[, matchDNAseqBeta]), colnames(mutRearrT1Robert101))

ClusteringSNFMethEventSilEnhRESET101 <- 
  SNFClustering(RNAseqCountMatrix = as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq]), 
                MvalueMatrix = EventMatrixNoZero[, matchDNAseqBeta],
                BetaMatrix = EventMatrixNoZero[, matchDNAseqBeta], 
                TargetedSeqMatrix = mutRearrT1Robert101,
                MethylationAnnotationFile = EPICAnnotationFile,
                RNAseqAnnotationFile = ENSMBLid, 
                QCMatrix = RNAseqQC18June2020T3Samples132$QCMatrixMatchedSampleFiltered[matchDNAseqRNAseq, ],
                TumorPurity = as.matrix(TumorPurity[matchDNAseqTumorPurity, ], ncol = 1),
                ClinicalFile = ClinicalFile_T1[matchDNAseqClinicalFile, ], 
                SurvivalFile = SurvivalFile,
                NumberOfClusters = 2,
                ImageName = paste0("SILENHFLmethTumorStatusAll", date()),
                PNGorPDF = "png",
                ProduceImages = "Yes") 



ClusteringSNFMethEventSilEnhRESET101 <- ClusterLabels

table(ClusteringSNFMethEventSilEnhRESET101)

InfiniumClustLabels <- readRDS(file = paste0(MethylationDirPath, "/InfiniumClustLabels2.RDS"))
InfSNFMethEvents <- match(colnames(EventMatrixNoZero[, matchDNAseqBeta]), InfiniumClustLabels$ID)
identical(colnames(EventMatrixNoZero[, matchDNAseqBeta]), unfactor(InfiniumClustLabels$ID[InfSNFMethEvents]))
# TRUE
table(ClusteringSNFMethEventSilEnhRESET101, InfiniumClustLabels$Cluster[InfSNFMethEvents])
chisq.test(table(ClusteringSNFMethEventSilEnhRESET101, InfiniumClustLabels$Cluster[InfSNF5000SD]))



ClinicalFile_T1SNFCluster <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventSilEnhRESET101)
matchDNAseqBetaMat_T1 <- match(colnames(as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq])), colnames(BetaMatrix_T1))
DifferentialMethylationProbesSNF <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1SNFCluster, 
                                                            MvalueMatrix = MvalueMatrix_T1[matchTSSProbes, matchDNAseqBetaMat_T1], 
                                                            BetaMatrix = BetaMatrix_T1[matchTSSProbes, matchDNAseqBetaMat_T1],
                                                            ProbeOrGene = "Probe", 
                                                            ContrastColumnName = "CLUSTER", 
                                                            RGChannelSet = NA, 
                                                            SampleSheet = NA, 
                                                            ProduceImages = "Yes", 
                                                            PNGorPDF = "png")


ClinicalFile_T1SNFCluster <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventSilEnhRESET101)
DifferentialMethylationProbesSNF <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1SNFCluster, 
                                                            MvalueMatrix = OutputSDeviation5000All$MvalueMatrix_SD_Filtered[, matchDNAseqBetaMat_T1], 
                                                            BetaMatrix = OutputSDeviation5000All$BetaMatrix_SD_Filtered[, matchDNAseqBetaMat_T1],
                                                            ProbeOrGene = "Probe", 
                                                            ContrastColumnName = "CLUSTER", 
                                                            RGChannelSet = NA, 
                                                            SampleSheet = NA, 
                                                            ProduceImages = "Yes", 
                                                            PNGorPDF = "png")

MethylationDensityOutput <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                   ClusterLabels = ClusteringSNFMethEventSilEnhRESET101, 
                                                   PlotWithinAnnotationCategories = "No", 
                                                   ClinicalCategoryToVisualize = "CLUSTER", 
                                                   BetaMatrix = BetaMatrix_T1[, matchDNAseqBetaMat_T1], 
                                                   AnnotationFile = EPICAnnotationFile, 
                                                   ClinicalFile = ClinicalFile_T1SNFCluster, 
                                                   SampleSheet = NA, 
                                                   FigureGenerate = "Yes", 
                                                   ImageName = paste0("MethylationPlot_SILFLmeth.tumor.status.all", ImageName), 
                                                   PNGorPDF = "png")

rownames(RNAseqQC18June2020T3Samples132Protein)


ClinicalFile_T1SNFCluster <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventSilEnhRESET101)
# identical(colnames(RNAseqT3ProteinCodeNormalized[, matchSamplesNorm]), colnames(BetaMatrix_T1[, matchDNAseqBetaMat_T1]))
matchSamplesNorm <- match(colnames(as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq])), colnames(RNAseqT3ProteinCodeNormalized))
BoxPlotsClusters <- BoxPlotsMethylation(BetaMatrix = BetaMatrix_T1[, matchDNAseqBetaMat_T1], 
                                        MvalueMatrix = MvalueMatrix_T1[, matchDNAseqBetaMat_T1],
                                        RNAseqMatrix = as.matrix(RNAseqT3ProteinCodeNormalized[, matchSamplesNorm]), 
                                        ClinicalFile = ClinicalFile_T1SNFCluster, 
                                        CategoryToVisualize = "CLUSTER", 
                                        ClusterLabels = ClinicalFile_T1SNFCluster2$CLUSTER, 
                                        TumorPurity = TumorPurity,
                                        PlotClustersWithinCategories = "Yes", 
                                        ProduceImages = "Yes",
                                        PNGorPDF = "png")





##### Connect SNF clusters (above) with Targeted Seq data: RESET_HypermethylatedProbesFL_(T3 n = 131-10 = 121)_HypomethylatedNormal(Silencers) ####
table(ClusteringSNFMethEventSilEnhRESET101)
#  1  2 
# 59 42


# calculating for all patients
DNAseqMutSignficance <- matrix(NA, 
                               ncol = 10, 
                               nrow = nrow(mutRearrT1Robert101))
colnames(DNAseqMutSignficance) <- c("Xsquared", "XsquaredP.Value", 
                                    "XsquaredRes1-C1", "XsquaredRes0-C1", 
                                    "XsquaredRes1-C2", "XsquaredRes0-C2", 
                                    "fisherTOddsRatio", "fisherTP.Value",
                                    "NumbMutations",
                                    "FrequencyMinOneMutation")
rownames(DNAseqMutSignficance) <- rownames(mutRearrT1Robert101)
for (k in 1:ncol(mutRearrT1Robert101)) { # for each gene/mutation
  cat("\nk is:", k)
  if(! all(mutRearrT1Robert101[k, ] == 0, na.rm = T)) { # only if all entries are not zero
    Table <- table(MutationPresence = mutRearrT1Robert101[k, ], 
                   Cluster = ClusteringSNFMethEventSilEnhRESET101)
    Table <- Table[c(2, 1), ]
    # The reference group will be the first cell \n")
    chisqTest  <- stats::chisq.test(Table)
    fisherTest <- stats::fisher.test(Table)
    
    DNAseqMutSignficance[k, 1] <- chisqTest$statistic
    DNAseqMutSignficance[k, 2] <- chisqTest$p.value
    DNAseqMutSignficance[k, 3] <- chisqTest$residuals[1, 1]
    DNAseqMutSignficance[k, 4] <- chisqTest$residuals[2, 1]
    DNAseqMutSignficance[k, 5] <- chisqTest$residuals[1, 2]
    DNAseqMutSignficance[k, 6] <- chisqTest$residuals[2, 2]
    DNAseqMutSignficance[k, 7] <- fisherTest$estimate
    DNAseqMutSignficance[k, 8] <- fisherTest$p.value
    DNAseqMutSignficance[k, 9] <- length(which(mutRearrT1Robert101[, k] == 1))
    DNAseqMutSignficance[k, 10] <- (length(which(mutRearrT1Robert101[, k] == 1)) / 
                                      length(mutRearrT1Robert101[, k]))
  }
}



DNAseqMutSignficancePlottingT1 <- data.frame(
  DNAseqMutSignficance, 
  NegLogOdds = -log2(DNAseqMutSignficance[, 7]),
  NegLogPval = -log10(DNAseqMutSignficance[, 8]))






# Basic dot plot
# plot P value against odds ratio for 60 matching genes of PM_FL_PANEL
DNAseqMutSignficanceGGplotT1 <- ggplot2::ggplot(
  data = DNAseqMutSignficancePlottingT1, 
  aes(x = NegLogOdds, 
      y = NegLogPval, 
      size = FrequencyMinOneMutation)) + 
  ggplot2::geom_point() + 
  scale_y_continuous(name = "-log10(P value)") +
  # geom_text(aes(label = rownames(RNAseqMutationCalls01ClusterLabsSignficancePMPanelFilteredPlotting)), 
  # hjust = 0, vjust = 0) + 
  ggrepel::geom_label_repel(aes(label = rownames(DNAseqMutSignficancePlottingT1)),
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            segment.color = 'grey50') +
  labs(x = "-log2(odds ratio)") +
  theme_bw() + 
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "red")

# CREBBP mutation
CREBBPMutLabs <- table(ClusteringSNFMethEventSilEnhRESET101, mutRearrT1Robert101[which(rownames(mutRearrT1Robert101) == "CREBBP"), ])[, c(2, 1)]
fisher.test(CREBBPMutLabs)




######################
#### Running RESET with Karin's normalized counts RNAseq dataset by Kallisto program ####

# Set pathways:
dim(RNAseqNormalizedCountDataTab) # 31600   132



#### RESET_HypomethylatedProbes_FL_(T3 n = 131-10 = 121)_HypermethylatedNormal(Enhancers) ####

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


# match methylation samples with RNAseq samples
matchBetaRNAseq <- match(colnames(RNAseqProtCodT3Mat[, c(11:131)]), 
                         colnames(BetaMatrix_T1))
# check colnames
identical(colnames(RNAseqProtCodT3Mat[, c(11:131)]), colnames(BetaMatrix_T1[, matchBetaRNAseq])) # TRUE

# Run reset function on all tumor samples - enhancers - FL only
set.seed(1234)
resetResultsENHFLonlyKallisto <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                                       meth.tumor = BetaMatrix_T1[, matchBetaRNAseq],
                                       transcriptome = RNAseqProtCodT3Mat[, c(11:131)],
                                       methylation.event = c('enh'),
                                       FDR.permutation.no = 100, 
                                       seed = 100)
saveRDS(resetResultsENHFLonlyKallisto, file = "40_resetResultsENHFLonlyKallisto_TumorHypo.rds")
resetResultsENHFLonlyKallisto <- readRDS("40_resetResultsENHFLonlyKallisto_TumorHypo.rds")

dim(resetResultsENHFLonlyKallisto$FDR.res) #  540   2
head(resetResultsENHFLonlyKallisto$meth.tumor.status.all)
ENHFLmeth.tumor.status.allKallisto <- data.frame(resetResultsENHFLonlyKallisto$meth.tumor.status.all, Gene = gsub(".*@", "", rownames(resetResultsENHFLonlyKallisto$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.allKallisto)
dim(ENHFLmeth.tumor.status.allKallisto) # 2798  122
dim(resetResultsENHFLonlyKallisto$meth.tumor.all) # 2798  121
dim(resetResultsENHFLonlyKallisto$transcriptome) # 2798  121
range(resetResultsENHFLonlyKallisto$FDR.res$FDR) # 0.1200000 0.9578519
base::range(resetResultsENHFLonlyKallisto$Score.report$Score, na.rm = TRUE) # 0.02489615 1.03836037

genelist <- c("SPO11", "FCRLB")
testingGene <- genelist[1]
ENHFLmeth.tumor.status.allKallisto[which(ENHFLmeth.tumor.status.allKallisto$Gene == testingGene), ]
rowSums(ENHFLmeth.tumor.status.allKallisto[which(ENHFLmeth.tumor.status.allKallisto$Gene == testingGene), c(1:121)])


# comparing events with stage
matchClinicalENHFLKallisto <- match(colnames(ENHFLmeth.tumor.status.allKallisto)[1:121], ClinicalFile_T1$SAMPLE_ID)
testingGene <- "MTSS1"
table(ClinicalFile_T1$STAGE[matchClinicalENHFLKallisto],
      ENHFLmeth.tumor.status.allKallisto[which(ENHFLmeth.tumor.status.allKallisto$Gene == testingGene)[1], c(1:121)])
chisq.test(table(ClinicalFile_T1$STAGE[matchClinicalENHFLKallisto],
                 ENHFLmeth.tumor.status.allKallisto[which(ENHFLmeth.tumor.status.allKallisto$Gene == testingGene)[1], c(1:121)]))


plotting <- function() {
  # If one probe
  if(length(which(gsub(".*@", "", rownames(resetResultsENHFLonlyKallisto$normal.meth)) == testingGene)) == 1) {
    
    
    testingDataFrame <- data.frame(Values = c(resetResultsENHFLonlyKallisto$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonlyKallisto$normal.meth)) == testingGene), ],
                                              resetResultsENHFLonlyKallisto$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonlyKallisto$meth.tumor.all)) == testingGene), ],
                                              resetResultsENHFLonlyKallisto$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonlyKallisto$transcriptome)) == testingGene), ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsENHFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsENHFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsENHFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    # only looking at samples with a status marked as 1
    hypoProbes <- which(ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.all$Gene == testingGene), ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ],
                                                    resetResultsENHFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), hypoProbes],
                                                    resetResultsENHFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) # 247 2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ]),
                                      MeanMethTumor = mean(resetResultsENHFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsENHFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ]),
                                         t(resetResultsENHFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsENHFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 250   4
    
    
    testingDataFrame %>%  # highest score: cg20802515.p1.SPO11 methylation only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsENHFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsENHFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsENHFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    # look at 2nd probe with highest score
    hypoProbes <- which(ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.all$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)[2], ],
                                                    resetResultsENHFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene)[2], hypoProbes],
                                                    resetResultsENHFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene)[2], hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) # 25  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsENHFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene)[2], hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene)[2], hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  }
}



resetScoreENHFLonlyKallisto <- data.frame(resetResultsENHFLonlyKallisto$Score.report[which(! is.na(resetResultsENHFLonlyKallisto$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFLonlyKallisto)) {
  resetScoreENHFLonlyKallisto[i, 6] <- resetResultsENHFLonlyKallisto$FDR.res$FDR[which(resetScoreENHFLonlyKallisto$Score[i] == resetResultsENHFLonlyKallisto$FDR.res$observed.scores)]
}
head(resetScoreENHFLonlyKallisto)
dim(resetScoreENHFLonlyKallisto) # 540 6
length(unique(resetScoreENHFLonlyKallisto$Gene)) # 452
range(resetScoreENHFLonlyKallisto$Score) # 0.02489615 1.03836037
range(resetScoreENHFLonlyKallisto$FDR) # 0.1200000 0.9578519
resetScoreENHFLonlyKallisto[which(resetScoreENHFLonlyKallisto$Gene == testingGene), ]


resetScoreENHFLonlyKallisto[order(resetScoreENHFLonlyKallisto$Score, decreasing = TRUE),]


# View(resetScoreENHFLonly)
# write.csv(resetScoreENHFLonly, file = "RESET_HypermethylatedProbes_FL_Enhancers.csv")

# Plot No.Methylation.Events vs Score
resetScoreENHFLonlyKallisto %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()

# Plot No.Methylation.Events vs Score
resetScoreENHFLonlyKallisto %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreENHFLonlyKallisto %>%  
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreENHFLonlyKallisto %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))





# plotting RNAseq vs methylation value - SPO11, FCRLB
# get the probe IDs
genelist <- c("SPO11", "FCRLB", "BCL2L10")
testingGene <- genelist[3]
probesToPlot <- resetScoreENHFLonlyKallisto[which(resetScoreENHFLonlyKallisto$Gene == testingGene), ]$ProbID
which(rownames(BetaMatrix_T1[, matchRNAseqBeta[11:131]]) == tolower(probesToPlot))
# 162297


cg09593462

# Get RNAseq counts
which(rownames(RNAseqNormalizedCountMatrixRESET) == testingGene) 
# 7077


# getting the events
which(ENHFLmeth.tumor.status.allKallisto$Gene == testingGene)
# 1586 1589 1592

GeneToPlot <- data.frame(Beta = (BetaMatrix_T1[which(rownames(BetaMatrix_T1[, matchRNAseqBeta[11:131]]) == tolower(probesToPlot)), matchRNAseqBeta[11:131]]),
                         Event = t(ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.allKallisto$Gene == testingGene), c(1:121)]),
                         RNAseqCounts = t(RNAseqNormalizedCountMatrixRESET[, c(1:121)][which(rownames(RNAseqNormalizedCountMatrixRESET[, c(1:121)]) == testingGene), c(11:131)]),
                         RNAseqNormalized = t(RNAseqT3ProteinCodeNormalized[which(rownames(RNAseqT3ProteinCode) == testingGene), ]))

ggplot2::ggplot(GeneToPlot, aes(x = Beta, y = FCRLB.1, 
                                color = factor(cg23641748.p17.FCRLB))) +
  geom_point(size = 2)+
  labs(title = "cg23641748 - FCRLB",
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





#### RESET_HypermethylatedProbesFL_(T3 n = 131-10 = 121)_HypomethylatedNormal(Silencers) ####

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


# match methylation samples with RNAseq samples
matchBetaRNAseq <- match(colnames(RNAseqProtCodT3Mat[, c(11:131)]), 
                         colnames(BetaMatrix_T1))
# check colnames
identical(colnames(RNAseqProtCodT3Mat[, c(11:131)]), colnames(BetaMatrix_T1[, matchBetaRNAseq])) # TRUE

# Run reset function on all tumor samples - silencers - FL only
set.seed(1234)
resetResultsSILFLonlyKallisto <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                                       meth.tumor = BetaMatrix_T1[, matchBetaRNAseq],
                                       transcriptome = RNAseqProtCodT3Mat[, c(11:131)],
                                       methylation.event = c('sil'),
                                       FDR.permutation.no = 100, 
                                       seed = 100)
# saveRDS(resetResultsSILFLonlyKallisto, file = "40_resetResultsSILFLonlyKallisto_TumorHyper.rds")
resetResultsSILFLonlyKallisto <- readRDS("40_resetResultsSILFLonlyKallisto_TumorHyper.rds")


dim(resetResultsSILFLonlyKallisto$FDR.res) #   5625   2
SILFLmeth.tumor.status.allKallisto <- data.frame(resetResultsSILFLonlyKallisto$meth.tumor.status.all, Gene = gsub(".*@", "", rownames(resetResultsSILFLonlyKallisto$meth.tumor.status.all)))
head(SILFLmeth.tumor.status.allKallisto)
dim(SILFLmeth.tumor.status.allKallisto) # 45229   122
dim(resetResultsSILFLonlyKallisto$meth.tumor.all) # 45229  121
dim(resetResultsSILFLonlyKallisto$transcriptome) # 45229  121
range(resetResultsSILFLonlyKallisto$FDR.res$FDR) # 0.3991304 1.1266667
range(resetResultsSILFLonlyKallisto$Score.report$Score,  na.rm = TRUE)  # 0.004880055 1.217787799


genelist <- c("BMP3", "MME", "HLTF", "CAMK1")
testingGene <- genelist[1]
SILFLmeth.tumor.status.allKallisto[which(SILFLmeth.tumor.status.allKallisto$Gene == testingGene), ]
rowSums(SILFLmeth.tumor.status.allKallisto[which(SILFLmeth.tumor.status.allKallisto$Gene == testingGene), c(1:121)])


# comparing events with stage
matchClinicalSILKallisto <- match(colnames(SILFLmeth.tumor.status.allKallisto)[1:121], ClinicalFile_T1$SAMPLE_ID)
testingGene <- "BMP7"
table(ClinicalFile_T1$STAGE[matchClinicalSILKallisto],
      SILFLmeth.tumor.status.allKallisto[which(SILFLmeth.tumor.status.allKallisto$Gene == testingGene)[1], c(1:121)])
chisq.test(table(ClinicalFile_T1$STAGE[matchClinicalSILKallisto],
                 SILFLmeth.tumor.status.allKallisto[which(SILFLmeth.tumor.status.allKallisto$Gene == testingGene)[1], c(1:121)]))



plotting <- function() {
  # If one probe
  if(length(which(gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene)) == 1) {
    testingDataFrame <- data.frame(Values = c(resetResultsSILFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene), ],
                                              resetResultsSILFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene), ],
                                              resetResultsSILFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene), ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsSILFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsSILFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsSILFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    hyperProbes <- which(SILFLmeth.tumor.status.all[which(SILFLmeth.tumor.status.all$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsSILFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene)[2], ],
                                                    resetResultsSILFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene)[2], hyperProbes],
                                                    resetResultsSILFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene)[2], hyperProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hyperProbes)), rep("TumorRNAseq", length(hyperProbes))))
    dim(testingDataFrameStatus) # 33  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsSILFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsSILFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene)[2], hyperProbes]),
                                      MeanTranscriptome = mean(resetResultsSILFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene)[2], hyperProbes]))
    
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsSILFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene), ]),
                                         t(resetResultsSILFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsSILFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247   4
    
    
    testingDataFrame %>%  # highest score: cg01941671@p1@BMP3 methylation only
      melt(id.vars = "Type")  %>%
      filter(variable == "cg01941671.p1.BMP3") %>% # highest score: cg01941671@p1@BMP3
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg01941671.p1.BMP3") %>% # highest score: cg01941671.p1.BMP3
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsSILFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsSILFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsSILFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
  }
}

resetScoreSILFLonlyKallisto <- data.frame(resetResultsSILFLonlyKallisto$Score.report[which(! is.na(resetResultsSILFLonlyKallisto$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSILFLonlyKallisto)) {
  resetScoreSILFLonlyKallisto[i, 6] <- resetResultsSILFLonlyKallisto$FDR.res$FDR[which(resetScoreSILFLonlyKallisto$Score[i] == resetResultsSILFLonlyKallisto$FDR.res$observed.scores)]
}
head(resetScoreSILFLonlyKallisto)
dim(resetScoreSILFLonlyKallisto) #  5621  6
length(unique(resetScoreSILFLonlyKallisto$Gene)) # 3425
range(resetScoreSILFLonlyKallisto$Score) # 0.004880055 1.217787799
range(resetScoreSILFLonlyKallisto$FDR) # 0.3904348 1.1033333
resetScoreSILFLonlyKallisto[which(resetScoreSILFLonlyKallisto$Gene == testingGene), ]



resetScoreSILFLonlyKallisto[order(resetScoreSILFLonlyKallisto$Score, decreasing = TRUE),]


# View(resetScoreSILFLonly)
# write.csv(resetScoreSILFLonly, file = "RESET_HypomethylatedProbes_FL_Silencers.csv")

# Plot No.Methylation.Events vs Score
resetScoreSILFLonlyKallisto %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()


# Plot No.Methylation.Events vs Score
resetScoreSILFLonlyKallisto %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreSILFLonlyKallisto %>% 
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreSILFLonlyKallisto %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# plotting RNAseq vs methylation value - BMP3, FCRLB
# get the probe IDs
genelist <- c("BMP3", "MME")
testingGene <- genelist[2]
probesToPlot <- resetScoreSILFLonlyKallisto[which(resetScoreSILFLonlyKallisto$Gene == testingGene), ]$ProbID
for (i in 1:length(probesToPlot)) {
  cat("\n", which(rownames(BetaMatrix_T1[, matchRNAseqBeta[11:131]]) == tolower(probesToPlot)[i]), "\n")
}
# 218098, 125779, 1060, 97335, 86432, 63883, 32900, 99749 


# Get RNAseq counts
which(rownames(RNAseqNormalizedCountMatrixRESET[, c(1:121)]) == testingGene) 
# 679


# getting the events
which(SILFLmeth.tumor.status.allKallisto$Gene == testingGene)
# 937 939 941

GeneToPlot <- data.frame("Beta" = t(BetaMatrix_T1[c(218098, 125779, 1060, 97335, 86432, 63883, 32900, 99749), matchRNAseqBeta[11:131]]),
                         "Event" = t(SILFLmeth.tumor.status.all[which(SILFLmeth.tumor.status.all$Gene == testingGene), c(1:121)]),
                         "RNAseqCounts" = t(RNAseqT3ProteinCode[which(rownames(RNAseqT3ProteinCode) == testingGene), c(11:131)]),
                         "RNAseqNormalized" = t(RNAseqNormalizedCountMatrixRESET[which(rownames(RNAseqNormalizedCountMatrixRESET) == testingGene), c(1:121)]))

ggplot2::ggplot(GeneToPlot, aes(x = Beta.cg07862554, y = MME.1, 
                                color = factor(Event.cg07862554.p8.MME))) +
  geom_point(size = 2)+
  labs(title = "cg07862554 - MME",
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



# Looking at events 



###########################################################################################
#### RESET_HypomethylatedProbes_FL_(T2 n = 104-10 = 94)_HypermethylatedNormal(Enhancers) ####


# Match samples between RNAseq and methylation - 104 tumor samples
matchRNAseqBeta104 <- match(colnames(RNAseqQC18June2020T2Samples104$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized),
                            colnames(BetaMatrix_T1))

# Run reset function on all tumor samples - enhancers - FL only
resetResultsENHFL104only <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                                  meth.tumor = BetaMatrix_T1[, matchRNAseqBeta104[11:104]],
                                  transcriptome = RNAseqT2ProteinCode[, c(11:104)],
                                  methylation.event = c('enh'),
                                  FDR.permutation.no = 100, 
                                  seed = 100)
dim(resetResultsENHFL104only$FDR.res) #  470   2
head(resetResultsENHFL104only$meth.tumor.status.all)
ENHFLmeth.tumor.status.all104 <- data.frame(resetResultsENHFL104only$meth.tumor.status.all, Gene = gsub(".*@", "", rownames(resetResultsENHFL104only$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.all104)
dim(ENHFLmeth.tumor.status.all104) # 2798  95
dim(resetResultsENHFL104only$meth.tumor.all) # 2798  94
dim(resetResultsENHFL104only$transcriptome) # 2798  94
resetResultsENHFL104only$FDR.res



genelist <- c("FCRLB")
testingGene <- genelist[1]
ENHFLmeth.tumor.status.all104[which(ENHFLmeth.tumor.status.all104$Gene == testingGene), ]
rowSums(ENHFLmeth.tumor.status.all104[which(ENHFLmeth.tumor.status.all104$Gene == testingGene), c(1:94)])


# If one probe
plotting <- function() {
  if(length(which(gsub(".*@", "", rownames(resetResultsENHFL104only$normal.meth)) == testingGene)) == 1) {
    
    
    testingDataFrame <- data.frame(Values = c(resetResultsENHFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFL104only$normal.meth)) == testingGene), ],
                                              resetResultsENHFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFL104only$meth.tumor.all)) == testingGene), ],
                                              resetResultsENHFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFL104only$transcriptome)) == testingGene), ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 94), rep("TumorRNAseq", 94)))
    dim(testingDataFrame) # 193 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsENHFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFL104only$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsENHFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFL104only$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsENHFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFL104only$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    # only looking at samples with a status marked as 1
    hypoProbes <- which(ENHFLmeth.tumor.status.all104[which(ENHFLmeth.tumor.status.all104$Gene == testingGene), ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFL104only$normal.meth)) == testingGene), ],
                                                    resetResultsENHFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFL104only$meth.tumor.all)) == testingGene), hypoProbes],
                                                    resetResultsENHFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFL104only$transcriptome)) == testingGene), hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) #  133   2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFL104only$normal.meth)) == testingGene), ]),
                                      MeanMethTumor = mean(resetResultsENHFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFL104only$meth.tumor.all)) == testingGene), hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFL104only$transcriptome)) == testingGene), hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsENHFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ]),
                                         t(resetResultsENHFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsENHFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 250   4
    
    
    testingDataFrame %>%  # highest score: cg20802515.p1.SPO11 methylation only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsENHFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsENHFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsENHFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    # look at 2nd probe with highest score
    hypoProbes <- which(ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.all$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)[2], ],
                                                    resetResultsENHFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene)[2], hypoProbes],
                                                    resetResultsENHFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene)[2], hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) # 25  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLonly$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsENHFLonly$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly$meth.tumor.all)) == testingGene)[2], hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLonly$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly$transcriptome)) == testingGene)[2], hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  }
}



resetScoreENHFL104only <- data.frame(resetResultsENHFL104only$Score.report[which(! is.na(resetResultsENHFL104only$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFL104only)) {
  resetScoreENHFL104only[i, 6] <- resetResultsENHFL104only$FDR.res$FDR[which(resetScoreENHFL104only$Score[i] == resetResultsENHFL104only$FDR.res$observed.scores)]
}
head(resetScoreENHFL104only)
dim(resetScoreENHFL104only) # 470   6
length(unique(resetScoreENHFL104only$Gene)) # 388
range(resetScoreENHFL104only$Score) # 0.02766526 1.05906088
resetScoreENHFL104only[which(resetScoreENHFL104only$Gene == testingGene), ]


# View(resetScoreENHFLonly)
# write.csv(resetScoreENHFLonly, file = "RESET_HypermethylatedProbes_FL_Enhancers.csv")

# Plot No.Methylation.Events vs Score
resetScoreENHFL104only %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()

# Plot No.Methylation.Events vs Score
resetScoreENHFL104only %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreENHFL104only %>%  
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreENHFL104only %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



#### RESET_HypermethylatedProbesFL_(T2 n = 104-10 = 94)_HypomethylatedNormal(Silencers) ####

# Match samples between RNAseq and methylation - 131 tumor samples
matchRNAseqBeta104 <- match(colnames(RNAseqQC18June2020T2Samples104$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized),
                            colnames(BetaMatrix_T1))

# Run reset function on all tumor samples - silencers - FL only
resetResultsSILFL104only <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                                  meth.tumor = BetaMatrix_T1[, matchRNAseqBeta104[11:104]],
                                  transcriptome = RNAseqT2ProteinCode[, c(11:104)],
                                  methylation.event = c('sil'),
                                  FDR.permutation.no = 100, 
                                  seed = 100)
dim(resetResultsSILFL104only$FDR.res) # 5257    2
SILFLmeth.tumor.status.all104 <- data.frame(resetResultsSILFL104only$meth.tumor.status.all, Gene = gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.status.all)))
head(SILFLmeth.tumor.status.all104)
dim(SILFLmeth.tumor.status.all104) # 45229   95
dim(resetResultsSILFL104only$meth.tumor.all) # 45229  94
dim(resetResultsSILFL104only$transcriptome) # 45229  94

genelist <- c("BMP7")
testingGene <- genelist[1]
SILFLmeth.tumor.status.all[which(SILFLmeth.tumor.status.all104$Gene == testingGene), ]
rowSums(SILFLmeth.tumor.status.all104[which(SILFLmeth.tumor.status.all104$Gene == testingGene), c(1:94)])



plotting <- function() {
  
  # If one probe
  if(length(which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene)) == 1) {
    testingDataFrame <- data.frame(Values = c(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ],
                                              resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ],
                                              resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    hyperProbes <- which(SILFLmeth.tumor.status.all104[which(SILFLmeth.tumor.status.all104$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene)[2], ],
                                                    resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene)[2], hyperProbes],
                                                    resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene)[2], hyperProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hyperProbes)), rep("TumorRNAseq", length(hyperProbes))))
    dim(testingDataFrameStatus) # 89  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene)[2], hyperProbes]),
                                      MeanTranscriptome = mean(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene)[2], hyperProbes]))
    
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ]),
                                         t(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247   4
    
    
    testingDataFrame %>%  # highest score: cg01941671@p1@BMP3 methylation only
      melt(id.vars = "Type")  %>%
      filter(variable == "cg01941671.p1.BMP3") %>% # highest score: cg01941671@p1@BMP3
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg01941671.p1.BMP3") %>% # highest score: cg01941671.p1.BMP3
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
  }
}


resetScoreSILFL104only <- data.frame(resetResultsSILFL104only$Score.report[which(! is.na(resetResultsSILFL104only$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSILFL104only)) {
  resetScoreSILFL104only[i, 6] <- resetResultsSILFL104only$FDR.res$FDR[which(resetScoreSILFL104only$Score[i] == resetResultsSILFL104only$FDR.res$observed.scores)]
}
head(resetScoreSILFL104only)
dim(resetScoreSILFL104only) #   5257    6
length(unique(resetScoreSILFL104only$Gene)) # 3120
range(resetScoreSILFL104only$Score) # 0.009003583 1.312216365
resetScoreSILFL104only[which(resetScoreSILFL104only$Gene == testingGene), ]

# View(resetScoreSILFLonly)
# write.csv(resetScoreSILFLonly, file = "RESET_HypomethylatedProbes_FL_Silencers.csv")

# Plot No.Methylation.Events vs Score
resetScoreSILFL104only %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()


# Plot No.Methylation.Events vs Score
resetScoreSILFL104only %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreSILFL104only %>% 
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreSILFL104only %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))




###########################################################################################
#### Comparison of T3 and T2 #####

resetScoreENHFLonlyGenes <- resetScoreENHFLonly$Gene
resetScoreSILFLonlyGenes <- resetScoreSILFLonly$Gene
resetScoreENHFL104onlyGenes <- resetScoreENHFL104only$Gene
resetScoreSILFL104onlyGenes <- resetScoreSILFL104only$Gene

n <- max(length(resetScoreENHFLonlyGenes), 
         length(resetScoreSILFLonlyGenes),
         length(resetScoreENHFL104onlyGenes),
         length(resetScoreSILFL104onlyGenes))

length(resetScoreENHFLonlyGenes) <- n                      
length(resetScoreSILFLonlyGenes) <- n
length(resetScoreENHFL104onlyGenes) <- n                      
length(resetScoreSILFL104onlyGenes) <- n

TableGenesENHSILT2T3 <- cbind(resetScoreENHFLonlyGenes,
                              resetScoreSILFLonlyGenes,
                              resetScoreENHFL104onlyGenes,
                              resetScoreSILFL104onlyGenes)

write.csv(TableGenesENHSILT2T3, "TableGenesENHSILT2T3_1Sept2020.csv")



resetScoreENHFLonlyProbes <- rownames(resetScoreENHFLonly)
resetScoreSILFLonlyProbes <- rownames(resetScoreSILFLonly)
resetScoreENHFL104onlyProbes <- rownames(resetScoreENHFL104only)
resetScoreSILFL104onlyProbes <- rownames(resetScoreSILFL104only)

n <- max(length(resetScoreENHFLonlyProbes), 
         length(resetScoreSILFLonlyProbes),
         length(resetScoreENHFL104onlyProbes),
         length(resetScoreSILFL104onlyProbes))

length(resetScoreENHFLonlyProbes) <- n                      
length(resetScoreSILFLonlyProbes) <- n
length(resetScoreENHFL104onlyProbes) <- n                      
length(resetScoreSILFL104onlyProbes) <- n

TableProbesENHSILT2T3 <- cbind(resetScoreENHFLonlyProbes,
                               resetScoreSILFLonlyProbes,
                               resetScoreENHFL104onlyProbes,
                               resetScoreSILFL104onlyProbes)

write.csv(TableProbesENHSILT2T3, "TableProbesENHSILT2T3_1Sept2020.csv")

###########################################################################################
#### RESET_HypomethylatedProbes_FL_(T3 with C1 samples as tumor)_HypermethylatedNormal(Enhancers) ####


# Get C1 FL only samples
matchRNAseqC1 <- match(InfiniumClustLabels$ID[which(InfiniumClustLabels$Cluster == 1)][-c(1:10)],
                       colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized))
matchRNAseqC1 <- matchRNAseqC1[! is.na(matchRNAseqC1)]
colnames(RNAseqT3ProteinCode[, matchRNAseqC1])# double check


# Match samples between RNAseq and methylation - 104 tumor samples
matchRNAseqBetaC1 <- match(colnames(RNAseqT3ProteinCode[, matchRNAseqC1]),
                           colnames(BetaMatrix_T1))


# ENSURE identical colnames
identical(colnames(BetaMatrix_T1[, matchRNAseqBetaC1]), colnames(RNAseqT3ProteinCode[, matchRNAseqC1]))
# TRUE


# Run reset function on C1 tumor samples - enhancers - FL only
resetResultsENHFLC1 <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                             meth.tumor = BetaMatrix_T1[, matchRNAseqBetaC1],
                             transcriptome = RNAseqT3ProteinCode[, matchRNAseqC1],
                             methylation.event = c('enh'),
                             FDR.permutation.no = 100, 
                             seed = 100)
dim(resetResultsENHFLC1$FDR.res) #  363   2
head(resetResultsENHFLC1$meth.tumor.status.all)
ENHFLmeth.tumor.status.allC1 <- data.frame(resetResultsENHFLC1$meth.tumor.status.all, Gene = gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.allC1)
dim(ENHFLmeth.tumor.status.allC1) # 2798  58
dim(resetResultsENHFLC1$meth.tumor.all) # 2798  57
dim(resetResultsENHFLC1$transcriptome) # 2798  57
resetResultsENHFLC1$FDR.res
range(resetResultsENHFLC1$Score.report$Score, na.rm = TRUE) # 0.03474982 0.96606185



genelist <- c("ATF5")
testingGene <- genelist[1]
ENHFLmeth.tumor.status.allC1[which(ENHFLmeth.tumor.status.allC1$Gene == testingGene), ]
rowSums(ENHFLmeth.tumor.status.allC1[which(ENHFLmeth.tumor.status.allC1$Gene == testingGene), c(1:57)])


# If one probe
plotting <- function() {
  if(length(which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene)) == 1) {
    
    
    testingDataFrame <- data.frame(Values = c(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ],
                                              resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), ],
                                              resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 94), rep("TumorRNAseq", 94)))
    dim(testingDataFrame) # 193 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    # only looking at samples with a status marked as 1
    hypoProbes <- which(ENHFLmeth.tumor.status.allC1[which(ENHFLmeth.tumor.status.allC1$Gene == testingGene), ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ],
                                                    resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), hypoProbes],
                                                    resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) #  133   2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ]),
                                      MeanMethTumor = mean(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ]),
                                         t(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 250   4
    
    
    testingDataFrame %>%  # highest score: cg20802515.p1.SPO11 methylation only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    # look at 2nd probe with highest score
    hypoProbes <- which(ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.all$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene)[2], ],
                                                    resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene)[2], hypoProbes],
                                                    resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene)[2], hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) # 25  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene)[2], hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene)[2], hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  }
}



resetScoreENHFLC1 <- data.frame(resetResultsENHFLC1$Score.report[which(! is.na(resetResultsENHFLC1$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFLC1)) {
  resetScoreENHFLC1[i, 6] <- resetResultsENHFLC1$FDR.res$FDR[which(resetScoreENHFLC1$Score[i] == resetResultsENHFLC1$FDR.res$observed.scores)]
}
head(resetScoreENHFLC1)
dim(resetScoreENHFLC1) # 363   6
length(unique(resetScoreENHFLC1$Gene)) # 312
range(resetScoreENHFLC1$Score) # 0.03474982 0.96606185
resetScoreENHFLC1[which(resetScoreENHFLC1$Gene == testingGene), ]


# View(resetScoreENHFLonly)
# write.csv(resetScoreENHFLonly, file = "RESET_HypermethylatedProbes_FL_Enhancers.csv")

# Plot No.Methylation.Events vs Score
resetScoreENHFLC1 %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()

# Plot No.Methylation.Events vs Score
resetScoreENHFLC1 %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreENHFLC1 %>%  
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL - C1")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreENHFLC1 %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


# GProfilerAnalysis
# get gene list of those with FDR < 0.5
resetScoreENHFLC1Genes <- resetScoreENHFLC1 %>%  
  dplyr::filter(FDR < 0.5) %>% 
  dplyr::select(Gene)

gprofilerResetScoreENHFLC1 <- GProfilerAnalysis(GeneIDs = list(resetScoreENHFLC1Genes$Gene),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = NA, # with respect to expression
                                                ConditionName = "All",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
gprofilerResetScoreENHFLC1$shortLink
# "https://biit.cs.ut.ee/gplink/l/PC-KEBe5Rj"

# run with all genes 
gprofilerResetScoreENHFLC1All <- GProfilerAnalysis(GeneIDs = list(resetScoreENHFLC1$Gene),
                                                   Organism = "hsapiens",
                                                   OrderedQuery = TRUE,
                                                   PvalAlphaLevel = 0.01,
                                                   PositiveorNegFC = NA, # with respect to expression
                                                   ConditionName = "All",
                                                   ProduceImages = "Yes", 
                                                   PNGorPDF = "png")
gprofilerResetScoreENHFLC1All$shortLink
# "https://biit.cs.ut.ee/gplink/l/MABDo3KOQG"

#### RESET_HypermethylatedProbesFL_(T3 with C1 samples as tumor)_HypomethylatedNormal(Silencers) ####

# Get C1 FL only samples
matchRNAseqC1 <- match(InfiniumClustLabels$ID[which(InfiniumClustLabels$Cluster == 1)][-c(1:10)],
                       colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized))
matchRNAseqC1 <- matchRNAseqC1[! is.na(matchRNAseqC1)]
colnames(RNAseqT3ProteinCode[, matchRNAseqC1])# double check


# Match samples between RNAseq and methylation - 104 tumor samples
matchRNAseqBetaC1 <- match(colnames(RNAseqT3ProteinCode[, matchRNAseqC1]),
                           colnames(BetaMatrix_T1))


# ENSURE identical colnames
identical(colnames(BetaMatrix_T1[, matchRNAseqBetaC1]), colnames(RNAseqT3ProteinCode[, matchRNAseqC1]))
# TRUE


# Run reset function on C1 tumor samples - enhancers - FL only
resetResultsSILFLC1 <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                             meth.tumor = BetaMatrix_T1[, matchRNAseqBetaC1],
                             transcriptome = RNAseqT3ProteinCode[, matchRNAseqC1],
                             methylation.event = c('sil'),
                             FDR.permutation.no = 100, 
                             seed = 100)
dim(resetResultsSILFLC1$FDR.res) # 4359    2
SILFLmeth.tumor.status.allC1 <- data.frame(resetResultsSILFLC1$meth.tumor.status.all, Gene = gsub(".*@", "", rownames(resetResultsSILFLC1$meth.tumor.status.all)))
head(SILFLmeth.tumor.status.allC1)
dim(SILFLmeth.tumor.status.allC1) #  45229    58
dim(resetResultsSILFLC1$meth.tumor.all) # 45229   57
dim(resetResultsSILFLC1$transcriptome) # 45229   57

genelist <- c("BMP7")
testingGene <- genelist[1]
SILFLmeth.tumor.status.allC1[which(SILFLmeth.tumor.status.allC1$Gene == testingGene), ]
rowSums(SILFLmeth.tumor.status.allC1[which(SILFLmeth.tumor.status.allC1$Gene == testingGene), c(1:57)])



plotting <- function() {
  
  # If one probe
  if(length(which(gsub(".*@", "", rownames(resetResultsSILFLC1$normal.meth)) == testingGene)) == 1) {
    testingDataFrame <- data.frame(Values = c(resetResultsSILFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC1$normal.meth)) == testingGene), ],
                                              resetResultsSILFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC1$meth.tumor.all)) == testingGene), ],
                                              resetResultsSILFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC1$transcriptome)) == testingGene), ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsSILFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC1$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsSILFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC1$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsSILFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC1$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    hyperProbes <- which(SILFLmeth.tumor.status.allC1[which(SILFLmeth.tumor.status.allC1$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsSILFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC1$normal.meth)) == testingGene)[2], ],
                                                    resetResultsSILFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC1$meth.tumor.all)) == testingGene)[2], hyperProbes],
                                                    resetResultsSILFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC1$transcriptome)) == testingGene)[2], hyperProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hyperProbes)), rep("TumorRNAseq", length(hyperProbes))))
    dim(testingDataFrameStatus) # 61  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsSILFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC1$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsSILFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC1$meth.tumor.all)) == testingGene)[2], hyperProbes]),
                                      MeanTranscriptome = mean(resetResultsSILFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC1$transcriptome)) == testingGene)[2], hyperProbes]))
    
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ]),
                                         t(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247   4
    
    
    testingDataFrame %>%  # highest score: cg01941671@p1@BMP3 methylation only
      melt(id.vars = "Type")  %>%
      filter(variable == "cg01941671.p1.BMP3") %>% # highest score: cg01941671@p1@BMP3
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg01941671.p1.BMP3") %>% # highest score: cg01941671.p1.BMP3
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
  }
}


resetScoreSILFLC1 <- data.frame(resetResultsSILFLC1$Score.report[which(! is.na(resetResultsSILFLC1$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSILFLC1)) {
  resetScoreSILFLC1[i, 6] <- resetResultsSILFLC1$FDR.res$FDR[which(resetScoreSILFLC1$Score[i] == resetResultsSILFLC1$FDR.res$observed.scores)]
}
head(resetScoreSILFLC1)
dim(resetScoreSILFLC1) #  4359    6
length(unique(resetScoreSILFLC1$Gene)) # 2493  
range(resetScoreSILFLC1$Score) # 0.007039039 1.722105323
resetScoreSILFLC1[which(resetScoreSILFLC1$Gene == testingGene), ]

# View(resetScoreSILFLonly)
# write.csv(resetScoreSILFLonly, file = "RESET_HypomethylatedProbes_FL_Silencers.csv")

# Plot No.Methylation.Events vs Score
resetScoreSILFLC1 %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()


# Plot No.Methylation.Events vs Score
resetScoreSILFLC1 %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreSILFLC1 %>% 
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreSILFLC1 %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


# GProfilerAnalysis
# get gene list of those with FDR < 0.5
resetScoreSILFLC1Genes <- resetScoreSILFLC1 %>%  
  dplyr::filter(FDR < 0.5) %>% 
  dplyr::select(Gene)

gprofilerResetScoreSILFLC1 <- GProfilerAnalysis(GeneIDs = list(resetScoreSILFLC1Genes$Gene),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = NA, # with respect to expression
                                                ConditionName = "All",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
gprofilerResetScoreSILFLC1$shortLink
# "https://biit.cs.ut.ee/gplink/l/Qn9NIgj9Qo"

# run with all genes 
gprofilerResetScoreSILFLC1All <- GProfilerAnalysis(GeneIDs = list(resetScoreSILFLC1$Gene),
                                                   Organism = "hsapiens",
                                                   OrderedQuery = TRUE,
                                                   PvalAlphaLevel = 0.01,
                                                   PositiveorNegFC = NA, # with respect to expression
                                                   ConditionName = "All",
                                                   ProduceImages = "Yes", 
                                                   PNGorPDF = "png")
gprofilerResetScoreSILFLC1All$shortLink
# "https://biit.cs.ut.ee/gplink/l/q8dXk_IkRa"


###########################################################################################
#### RESET_HypermethylatedProbesFL_(T3 with C2 samples as tumor)_HypomethylatedNormal(Silencers) ####

# Get C2 FL only samples
matchRNAseqC2 <- match(InfiniumClustLabels$ID[which(InfiniumClustLabels$Cluster == 2)][- c(92:96)],
                       colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized))
matchRNAseqC2 <- matchRNAseqC2[! is.na(matchRNAseqC2)]
colnames(RNAseqT3ProteinCode[, matchRNAseqC2])# double check


# Match samples between RNAseq and methylation - 104 tumor samples
matchRNAseqBetaC2 <- match(colnames(RNAseqT3ProteinCode[, matchRNAseqC2]),
                           colnames(BetaMatrix_T1))


# ENSURE identical colnames
identical(colnames(BetaMatrix_T1[, matchRNAseqBetaC2]), colnames(RNAseqT3ProteinCode[, matchRNAseqC2]))
# TRUE


# Run reset function on C2 tumor samples - enhancers - FL only
resetResultsSILFLC2 <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                             meth.tumor = BetaMatrix_T1[, matchRNAseqBetaC2],
                             transcriptome = RNAseqT3ProteinCode[, matchRNAseqC2],
                             methylation.event = c('sil'),
                             FDR.permutation.no = 100, 
                             seed = 100)
dim(resetResultsSILFLC2$FDR.res) # 3808    2
range(resetResultsSILFLC2$Score.report$Score, na.rm = TRUE)

SILFLmeth.tumor.status.allC2 <- data.frame(resetResultsSILFLC2$meth.tumor.status.all, Gene = gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.status.all)))
head(SILFLmeth.tumor.status.allC2)
dim(SILFLmeth.tumor.status.allC2) # 45229    65
dim(resetResultsSILFLC2$meth.tumor.all) # 45229    64
dim(resetResultsSILFLC2$transcriptome) # 45229    64
dim(resetResultsSILFLC2$score.cutoff)  # NULL

genelist <- c("WNT9A")
testingGene <- genelist[1]
SILFLmeth.tumor.status.allC2[which(SILFLmeth.tumor.status.allC2$Gene == testingGene), ]
rowSums(SILFLmeth.tumor.status.allC2[which(SILFLmeth.tumor.status.allC2$Gene == testingGene), c(1:57)])



plotting <- function() {
  
  # If one probe
  if(length(which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene)) == 1) {
    testingDataFrame <- data.frame(Values = c(resetResultsSILFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene), ],
                                              resetResultsSILFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.all)) == testingGene), ],
                                              resetResultsSILFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC2$transcriptome)) == testingGene), ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsSILFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsSILFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsSILFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC2$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    hyperProbes <- which(SILFLmeth.tumor.status.allC2[which(SILFLmeth.tumor.status.allC2$Gene == testingGene), ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsSILFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene), ],
                                                    resetResultsSILFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.all)) == testingGene), hyperProbes],
                                                    resetResultsSILFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC2$transcriptome)) == testingGene), hyperProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hyperProbes)), rep("TumorRNAseq", length(hyperProbes))))
    dim(testingDataFrameStatus) # 9  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsSILFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene), ]),
                                      MeanMethTumor = mean(resetResultsSILFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.all)) == testingGene), hyperProbes]),
                                      MeanTranscriptome = mean(resetResultsSILFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC2$transcriptome)) == testingGene), hyperProbes]))
    
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ]),
                                         t(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247   4
    
    
    testingDataFrame %>%  # highest score: cg01941671@p1@BMP3 methylation only
      melt(id.vars = "Type")  %>%
      filter(variable == "cg01941671.p1.BMP3") %>% # highest score: cg01941671@p1@BMP3
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg01941671.p1.BMP3") %>% # highest score: cg01941671.p1.BMP3
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
  }
}


resetScoreSILFLC2 <- data.frame(resetResultsSILFLC2$Score.report[which(! is.na(resetResultsSILFLC2$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSILFLC2)) {
  resetScoreSILFLC2[i, 6] <- resetResultsSILFLC2$FDR.res$FDR[which(resetScoreSILFLC2$Score[i] == resetResultsSILFLC2$FDR.res$observed.scores)]
}
head(resetScoreSILFLC2)
dim(resetScoreSILFLC2) #  3808    6
range(resetScoreSILFLC2$FDR) # 0.0900000 0.9553846
length(unique(resetScoreSILFLC2$Gene)) # 2386
range(resetScoreSILFLC2$Score) # 0.001555952 1.583194347
resetScoreSILFLC2[which(resetScoreSILFLC2$Gene == testingGene), ]

# View(resetScoreSILFLonly)
# write.csv(resetScoreSILFLonly, file = "RESET_HypomethylatedProbes_FL_Silencers.csv")

# Plot No.Methylation.Events vs Score
resetScoreSILFLC2 %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()


# Plot No.Methylation.Events vs Score
resetScoreSILFLC2 %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreSILFLC2 %>% 
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreSILFLC2 %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


# run with all genes 
gprofilerResetScoreSILFLC2All <- GProfilerAnalysis(GeneIDs = list(resetScoreSILFLC2$Gene),
                                                   Organism = "hsapiens",
                                                   OrderedQuery = TRUE,
                                                   PvalAlphaLevel = 0.01,
                                                   PositiveorNegFC = NA, # with respect to expression
                                                   ConditionName = "All",
                                                   ProduceImages = "Yes", 
                                                   PNGorPDF = "png")
gprofilerResetScoreSILFLC2All$shortLink
# "https://biit.cs.ut.ee/gplink/l/KBybBHDXSV"


#### RESET_HypomethylatedProbes_FL_(T3 with C2 samples as tumor)_HypermethylatedNormal(Enhancers) ####


# Get C2 FL only samples
matchRNAseqC2 <- match(InfiniumClustLabels$ID[which(InfiniumClustLabels$Cluster == 2)][- c(92:96)],
                       colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized))
matchRNAseqC2 <- matchRNAseqC2[! is.na(matchRNAseqC2)]
colnames(RNAseqT3ProteinCode[, matchRNAseqC2])# double check


# Match samples between RNAseq and methylation - 104 tumor samples
matchRNAseqBetaC2 <- match(colnames(RNAseqT3ProteinCode[, matchRNAseqC2]),
                           colnames(BetaMatrix_T1))


# ENSURE identical colnames
identical(colnames(BetaMatrix_T1[, matchRNAseqBetaC2]), colnames(RNAseqT3ProteinCode[, matchRNAseqC2]))


# Run reset function on C2 tumor samples - enhancers - FL only
resetResultsENHFLC2 <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                             meth.tumor = BetaMatrix_T1[, matchRNAseqBetaC2],
                             transcriptome = RNAseqT3ProteinCode[, matchRNAseqC2],
                             methylation.event = c('enh'),
                             FDR.permutation.no = 100, 
                             seed = 100)

dim(resetResultsENHFLC2$FDR.res) #  494   2
head(resetResultsENHFLC2$meth.tumor.status.all)
ENHFLmeth.tumor.status.allC2 <- data.frame(resetResultsENHFLC2$meth.tumor.status.all, Gene = gsub(".*@", "", rownames(resetResultsENHFLC2$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.allC2)
dim(ENHFLmeth.tumor.status.allC2) # 2798  65
dim(resetResultsENHFLC2$meth.tumor.all) # 2798  64
dim(resetResultsENHFLC2$transcriptome) # 2798  64
resetResultsENHFLC2$FDR.res
range(resetResultsENHFLC2$Score.report$Score, na.rm = TRUE) # 0.009277029 1.168854874



genelist <- c("TBC1D12")
testingGene <- genelist[1]
ENHFLmeth.tumor.status.allC2[which(ENHFLmeth.tumor.status.allC2$Gene == testingGene), ]
rowSums(ENHFLmeth.tumor.status.allC2[which(ENHFLmeth.tumor.status.allC2$Gene == testingGene), c(1:57)])


# If one probe
plotting <- function() {
  if(length(which(gsub(".*@", "", rownames(resetResultsENHFLC2$normal.meth)) == testingGene)) == 1) {
    
    
    testingDataFrame <- data.frame(Values = c(resetResultsENHFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC2$normal.meth)) == testingGene), ],
                                              resetResultsENHFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC2$meth.tumor.all)) == testingGene), ],
                                              resetResultsENHFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC2$transcriptome)) == testingGene), ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 94), rep("TumorRNAseq", 94)))
    dim(testingDataFrame) # 193 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsENHFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC2$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsENHFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC2$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsENHFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC2$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    # only looking at samples with a status marked as 1
    hypoProbes <- which(ENHFLmeth.tumor.status.allC2[which(ENHFLmeth.tumor.status.allC2$Gene == testingGene), ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC2$normal.meth)) == testingGene), ],
                                                    resetResultsENHFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC2$meth.tumor.all)) == testingGene), hypoProbes],
                                                    resetResultsENHFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC2$transcriptome)) == testingGene), hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) #  133   2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC2$normal.meth)) == testingGene), ]),
                                      MeanMethTumor = mean(resetResultsENHFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC2$meth.tumor.all)) == testingGene), hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC2$transcriptome)) == testingGene), hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsENHFLC2$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsENHFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC2$normal.meth)) == testingGene), ]),
                                         t(resetResultsENHFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC2$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsENHFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC2$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 250   4
    
    
    testingDataFrame %>%  # highest score: cg20802515.p1.SPO11 methylation only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsENHFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC2$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsENHFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC2$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsENHFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC2$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    # look at 2nd probe with highest score
    hypoProbes <- which(ENHFLmeth.tumor.status.allC2[which(ENHFLmeth.tumor.status.allC2$Gene == testingGene), ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC2$normal.meth)) == testingGene), ],
                                                    resetResultsENHFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC2$meth.tumor.all)) == testingGene), hypoProbes],
                                                    resetResultsENHFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC2$transcriptome)) == testingGene), hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) # 25  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC2$normal.meth)) == testingGene), ]),
                                      MeanMethTumor = mean(resetResultsENHFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC2$meth.tumor.all)) == testingGene), hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC2$transcriptome)) == testingGene), hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  }
}



resetScoreENHFLC2 <- data.frame(resetResultsENHFLC2$Score.report[which(! is.na(resetResultsENHFLC2$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFLC2)) {
  resetScoreENHFLC2[i, 6] <- resetResultsENHFLC2$FDR.res$FDR[which(resetScoreENHFLC2$Score[i] == resetResultsENHFLC2$FDR.res$observed.scores)]
}
head(resetScoreENHFLC2)
dim(resetScoreENHFLC2) # 494    6
length(unique(resetScoreENHFLC2$Gene)) # 424
range(resetScoreENHFLC2$Score) # 0.009277029 1.168854874
resetScoreENHFLC2[which(resetScoreENHFLC2$Gene == testingGene), ]


# View(resetScoreENHFLonly)
# write.csv(resetScoreENHFLonly, file = "RESET_HypermethylatedProbes_FL_Enhancers.csv")

# Plot No.Methylation.Events vs Score
resetScoreENHFLC2 %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()

# Plot No.Methylation.Events vs Score
resetScoreENHFLC2 %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreENHFLC2 %>%  
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL - C2")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreENHFLC2 %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))

# GProfilerAnalysis
# get gene list of those with FDR < 0.5
resetScoreENHFLC2Genes <- resetScoreENHFLC2 %>%  
  dplyr::filter(FDR < 0.5) %>% 
  dplyr::select(Gene)

gprofilerResetScoreENHFLC2 <- GProfilerAnalysis(GeneIDs = list(resetScoreENHFLC2Genes$Gene),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = NA, # with respect to expression
                                                ConditionName = "All",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
# No results to show
# Please make sure that the organism is correct or set significant = FALSE


# run with all genes 
gprofilerResetScoreENHFLC2All <- GProfilerAnalysis(GeneIDs = list(resetScoreENHFLC2$Gene),
                                                   Organism = "hsapiens",
                                                   OrderedQuery = TRUE,
                                                   PvalAlphaLevel = 0.01,
                                                   PositiveorNegFC = NA, # with respect to expression
                                                   ConditionName = "All",
                                                   ProduceImages = "Yes", 
                                                   PNGorPDF = "png")
gprofilerResetScoreENHFLC2All$shortLink

###########################################################################################
#### Run first function from RESET with C2 as normal####

# Get C2 FL only samples
matchRNAseqC2 <- match(InfiniumClustLabels$ID[which(InfiniumClustLabels$Cluster == 2)][- c(92:96)],
                       colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized))
matchRNAseqC2 <- matchRNAseqC2[! is.na(matchRNAseqC2)]
colnames(RNAseqT3ProteinCode[, matchRNAseqC2])# double check


# Match samples between RNAseq and methylation - 104 tumor samples
matchRNAseqBetaC2 <- match(colnames(RNAseqT3ProteinCode[, matchRNAseqC2]),
                           colnames(BetaMatrix_T1))

# ENSURE identical colnames
identical(colnames(BetaMatrix_T1[, matchRNAseqBetaC2]), colnames(RNAseqT3ProteinCode[, matchRNAseqC2]))
# TRUE

methNorSelOutputC2 <- methNorSel(normal.mtx = BetaMatrix_T1[matchTSSProbes, matchRNAseqBetaC2],
                                 probe.list = promoter.probes.list)
# this function first extract the data regarding the selected probes, and later based on the beta-values prepare two seperate normal-sample data-sets:
## 1. hypermethylated probes in normal condition: Enh.probes
## 2. hypomethylated probes in normal condition: Sil.probes
names(methNorSelOutputC2) # "normal.sil.probes" "normal.enh.probes"
typeof(methNorSelOutputC2)
dim(methNorSelOutputC2$normal.sil.probes) # 40718     64
dim(methNorSelOutputC2$normal.enh.probes) # 1715  64
class(methNorSelOutputC2$normal.sil.probes)
normal.sil.probesC2 <- data.frame(methNorSelOutputC2$normal.sil.probes, rowMedians = rowMedians(methNorSelOutputC2$normal.sil.probes))

ggplot(normal.sil.probesC2, aes(x = c(1:40718), y = rowMedians)) + # plot rowMeans of sil.probes
  geom_point(size = 2) +
  labs(x = "probe", y = "Median Beta Value") +
  theme_bw() + 
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  scale_color_manual(values = c("#d6604d", "#66bd63", "#4575b4")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))

silProbesFLC2 <- gsub(".*@", "", rownames(methNorSelOutputC2$normal.sil.probes))
length(silProbesFLC2) # 40718
# write.csv(silProbesFLC2, file = "RESET_HypermethylatedProbesNormalCondition_Enhancers.csv")
source("26_Gprofiler.R")
gprofilersilProbesFLC2 <- GProfilerAnalysis(GeneIDs = list(silProbesFLC2),
                                            Organism = "hsapiens",
                                            OrderedQuery = TRUE,
                                            PvalAlphaLevel = 0.01,
                                            PositiveorNegFC = NA, # with respect to expression
                                            ConditionName = "All",
                                            ProduceImages = "Yes", 
                                            PNGorPDF = "png")
gprofilersilProbesFLC2$shortLink  # "https://biit.cs.ut.ee/gplink/l/ZLxq3tvKRC"



enhProbesFLC2 <- gsub(".*@", "", rownames(methNorSelOutputC2$normal.enh.probes))
length(enhProbesFLC2) # 1715
# write.csv(enhProbesFL, file = "RESET_HypomethylatedProbesNormalCondition_Silencers.csv")
normal.enh.probesC2 <- data.frame(methNorSelOutputC2$normal.enh.probes, rowMedians = rowMedians(methNorSelOutputC2$normal.enh.probes))

ggplot(normal.enh.probesC2, aes(x = c(1:1715), y = rowMedians)) + # plot rowMeans of enh.probes
  geom_point(size = 2) +
  labs(x = "probe", y = "Median Beta Value") +
  theme_bw() + 
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  scale_color_manual(values = c("#d6604d", "#66bd63", "#4575b4")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))
silProbesFLC2 <- gsub(".*@", "", rownames(methNorSelOutput$normal.sil.probes))

gprofilerenhProbesFLC2 <- GProfilerAnalysis(GeneIDs = list(silProbesFLC2),
                                            Organism = "hsapiens",
                                            OrderedQuery = TRUE,
                                            PvalAlphaLevel = 0.01,
                                            PositiveorNegFC = NA, # with respect to expression
                                            ConditionName = "All",
                                            ProduceImages = "Yes", 
                                            PNGorPDF = "png")
gprofilerenhProbesFLC2$shortLink # "https://biit.cs.ut.ee/gplink/l/Aw1-UW-RS_"




#### RESET_HypermethylatedProbesFL_(T3 with C1 samples as tumor)_HypomethylatedNormalC2(Silencers) ####

# Get C1 FL only samples
matchRNAseqC1 <- match(InfiniumClustLabels$ID[which(InfiniumClustLabels$Cluster == 1)][-c(1:10)],
                       colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized))
matchRNAseqC1 <- matchRNAseqC1[! is.na(matchRNAseqC1)]
colnames(RNAseqT3ProteinCode[, matchRNAseqC1])# double check


# Match samples between RNAseq and methylation - 104 tumor samples
matchRNAseqBetaC1 <- match(colnames(RNAseqT3ProteinCode[, matchRNAseqC1]),
                           colnames(BetaMatrix_T1))


# ENSURE identical colnames
identical(colnames(BetaMatrix_T1[, matchRNAseqBetaC1]), colnames(RNAseqT3ProteinCode[, matchRNAseqC1]))
# TRUE


# Run reset function on C1 tumor samples - enhancers - FL only
resetResultsSILFLC1againstC2 <- reset(normal.db = methNorSelOutputC2$normal.sil.probes,
                                      meth.tumor = BetaMatrix_T1[, matchRNAseqBetaC1],
                                      transcriptome = RNAseqT3ProteinCode[, matchRNAseqC1],
                                      methylation.event = c('sil'),
                                      FDR.permutation.no = 100, 
                                      seed = 100)
dim(resetResultsSILFLC1againstC2$FDR.res) # 1103    2
range(resetResultsSILFLC1againstC2$Score.report$Score, na.rm = TRUE) # 0.00333016 1.70897569

SILFLmeth.tumor.status.allC1againstC2 <- data.frame(resetResultsSILFLC1againstC2$meth.tumor.status.all, 
                                                    Gene = gsub(".*@", "", rownames(resetResultsSILFLC1againstC2$meth.tumor.status.all)),
                                                    Score = )
head(SILFLmeth.tumor.status.allC1againstC2)
dim(SILFLmeth.tumor.status.allC1againstC2) # 40718    58
dim(resetResultsSILFLC1againstC2$meth.tumor.all) # 40718    57
dim(resetResultsSILFLC1againstC2$transcriptome) # 40718    57
resetResultsSILFLC1againstC2$score.cutoff # 1.708976

genelist <- c("MME", "PPP1R14C")
testingGene <- genelist[1]
SILFLmeth.tumor.status.allC1againstC2[which(SILFLmeth.tumor.status.allC1againstC2$Gene == testingGene), ]
rowSums(SILFLmeth.tumor.status.allC1againstC2[which(SILFLmeth.tumor.status.allC1againstC2$Gene == testingGene), c(1:57)])



plotting <- function() {
  
  which(rownames(resetResultsSILFLC1againstC2$normal.meth) == "cg22001630@p1@MME") # probe with highest score of 1.7
  # If one probe
  if(length(which(gsub(".*@", "", rownames(resetResultsSILFLC1againstC2$normal.meth)) == testingGene)[5]) == 1) {
    testingDataFrame <- data.frame(Values = c(resetResultsSILFLC1againstC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC1againstC2$normal.meth)) == testingGene)[5], ],
                                              resetResultsSILFLC1againstC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC1againstC2$meth.tumor.all)) == testingGene)[5], ],
                                              resetResultsSILFLC1againstC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC1againstC2$transcriptome)) == testingGene)[5], ]),
                                   Type = c(rep("C2Methylation", 64), rep("C1Methylation", 57), rep("C1RNAseq", 57)))
    dim(testingDataFrame) # 247 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "C1RNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "C1RNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsSILFLC1againstC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC1againstC2$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsSILFLC1againstC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC1againstC2$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsSILFLC1againstC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC1againstC2$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    hyperProbes <- which(SILFLmeth.tumor.status.allC2[which(SILFLmeth.tumor.status.allC2$Gene == testingGene), ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsSILFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene), ],
                                                    resetResultsSILFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.all)) == testingGene), hyperProbes],
                                                    resetResultsSILFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC2$transcriptome)) == testingGene), hyperProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hyperProbes)), rep("TumorRNAseq", length(hyperProbes))))
    dim(testingDataFrameStatus) # 9  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsSILFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene), ]),
                                      MeanMethTumor = mean(resetResultsSILFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.all)) == testingGene), hyperProbes]),
                                      MeanTranscriptome = mean(resetResultsSILFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC2$transcriptome)) == testingGene), hyperProbes]))
    
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ]),
                                         t(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247   4
    
    
    testingDataFrame %>%  # highest score: cg01941671@p1@BMP3 methylation only
      melt(id.vars = "Type")  %>%
      filter(variable == "cg01941671.p1.BMP3") %>% # highest score: cg01941671@p1@BMP3
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg01941671.p1.BMP3") %>% # highest score: cg01941671.p1.BMP3
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
  }
}


resetScoreSILFLC1againstC2 <- data.frame(resetResultsSILFLC1againstC2$Score.report[which(! is.na(resetResultsSILFLC1againstC2$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSILFLC1againstC2)) {
  resetScoreSILFLC1againstC2[i, 6] <- resetResultsSILFLC1againstC2$FDR.res$FDR[which(resetScoreSILFLC1againstC2$Score[i] == resetResultsSILFLC1againstC2$FDR.res$observed.scores)]
}
head(resetScoreSILFLC1againstC2)
dim(resetScoreSILFLC1againstC2) #  1103    6
range(resetScoreSILFLC1againstC2$FDR, na.rm = TRUE) # 0.0700000 0.9888551
length(unique(resetScoreSILFLC1againstC2$Gene)) # 647
range(resetScoreSILFLC1againstC2$Score) # 0.00333016 1.70897569
resetScoreSILFLC1againstC2[which(resetScoreSILFLC1againstC2$Gene == testingGene), ]

# View(resetScoreSILFLC1againstC2)
# write.csv(resetScoreSILFLC1againstC2, file = "RESET_HypomethylatedProbes_FL_Silencers.csv")

# Plot No.Methylation.Events vs Score
resetScoreSILFLC1againstC2 %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot No.Methylation.Events vs Score
resetScoreSILFLC1againstC2 %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreSILFLC1againstC2 %>% 
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in C1")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreSILFLC1againstC2 %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in C1")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


# run with all genes 
gprofilerResetScoreSILFLC2All <- GProfilerAnalysis(GeneIDs = list(resetScoreSILFLC1againstC2$Gene),
                                                   Organism = "hsapiens",
                                                   OrderedQuery = TRUE,
                                                   PvalAlphaLevel = 0.01,
                                                   PositiveorNegFC = NA, # with respect to expression
                                                   ConditionName = "All",
                                                   ProduceImages = "Yes", 
                                                   PNGorPDF = "png")
gprofilerResetScoreSILFLC2All$shortLink
# "https://biit.cs.ut.ee/gplink/l/KBybBHDXSV"



#### RESET_HypomethylatedProbes_FL_(T3 with C1 samples as tumor)_HypermethylatedNormalC2(Enhancers) ####


# Get C1 FL only samples
matchRNAseqC1 <- match(InfiniumClustLabels$ID[which(InfiniumClustLabels$Cluster == 1)][-c(1:10)],
                       colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized))
matchRNAseqC1 <- matchRNAseqC1[! is.na(matchRNAseqC1)]
colnames(RNAseqT3ProteinCode[, matchRNAseqC1])# double check


# Match samples between RNAseq and methylation - 104 tumor samples
matchRNAseqBetaC1 <- match(colnames(RNAseqT3ProteinCode[, matchRNAseqC1]),
                           colnames(BetaMatrix_T1))


# ENSURE identical colnames
identical(colnames(BetaMatrix_T1[, matchRNAseqBetaC1]), colnames(RNAseqT3ProteinCode[, matchRNAseqC1]))
# TRUE


# Run reset function on C1 tumor samples - enhancers - FL only
resetResultsENHFLC1againstC2 <- reset(normal.db = methNorSelOutputC2$normal.enh.probes,
                                      meth.tumor = BetaMatrix_T1[, matchRNAseqBetaC1],
                                      transcriptome = RNAseqT3ProteinCode[, matchRNAseqC1],
                                      methylation.event = c('enh'),
                                      FDR.permutation.no = 100, 
                                      seed = 100)
dim(resetResultsENHFLC1againstC2$FDR.res) #  12   2
head(resetResultsENHFLC1againstC2$meth.tumor.status.all)
ENHFLmeth.tumor.status.allC1againstC2 <- data.frame(resetResultsENHFLC1againstC2$meth.tumor.status.all, Gene = gsub(".*@", "", rownames(resetResultsENHFLC1againstC2$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.allC1againstC2)
dim(ENHFLmeth.tumor.status.allC1againstC2) # 1715  58
dim(resetResultsENHFLC1againstC2$meth.tumor.all) # 1715  57
dim(resetResultsENHFLC1againstC2$transcriptome) # 1715  57
resetResultsENHFLC1againstC2$FDR.res
range(resetResultsENHFLC1againstC2$Score.report$Score, na.rm = TRUE) # 0.1814491 0.5201030



genelist <- c("AQP8")
testingGene <- genelist[1]
ENHFLmeth.tumor.status.allC1againstC2[which(ENHFLmeth.tumor.status.allC1againstC2$Gene == testingGene), ]
rowSums(ENHFLmeth.tumor.status.allC1againstC2[which(ENHFLmeth.tumor.status.allC1againstC2$Gene == testingGene), c(1:57)])


# If one probe
plotting <- function() {
  if(length(which(gsub(".*@", "", rownames(resetResultsENHFLC1againstC2$normal.meth)) == testingGene)) == 1) {
    
    
    testingDataFrame <- data.frame(Values = c(resetResultsENHFLC1againstC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1againstC2$normal.meth)) == testingGene), ],
                                              resetResultsENHFLC1againstC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1againstC2$meth.tumor.all)) == testingGene), ],
                                              resetResultsENHFLC1againstC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1againstC2$transcriptome)) == testingGene), ]),
                                   Type = c(rep("C2Methylation", 64), rep("C1Methylation", 57), rep("C1RNAseq", 57)))
    dim(testingDataFrame) # 178 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "C1RNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "C1RNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsENHFLC1againstC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1againstC2$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsENHFLC1againstC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1againstC2$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsENHFLC1againstC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1againstC2$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    # only looking at samples with a status marked as 1
    hypoProbes <- which(ENHFLmeth.tumor.status.allC1[which(ENHFLmeth.tumor.status.allC1$Gene == testingGene), ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ],
                                                    resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), hypoProbes],
                                                    resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) #  133   2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ]),
                                      MeanMethTumor = mean(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ]),
                                         t(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 250   4
    
    
    testingDataFrame %>%  # highest score: cg20802515.p1.SPO11 methylation only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    # look at 2nd probe with highest score
    hypoProbes <- which(ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.all$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene)[2], ],
                                                    resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene)[2], hypoProbes],
                                                    resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene)[2], hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) # 25  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene)[2], hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene)[2], hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  }
}

plottingRecRobert <- function() {
  
  RNAseqT3ProteinCodeNorm <- RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations[which(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations$type == "protein_coding"), ]
  dim(RNAseqT3ProteinCodeNorm) # 17190   136
  
  
  identical(rownames(resetResultsENHFLC1againstC2$meth.tumor.all),
            rownames(resetResultsENHFLC1againstC2$transcriptome))
  
  
  dim(resetResultsENHFLC1againstC2$meth.tumor.all) # 1715  57
  dim(resetResultsENHFLC1againstC2$transcriptome) # 1715  57
  
  
  MethDataUpdate <- data.frame(resetResultsENHFLC1againstC2$meth.tumor.all[- which(is.na(rowSums(resetResultsENHFLC1againstC2$transcriptome)) == TRUE), ])
  TranscriptDataUpdate <- data.frame(resetResultsENHFLC1againstC2$transcriptome[- which(is.na(rowSums(resetResultsENHFLC1againstC2$transcriptome)) == TRUE), ])
  
  # getting normalized counts
  matchNames <- match(sub(".*@", "", rownames(TranscriptDataUpdate)), RNAseqT3ProteinCodeNorm$name)
  matchColNames <- match(colnames(TranscriptDataUpdate), colnames(RNAseqT3ProteinCodeNorm))
  RNAseqT3ProteinCodeNormMatched <- RNAseqT3ProteinCodeNorm[matchNames, c(1:4, matchColNames) ]
  
  # dataframes - tumor
  rawDataFrame <- data.frame(Meth = MethDataUpdate, RNAseq = TranscriptDataUpdate)
  normDataFrame <- data.frame(Meth = MethDataUpdate,  RNAseq = RNAseqT3ProteinCodeNormMatched)
  
  # plot(MethDataUpdate$LY_FL_013_T1, TranscriptDataUpdate$LY_FL_013_T1)
  # plot(MethDataUpdate$LY_FL_013_T1, RNAseqT3ProteinCodeNormMatched$LY_FL_013_T1)
  
  # non-normalized counts
  rawPlot <- ggplot2::ggplot(rawDataFrame, aes(x = Meth.LY_FL_013_T1, y = RNAseq.LY_FL_013_T1)) + 
    geom_point() +
    theme_bw() + 
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"))
  
  # normalized counts
  normPlotLY_FL_013 <- ggplot2::ggplot(normDataFrame, aes(x = Meth.LY_FL_013_T1, y = RNAseq.LY_FL_013_T1)) + 
    geom_point() +
    theme_bw() + 
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"))
  
  normPlotLY_FL_032 <- ggplot2::ggplot(normDataFrame, aes(x = Meth.LY_FL_032_T1, y = RNAseq.LY_FL_032_T1)) + 
    geom_point() +
    theme_bw() + 
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"))
  
  normPlotLY_FL_038 <- ggplot2::ggplot(normDataFrame, aes(x = Meth.LY_FL_038_T1, y = RNAseq.LY_FL_038_T1)) + 
    geom_point() +
    theme_bw() + 
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"))
  
  normPlotLY_FL_046 <- ggplot2::ggplot(normDataFrame, aes(x = Meth.LY_FL_046_T1, y = RNAseq.LY_FL_046_T1)) + 
    geom_point() +
    theme_bw() + 
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"))
  
  library(magrittr)
  library(multipanelfigure)
  figure1 <- multipanelfigure::multi_panel_figure(columns = 2, rows = 2, panel_label_type = "upper-alpha")
  figure1 %<>%
    fill_panel(normPlotLY_FL_013, column = 1, row = 1) %<>%
    fill_panel(normPlotLY_FL_032, column = 2, row = 1) %<>%
    fill_panel(normPlotLY_FL_038, column = 1, row = 2) %<>%
    fill_panel(normPlotLY_FL_046, column = 2, row = 2)
  figure1
  
  
  # dataframes - Normal
  
  matchNames <- match(sub(".*@", "", rownames((methNorSelOutputC2$normal.enh.probes))), 
                          RNAseqT3ProteinCodeNorm$name)
  matchColNames <- match(colnames(methNorSelOutputC2$normal.enh.probes), colnames(RNAseqT3ProteinCodeNorm))
  RNAseqT3ProteinCodeNormMatchedNorm <- RNAseqT3ProteinCodeNorm[matchNames, c(1:4, matchColNames) ]
  
  normDataFrameC2normal <- data.frame(Meth = methNorSelOutputC2$normal.enh.probes,  RNAseq = RNAseqT3ProteinCodeNormMatchedNorm)
  
  
  # normalized counts
  normPlotLY_FL_008_T1 <- ggplot2::ggplot(normDataFrameC2normal, aes(x = Meth.LY_FL_008_T1, y = RNAseq.LY_FL_008_T1)) + 
    geom_point() +
    theme_bw() + 
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"))
  
  normPlotLY_FL_012_T1 <- ggplot2::ggplot(normDataFrameC2normal, aes(x = Meth.LY_FL_012_T1, y = RNAseq.LY_FL_012_T1)) + 
    geom_point() +
    theme_bw() + 
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"))
  
  normPlotLY_FL_018_T1 <- ggplot2::ggplot(normDataFrameC2normal, aes(x = Meth.LY_FL_018_T1, y = RNAseq.LY_FL_018_T1)) + 
    geom_point() +
    theme_bw() + 
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"))
  
  normPlotLY_FL_020_T1 <- ggplot2::ggplot(normDataFrameC2normal, aes(x = Meth.LY_FL_020_T1, y = RNAseq.LY_FL_020_T1)) + 
    geom_point() +
    theme_bw() + 
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"))
  
}


resetScoreENHFLC2againstC2 <- data.frame(resetResultsENHFLC1againstC2$Score.report[which(! is.na(resetResultsENHFLC1againstC2$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSILFLC2againstC2)) {
  resetScoreENHFLC2againstC2[i, 6] <- resetResultsENHFLC1againstC2$FDR.res$FDR[which(resetScoreENHFLC2againstC2$Score[i] == resetResultsENHFLC1againstC2$FDR.res$observed.scores)]
}
head(resetScoreENHFLC2againstC2)
dim(resetScoreENHFLC2againstC2) #  12  6
range(resetScoreENHFLC2againstC2$FDR) # 0.659 1.920
length(unique(resetScoreENHFLC2againstC2$Gene)) # 12
range(resetScoreENHFLC2againstC2$Score) # 0.1814491 0.5201030
resetScoreENHFLC2againstC2[which(resetScoreENHFLC2againstC2$Gene == testingGene), ]
# ProbID Gene Promoter.no    Score No.Methylation.Events  FDR
# cg19984833@p1@AQP8 CG19984833 AQP8          P1 0.520103                     8 1.92


# View(resetScoreENHFLC2againstC2)
# write.csv(resetScoreENHFLC2againstC2, file = "RESET_HypermethylatedProbes_FL_Enhancers.csv")

# Plot No.Methylation.Events vs Score
resetScoreENHFLC2againstC2 %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()

# Plot No.Methylation.Events vs Score
resetScoreENHFLC2againstC2 %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreENHFLC2againstC2 %>%  
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL - C1")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreENHFLC2againstC2 %>%  
  #filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


# GProfilerAnalysis
# get gene list of those with FDR < 0.5
resetScoreENHFLC1Genes <- resetScoreENHFLC2againstC2 %>%  
  dplyr::filter(FDR < 0.5) %>% 
  dplyr::select(Gene)

gprofilerResetScoreENHFLC1 <- GProfilerAnalysis(GeneIDs = list(resetScoreENHFLC1Genes$Gene),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = NA, # with respect to expression
                                                ConditionName = "All",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
gprofilerResetScoreENHFLC1$shortLink
# "https://biit.cs.ut.ee/gplink/l/PC-KEBe5Rj"

# run with all genes 
gprofilerResetScoreENHFLC1All <- GProfilerAnalysis(GeneIDs = list(resetScoreENHFLC1$Gene),
                                                   Organism = "hsapiens",
                                                   OrderedQuery = TRUE,
                                                   PvalAlphaLevel = 0.01,
                                                   PositiveorNegFC = NA, # with respect to expression
                                                   ConditionName = "All",
                                                   ProduceImages = "Yes", 
                                                   PNGorPDF = "png")
gprofilerResetScoreENHFLC1All$shortLink
# "https://biit.cs.ut.ee/gplink/l/MABDo3KOQG"






##### SNF: RESET_HypomethylatedProbes_FL_(T3 with C1 samples as tumor)_HypermethylatedNormalC2(Enhancers) ####
# bind tumor + normal samples
identical(rownames(ENHFLmeth.tumor.status.allC1againstC2), 
          rownames(methNorSelOutputC2$normal.enh.probes)) # TRUE

dim(methNorSelOutputC2$normal.enh.probes) # 1715   64
 
EventMatrix <- cbind(ENHFLmeth.tumor.status.allC1againstC2[, 1:57], methNorSelOutputC2$normal.enh.probes)
dim(EventMatrix) # 1715  121
EventMatrix[, c(58:121)] <- 0




TumorPurity <- readRDS(file = paste0(MethylationDirPath, "/Purity_281probes_10Jan2020.rds"))
dim(TumorPurity) # 170   1

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

# cluster only RNAseq protein coding genes
RNAseqQC18June2020T3Samples132Protein <- RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations[which(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsWithAnnotations$type == "protein_coding"), ]
dim(RNAseqQC18June2020T3Samples132Protein) # 17190   135
RNAseqQC18June2020T3Samples132Protein <- as.matrix(RNAseqQC18June2020T3Samples132Protein[, c(4:135)])
dim(RNAseqQC18June2020T3Samples132Protein) # 17190   132


TargetedDNAseqDirPath <- "~/Desktop/UHN/FLOMICS/TargetedSequencing"
mutMergedT1Robert <- read.csv(paste0(TargetedDNAseqDirPath, "/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl_25Aug2020.csv"), row.names = 1)
dim(mutMergedT1Robert) # 55 138
colnames(mutMergedT1Robert)
class(mutMergedT1Robert)

# Fixed in mut.merged.df.T1.poor.cov.excl_25Aug2020 file 
# mutMergedT1Robert[, which(colnames(mutMergedT1Robert) == "LY_FL_179_T1")] # NA values
# mutMergedT1Robert[, which(colnames(mutMergedT1Robert) == "LY_FL_181_T1")] # NA values
# which(is.na(colSums(mutMergedT1Robert)) == TRUE) # LY_FL_179_T1 50; LY_FL_181_T1  51 
# mutMergedT1Robert[, c(50, 51)] <- 0
which(is.na(colSums(mutMergedT1Robert)) == TRUE) # named integer(0)




matchRNAseqMutation <- match(colnames(RNAseqQC18June2020T3Samples132Protein), colnames(mutMergedT1Robert))
dim(mutMergedT1Robert[, matchRNAseqMutation[! is.na(matchRNAseqMutation)]]) # 55 101

mutMergedT1Robert101 <- mutMergedT1Robert[, matchRNAseqMutation[! is.na(matchRNAseqMutation)]]
dim(mutMergedT1Robert101) # 55 101

matchDNAseqBeta <- match(colnames(mutMergedT1Robert101), colnames(EventMatrix))

matchDNAseqRNAseq <- match(colnames(mutMergedT1Robert101), colnames(RNAseqQC18June2020T3Samples132Protein))

matchDNAseqClinicalFile <- match(colnames(mutMergedT1Robert101), ClinicalFile_T1$SAMPLE_ID)


# Combine mutations with BCL2 and BCL6 translocations
BCL2BCL6BreakapartPredictions <- data.table::fread(file = paste0(MethylationDirPath, "/2020-07-05_BCL2_BCL6_rearrangements.txt"))
BCL2BCL6BreakapartPredictions <- data.frame(BCL2BCL6BreakapartPredictions)
dim(BCL2BCL6BreakapartPredictions) # 183   3
rownames(BCL2BCL6BreakapartPredictions) <- BCL2BCL6BreakapartPredictions$External_ID
BCL2BCL6BreakapartPredictions <- BCL2BCL6BreakapartPredictions[, -1]
matchMutationRearrangements <- match(colnames(mutMergedT1Robert101), rownames(BCL2BCL6BreakapartPredictions))
rearrangementsT1Robert108 <- BCL2BCL6BreakapartPredictions[matchMutationRearrangements, ]
identical(colnames(t(rearrangementsT1Robert108)), colnames(mutMergedT1Robert101)) #TRUE
mutRearrT1Robert101 <- rbind(mutMergedT1Robert101, t(rearrangementsT1Robert108))
dim(mutRearrT1Robert101) # 57 101

which(is.na(colSums(mutRearrT1Robert101)) == TRUE)
# LY_FL_244_T1 LY_FL_253_T1 LY_FL_264_T1 
# 62           71           82
mutRearrT1Robert101[56:57, which(is.na(colSums(mutRearrT1Robert101)) == TRUE)]
# manual adjustment
mutRearrT1Robert101[56:57, 62] <- c(0, 0)
mutRearrT1Robert101[56:57, 71] <- c(1, 0)
mutRearrT1Robert101[56:57, 82] <- c(0, 0)

# check
mutRearrT1Robert101[56:57, c(62, 71, 82)] 





ClusteringSNFMethEventRESET101 <- SNFClustering(RNAseqCountMatrix = as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq]), 
                                                 MvalueMatrix = EventMatrix[, matchDNAseqBeta],
                                                 BetaMatrix = EventMatrix[, matchDNAseqBeta], 
                                                 TargetedSeqMatrix = mutRearrT1Robert101,
                                                 MethylationAnnotationFile = EPICAnnotationFile,
                                                 RNAseqAnnotationFile = ENSMBLid, 
                                                 QCMatrix = RNAseqQC18June2020T3Samples132$QCMatrixMatchedSampleFiltered[matchDNAseqRNAseq, ],
                                                 TumorPurity = as.matrix(TumorPurity[matchDNAseqTumorPurity, ], ncol = 1),
                                                 ClinicalFile = ClinicalFile_T1[matchDNAseqClinicalFile, ], 
                                                 SurvivalFile = SurvivalFile,
                                                 NumberOfClusters = 2,
                                                 ImageName = paste0("Image_MethEvents", date()),
                                                 PNGorPDF = "png",
                                                 ProduceImages = "Yes") 

SNFClusterLabels$SampleID[SNFClusterLabels$SampleID %in% colnames(ENHFLmeth.tumor.status.allC1againstC2[, 1:57])]
table(SNFClusterLabels$Cluster[SNFClusterLabels$SampleID %in% colnames(ENHFLmeth.tumor.status.allC1againstC2[, 1:57])])
#  1  2 
# 5 44 


ClusteringSNFMethEventRESET101 <- SNFClusterLabels$Cluster

table(ClusteringSNFMethEventRESET101)

InfiniumClustLabels <- readRDS(file = paste0(MethylationDirPath, "/InfiniumClustLabels2.RDS"))
InfSNFMethEvents <- match(colnames(EventMatrix[, matchDNAseqBeta]), InfiniumClustLabels$ID)
identical(colnames(EventMatrix[, matchDNAseqBeta]), unfactor(InfiniumClustLabels$ID[InfSNFMethEvents]))
# TRUE
table(ClusteringSNFMethEventRESET101, InfiniumClustLabels$Cluster[InfSNFMethEvents])
chisq.test(table(ClusteringSNFMethEventRESET101, InfiniumClustLabels$Cluster[InfSNF5000SD]))



ClinicalFile_T1SNFCluster <- data.frame(ClinicalFile_T1[matchDNAseqClinicalFile, ], "CLUSTER" = ClusteringSNFMethEventRESET101)
matchDNAseqBetaMat_T1 <- match(colnames(as.matrix(RNAseqQC18June2020T3Samples132Protein[, matchDNAseqRNAseq])), colnames(BetaMatrix_T1))

DifferentialMethylationProbesSNF <- DifferentialMethylation(ClinicalFile = ClinicalFile_T1SNFCluster, 
                                                            MvalueMatrix = MvalueMatrix_T1[, matchDNAseqBetaMat_T1], 
                                                            BetaMatrix = BetaMatrix_T1[, matchDNAseqBetaMat_T1],
                                                            ProbeOrGene = "Probe", 
                                                            ContrastColumnName = "CLUSTER", 
                                                            RGChannelSet = NA, 
                                                            SampleSheet = NA, 
                                                            ProduceImages = "Yes", 
                                                            PNGorPDF = "png")

MethylationDensityOutput <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                   ClusterLabels = ClusteringSNFMethEventRESET101, 
                                                   PlotWithinAnnotationCategories = "No", 
                                                   ClinicalCategoryToVisualize = "CLUSTER", 
                                                   BetaMatrix = BetaMatrix_T1[, matchDNAseqBetaMat_T1], 
                                                   AnnotationFile = EPICAnnotationFile, 
                                                   ClinicalFile = ClinicalFile_T1SNFCluster, 
                                                   SampleSheet = NA, 
                                                   FigureGenerate = "Yes", 
                                                   ImageName = paste0("MethylationPlot_SNFNoMeth", ImageName), 
                                                   PNGorPDF = "png")


###########################################################################################
#### REDO Run first function from RESET with RLN samples as (normal x 10) samples ####
methNorSelOutput10 <- methNorSel(normal.mtx = BetaMatrix_T1[matchTSSProbes, rep((166:170), 2)],
                               probe.list = promoter.probes.list)
# this function first extract the data regarding the selected probes, and later based on the beta-values prepare two seperate normal-sample data-sets:
## 1. hypermethylated probes in normal condition: Enh.probes
## 2. hypomethylated probes in normal condition: Sil.probes
names(methNorSelOutput10) # "normal.sil.probes" "normal.enh.probes"
typeof(methNorSelOutput10)
dim(methNorSelOutput10$normal.sil.probes) # 45250 (45229; difference of 21)     5
class(methNorSelOutput10$normal.sil.probes)
normal.sil.probes10 <- data.frame(methNorSelOutput10$normal.sil.probes, 
                                  rowMedians = rowMedians(methNorSelOutput10$normal.sil.probes))
ggplot(normal.sil.probes10, aes(x = c(1:45250), y = rowMedians)) + # plot rowMeans of sil.probes
  geom_point(size = 2) +
  labs(x = "probe", y = "Median Beta Value") +
  theme_bw() + 
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  scale_color_manual(values = c("#d6604d", "#66bd63", "#4575b4")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))
silProbesFL10 <- gsub(".*@", "", rownames(methNorSelOutput10$normal.sil.probes))
length(silProbesFL10) # 45250 (45229)
# write.csv(silProbesFL, file = "RESET_HypermethylatedProbesNormalCondition_Enhancers.csv")
source("26_Gprofiler.R")
gprofilersilProbesFL <- GProfilerAnalysis(GeneIDs = list(silProbesFL10),
                                          Organism = "hsapiens",
                                          OrderedQuery = TRUE,
                                          PvalAlphaLevel = 0.01,
                                          PositiveorNegFC = NA, # with respect to expression
                                          ConditionName = "All",
                                          ProduceImages = "Yes", 
                                          PNGorPDF = "png")
gprofilersilProbesFL$shortLink  # "https://biit.cs.ut.ee/gplink/l/ZLxq3tvKRC"



enhProbesFL10 <- gsub(".*@", "", rownames(methNorSelOutput10$normal.enh.probes))
length(enhProbesFL10) # 2838 (2798; difference of 40) 
# write.csv(enhProbesFL, file = "RESET_HypomethylatedProbesNormalCondition_Silencers.csv")
normal.enh.probes10 <- data.frame(methNorSelOutput10$normal.enh.probes, 
                                  rowMedians = rowMedians(methNorSelOutput10$normal.enh.probes))
ggplot(normal.enh.probes10, aes(x = c(1:2838), y = rowMedians)) + # plot rowMeans of enh.probes
  geom_point(size = 2) +
  labs(x = "probe", y = "Median Beta Value") +
  theme_bw() + 
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  scale_color_manual(values = c("#d6604d", "#66bd63", "#4575b4")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))
enhProbesFL10 <- gsub(".*@", "", rownames(methNorSelOutput10$normal.enh.probes))

gprofilerenhProbesFL <- GProfilerAnalysis(GeneIDs = list(enhProbesFL10),
                                          Organism = "hsapiens",
                                          OrderedQuery = TRUE,
                                          PvalAlphaLevel = 0.01,
                                          PositiveorNegFC = NA, # with respect to expression
                                          ConditionName = "All",
                                          ProduceImages = "Yes", 
                                          PNGorPDF = "png")
gprofilerenhProbesFL$shortLink # "https://biit.cs.ut.ee/gplink/l/Aw1-UW-RS_"




#### RESET_HypomethylatedProbes_FL_(T3 n = 131-10 = 121)_HypermethylatedNormal(Enhancers; normal x 10) ####

# Run reset function on all tumor samples - enhancers - FL only
resetResultsENHFLonly10 <- reset(normal.db = methNorSelOutput10$normal.enh.probes,
                               meth.tumor = BetaMatrix_T1[, matchRNAseqBeta[11:131]],
                               transcriptome = RNAseqT3ProteinCode[, c(11:131)],
                               methylation.event = c('enh'),
                               FDR.permutation.no = 100, 
                               seed = 100)

dim(resetResultsENHFLonly10$FDR.res) #  591 (588)   2
head(resetResultsENHFLonly10$meth.tumor.status.all)
ENHFLmeth.tumor.status.all10 <- data.frame(resetResultsENHFLonly10$meth.tumor.status.all, Gene = gsub(".*@", "", rownames(resetResultsENHFLonly10$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.all10)
dim(ENHFLmeth.tumor.status.all10) # 2838 (2798)  122
dim(resetResultsENHFLonly10$meth.tumor.all) # 2838 (2798)  121
dim(resetResultsENHFLonly10$transcriptome) # 2838 (2798)  121

genelist <- c("SPO11", "FCRLB")
testingGene <- genelist[1]
ENHFLmeth.tumor.status.all10[which(ENHFLmeth.tumor.status.all10$Gene == testingGene), ]
rowSums(ENHFLmeth.tumor.status.all10[which(ENHFLmeth.tumor.status.all10$Gene == testingGene), c(1:121)])



plotting <- function() {
  # If one probe
  if(length(which(gsub(".*@", "", rownames(resetResultsENHFLonly10$normal.meth)) == testingGene)) == 1) {
    
    
    testingDataFrame <- data.frame(Values = c(resetResultsENHFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$normal.meth)) == testingGene), ],
                                              resetResultsENHFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$meth.tumor.all)) == testingGene), ],
                                              resetResultsENHFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$transcriptome)) == testingGene), ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsENHFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsENHFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsENHFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    # only looking at samples with a status marked as 1
    hypoProbes <- which(ENHFLmeth.tumor.status.all10[which(ENHFLmeth.tumor.status.all10$Gene == testingGene), ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$normal.meth)) == testingGene), ],
                                                    resetResultsENHFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$meth.tumor.all)) == testingGene), hypoProbes],
                                                    resetResultsENHFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$transcriptome)) == testingGene), hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) # 247 2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$normal.meth)) == testingGene), ]),
                                      MeanMethTumor = mean(resetResultsENHFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$meth.tumor.all)) == testingGene), hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$transcriptome)) == testingGene), hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsENHFLonly10$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsENHFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$normal.meth)) == testingGene), ]),
                                         t(resetResultsENHFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsENHFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 250   4
    
    
    testingDataFrame %>%  # highest score: cg20802515.p1.SPO11 methylation only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsENHFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsENHFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsENHFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    # look at 2nd probe with highest score
    hypoProbes <- which(ENHFLmeth.tumor.status.all10[which(ENHFLmeth.tumor.status.all10$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$normal.meth)) == testingGene)[2], ],
                                                    resetResultsENHFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$meth.tumor.all)) == testingGene)[2], hypoProbes],
                                                    resetResultsENHFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$transcriptome)) == testingGene)[2], hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) # 25  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsENHFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$meth.tumor.all)) == testingGene)[2], hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLonly10$transcriptome)) == testingGene)[2], hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  }
}



resetScoreENHFLonly10 <- data.frame(resetResultsENHFLonly10$Score.report[which(! is.na(resetResultsENHFLonly10$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFLonly10)) {
  resetScoreENHFLonly10[i, 6] <- resetResultsENHFLonly10$FDR.res$FDR[which(resetScoreENHFLonly10$Score[i] == resetResultsENHFLonly10$FDR.res$observed.scores)]
}
head(resetScoreENHFLonly10)
dim(resetScoreENHFLonly10) # 591 (588) 6
length(unique(resetScoreENHFLonly10$Gene)) # 489
range(resetScoreENHFLonly10$Score) # 0.004389081 1.150200828 (before 0.004389081 1.150200828)
range(resetScoreENHFLonly10$FDR) # 0.0450000 0.9588983 (0.0540000 0.9513823)
resetScoreENHFLonly10[which(resetScoreENHFLonly10$Gene == testingGene), ]


# View(resetScoreENHFLonly10)
# write.csv(resetScoreENHFLonly10, file = "RESET_HypermethylatedProbes_FL_Enhancers.csv")

# Plot No.Methylation.Events vs Score
resetScoreENHFLonly10 %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()

# Plot No.Methylation.Events vs Score
resetScoreENHFLonly10 %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreENHFLonly10 %>%  
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreENHFLonly10 %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))





# plotting RNAseq vs methylation value - SPO11, FCRLB
# get the probe IDs
genelist <- c("SPO11", "FCRLB")
testingGene <- genelist[2]
probesToPlot <- resetScoreENHFLonly10[which(resetScoreENHFLonly10$Gene == testingGene), ]$ProbID
which(rownames(BetaMatrix_T1[, matchRNAseqBeta[11:131]]) == tolower(probesToPlot))
# 35143 54921 78854

# Get RNAseq counts
which(rownames(RNAseqT3ProteinCode) == testingGene) 
# 679

# Get normalized RNAseq values
RNAseqT3ProteinCodeNormalized <- RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations[which(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations$type == "protein_coding"), ]
dim(RNAseqT3ProteinCodeNormalized) # 17190   136
rownames(RNAseqT3ProteinCodeNormalized) <- RNAseqT3ProteinCodeNormalized$name
RNAseqT3ProteinCodeNormalized <- RNAseqT3ProteinCodeNormalized[, c(15:135)]
dim(RNAseqT3ProteinCodeNormalized) # 17190   121
range(RNAseqT3ProteinCodeNormalized) # -2.72230 14.11687

# getting the events
which(ENHFLmeth.tumor.status.all10$Gene == testingGene)
# 937 939 941

GeneToPlot <- data.frame(Beta = (BetaMatrix_T1[which(rownames(BetaMatrix_T1[, matchRNAseqBeta[11:131]]) == tolower(probesToPlot)), matchRNAseqBeta[11:131]]),
                         Event = t(ENHFLmeth.tumor.status.all10[which(ENHFLmeth.tumor.status.all10$Gene == testingGene), c(1:121)]),
                         RNAseqCounts = t(RNAseqT3ProteinCode[which(rownames(RNAseqT3ProteinCode) == testingGene), c(11:131)]),
                         RNAseqNormalized = t(RNAseqT3ProteinCodeNormalized[which(rownames(RNAseqT3ProteinCode) == testingGene), ]))

ggplot2::ggplot(GeneToPlot, aes(x = Beta, y = FCRLB.1, 
                                color = factor(cg23641748.p17.FCRLB))) +
  geom_point(size = 2)+
  labs(title = "cg23641748 - FCRLB",
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







#### RESET_HypermethylatedProbesFL_(T3 n = 131-10 = 121)_HypomethylatedNormal(Silencers; normal x 10) ####

# Match samples between RNAseq and methylation - 131 tumor samples
matchRNAseqBeta <- match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized),
                         colnames(BetaMatrix_T1))

# Run reset function on all tumor samples - silencers - FL only
resetResultsSILFLonly10 <- reset(normal.db = methNorSelOutput10$normal.sil.probes,
                               meth.tumor = BetaMatrix_T1[, matchRNAseqBeta[11:131]],
                               transcriptome = RNAseqT3ProteinCode[, c(11:131)],
                               methylation.event = c('sil'),
                               FDR.permutation.no = 100, 
                               seed = 100)
dim(resetResultsSILFLonly10$FDR.res) #  6645 (6452)   2
SILFLmeth.tumor.status.all10 <- data.frame(resetResultsSILFLonly10$meth.tumor.status.all, Gene = gsub(".*@", "", rownames(resetResultsSILFLonly10$meth.tumor.status.all)))
head(SILFLmeth.tumor.status.all10)
dim(SILFLmeth.tumor.status.all10) #  45250 (45229)   122
dim(resetResultsSILFLonly10$meth.tumor.all) # 45250 (45229)  121
dim(resetResultsSILFLonly10$transcriptome) # 45250 (45229)  121

genelist <- c("BMP3", "MME", "HLTF", "CAMK1")
testingGene <- genelist[1]
SILFLmeth.tumor.status.all10[which(SILFLmeth.tumor.status.all10$Gene == testingGene), ]
rowSums(SILFLmeth.tumor.status.all10[which(SILFLmeth.tumor.status.all10$Gene == testingGene), c(1:121)])




plotting <- function() {
  # If one probe
  if(length(which(gsub(".*@", "", rownames(resetResultsSILFLonly10$normal.meth)) == testingGene)) == 1) {
    testingDataFrame <- data.frame(Values = c(resetResultsSILFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$normal.meth)) == testingGene), ],
                                              resetResultsSILFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$meth.tumor.all)) == testingGene), ],
                                              resetResultsSILFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$transcriptome)) == testingGene), ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsSILFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsSILFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsSILFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    hyperProbes <- which(SILFLmeth.tumor.status.all10[which(SILFLmeth.tumor.status.all10$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsSILFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$normal.meth)) == testingGene)[2], ],
                                                    resetResultsSILFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$meth.tumor.all)) == testingGene)[2], hyperProbes],
                                                    resetResultsSILFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$transcriptome)) == testingGene)[2], hyperProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hyperProbes)), rep("TumorRNAseq", length(hyperProbes))))
    dim(testingDataFrameStatus) # 33  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsSILFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsSILFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$meth.tumor.all)) == testingGene)[2], hyperProbes]),
                                      MeanTranscriptome = mean(resetResultsSILFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$transcriptome)) == testingGene)[2], hyperProbes]))
    
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsENHFLonly$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsSILFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$normal.meth)) == testingGene), ]),
                                         t(resetResultsSILFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsSILFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247   4
    
    
    testingDataFrame %>%  # highest score: cg01941671@p1@BMP3 methylation only
      melt(id.vars = "Type")  %>%
      filter(variable == "cg01941671.p1.BMP3") %>% # highest score: cg01941671@p1@BMP3
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg01941671.p1.BMP3") %>% # highest score: cg01941671.p1.BMP3
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsSILFLonly10$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsSILFLonly10$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsSILFLonly10$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLonly10$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
  }
}

resetScoreSILFLonly10 <- data.frame(resetResultsSILFLonly10$Score.report[which(! is.na(resetResultsSILFLonly10$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFLonly)) {
  resetScoreSILFLonly10[i, 6] <- resetResultsSILFLonly10$FDR.res$FDR[which(resetScoreSILFLonly10$Score[i] == resetResultsSILFLonly10$FDR.res$observed.scores)]
}
head(resetScoreSILFLonly10)
dim(resetScoreSILFLonly10) #  6645 (6452)  6
length(unique(resetScoreSILFLonly10$Gene)) # 3994 (3848)
range(resetScoreSILFLonly10$Score) # 0.003460362 1.267243664 (0.003460362 1.282929804)
range(resetScoreSILFLonly10$FDR, na.rm = TRUE) # 0.3637097 0.7106357 (0.3518487 0.7150675)
resetScoreSILFLonly10[which(resetScoreSILFLonly10$Gene == testingGene), ]

# View(resetScoreSILFLonly10)
# write.csv(resetScoreSILFLonly10, file = "RESET_HypomethylatedProbes_FL_Silencers.csv")

# Plot No.Methylation.Events vs Score
resetScoreSILFLonly10 %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()


# Plot No.Methylation.Events vs Score
resetScoreSILFLonly10 %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreSILFLonly10 %>% 
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreSILFLonly10 %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# plotting RNAseq vs methylation value - BMP3, FCRLB
# get the probe IDs
genelist <- c("BMP3", "MME")
testingGene <- genelist[2]
probesToPlot <- resetScoreSILFLonly10[which(resetScoreSILFLonly10$Gene == testingGene), ]$ProbID
which(rownames(BetaMatrix_T1[, matchRNAseqBeta[11:131]]) == tolower(probesToPlot)[8])
# 218098 125779 1060 97335 86432 63883 32900 99749

# Get RNAseq counts
which(rownames(RNAseqT3ProteinCode) == testingGene) 
# 679

# Get normalized RNAseq values
RNAseqT3ProteinCodeNormalized <- RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations[which(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalizedWithAnnotations$type == "protein_coding"), ]
dim(RNAseqT3ProteinCodeNormalized) # 17190   136
rownames(RNAseqT3ProteinCodeNormalized) <- RNAseqT3ProteinCodeNormalized$name
RNAseqT3ProteinCodeNormalized <- RNAseqT3ProteinCodeNormalized[, c(15:135)]
dim(RNAseqT3ProteinCodeNormalized) # 17190   121
range(RNAseqT3ProteinCodeNormalized) # -2.72230 14.11687

# getting the events
which(SILFLmeth.tumor.status.all10$Gene == testingGene)
# 937 939 941

GeneToPlot <- data.frame(Beta = t(BetaMatrix_T1[c(218098, 125779, 1060, 97335, 86432, 63883, 32900, 99749), matchRNAseqBeta[11:131]]),
                         Event = t(SILFLmeth.tumor.status.all10[which(SILFLmeth.tumor.status.all10$Gene == testingGene), c(1:121)]),
                         RNAseqCounts = t(RNAseqT3ProteinCode[which(rownames(RNAseqT3ProteinCode) == testingGene), c(11:131)]),
                         RNAseqNormalized = t(RNAseqT3ProteinCodeNormalized[which(rownames(RNAseqT3ProteinCode) == testingGene), ]))

ggplot2::ggplot(GeneToPlot, aes(x = Beta.cg07862554, y = MME.1, 
                                color = factor(Event.cg07862554.p8.MME))) +
  geom_point(size = 2)+
  labs(title = "cg07862554 - MME",
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



###########################################################################################
##### Back to edgeR (non-log transformed) with SNF clusters ################################################
#### select TSS probes ####

AnnotationFileEdited <- EPICAnnotationFile
matchBetaProbes <- match(rownames(BetaMatrix_T1), AnnotationFileEdited$V1)
AnnotationFileEditedBetaMatrix <- AnnotationFileEdited[matchBetaProbes, ]
dim(AnnotationFileEditedBetaMatrix) # 595564     47
AnnotationFileEditedBetaMatrix$UCSC_RefGene_Group <- sub("\\;.*", "", AnnotationFileEditedBetaMatrix$UCSC_RefGene_Group)
promoterProbes <- AnnotationFileEditedBetaMatrix %>% filter(substr(UCSC_RefGene_Group, 1, 3) == "TSS") %>% pull("V1")
matchTSSProbes <- match(promoterProbes, rownames(BetaMatrix_T1))
length(matchTSSProbes) # 114897


#### select RNAseq samples (T3) ####

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

RNAseqProtCodT3MatName <- RNAseqProtCodT3Mat
rownames(RNAseqProtCodT3MatName) <- RNAseqProtCodT3Mat$name
rownames(RNAseqProtCodT3Mat) <- RNAseqProtCodT3Mat$ENSEMBLID


# match to 101 samples from SNF
matchtoSNF101 <- match(AllClustersLab$ID[- which(is.na(AllClustersLab$SNFClust) == TRUE)], colnames(RNAseqProtCodT3MatName))
RNAseqProtCodT3MatNameSNF <- RNAseqProtCodT3MatName[, matchtoSNF101]
dim(RNAseqProtCodT3MatNameSNF) # 17190   101
range(RNAseqProtCodT3MatNameSNF) # 0.00 17766.33

matchtoSNF101 <- match(AllClustersLab$ID[- which(is.na(AllClustersLab$SNFClust) == TRUE)], colnames(BetaMatrix_T1))
BetaMatrixT1SNF <- BetaMatrix_T1[, matchtoSNF101]
dim(BetaMatrixT1SNF) # 595564    101



#### Run first function from RESET with RLN samples as normal ####
set.seed(1234)
methNorSelOutput <- methNorSel(normal.mtx = BetaMatrix_T1[matchTSSProbes, c(166:170)],
                               probe.list = promoter.probes.list)
# this function first extract the data regarding the selected probes, and later based on the beta-values prepare two seperate normal-sample data-sets:
## 1. hypermethylated probes in normal condition: Enh.probes
## 2. hypomethylated probes in normal condition: Sil.probes
names(methNorSelOutput) # "normal.sil.probes" "normal.enh.probes"
typeof(methNorSelOutput)
dim(methNorSelOutput$normal.sil.probes) # 45229     5
dim(methNorSelOutput$normal.enh.probes) # 2798    5
class(methNorSelOutput$normal.sil.probes)
normal.sil.probes <- data.frame(methNorSelOutput$normal.sil.probes, rowMedians = rowMedians(methNorSelOutput$normal.sil.probes))

# isolate genes
length(unique(gsub(".*@", "", rownames(methNorSelOutput$normal.sil.probes)))) # 12704
length(unique(gsub(".*@", "", rownames(methNorSelOutput$normal.enh.probes)))) # 1729

ggplot(normal.sil.probes, aes(x = c(1:45229), y = rowMedians)) + # plot rowMeans of sil.probes
  geom_point(size = 2) +
  labs(x = "probe", y = "Median Beta Value") +
  theme_bw() + 
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  scale_color_manual(values = c("#d6604d", "#66bd63", "#4575b4")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))
silProbesFL <- gsub(".*@", "", rownames(methNorSelOutput$normal.sil.probes))
length(silProbesFL) # 45229
# write.csv(silProbesFL, file = "RESET_HypermethylatedProbesNormalCondition_Enhancers.csv")
source("26_Gprofiler.R")
gprofilersilProbesFL <- GProfilerAnalysis(GeneIDs = list(silProbesFL),
                                          Organism = "hsapiens",
                                          OrderedQuery = TRUE,
                                          PvalAlphaLevel = 0.01,
                                          PositiveorNegFC = NA, # with respect to expression
                                          ConditionName = "All",
                                          ProduceImages = "Yes", 
                                          PNGorPDF = "png")
gprofilersilProbesFL$shortLink  # "https://biit.cs.ut.ee/gplink/l/ZLxq3tvKRC"



enhProbesFL <- gsub(".*@", "", rownames(methNorSelOutput$normal.enh.probes))
length(enhProbesFL) # 2798
# write.csv(enhProbesFL, file = "RESET_HypomethylatedProbesNormalCondition_Silencers.csv")
normal.enh.probes <- data.frame(methNorSelOutput$normal.enh.probes, rowMedians = rowMedians(methNorSelOutput$normal.enh.probes))
ggplot(normal.enh.probes, aes(x = c(1:2798), y = rowMedians)) + # plot rowMeans of enh.probes
  geom_point(size = 2) +
  labs(x = "probe", y = "Median Beta Value") +
  theme_bw() + 
  theme(text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  scale_color_manual(values = c("#d6604d", "#66bd63", "#4575b4")) +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.background = element_rect(colour = "black", size=1.5),  
        axis.title =  element_text(face = "bold"))
silProbesFL <- gsub(".*@", "", rownames(methNorSelOutput$normal.sil.probes))

gprofilerenhProbesFL <- GProfilerAnalysis(GeneIDs = list(enhProbesFL),
                                          Organism = "hsapiens",
                                          OrderedQuery = TRUE,
                                          PvalAlphaLevel = 0.01,
                                          PositiveorNegFC = NA, # with respect to expression
                                          ConditionName = "All",
                                          ProduceImages = "Yes", 
                                          PNGorPDF = "png")
gprofilerenhProbesFL$shortLink # "https://biit.cs.ut.ee/gplink/l/Aw1-UW-RS_"





#### RESET_HypomethylatedProbes_FL_(T3 with C1 and C2 SNF samples as tumor)_HypermethylatedNormal(Enhancers) - SNF clusters + n=121 (edgeR normalized counts) ####


# Run reset function on 101 tumor samples - enhancers - FL only
set.seed(1234)
resetResultsENHFLAll101 <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                             meth.tumor = BetaMatrixT1SNF,
                             transcriptome = RNAseqProtCodT3MatNameSNF,
                             methylation.event = c('enh'),
                             FDR.permutation.no = 1000, 
                             seed = 100)
dim(resetResultsENHFLAll101$FDR.res) #  493   2
range(resetResultsENHFLAll101$FDR.res) # 0.0133416 1.1358012
head(resetResultsENHFLAll101$meth.tumor.status.all)
ENHFLmeth.tumor.status.all <- data.frame(resetResultsENHFLAll101$meth.tumor.status.all, 
                                           Gene = gsub(".*@", "", rownames(resetResultsENHFLAll101$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.all)
dim(ENHFLmeth.tumor.status.all) # 2798  102
dim(resetResultsENHFLAll101$meth.tumor.all) # 2798  101
dim(resetResultsENHFLAll101$transcriptome) # 2798  101
resetResultsENHFLAll101$FDR.res
range(resetResultsENHFLAll101$Score.report$Score, na.rm = TRUE) # 0.0133416 1.1358012
# saveRDS(resetResultsENHFLAll101, file = "40_resetResultsENHFLonlyedgeRnonLog_TumorHypo.rds")
resetResultsENHFLAll101 <- readRDS("40_resetResultsENHFLonlyedgeRnonLog_TumorHypo.rds")

# Run reset function on C1 tumor samples - enhancers - FL only
SNFClusterLabsC1 <- which(SNFClusterLabs == 1)
identical(colnames(BetaMatrixT1SNF[, SNFClusterLabsC1]), 
          colnames(RNAseqProtCodT3MatNameSNF[, SNFClusterLabsC1])) # TRUE

# Extract genes
resetScoreENHSNFAll <- data.frame(resetResultsENHFLAll101$Score.report[which(! is.na(resetResultsENHFLAll101$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHSNFAll)) {
  resetScoreENHSNFAll[i, 6] <- resetResultsENHFLAll101$FDR.res$FDR[which(resetScoreENHSNFAll$Score[i] == resetResultsENHFLAll101$FDR.res$observed.scores)]
}
head(resetScoreENHSNFAll)
dim(resetScoreENHSNFAll) # 493 6
length(unique(resetScoreENHSNFAll$Gene)) # 417
round(range(resetScoreENHSNFAll$Score), 3) # 0.0133416 1.1358012
round(range(resetScoreENHSNFAll$FDR), 3) # 0.05933333 1.02532507
resetScoreENHSNFAll[order(resetScoreENHSNFAll$FDR, decreasing = FALSE),]
testingGene <-
resetScoreENHSNFAll[which(resetScoreENHSNFAll$Gene == testingGene), ]




# looking at cluster 1
set.seed(1234)
resetResultsENHFLC1SNF <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                                  meth.tumor = BetaMatrixT1SNF[, SNFClusterLabsC1],
                                  transcriptome = RNAseqProtCodT3MatNameSNF[, SNFClusterLabsC1],
                                  methylation.event = c('enh'),
                                  FDR.permutation.no = 1000, 
                                  seed = 100)
# saveRDS(resetResultsENHFLC1SNF, file = "40_resetResultsENHFLonlyedgeRnonLog_TumorHypoC1.rds")
resetResultsENHFLC1SNF <- readRDS("40_resetResultsENHFLonlyedgeRnonLog_TumorHypoC1.rds")


dim(resetResultsENHFLC1SNF$FDR.res) #  266   2
range(resetResultsENHFLC1SNF$FDR.res) # 0.01682077 1.17079453
head(resetResultsENHFLC1SNF$meth.tumor.status.all)
ENHFLmeth.tumor.status.allC1SNF <- data.frame(resetResultsENHFLC1SNF$meth.tumor.status.all, 
                                           Gene = gsub(".*@", "", rownames(resetResultsENHFLC1SNF$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.allC1SNF)
dim(ENHFLmeth.tumor.status.allC1SNF) # 2798  54
dim(resetResultsENHFLC1SNF$meth.tumor.all) # 2798  53
dim(resetResultsENHFLC1SNF$transcriptome) # 2798  53
resetResultsENHFLC1SNF$FDR.res
range(resetResultsENHFLC1SNF$Score.report$Score, na.rm = TRUE) # 0.01682077 1.17079453

# Extract genes
resetScoreENHSNFC1 <- data.frame(resetResultsENHFLC1SNF$Score.report[which(! is.na(resetResultsENHFLC1SNF$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHSNFC1)) {
  resetScoreENHSNFC1[i, 6] <- resetResultsENHFLC1SNF$FDR.res$FDR[which(resetScoreENHSNFC1$Score[i] == resetResultsENHFLC1SNF$FDR.res$observed.scores)]
}
head(resetScoreENHSNFC1)
dim(resetScoreENHSNFC1) # 266 6
length(unique(resetScoreENHSNFC1$Gene)) # 240
round(range(resetScoreENHSNFC1$Score), 3) #  0.017 1.171
round(range(resetScoreENHSNFC1$FDR), 3) # 0.102 0.997
resetScoreENHSNFC1[order(resetScoreENHSNFC1$FDR, decreasing = FALSE), ]
testingGene <-
  resetScoreENHSNFC1[which(resetScoreENHSNFC1$Gene == testingGene), ]





# looking at cluster 2
SNFClusterLabsC2 <- which(SNFClusterLabs == 2)
identical(colnames(BetaMatrixT1SNF[, SNFClusterLabsC2]), 
          colnames(RNAseqProtCodT3MatNameSNF[, SNFClusterLabsC2])) # TRUE

set.seed(1234)
resetResultsENHFLC2SNF <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                                meth.tumor = BetaMatrixT1SNF[, SNFClusterLabsC2],
                                transcriptome = RNAseqProtCodT3MatNameSNF[, SNFClusterLabsC2],
                                methylation.event = c('enh'),
                                FDR.permutation.no = 1000, 
                                seed = 100)
# saveRDS(resetResultsENHFLC2SNF, file = "40_resetResultsENHFLonlyedgeRnonLog_TumorHypoC2.rds")
resetResultsENHFLC2SNF <- readRDS(file = "40_resetResultsENHFLonlyedgeRnonLog_TumorHypoC2.rds")


dim(resetResultsENHFLC2SNF$FDR.res) #  355   2
range(resetResultsENHFLC2SNF$FDR.res) # 0.02886528 1.17079453
head(resetResultsENHFLC2SNF$meth.tumor.status.all)
ENHFLmeth.tumor.status.allC2SNF <- data.frame(resetResultsENHFLC2SNF$meth.tumor.status.all, 
                                              Gene = gsub(".*@", "", rownames(resetResultsENHFLC2SNF$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.allC2SNF)
dim(ENHFLmeth.tumor.status.allC2SNF) # 2798  54
dim(resetResultsENHFLC2SNF$meth.tumor.all) # 2798  53
dim(resetResultsENHFLC2SNF$transcriptome) # 2798  53
resetResultsENHFLC2SNF$FDR.res
range(resetResultsENHFLC2SNF$Score.report$Score, na.rm = TRUE) # 0.03082655 1.43586074

# Extract genes
resetScoreENHSNFC2 <- data.frame(resetResultsENHFLC2SNF$Score.report[which(! is.na(resetResultsENHFLC2SNF$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHSNFC2)) {
  resetScoreENHSNFC2[i, 6] <- resetResultsENHFLC2SNF$FDR.res$FDR[which(resetScoreENHSNFC2$Score[i] == resetResultsENHFLC2SNF$FDR.res$observed.scores)]
}
head(resetScoreENHSNFC2)
dim(resetScoreENHSNFC2) # 355 6
length(unique(resetScoreENHSNFC2$Gene)) # 314
round(range(resetScoreENHSNFC2$Score), 3) # 0.031 1.436
round(range(resetScoreENHSNFC2$FDR), 3) # 0.016 1.085
resetScoreENHSNFC2[order(resetScoreENHSNFC2$FDR, decreasing = FALSE),]
which(resetScoreENHSNFC2[order(resetScoreENHSNFC2$FDR, decreasing = FALSE), ]$FDR < 0.05)
testingGene <-
  resetScoreENHSNFC2[which(resetScoreENHSNFC2$Gene == testingGene), ]




# looking at n = 121
# match methylation samples with RNAseq samples
matchBetaRNAseq <- match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized[, c(11:131)]), 
                         colnames(BetaMatrix_T1))


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
BetaMatrix_T1MatMatched121 <- match(colnames(RNAseqProtCodT3Mat), colnames(BetaMatrix_T1[, c(11:165)]))
BetaMatrix_T1MatMatched121 <- BetaMatrix_T1MatMatched121[! is.na(BetaMatrix_T1MatMatched121)]
length(BetaMatrix_T1MatMatched121) # 121
dim(BetaMatrix_T1[, c(11:165)][, BetaMatrix_T1MatMatched121]) # 595564    121
BetaMatrix_T1121 <- BetaMatrix_T1[, c(11:165)][, BetaMatrix_T1MatMatched121]
dim(BetaMatrix_T1121) # 595564    121


RNAseqProtCodT3MatMatMatched121 <- match(colnames(BetaMatrix_T1121), 
                                         colnames(RNAseqProtCodT3Mat))
RNAseqProtCod121 <- RNAseqProtCodT3Mat[, RNAseqProtCodT3MatMatMatched121]
dim(RNAseqProtCod121) #  17190   121
rownames(RNAseqProtCod121) <- RNAseqProtCodT3Mat$name
range(RNAseqProtCod121) # 0.00 17766.33
dim(RNAseqProtCod121) # 17190   121


# check colnames
identical(colnames(RNAseqProtCod121), colnames(BetaMatrix_T1121)) # TRUE

set.seed(1234)
resetResultsENHFL121 <- reset(normal.db = methNorSelOutput$normal.enh.probes,
                                meth.tumor = as.matrix(BetaMatrix_T1121),
                                transcriptome = as.matrix(RNAseqProtCod121),
                                methylation.event = c('enh'),
                                FDR.permutation.no = 1000, 
                                seed = 100)
# saveRDS(resetResultsENHFL121, file = "40_resetResultsENHFLonlyedgeRnonLog_TumorHypo121.rds")
resetResultsENHFL121 <- readRDS(file = "40_resetResultsENHFLonlyedgeRnonLog_TumorHypo121.rds")

dim(resetResultsENHFL121$FDR.res) #  552   2
range(resetResultsENHFL121$FDR.res) # 0.02886528 1.17079453
head(resetResultsENHFL121$meth.tumor.status.all)
ENHFLmeth.tumor.status.all121 <- data.frame(resetResultsENHFL121$meth.tumor.status.all, 
                                              Gene = gsub(".*@", "", rownames(resetResultsENHFL121$meth.tumor.status.all)))
head(ENHFLmeth.tumor.status.all121)
dim(ENHFLmeth.tumor.status.all121) # 2798  54
dim(resetResultsENHFL121$meth.tumor.all) # 2798  53
dim(resetResultsENHFL121$transcriptome) # 2798  53
resetResultsENHFL121$FDR.res
range(resetResultsENHFL121$Score.report$Score, na.rm = TRUE) # 0.02683298 1.08602792

# Extract genes
resetScoreENHFL121 <- data.frame(resetResultsENHFL121$Score.report[which(! is.na(resetResultsENHFL121$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreENHFL121)) {
  resetScoreENHFL121[i, 6] <- resetResultsENHFL121$FDR.res$FDR[which(resetScoreENHFL121$Score[i] == resetResultsENHFL121$FDR.res$observed.scores)]
}
head(resetScoreENHFL121)
dim(resetScoreENHFL121) # 552   6
length(unique(resetScoreENHFL121$Gene)) # 465
round(range(resetScoreENHFL121$Score), 3) # 0.027 1.086
round(range(resetScoreENHFL121$FDR), 3) # 0.044 0.994
resetScoreENHFL121[order(resetScoreENHFL121$FDR, decreasing = FALSE),]
testingGene <-
  resetScoreENHFL121[which(resetScoreENHFL121$Gene == testingGene), ]
resetScoreENHFL121[order(resetScoreENHFL121$FDR, decreasing = FALSE),]






genelist <- c("MKRN3")
testingGene <- genelist[1]
ENHFLmeth.tumor.status.all121[which(ENHFLmeth.tumor.status.all121$Gene == testingGene), ]
rowSums(ENHFLmeth.tumor.status.all121[which(ENHFLmeth.tumor.status.all121$Gene == testingGene), c(1:57)])


# If one probe
plotting <- function() {
  if(length(which(gsub(".*@", "", rownames(resetResultsENHFL121$normal.meth)) == testingGene)) == 1) {
    
    
    testingDataFrame <- data.frame(Values = c(resetResultsENHFL121$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFL121$normal.meth)) == testingGene)[1], ],
                                              resetResultsENHFL121$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFL121$meth.tumor.all)) == testingGene)[1], ],
                                              resetResultsENHFL121$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFL121$transcriptome)) == testingGene)[1], ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 193 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ]),
                                MeanMethTumor = mean(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = mean(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    # only looking at samples with a status marked as 1
    hypoProbes <- which(ENHFLmeth.tumor.status.allC1[which(ENHFLmeth.tumor.status.allC1$Gene == testingGene), ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ],
                                                    resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), hypoProbes],
                                                    resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) #  133   2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ]),
                                      MeanMethTumor = mean(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ]),
                                         t(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 250   4
    
    
    testingDataFrame %>%  # highest score: cg20802515.p1.SPO11 methylation only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg20802515.p1.SPO11") %>% # highest score: cg20802515.p1.SPO11
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg20802515.p1.SPO11") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg20802515.p1.SPO11 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    # look at 2nd probe with highest score
    hypoProbes <- which(ENHFLmeth.tumor.status.all[which(ENHFLmeth.tumor.status.all$Gene == testingGene)[2], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene)[2], ],
                                                    resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene)[2], hypoProbes],
                                                    resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene)[2], hypoProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hypoProbes)), rep("TumorRNAseq", length(hypoProbes))))
    dim(testingDataFrameStatus) # 25  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsENHFLC1$normal.meth[which(gsub(".*@", "", rownames(resetResultsENHFLC1$normal.meth)) == testingGene)[2], ]),
                                      MeanMethTumor = mean(resetResultsENHFLC1$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsENHFLC1$meth.tumor.all)) == testingGene)[2], hypoProbes]),
                                      MeanTranscriptome = mean(resetResultsENHFLC1$transcriptome[which(gsub(".*@", "", rownames(resetResultsENHFLC1$transcriptome)) == testingGene)[2], hypoProbes]))
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
  }
}

# View(resetScoreENHFLonly)
# write.csv(resetScoreENHFLonly, file = "RESET_HypermethylatedProbes_FL_Enhancers.csv")

# Plot No.Methylation.Events vs Score
resetScoreENHFLC1 %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()

# Plot No.Methylation.Events vs Score
resetScoreENHFLC1 %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreENHFLC1 %>%  
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL - C1")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreENHFLC1 %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypomethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


# GProfilerAnalysis
# get gene list of those with FDR < 0.5
resetScoreENHFLC1Genes <- resetScoreENHFLC1 %>%  
  dplyr::filter(FDR < 0.5) %>% 
  dplyr::select(Gene)

gprofilerResetScoreENHFLC1 <- GProfilerAnalysis(GeneIDs = list(resetScoreENHFLC1Genes$Gene),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = NA, # with respect to expression
                                                ConditionName = "All",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
gprofilerResetScoreENHFLC1$shortLink
# "https://biit.cs.ut.ee/gplink/l/PC-KEBe5Rj"

# run with all genes 
gprofilerResetScoreENHFLC1All <- GProfilerAnalysis(GeneIDs = list(resetScoreENHFLC1$Gene),
                                                   Organism = "hsapiens",
                                                   OrderedQuery = TRUE,
                                                   PvalAlphaLevel = 0.01,
                                                   PositiveorNegFC = NA, # with respect to expression
                                                   ConditionName = "All",
                                                   ProduceImages = "Yes", 
                                                   PNGorPDF = "png")
gprofilerResetScoreENHFLC1All$shortLink
# "https://biit.cs.ut.ee/gplink/l/MABDo3KOQG"

#### RESET_HypermethylatedProbesFL_(T3 with C1 and C2 SNF samples as tumor)_HypomethylatedNormal(Silencers) - SNF clusters (edgeR normalized counts) ####


# Run reset function on 101 tumor samples - enhancers - FL only
set.seed(123)
resetResultsSILFLAll101 <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                                  meth.tumor = BetaMatrixT1SNF,
                                  transcriptome = RNAseqProtCodT3MatNameSNF,
                                  methylation.event = c('sil'),
                                  FDR.permutation.no = 1000, 
                                  seed = 100)
dim(resetResultsSILFLAll101$FDR.res) #  4803    2
range(resetResultsSILFLAll101$FDR.res) # 0.00444775 1.47456633
head(resetResultsSILFLAll101$meth.tumor.status.all)
SILFLmeth.tumor.status.all <- data.frame(resetResultsSILFLAll101$meth.tumor.status.all, 
                                           Gene = gsub(".*@", "", rownames(resetResultsSILFLAll101$meth.tumor.status.all)))
head(SILFLmeth.tumor.status.all)
dim(SILFLmeth.tumor.status.all) # 45229   102
dim(resetResultsSILFLAll101$meth.tumor.all) # 45229   101
dim(resetResultsSILFLAll101$transcriptome) # 45229   101
resetResultsSILFLAll101$FDR.res
range(resetResultsSILFLAll101$Score.report$Score, na.rm = TRUE) # 0.00444775 1.47456633
# saveRDS(resetResultsSILFLAll101, file = "40_resetResultsSILFLonlyedgeRnonLog_TumorHyper.rds")
resetResultsSILFLAll101 <- readRDS("40_resetResultsSILFLonlyedgeRnonLog_TumorHyper.rds")


# Extract genes
resetScoreSILSNFAll <- data.frame(resetResultsSILFLAll101$Score.report[which(! is.na(resetResultsSILFLAll101$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSILSNFAll)) {
  resetScoreSILSNFAll[i, 6] <- resetResultsSILFLAll101$FDR.res$FDR[which(resetScoreSILSNFAll$Score[i] == resetResultsSILFLAll101$FDR.res$observed.scores)]
}
head(resetScoreSILSNFAll)
dim(resetScoreSILSNFAll) # 4803 6
length(unique(resetScoreSILSNFAll$Gene)) # 2752
round(range(resetScoreSILSNFAll$Score), 3) # 0.004 1.475
round(range(resetScoreSILSNFAll$FDR), 3) #  0.101 0.851
resetScoreSILSNFAll[order(resetScoreSILSNFAll$FDR, decreasing = FALSE),]
testingGene <-
  resetScoreSILSNFAll[which(resetScoreSILSNFAll$Gene == testingGene), ]






# Run reset function on C1 tumor samples - enhancers - FL only
resetResultsSILFLC1 <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                             meth.tumor = BetaMatrixT1SNF[, SNFClusterLabsC1],
                             transcriptome = RNAseqProtCodT3MatNameSNF[, SNFClusterLabsC1],
                             methylation.event = c('sil'),
                             FDR.permutation.no = 100, 
                             seed = 100)
dim(resetResultsSILFLC1$FDR.res) # 3139    2
SILFLmeth.tumor.status.allC1 <- data.frame(resetResultsSILFLC1$meth.tumor.status.all, 
                                           Gene = gsub(".*@", "", rownames(resetResultsSILFLC1$meth.tumor.status.all)))
head(SILFLmeth.tumor.status.allC1)
dim(SILFLmeth.tumor.status.allC1) #  45229    58
dim(resetResultsSILFLC1$meth.tumor.all) # 45229   54
dim(resetResultsSILFLC1$transcriptome) # 45229   53
# saveRDS(resetResultsSILFLC1, file = "40_resetResultsSILFLonlyedgeRnonLog_TumorHyperC1.rds")
resetResultsSILFLC1 <- readRDS("40_resetResultsSILFLonlyedgeRnonLog_TumorHyperC1.rds")


# Extract genes
resetScoreSILSNFC1 <- data.frame(resetResultsSILFLC1$Score.report[which(! is.na(resetResultsSILFLC1$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSILSNFC1)) {
  resetScoreSILSNFC1[i, 6] <- resetResultsSILFLC1$FDR.res$FDR[which(resetScoreSILSNFC1$Score[i] == resetResultsSILFLC1$FDR.res$observed.scores)]
}
head(resetScoreSILSNFC1)
dim(resetScoreSILSNFC1) # 3139 6
length(unique(resetScoreSILSNFC1$Gene)) # 1866
round(range(resetScoreSILSNFC1$Score), 3) # 0.000 1.332
round(range(resetScoreSILSNFC1$FDR), 3) #  0.579 1.680
resetScoreSILSNFC1[order(resetScoreSILSNFC1$FDR, decreasing = FALSE),]






# Run reset function on C2 tumor samples - enhancers - FL only
set.seed(123)
resetResultsSILFLC2 <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                             meth.tumor = BetaMatrixT1SNF[, SNFClusterLabsC2],
                             transcriptome = RNAseqProtCodT3MatNameSNF[, SNFClusterLabsC2],
                             methylation.event = c('sil'),
                             FDR.permutation.no = 100, 
                             seed = 100)
dim(resetResultsSILFLC2$FDR.res) # 3320    2
SILFLmeth.tumor.status.allC2 <- data.frame(resetResultsSILFLC2$meth.tumor.status.all, 
                                           Gene = gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.status.all)))
head(SILFLmeth.tumor.status.allC2)
dim(SILFLmeth.tumor.status.allC2) #  45229    49
dim(resetResultsSILFLC2$meth.tumor.all) # 45229   48
dim(resetResultsSILFLC2$transcriptome) # 45229   48
# saveRDS(resetResultsSILFLC2, file = "40_resetResultsSILFLonlyedgeRnonLog_TumorHyperC2.rds")
resetResultsSILFLC2 <- readRDS(file = "40_resetResultsSILFLonlyedgeRnonLog_TumorHyperC2.rds")


# Extract genes
resetScoreSILSNFC2 <- data.frame(resetResultsSILFLC2$Score.report[which(! is.na(resetResultsSILFLC2$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSILSNFC2)) {
  resetScoreSILSNFC2[i, 6] <- resetResultsSILFLC2$FDR.res$FDR[which(resetScoreSILSNFC2$Score[i] == resetResultsSILFLC2$FDR.res$observed.scores)]
}
head(resetScoreSILSNFC2)
dim(resetScoreSILSNFC2) # 3320 6
length(unique(resetScoreSILSNFC2$Gene)) # 1800
round(range(resetScoreSILSNFC2$Score), 3) # 0.004 2.158
round(range(resetScoreSILSNFC2$FDR), 3) #  0.579 1.680
resetScoreSILSNFC2[order(resetScoreSILSNFC2$FDR, decreasing = FALSE), ]
# get how many signficant entries
length(which(resetScoreSILSNFC2[order(resetScoreSILSNFC2$FDR, decreasing = FALSE), ]$FDR <= 0.05)) # 61
resetScoreSILSNFC2TabOrdered <- resetScoreSILSNFC2[order(resetScoreSILSNFC2$FDR, decreasing = FALSE), ]
resetScoreSILSNFC2TabOrdered[which(resetScoreSILSNFC2TabOrdered$FDR <= 0.05), ]
# run gProfiler analysis
length(unique(resetScoreSILSNFC2TabOrdered[which(resetScoreSILSNFC2TabOrdered$FDR <= 0.05), ]$Gene)) # 35 genes


gprofilerReset61 <- GProfilerAnalysis(GeneIDs = list(unique(resetScoreSILSNFC2TabOrdered[which(resetScoreSILSNFC2TabOrdered$FDR <= 0.05), ]$Gene)),
                                                 Organism = "hsapiens",
                                                 OrderedQuery = TRUE,
                                                 PvalAlphaLevel = 0.01,
                                                 PositiveorNegFC = NA, # with respect to expression
                                                 ConditionName = "61RESETgenes",
                                                 ProduceImages = "Yes", 
                                                 PNGorPDF = "png")
gprofilerReset61$shortLink
# https://biit.cs.ut.ee/gplink/l/KMfZtbgSSV
dim(gprofilerReset61$TableAllValues) # 44  6
 

gprofilerReset61FALSE <- GProfilerAnalysis(GeneIDs = list(unique(resetScoreSILSNFC2TabOrdered[which(resetScoreSILSNFC2TabOrdered$FDR <= 0.05), ]$Gene)),
                                      Organism = "hsapiens",
                                      OrderedQuery = FALSE,
                                      PvalAlphaLevel = 0.01,
                                      PositiveorNegFC = NA, # with respect to expression
                                      ConditionName = "61RESETgenes",
                                      ProduceImages = "Yes", 
                                      PNGorPDF = "png")
gprofilerReset61FALSE$shortLink
# https://biit.cs.ut.ee/gplink/l/KMfZtbgSSV
dim(gprofilerReset61FALSE$TableAllValues) # 28  6




# looking at n = 121
# match methylation samples with RNAseq samples
matchBetaRNAseq <- match(colnames(RNAseqQC18June2020T3Samples132$RNASeqCountMatrixMatchedSampleFilteredFeatureFilteredNoSexENSEMBLIDsNormalized[, c(11:131)]), 
                         colnames(BetaMatrix_T1))


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
BetaMatrix_T1MatMatched121 <- match(colnames(RNAseqProtCodT3Mat), colnames(BetaMatrix_T1[, c(11:165)]))
BetaMatrix_T1MatMatched121 <- BetaMatrix_T1MatMatched121[! is.na(BetaMatrix_T1MatMatched121)]
length(BetaMatrix_T1MatMatched121) # 121
dim(BetaMatrix_T1[, c(11:165)][, BetaMatrix_T1MatMatched121]) # 595564    121
BetaMatrix_T1121 <- BetaMatrix_T1[, c(11:165)][, BetaMatrix_T1MatMatched121]
dim(BetaMatrix_T1121) # 595564    121


RNAseqProtCodT3MatMatMatched121 <- match(colnames(BetaMatrix_T1121), 
                                         colnames(RNAseqProtCodT3Mat))
RNAseqProtCod121 <- RNAseqProtCodT3Mat[, RNAseqProtCodT3MatMatMatched121]
dim(RNAseqProtCod121) #  17190   121
rownames(RNAseqProtCod121) <- RNAseqProtCodT3Mat$name
range(RNAseqProtCod121) # 0.00 17766.33
dim(RNAseqProtCod121) # 17190   121


# check colnames
identical(colnames(RNAseqProtCod121), colnames(BetaMatrix_T1121)) # TRUE

set.seed(1234)
resetResultsSILFL121 <- reset(normal.db = methNorSelOutput$normal.sil.probes,
                              meth.tumor = as.matrix(BetaMatrix_T1121),
                              transcriptome = as.matrix(RNAseqProtCod121),
                              methylation.event = c('sil'),
                              FDR.permutation.no = 100, 
                              seed = 100)
# saveRDS(resetResultsSILFL121, file = "40_resetResultsSILFLonlyedgeRnonLog_TumorHypo121.rds")
resetResultsSILFL121 <- readRDS("40_resetResultsSILFLonlyedgeRnonLog_TumorHypo121.rds")


dim(resetResultsSILFL121$FDR.res) #  5814   2
range(resetResultsSILFL121$FDR.res) # 0.0005516156 1.4074279200
head(resetResultsSILFL121$meth.tumor.status.all)
SILFLmeth.tumor.status.all121 <- data.frame(resetResultsSILFL121$meth.tumor.status.all, 
                                            Gene = gsub(".*@", "", rownames(resetResultsSILFL121$meth.tumor.status.all)))
head(SILFLmeth.tumor.status.all121)
dim(SILFLmeth.tumor.status.all121) # 45229   122
dim(resetResultsSILFL121$meth.tumor.all) # 45229   121
dim(resetResultsSILFL121$transcriptome) # 45229   122
resetResultsSILFL121$FDR.res
round(range(resetResultsSILFL121$Score.report$Score, na.rm = TRUE), 2) #  0.00 1.41

# Extract genes
resetScoreSILFL121 <- data.frame(resetResultsSILFL121$Score.report[which(! is.na(resetResultsSILFL121$Score.report$Score) == TRUE), ], FDR = NA)
for (i in 1:nrow(resetScoreSILFL121)) {
  resetScoreSILFL121[i, 6] <- resetResultsSILFL121$FDR.res$FDR[which(resetScoreSILFL121$Score[i] == resetResultsSILFL121$FDR.res$observed.scores)]
}
head(resetScoreSILFL121)
dim(resetScoreSILFL121) # 5814    6
length(unique(resetScoreSILFL121$Gene)) # 3424
round(range(resetScoreSILFL121$Score), 3) # 0.001 1.407
round(range(resetScoreSILFL121$FDR), 3) # 0.107 0.808
resetScoreSILFL121[order(resetScoreSILFL121$FDR, decreasing = FALSE),]
testingGene <-
  resetScoreSILFL121[which(resetScoreSILFL121$Gene == testingGene), ]





gprofilerResetScoreSILSNFC2 <- GProfilerAnalysis(GeneIDs = list(resetScoreSILFL121$Gene),
                                                Organism = "hsapiens",
                                                OrderedQuery = TRUE,
                                                PvalAlphaLevel = 0.01,
                                                PositiveorNegFC = NA, # with respect to expression
                                                ConditionName = "ResetScoreSILSNFC2",
                                                ProduceImages = "Yes", 
                                                PNGorPDF = "png")
gprofilerResetScoreSILSNFC2$shortLink
# "https://biit.cs.ut.ee/gplink/l/cxBTXCoGS-"
dim(gprofilerResetScoreSILSNFC2$)






testingGene <- resetScoreSILSNFC2TabOrdered$Gene[which(resetScoreSILSNFC2TabOrdered$FDR <= 0.05)][1]

SILFLmeth.tumor.status.allC2[which(SILFLmeth.tumor.status.allC2$Gene == testingGene), ]
rowSums(SILFLmeth.tumor.status.allC2[which(SILFLmeth.tumor.status.allC2$Gene == testingGene), c(1:48)])





plotting <- function() {
  
  # If one probe
  if(length(which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene)) == 1) {
    testingDataFrame <- data.frame(Values = c(resetResultsSILFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene)[1], ],
                                              resetResultsSILFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.all)) == testingGene)[1], ],
                                              resetResultsSILFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC2$transcriptome)) == testingGene)[1], ]),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 48), rep("TumorRNAseq", 48)))
    dim(testingDataFrame) # 101 2
    head(testingDataFrame)
    
    testingDataFrame %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrame %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = mean(resetResultsSILFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene)[1], ]),
                                MeanMethTumor = mean(resetResultsSILFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.all)) == testingGene)[1], ]),
                                MeanTranscriptome = mean(resetResultsSILFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC2$transcriptome)) == testingGene)[1], ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
    
    
    # only looking at samples with a status marked as 1
    hyperProbes <- which(SILFLmeth.tumor.status.allC2[which(SILFLmeth.tumor.status.allC2$Gene == testingGene)[1], ] == 1)
    testingDataFrameStatus <- data.frame(Values = c(resetResultsSILFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene)[1], ],
                                                    resetResultsSILFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.all)) == testingGene)[1], hyperProbes],
                                                    resetResultsSILFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC2$transcriptome)) == testingGene)[1], hyperProbes]),
                                         Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", length(hyperProbes)), rep("TumorRNAseq", length(hyperProbes))))
    dim(testingDataFrameStatus) # 61  2
    head(testingDataFrameStatus)
    
    
    testingDataFrameStatus %>%  # methylation only
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    testingDataFrameStatus %>%  # RNAseq only
      filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrameStatus %>%  # All
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "Values", fill = "Type", add = "jitter",
                        ylab = "Values", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrameStatus <- data.frame(MeanNormal = mean(resetResultsSILFLC2$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFLC2$normal.meth)) == testingGene)[1], ]),
                                      MeanMethTumor = mean(resetResultsSILFLC2$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFLC2$meth.tumor.all)) == testingGene)[1], hyperProbes]),
                                      MeanTranscriptome = mean(resetResultsSILFLC2$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFLC2$transcriptome)) == testingGene)[1], hyperProbes]))
    
    
    MeanDataFrameStatus %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
    
    
    
  } else if(length(which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene)) > 1) {
    
    # IF MULTIPLE probes
    testingDataFrame <- data.frame(rbind(t(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ]),
                                         t(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ]),
                                         t(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ])),
                                   Type = c(rep("RLNMethylation", 5), rep("TumorMethylation", 121), rep("TumorRNAseq", 121)))
    dim(testingDataFrame) # 247   4
    
    
    testingDataFrame %>%  # highest score: cg01941671@p1@BMP3 methylation only
      melt(id.vars = "Type")  %>%
      filter(variable == "cg01941671.p1.BMP3") %>% # highest score: cg01941671@p1@BMP3
      filter(Type != "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 RNAseq only
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      filter(variable  == "cg01941671.p1.BMP3") %>% # highest score: cg01941671.p1.BMP3
      #filter(Type   == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Beta Value: cg01941671.p1.BMP3") +
      stat_compare_means() + stat_n_text()
    
    
    testingDataFrame %>% # highest score: cg01941671.p1.BMP3 methylation + RNAseq
      melt(id.vars = "Type")  %>%
      data.frame() %>%
      #filter(Type == "TumorRNAseq") %>%
      ggpubr::ggboxplot(x = "Type", y = "value", fill = "Type", add = "jitter",
                        ylab = "Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Value: all 3 probes") +
      stat_compare_means() + stat_n_text()
    
    
    MeanDataFrame <- data.frame(MeanNormal = rowMeans(resetResultsSILFL104only$normal.meth[which(gsub(".*@", "", rownames(resetResultsSILFL104only$normal.meth)) == testingGene), ]),
                                MeanMethTumor = rowMeans(resetResultsSILFL104only$meth.tumor.all[which(gsub(".*@", "", rownames(resetResultsSILFL104only$meth.tumor.all)) == testingGene), ]),
                                MeanTranscriptome = rowMeans(resetResultsSILFL104only$transcriptome[which(gsub(".*@", "", rownames(resetResultsSILFL104only$transcriptome)) == testingGene), ]))
    
    MeanDataFrame %>%
      reshape2::melt() %>%
      ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable", add = "jitter",
                        ylab = "Mean Value", 
                        font.label = list(size = 20, color = "black")) +
      ggtitle("Type vs. Mean Value") +
      stat_compare_means() + stat_n_text()
  }
}

# View(resetScoreSILFLonly)
# write.csv(resetScoreSILFLonly, file = "RESET_HypomethylatedProbes_FL_Silencers.csv")

# Plot No.Methylation.Events vs Score
resetScoreSILSNFC2 %>%
  ggplot(aes(x = No.Methylation.Events, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor()


# Plot No.Methylation.Events vs Score
resetScoreSILSNFC2 %>%
  ggplot(aes(x = FDR, y = Score)) + geom_point() + 
  geom_smooth(method = lm) +
  stat_cor() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



# Plot top 100 Genes vs Score
resetScoreSILSNFC2 %>% 
  #filter(FDR < 0.5) %>% # FDR < 0.5
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))



resetScoreSILSNFC2 %>%  
  filter(FDR < 0.5) %>% # FDR < 0.1
  arrange(desc(Score)) %>% # arrange highest to lowest by score value
  dplyr::distinct(Gene, .keep_all = TRUE) %>% # keep only distinct entries (unique genes)
  top_n(- 100) %>%  # select top 100 genes based on score
  ggplot(aes(x = reorder(Gene, -Score), y = Score, fill = FDR)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene") +
  ggtitle(paste0("Plot of Hypermethylated Probes in FL")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15))


# GProfilerAnalysis
# get gene list of those with FDR < 0.5
resetScoreSILFLC1Genes <- resetScoreSILSNFC2 %>%  
  dplyr::filter(FDR < 0.5) %>% 
  dplyr::select(Gene)





#### Looking at overalp - SNF clusters (edgeR normalized counts) #####

# compare RESET_HypomethylatedProbes_FL_(T3 with C1 and C2 SNF samples as tumor)_HypermethylatedNormal(Enhancers)
mainlistRESETenh <- list(All = rownames(resetScoreENHSNFAll),
                         C1 = rownames(resetScoreENHSNFC1),
                         C2 = rownames(resetScoreENHSNFC2))

OutputSNFClustersENH <- VennDiagramAnalysis(MainList = mainlistRESETenh, 
                                         Labels = c("All","C1","C2"),
                                         FigureGenerate = "Yes", 
                                         ImageName = "MethylMixResultsENH", 
                                         PNGorPDF = "png")

mainlistRESETsil <- list(All = rownames(resetScoreSILSNFAll),
                         C1 = rownames(resetScoreSILSNFC1),
                         C2 = rownames(resetScoreSILSNFC2))

OutputSNFClustersSIL <- VennDiagramAnalysis(MainList = mainlistRESETsil, 
                                            Labels = c("All","C1","C2"),
                                            FigureGenerate = "Yes", 
                                            ImageName = "MethylMixResultsSIL", 
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


# RESET_HypermethylatedProbesFL_(T3 with C1 and C2 SNF samples as tumor)_HypomethylatedNormal(Silencers)


#### Compared results with MethylMix results ####

MethylMixResultsTSSSNFC2 <- readRDS(file = "41_MethylMixResultsTSSC2_SNF_Nov2020.rds")

resetScoreENHSNFC2TabOrdered <- resetScoreENHSNFC2[order(resetScoreENHSNFC2$FDR, decreasing = FALSE), ]
C2MethylMixRESEToverlapENH <- intersect(resetScoreENHSNFC2TabOrdered$Gene[which(resetScoreENHSNFC2TabOrdered$FDR <= 0.05)],
                                     MethylMixResultsTSSSNFC2$MethylationDrivers)
length(C2MethylMixRESEToverlapENH) # 1; MYLK3

C2MethylMixRESEToverlap <- intersect(resetScoreSILSNFC2TabOrdered$Gene[which(resetScoreSILSNFC2TabOrdered$FDR <= 0.05)],
          MethylMixResultsTSSSNFC2$MethylationDrivers)
length(C2MethylMixRESEToverlap) # 23
noquote(C2MethylMixRESEToverlap)

gprofilerC2MethylMixRESEToverlap <- GProfilerAnalysis(GeneIDs = list(C2MethylMixRESEToverlap),
                                                      Organism = "hsapiens",
                                                      OrderedQuery = TRUE,
                                                      PvalAlphaLevel = 0.01,
                                                      PositiveorNegFC = NA, # with respect to expression
                                                      ConditionName = "C2MethylMixRESEToverlap",
                                                      ProduceImages = "Yes", 
                                                      PNGorPDF = "png")
gprofilerC2MethylMixRESEToverlap$shortLink
# "https://biit.cs.ut.ee/gplink/l/4QXpB9SWRl"
dim(gprofilerC2MethylMixRESEToverlap$TableAllValues) # 26  6



gprofilerC2MethylMixRESEToverlapENH <- GProfilerAnalysis(GeneIDs = list(c(C2MethylMixRESEToverlap, C2MethylMixRESEToverlapENH)),
                                       Organism = "hsapiens",
                                       OrderedQuery = TRUE,
                                       PvalAlphaLevel = 0.01,
                                       PositiveorNegFC = NA, # with respect to expression
                                       ConditionName = "C2MethylMixRESEToverlapENH",
                                       ProduceImages = "Yes", 
                                       PNGorPDF = "png")
gprofilerC2MethylMixRESEToverlapENH$shortLink
# "https://biit.cs.ut.ee/gplink/l/FcC35UQTQY"
dim(gprofilerC2MethylMixRESEToverlapENH$TableAllValues) # 26  6


MethylMixResultsTSS121 <- readRDS(file = "41_MethylMixResultsTSS_121_Dec2020.rds")
resetScoreENHFL121Ordered <- resetScoreENHFL121[order(resetScoreENHFL121$FDR, decreasing = FALSE), ]
resetScoreENHFL121Ordered[which(resetScoreENHFL121Ordered$FDR <= 0.05), ]
MethylMix121RESEToverlap <- intersect(resetScoreENHFL121Ordered$Gene[which(resetScoreENHFL121Ordered$FDR <= 0.05)],
                                      MethylMixResultsTSS121$MethylationDrivers)
MethylMix121RESEToverlap # "ENPP3" "MKRN3"


resetScoreSILFL121Ordered <- resetScoreSILFL121[order(resetScoreSILFL121$FDR, decreasing = FALSE), ]
resetScoreENHFL121Ordered[which(resetScoreENHFL121Ordered$FDR <= 0.05), ]
MethylMix121RESEToverlap <- intersect(resetScoreENHFL121Ordered$Gene[which(resetScoreENHFL121Ordered$FDR <= 0.05)],
                                      MethylMixResultsTSS121$MethylationDrivers)




#### save image ####
# save.image("40_RESET_RNAseqVsMethFC_script.RData")
# [END]
