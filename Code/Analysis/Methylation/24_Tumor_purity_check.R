# Date: 29 September 2019
# Function: Tumor purity check with purity provided by TGL. Not a function.
# Author: Anjali Silva

# install.packages("InfiniumPurify")
library(InfiniumPurify)
library(mclust)
library(minfi) # /Users/anjalisilva/Library/R/3.6/library/
# install.packages("factoextra")
# library(factoextra)
library(devtools)
library(gridExtra)
library(ggpubr)
library(dplyr)
# devtools::install_github("brentp/celltypes450")
# install.packages("FactoMineR", force =TRUE)
# install.packages(c("Factoshiny","missMDA","FactoInvestigate"))
# library(FactoMineR)
# library(Factoshiny)
# library(missMDA)
# library(FactoInvestigate)
# library(gridExtra)
# library(reshape)

# Check tumor purity provided by TGL 
TumorPurity <- function() {
purity <- read.delim(file = "purity_Kridel_ifpf-DLBC.txt")
dim(purity) # 177   2
length(unique(purity$SampleID)) # 172

# Loading needed packages
# LoadCheckPkg(RegularPckgs=c("tidyr","ggplot2","plotly"))
library(tidyr)
library(ggplot2)
library(plotly)
library(ggpubr)

purity_names <- c(substr(purity$SampleID[1:78], 5, 13), 
                  substr(purity$SampleID[79:92], 1, 9),
                  substr(purity$SampleID[93:107], 1, 10),
                  substr(purity$SampleID[108:177], 1, 9))
purity[, 1] <- purity_names

# beta_names <- c(substr(colnames(Output_Clustering_0.25SD_noSexChsomes_AllSamples_14$BetaMatrixProbesUsedClustering)[1:90], 1, 9),
#                colnames(Output_Clustering_0.25SD_noSexChsomes_AllSamples_14$BetaMatrixProbesUsedClustering)[91:105],
#                substr(colnames(Output_Clustering_0.25SD_noSexChsomes_AllSamples_14$BetaMatrixProbesUsedClustering)[106:171], 1, 9))

length(which(substr(colnames(Output_Remove_2_BetaMatrix), 4, 5) == "FL"))
Output_Remove_2_BetaMatrix_editnames <- Output_Remove_2_BetaMatrix

colnames(Output_Remove_2_BetaMatrix_editnames)[which(substr(colnames(Output_Remove_2_BetaMatrix_editnames), 4, 5) == "FL")] <-
substr(colnames(Output_Remove_2_BetaMatrix_editnames)[which(substr(colnames(Output_Remove_2_BetaMatrix_editnames), 4, 5) == "FL")], 1, 9)

match_beta_bypurity <- match(colnames(Output_Remove_2_BetaMatrix_editnames), purity[, 1])
purity_edited <- purity[match_beta_bypurity[!is.na(match_beta_bypurity)], ]
dim(purity_edited) # 170   2
length(unique(purity_edited$SampleID)) #170


purity_names_ordered <- c(substr(purity$SampleID[match_beta_bypurity][1:30], 5, 13),
                          substr(purity$SampleID[match_beta_bypurity][31:77], 5, 13),
                          substr(purity$SampleID[match_beta_bypurity][78:90], 1, 9),
                          substr(purity$SampleID[match_beta_bypurity][91:105], 1, 10),
                          substr(purity$SampleID[match_beta_bypurity][106:171], 1, 9))

# corrected purity table for 171 patients
testdata2 <- data.frame(fullname = purity_names_ordered,
                        names = c(substr(factor(purity_names_ordered), 4, 5)),
                        values = purity[match_beta_bypurity, 2])


p2 <- ggplot(testdata2, aes(x = names, y = values, color = factor(names))) +
  geom_boxplot(size = 2, outlier.size = 3) +
  scale_y_continuous(name = "Purity") +
  geom_jitter(aes(color = names), size = 3, alpha = 1, width = 0.2) +
  labs(x = "Type") +
  theme_bw() + theme(text = element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face = "bold")) +
  scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Type") +
  theme(aspect.ratio = 1, legend.position = "right", panel.background = element_rect(colour = "black", size = 1.5),  axis.title =  element_text(face = "bold"))


# If more than two groups to compare, use One-Way ANOVA Test
# http://www.sthda.com/english/wiki/one-way-anova-test-in-r
# This is an extension of independent two-samples t-test for comparing means in a situation where there are more than two groups

# compute One-Way ANOVA 
test_results <- summary(aov(values ~ names, data = testdata2))
test_pvalue <- summary(aov(values ~ names, data = testdata2))[[1]][5][[1]][1]
cat("\n One-Way ANOVA performed using proportions calculated based on Beta values.\n")

# In one-way ANOVA test, a significant p-value indicates that some of the group means are different
# To determine if the mean difference between specific pairs of group are statistically significant, 
# perform multiple pairwise-comparison between the means of groups.
if (test_pvalue < 0.05) {
  test_results <- TukeyHSD(aov(values ~ names, data = testdata2))
}

ggplot(testdata2, aes(x = names, y = values, color=factor(names))) +
  #geom_point(size=3)+
  theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "Purity") +
  geom_jitter(aes(color = names), size = 3, alpha = 1, width = 0.2) +
  labs(x = "Type") +
  annotate("text", x = 1.5, y = max(testdata2$values) + 0.1, label = paste("P.Value = ", round(test_pvalue,5)), size = 7) +
  theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
  scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
  theme(aspect.ratio = 1, legend.position = "right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))

ggplot(testdata2, aes(x = names, y = values, fill = names)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width=0.1) +
  scale_y_continuous("Purity") +
  scale_x_discrete(name="Type") +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(aspect.ratio = 1)

# Setting the number of comparisons for ggpubr
my_comparisons <- list(c("FL","DL"),c("FL","RL"),c("DL","RL"))

p4 <- ggboxplot(testdata2, x = "names", y = "values",
                color = "names", palette = "jco") + 
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.5, paired = FALSE)     # Add global p-value

ggviolin(testdata2, x = "names", y = "values", fill = "names",
         add = "boxplot") + 
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.5, paired = FALSE)


# remove patients of FL with < 0.3 purity
testing <- testdata2[testdata2$values < 0.3, ] 
testing_RemoveBelow0.3 <- testdata2

for(i in 1:nrow(testing)) {
  if(testing$names[i] == "FL") {
    testing_RemoveBelow0.3[as.numeric(rownames(testing[i,])), 3] <- NA
  }
}

# only 156-9 = 147 patients have purity above 0.3, total of 162 patients with DLBCL and RLN
dim(Output_Clustering_0.25SD_noSexChsomes_AllSamples_14$BetaMatrixProbesUsedClustering[, ! is.na(testing_RemoveBelow0.3$values)])
ClinicalFile_updSamples[! is.na(testing_RemoveBelow0.3$values), ]

my_comparisons <- list(c("FL","DL"),c("FL","RL"),c("DL","RL"))

ggviolin(testing_RemoveBelow0.3, x = "names", y = "values", fill = "names",
         add = "boxplot") + 
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.5, paired = FALSE)

# 585 probes
Output_Clustering_0.25SD_585probes_PurityRemoved <- Clustering(BetaMatrix = as.matrix(Output_Clustering_0.25SD_noSexChsomes_AllSamples_14$BetaMatrixProbesUsedClustering[, ! is.na(testing_RemoveBelow0.3$values)]), 
                                                                     MvalueMatrix = NA, 
                                                                     ListofProbes = rownames(Output_Clustering_0.25SD_noSexChsomes_AllSamples_14$BetaMatrixProbesUsedClustering[, ! is.na(testing_RemoveBelow0.3$values)]), 
                                                                     ClinicalFile = ClinicalFile_updSamples[! is.na(testing_RemoveBelow0.3$values), ],
                                                                     FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.25SDNoSexChsomes_Allsamples_585probes_PurityRemoved")
# 607 probes
Output_Clustering_0.25SD_607probes_PurityRemoved <- Clustering(BetaMatrix = as.matrix(Output_Clustering_0.25SD_noSexChsomes_AllSamples_14_607probes$BetaMatrixProbesUsedClustering[, ! is.na(testing_RemoveBelow0.3$values)]), 
                                                               MvalueMatrix = NA, 
                                                               ListofProbes = rownames(Output_Clustering_0.25SD_noSexChsomes_AllSamples_14_607probes$BetaMatrixProbesUsedClustering[, ! is.na(testing_RemoveBelow0.3$values)]), 
                                                               ClinicalFile = ClinicalFile_updSamples[! is.na(testing_RemoveBelow0.3$values), ],
                                                               FigureGenerate = "Yes", PNGorPDF = "png", ImageName = "0.25SDNoSexChsomes_Allsamples_607probes_PurityRemoved")


############################### #
# Divide purity based on cluster - 3 RPMM cluster model
purityIDs <- c(substr(as.character(factor(purity[c(1:78),1])), 5, 16), as.character(factor(purity[c(79:177),1])))
purityIDs <- c(paste0(purityIDs[c(1:30)], "_T1"), purityIDs[31:177])
match_purity_byBetaOrder <- match(colnames(Output_Clustering_0.25SD_noSexChsomes_AllSamples_14$BetaMatrixProbesUsedClustering), purityIDs)

purity_3clusterlabs <- data.frame(names = Output_Clustering_0.25SD_noSexChsomes_AllSamples_14$RPMM$RPMMoutputLabels,
                                  values = purity$purity[match_purity_byBetaOrder])

my_comparisons <- list(c("1", "2"), c("1", "3"), c("2", "3"))
ggviolin(purity_3clusterlabs, x = "names", y = "values", fill = "names",
         add = "boxplot") + 
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.5, paired = FALSE)

############################### #
# Divide purity based on cluster - 2 RPMM cluster model

purityIDs <- c(substr(as.character(factor(purity[c(1:78),1])), 5, 16), as.character(factor(purity[c(79:177),1])))
purityIDs <- c(paste0(purityIDs[c(1:30)],"_T1"), purityIDs[31:177])
match_purity_byBetaOrder_FLonly<- match(colnames(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$BetaMatrixProbesUsedClustering), purityIDs)

purity_2clusterlabs <- data.frame(names = Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$RPMM$RPMMoutputLabels,
                                  values = purity$purity[match_purity_byBetaOrder_FLonly])

my_comparisons_FLonly<-list(c("1","2"))
ggviolin(purity_2clusterlabs, x = "names", y = "values", fill = "names",
         add = "boxplot")+ 
  stat_compare_means(comparisons = my_comparisons_FLonly)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.5, paired = FALSE)

############################### #
# Divide purity based on cluster - 2 SNF cluster model

purityIDs <- c(substr(as.character(factor(purity[c(1:78),1])), 5, 16), as.character(factor(purity[c(79:177), 1])))
purityIDs <- c(paste0(purityIDs[c(1:30)],"_T1"), purityIDs[31:177])
match_purity_byBetaOrder_FLonly <- match(colnames(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$BetaMatrixProbesUsedClustering), purityIDs)
match_clinical_byBetaOrder_FLonly <- match(colnames(Output_Clustering_0.25SD_noSexChsomes_FLonlysamples_14$BetaMatrixProbesUsedClustering), Output_Remove_2$ClinicalFile_updSamples$SAMPLE_ID)


purity_2clusterlabs <- data.frame(names = Output_Remove_2$ClinicalFile_updSamples$STAGE[match_clinical_byBetaOrder_FLonly],
                                  values = purity$purity[match_purity_byBetaOrder_FLonly])

my_comparisons_FLonly<-list(c("ADVANCED","LIMITED"))
ggviolin(purity_2clusterlabs, x = "names", y = "values", fill = "names",
         add = "boxplot")+ 
  stat_compare_means(comparisons = my_comparisons_FLonly)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.5, paired = FALSE)



################################ #
# 1 August 2019
# Working with 24 cases with purified FL cells and microenvironment
TumEnvGeneExp <- read.delim(file = "differential_gene_expression_tum_env_rcd1Aug2019.txt")
IndicatorVariable <- vector(mode="character",length = length(TumEnvGeneExp$GENE))
for (i in 1:length(TumEnvGeneExp$GENE)) {
  if (TumEnvGeneExp$logFC[i] > 0) {
    print(i)
    IndicatorVariable[i] <- "tumour"
  } else{
    IndicatorVariable[i] <- "microenvironment"
  }
}

TumEnvGeneExp2 <- data.frame(gene = TumEnvGeneExp$GENE, tvalue = TumEnvGeneExp$t, purity = IndicatorVariable)
}

# Test InfiniumPurify purity 
RunningExample <- function() {
  # This data set lists abbreviations for all TCGA cancer types.
  data(abbr)
  dim(abbr)
  head(abbr)
  
  # An example data set for InfiniumClust and InfiniumPurify
  data(beta.emp)
  dim(beta.emp) # 1417   61
  head(beta.emp)
  
  data(abbr)
  # Print abbreviations of cancer types with known iDMCs.
  CancerTypeAbbr()
  TypesExample <- CancerTypeAbbr()
  length(TypesExample[[1]]) #32
  # DLBC                   Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
  
  # Estimate the tumor purity for 450K methylation data
  # getPurity(tumor.data,normal.data = NULL,tumor.type = NULL)
  
  ## load example data
  data(beta.emp)
  colnames(beta.emp)
  normal.data <- beta.emp[,1:21]
  tumor.data <- beta.emp[,22:61]
  ## call purity for single tumor sample
  purity <- getPurity(tumor.data = tumor.data[,1],normal.data = NULL,tumor.type= "LUAD")
  # purity # 0.6578583
  ## call purity for less than 20 tumor samples
  purity <- getPurity(tumor.data = tumor.data[,1:10],normal.data = NULL,tumor.type= "LUAD")
  # purity
  ## call purity for more than 20 tumor samples with matched normal samples
  # ncol(normal.data) # 21
  # ncol(tumor.data[,1:40]) # 40
  purity <- getPurity(tumor.data = tumor.data[,1:40],normal.data = normal.data)
  
  ################################### #
  # Clustering
  # InfiniumClust Tumor sample clustering from Infinium 450k array data
  # InfiniumClust(tumor.data, purity, K, maxiter = 100, tol = 0.001)
  
  ## load example data
  data(beta.emp)
  normal.data <- beta.emp[,1:21]
  tumor.data <- beta.emp[,22:31]
  ## estimate tumor purity 
  ncol(normal.data)# 21
  ncol(tumor.data) # 10
  purity <- getPurity(tumor.data = tumor.data, tumor.type= "LUAD")
  ## cluster tumor samples accounting for tumor purity
  out <- InfiniumClust(tumor.data, purity, K=3, maxiter=5, tol=0.001)
  mclust::map(out$Z)
  
  
  ################################### #
  # InfiniumDMC Differentially Methylation Calling accounting for tumor purity
  # Infer differentially methylated CpG sites with the consideration of tumor purities.
  ## load example data
  data(beta.emp)
  normal.data <- beta.emp[,1:21]
  tumor.data <- beta.emp[,22:61]
  ## estimate tumor purity
  purity <- getPurity(tumor.data = tumor.data,normal.data = normal.data)
  ## DM calling with normal controls
  DMC = InfiniumDMC(tumor.data = tumor.data,normal.data = normal.data,purity = purity)
  ## DM calling without normal control
  DMC_ctlFree = InfiniumDMC(tumor.data = tumor.data,purity = purity)
  
  
  ################################### #  
  # InfiniumPurify Purify tumor methylomes caused by normal cell contamination.
  # The function deconvolutes purified tumor methylomes by a linear regression model
  # Returns a matrix of purified beta values for all CpG sites (row) and tumor samples (column).
  ## load example data
  data(beta.emp)
  normal.data <- beta.emp[,1:21]
  tumor.data <- beta.emp[,22:61]
  ## estimate tumor purity
  purity <- getPurity(tumor.data = tumor.data,normal.data = NULL,tumor.type= "LUAD")
  length(purity) # 40
  ## correct tumor methylome by tumor purity
  # ncol(normal.data)
  # ncol(tumor.data)
  tumor.purified = InfiniumPurify(tumor.data = tumor.data[1:100,],
                                  normal.data = normal.data[1:100,],
                                  purity = purity)
  head(tumor.purified)
  dim(tumor.purified) # 100  40
  minfi::densityPlot(tumor.purified)
}

# Apply InfiniumPurify to our data
ApplyingData <- function() {
  ## load example data
  
  # A list containing informative Differential methylation CpG sites 
  # (iDMC) and their average methylation levels in tumor and normal samples.
  data(iDMC) 
  head(iDMC)
  typeof(iDMC) # list
  length(iDMC$DLBC) #1000
  tumor.type <- "DLBC"
  
  normal.data <- BetaMatrix_T1[, 166:170]
  tumor.data <- BetaMatrix_T1[, 1:165]
  tumor.data_noDLBCL <- BetaMatrix_T1[, 11:165]
  # with all probes
  # purity <- getPurity(tumor.data = tumor.data, normal.data = normal.data, tumor.type = "DLBC")
  # purity_noDLBCL <- getPurity(tumor.data = tumor.data_noDLBCL, normal.data = normal.data, tumor.type = "DLBC")
  # purity_withNormal <- data.frame(purity = getPurity(tumor.data = BetaMatrix_T1, tumor.type = "DLBC"))
  
  
  
  # with selected probes, as not enough controls = RECOMMENDED
  probes <- iDMC[[tumor.type]]
  probes.true <- names(probes[probes == T])
  beta.sel <- BetaMatrix_T1[row.names(BetaMatrix_T1) %in% probes.true,]
  # data(beta.emp) # An example data set for InfiniumClust and InfiniumPurify.
  # purity_OICR_withNormal <- data.frame(purity = getPurity(tumor.data = BetaMatrix_T1, tumor.type = tumor.type))
  purity_OICR_withNormal <- getPurity(tumor.data = beta.sel, tumor.type = tumor.type)
  

  # saveRDS(purity_OICR_withNormal, file = "Purity_281probes_19Dec2019.rds")
  # write.csv(purity_OICR_withNormal, file ="Purity_281probes_19Dec2019.csv")
  # purity_OICR_withNormal <- readRDS(file = "Purity_281probes_10Jan2020.rds")
  
  purity_OICR <- data.frame(purity = getPurity(tumor.data = beta.sel [, 1:165], normal.data = beta.sel [, 166:170], tumor.type = tumor.type))
  purity_OICR_noDLBCL <- data.frame(purity = getPurity(tumor.data = beta.sel [, 11:165], normal.data = beta.sel [, 166:170], tumor.type = tumor.type))
  
  plot(purity, purity_OICR$purity, col = c(rep(1, 10), rep(3, 155)), ylab = "OICRmethod", xlab = "RegularMethod")
  plot(purity_noDLBCL, purity_OICR_noDLBCL$purity, col = 3, ylab = "OICRmethod_noDLBCL", xlab = "RegularMethod_noDLBCL")
  plot(cbind(c(1:165), purity_OICR))
  plot(cbind(c(1:170), purity_OICR_withNormal))
  plot(purity_withNormal$purity, purity_OICR_withNormal$purity)
  plot(cbind(c(1:170), purity_withNormal))
  
  # reading data from TGL
  TGL.purity <- data.table::fread(file = "purity_Kridel_ifpf-DLBC.txt",
                                  sep = "\t", header = TRUE)
  TGL.purity$SampleID[which(substr(TGL.purity$SampleID, 1, 4) == "KRI_")] <- 
    substr(TGL.purity$SampleID[which(substr(TGL.purity$SampleID, 1, 4) == "KRI_")], 5, 16)
  TGL.purity$SampleID[which(substr(TGL.purity$SampleID, 10, 12) == "")] <- 
    paste0(TGL.purity$SampleID[which(substr(TGL.purity$SampleID, 10, 12) == "")], "_T1")
  TGL.purity <- TGL.purity[- which(TGL.purity$SampleID == "LY_FL_159_T1"), ]
  TGL.purity$SampleID[which(TGL.purity$SampleID == "LY_FL_159_T1_rep")] <- "LY_FL_159_T1"
  MatchingSamples <- match(colnames(BetaMatrix_T1), TGL.purity$SampleID)
  TGL.purity <- TGL.purity[MatchingSamples, ]
  plot(TGL.purity$purity, purity_OICR_withNormal, 
       ylab = "Anjali_InfiniumPurify_Method", xlab = "TGLPurity",
       col= c(rep(1,10), rep(2, 155)))
  
  
  PurityWithType <- data.frame(Purity = purity_OICR_withNormal, Type = ClinicalFile_T1$TYPE, Stage =  ClinicalFile_T1$STAGE)
  # plotting high purity samples by cluster and type
  p <- ggplot2::ggplot(PurityWithType, aes(x = factor(Type), y = Purity, color=factor(Type))) +
    geom_boxplot(size = 2, outlier.size = 3) +
    ggtitle("Purity by type") +
    scale_y_continuous(name = "Purity") +
    geom_jitter(aes(color = factor(Type)), size = 2, alpha = 1, width = 0.2) +
    labs(x = "Type") +
    theme_bw() + theme(text = element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face = "bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Type") +
    theme(aspect.ratio = 1, legend.position = "right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  my_comparisons <- list(c("FL","RLN"), c("RLN","DLBCL"), c("FL","DLBCL"))
  p_type <- ggpubr::ggviolin(PurityWithType, x = "Type", y = "Purity", fill = "Type",
                             add = "boxplot", ylab=" Purity") +
    ggtitle("Purity by type") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  my_comparisons <- list(c("LIMITED","ADVANCED"))
  p_stage <- ggpubr::ggviolin(PurityWithType[which(PurityWithType$Type == "FL"),], x = "Stage", y = "Purity", fill = "Stage",
                              add = "boxplot", ylab=" Purity") +
    ggtitle("Purity by stage") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  # ############################## #  
  # Analyze the 201 probes
  
  table(AnnotationFile$chr[AnnotationFile$V1 %in% rownames(beta.sel)])
  table(AnnotationFile$Relation_to_Island[AnnotationFile$V1 %in% rownames(beta.sel)])
  table(AnnotationFile$UCSC_RefGene_Name[AnnotationFile$V1 %in% rownames(beta.sel)])
  table(AnnotationFile$SBE_maf[AnnotationFile$V1 %in% rownames(beta.sel)])
  PieCharts(BetaMatrix = BetaMatrix_T1[rownames(BetaMatrix_T1) %in%  rownames(beta.sel), ], 
            AnnotationFile = AnnotationFile, 
            ClinicalFile = ClinicalFile_T1,
            ProduceImages = "Yes", 
            PNGorPDF = "png",  
            ImageName = "InfiniumPurify_201Probes")

  
  # ############################## #
  # Clustering
  # running for all samples, DLBCL, FL, RLN
  InfiniumClust_AllSamples_100iterations <- InfiniumPurify::InfiniumClust(tumor.data = beta.sel, purity = purity_OICR_withNormal, K = 2, maxiter = 100, tol = 0.001)
  table(mclust::map(InfiniumClust_AllSamples_100iterations$Z))
  # 1  2 
  # 86 84 
  
  InfiniumClust_AllSamples_100iterations_Tol0.001_WithPurity <- data.frame(cluster = mclust::map(InfiniumClust_AllSamples_100iterations$Z), 
                                                                           purity = purity_OICR_withNormal,
                                                                           type = ClinicalFile_T1$TYPE,
                                                                           stage = ClinicalFile_T1$STAGE )
  p1 <- ggplot2::ggplot(InfiniumClust_AllSamples_100iterations_Tol0.001_WithPurity, aes(x = factor(cluster), y = purity, color=factor(cluster))) +
    geom_boxplot(size = 2, outlier.size = 3) +
    ggtitle("All samples (100 iterations, 0.001 tolerance)") +
    scale_y_continuous(name = "Purity (170 samples, 201 probes)") +
    geom_jitter(aes(color = factor(cluster)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  my_comparisons <- list(c("1", "2"))
  p2 <- ggpubr::ggviolin(InfiniumClust_AllSamples_100iterations_Tol0.001_WithPurity, x = "cluster", y = "purity", fill = "cluster",
                         add = "boxplot", ylab=" Purity") +
    ggtitle("Plot - all samples with 20 iterations, 0.001 tolerance") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  InfiniumClust_AllSamples_50iterations <- InfiniumPurify::InfiniumClust(tumor.data = beta.sel, purity = purity_OICR_withNormal, K = 2, maxiter = 50, tol = 0.001)
  table(mclust::map(InfiniumClust_AllSamples_50iterations$Z))
  #  1  2 
  # 84 86
  
  # 50 to 100 iterations give consistent results
  InfiniumClust_AllSamples_100iterations <- InfiniumPurify::InfiniumClust(tumor.data = beta.sel, purity = purity_OICR_withNormal, K = 2, maxiter = 50, tol = 0.001)
  table(mclust::map(InfiniumClust_AllSamples_100iterations$Z))
  # 1  2 
  # 84 86 
  
  set.seed(1234)
  InfiniumClust_AllSamples_100iterations_Tol0.0001 <- InfiniumPurify::InfiniumClust(tumor.data = beta.sel, purity = purity_OICR_withNormal, K = 2, maxiter = 100, tol = 0.0001)
  table(mclust::map(InfiniumClust_AllSamples_100iterations_Tol0.0001$Z))
  # 1  2 
  # 87 83 
  
  table(mclust::map(InfiniumClust_AllSamples_100iterations_Tol0.0001$Z), ClinicalFile_T1$STAGE)
  # ADVANCED LIMITED
  # 1        49      24
  # 2       31      51
  chisq.test(table(mclust::map(InfiniumClust_AllSamples_100iterations_Tol0.0001$Z), ClinicalFile_T1$STAGE))
  # X-squared = 12.144, df = 1, p-value = 0.0004924
  
  
  table(mclust::map(InfiniumClust_AllSamples_100iterations_Tol0.0001$Z), ClinicalFile_T1$TRANSLOC_14_18)
  #  0  1
  # 1 14 29
  # 2 17 10
  
  chisq.test(table(mclust::map(InfiniumClust_AllSamples_100iterations_Tol0.0001$Z), ClinicalFile_T1$TRANSLOC_14_18))
  # X-squared = 5.0431, df = 1, p-value = 0.02472
  
  table(mclust::map(InfiniumClust_AllSamples_100iterations_Tol0.0001$Z), ClinicalFile_T1$TYPE)
  #    DLBCL FL RLN
  # 1     9 73   5
  # 2      1 82   0
  chisq.test(table(mclust::map(InfiniumClust_AllSamples_100iterations_Tol0.0001$Z), ClinicalFile_T1$TYPE))
  # X-squared = 11.835, df = 2, p-value = 0.002692
  
  InfiniumClust_AllSamples_201probes_100iterations_100iterations_Tol0.0001_WithPurity <- data.frame(cluster = mclust::map(InfiniumClust_AllSamples_100iterations_Tol0.0001$Z), 
                                                                                                    purity = purity_OICR_withNormal,
                                                                                                    stage = ClinicalFile_T1$STAGE,
                                                                                                    translocation = ClinicalFile_T1$TRANSLOC_14_18,
                                                                                                    type = ClinicalFile_T1$TYPE)
  
  my_comparisons <- list(c("1", "2"))
  p2 <- ggpubr::ggviolin(InfiniumClust_AllSamples_201probes_100iterations_100iterations_Tol0.0001_WithPurity, x = "cluster", y = "purity", fill = "cluster",
                         add = "boxplot", ylab=" Purity") +
    ggtitle("All samples") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  p1_201probes <- ggplot2::ggplot(InfiniumClust_AllSamples_201probes_100iterations_100iterations_Tol0.0001_WithPurity, 
                                  aes(x = factor(cluster), y = purity, fill=stage)) +
    geom_boxplot(size = 1, outlier.size = 3) +
    ggtitle("All samples") +
    scale_y_continuous(name = "Purity") +
    # geom_jitter(aes(color = factor(stage)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  
  
  p3_201probes <- ggplot2::ggplot(InfiniumClust_AllSamples_201probes_100iterations_100iterations_Tol0.0001_WithPurity, 
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
  
  
  # Cross tabulation between 100 iterations, with 0.001 and 0.0001 tolerance. 
  table(mclust::map(InfiniumClust_AllSamples_100iterations_Tol0.0001$Z), mclust::map(InfiniumClust_AllSamples_100iterations$Z))
  # 1  2
  # 1 83  0
  # 2  1 86
  
  # Cross tabulation between 20 and 100 iterations
  table(mclust::map(InfiniumClust_AllSamples_100iterations_Tol0.0001$Z), mclust::map(InfiniumClust_AllSamples_20iterations$Z))
  # 1  2
  # 1 75  8
  # 2 12 75
  
  # Outcome analysis
  SurvivalFile_All <- SurvivalFile[match(substr(colnames(BetaMatrix_T1), 1, 9), SurvivalFile$LY_FL_ID), ]
  BetaMatrix = BetaMatrix_T1
  ClinicalFile = ClinicalFile_T1
  SurvivalFile = SurvivalFile_All
  ClusterLabels = mclust::map(InfiniumClust_AllSamples_100iterations_Tol0.0001$Z)
  # Then see 15_SurivivalAnalysis.R
  
  
  
  # ############################## #
  # running for samples, DLBCL + FL
  set.seed(1234)
  InfiniumClust_DLBCL_FL <- InfiniumPurify::InfiniumClust(tumor.data = beta.sel[, c(1:165)], purity = purity_OICR_withNormal[c(1:165)], K=2, maxiter = 100, tol = 0.0001)
  table(mclust::map(InfiniumClust_DLBCL_FL$Z))
  # 1  2 
  # 88 77
  
  table(mclust::map(InfiniumClust_DLBCL_FL$Z), ClinicalFile_T1$STAGE[c(1:165)])
  # ADVANCED LIMITED
  # 1        52      27
  # 2       28      48
  chisq.test(table(mclust::map(InfiniumClust_DLBCL_FL$Z), ClinicalFile_T1$STAGE[c(1:165)]))
  # X-squared = 11.892, df = 1, p-value = 0.0005637
  
  
  table(mclust::map(InfiniumClust_DLBCL_FL$Z), ClinicalFile_T1$TRANSLOC_14_18[c(1:165)])
  #  0  1
  # 1 16 30
  # 2 15  9
  
  chisq.test(table(mclust::map(InfiniumClust_DLBCL_FL$Z), ClinicalFile_T1$TRANSLOC_14_18[c(1:165)]))
  # X-squared = 3.8516, df = 1, p-value = 0.0497
  
  InfiniumClust_DLBCL_FL_WithPurity <- data.frame(cluster = mclust::map(InfiniumClust_DLBCL_FL$Z), 
                                                  purity = purity_OICR_withNormal[c(1:165)],
                                                  stage = ClinicalFile_T1$STAGE[c(1:165)],
                                                  translocation = ClinicalFile_T1$TRANSLOC_14_18[c(1:165)],
                                                  type = ClinicalFile_T1$TYPE[c(1:165)])
  
  
  p5 <- ggplot2::ggplot(InfiniumClust_DLBCL_FL_WithPurity, aes(x = factor(cluster), y = purity, color=factor(cluster))) +
    geom_boxplot(size = 2, outlier.size = 3) +
    ggtitle("Plot - DLBCL+FL samples with 100 iterations, 0.0001 tolerance") +
    scale_y_continuous(name = "Purity (170 samples and 201 probes)") +
    geom_jitter(aes(color = factor(cluster)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  my_comparisons <- list(c("1", "2"))
  p6 <- ggpubr::ggviolin(InfiniumClust_DLBCL_FL_WithPurity, x = "cluster", y = "purity", fill = "cluster",
                         add = "boxplot", ylab=" Purity") +
    ggtitle("DLBCL+FL") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  
  DLBCL_FL_201probes1 <- ggplot2::ggplot(InfiniumClust_DLBCL_FL_WithPurity, 
                                         aes(x = factor(cluster), y = purity, fill=stage)) +
    geom_boxplot(size = 1, outlier.size = 3) +
    ggtitle("DLBCL+FL samples") +
    scale_y_continuous(name = "Purity") +
    # geom_jitter(aes(color = factor(stage)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  
  
  DLBCL_FL_201probes2 <- ggplot2::ggplot(InfiniumClust_DLBCL_FL_WithPurity, 
                                         aes(x = factor(cluster), y = purity, fill=factor(translocation) )) +
    geom_boxplot(size = 1, outlier.size = 3) +
    ggtitle("DLBCL+FL samples") +
    scale_y_continuous(name = "Purity") +
    # geom_jitter(aes(color = factor(stage)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    labs(fill = "translocation") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  # Outcome analysis
  SurvivalFile_DLBCL_FLonly <- SurvivalFile[match(substr(colnames(BetaMatrix_T1[, c(1:165)]), 1, 9), SurvivalFile$LY_FL_ID), ]
  BetaMatrix = BetaMatrix_T1[, c(1:165)]
  ClinicalFile = ClinicalFile_T1[c(1:165), ]
  SurvivalFile = SurvivalFile_DLBCL_FLonly
  ClusterLabels = mclust::map(InfiniumClust_DLBCL_FL$Z)
  
  # Then see 15_SurivivalAnalysis.R
  
  
  # ############################## #
  # running for FL only samples
  set.seed(1234)
  InfiniumClust_FLonlySamples <- InfiniumPurify::InfiniumClust(tumor.data = beta.sel[, c(11:165)], purity = purity_OICR_withNormal[c(11:165)], K=2, maxiter = 100, tol = 0.0001)
  table(mclust::map(InfiniumClust_FLonlySamples$Z))
  # 1  2 
  # 87 68 
  
  # Checking posterior probability
  InfiniumClust_FLonlySamples$Z
  pairs(InfiniumClust_FLonlySamples$Z)
  ggplot2::ggplot(data.frame(InfiniumClust_FLonlySamples$Z), aes(x=X1, y=X2)) + geom_point()
  hist(InfiniumClust_FLonlySamples$Z)
  PosteriorProbWithCluster <- data.frame(ID = c(1:155), InfiniumClust_FLonlySamples$Z, Cluster = mclust::map(InfiniumClust_FLonlySamples$Z))
  
  # scatter plot
  ggplot(PosteriorProbWithCluster, aes(ID)) + 
         geom_point(aes(y = X1, colour = "one")) + 
         geom_point(aes(y = X2, colour = "two")) + 
         #ggtitle("Posterior probabilities of FL only 2Cluster model by InfiniumClust") +
         scale_y_continuous(name = "Posterior probability") +
         scale_x_continuous(name = "Sample" ) +
         theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
         scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") #+
         #theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  # line plot
  ggplot(PosteriorProbWithCluster, aes(ID)) + 
         geom_line(aes(y = X1, colour = "one")) + 
         geom_line(aes(y = X2, colour = "two")) + 
         ggtitle("Posterior probabilities") +
         scale_y_continuous(name = "Posterior probability") +
         labs(x = "Sample") +
         theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold"))# +
    # scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    # theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  # density plot - under construction 
  ggplot(PosteriorProbWithCluster, aes(x = ID)) +
         geom_density()
        
  # bar plot
  tableProp3 <- as.data.frame(cbind(colnames(beta.sel[, c(11:165)]),
                                    mclust::map(InfiniumClust_FLonlySamples$Z), 
                                    InfiniumClust_FLonlySamples$Z[, 1], 
                                    InfiniumClust_FLonlySamples$Z[, 2]))
  names(tableProp3) <- c("Sample", "Cluster", "Probability 1", "Probability 2")
  tableProp3$Cluster <- factor(tableProp3$Cluster, levels = c("1", "2"), labels = c("C1", "C2"))
  tableProp3Melt <- reshape::melt(tableProp3, id.vars = c("Sample", "Cluster"), 
                                  measure.vars = c("Probability 1", "Probability 2"))
  tableProp3Melt$value <- as.numeric(levels(tableProp3Melt$value))[tableProp3Melt$value]
  ggplot2::ggplot(tableProp3Melt, aes(fill = variable, y = value, x = Sample)) +
                  geom_bar(position = "fill", stat = "identity") +
                  theme(axis.text.x = element_text(angle = 90)) + 
                  coord_cartesian(ylim = c(0, 1)) +
                  # labs(y = "Posterior probability") +
                  # y axis tick mark lables
                  scale_y_continuous(name = "Posterior probability", limits = c(0: 1)) +
                  scale_fill_discrete(name = "Cluster")
                                
                                  
  # comparing labels with clinical data
  table(mclust::map(InfiniumClust_FLonlySamples$Z), ClinicalFile_T1$STAGE[c(11:165)])
  # ADVANCED LIMITED
  # 1       60      27
  # 2       20      48
  chisq.test(table(mclust::map(InfiniumClust_FLonlySamples$Z), ClinicalFile_T1$STAGE[c(11:165)]))
  # X-squared = 22.353, df = 1, p-value = 2.269e-06
  
  
  table(mclust::map(InfiniumClust_FLonlySamples$Z), ClinicalFile_T1$TRANSLOC_14_18[c(11:165)])
  #   0  1
  #1 17 31
  #2 14  8
  chisq.test(table(mclust::map(InfiniumClust_FLonlySamples$Z), ClinicalFile_T1$TRANSLOC_14_18[c(11:165)]))
  # X-squared = 3.7924, df = 1, p-value = 0.05148
  
  # based on institute
  table(mclust::map(InfiniumClust_FLonlySamples$Z), ClinicalFile_T1$INSTITUTION[c(11:165)])
  # ADVANCED LIMITED
  # 1       60      27
  # 2       20      48
  
  chisq.test(  table(mclust::map(InfiniumClust_FLonlySamples$Z), ClinicalFile_T1$INSTITUTION[c(11:165)]))
  
  InfiniumClust_FLonlySamples_WithPurity <- data.frame(cluster = mclust::map(InfiniumClust_FLonlySamples$Z), 
                                                       purity = purity_OICR_withNormal[c(11:165)], 
                                                       stage = ClinicalFile_T1$STAGE[c(11:165)],
                                                       translocation = ClinicalFile_T1$TRANSLOC_14_18[c(11:165)],
                                                       type = ClinicalFile_T1$TYPE[c(11:165)])
  
  p7 <- ggplot2::ggplot(InfiniumClust_FLonlySamples_WithPurity, aes(x = factor(cluster), y = purity, color=factor(cluster))) +
    geom_boxplot(size=2, outlier.size = 3) +
    ggtitle("FL samples ") +
    scale_y_continuous(name = "Purity") +
    geom_jitter(aes(color = factor(cluster)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  my_comparisons <- list(c("1","2"))
  p8 <- ggpubr::ggviolin(InfiniumClust_FLonlySamples_WithPurity, x = "cluster", y = "purity", fill = "cluster",
                         add = "boxplot", ylab=" Purity") +
    ggtitle("FL samples ") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  
  FLonly_201probes1 <- ggplot2::ggplot(InfiniumClust_FLonlySamples_WithPurity, 
                                       aes(x = factor(cluster), y = purity, fill=stage)) +
    geom_boxplot(size = 1, outlier.size = 3) +
    ggtitle("FL samples") +
    scale_y_continuous(name = "Purity") +
    # geom_jitter(aes(color = factor(stage)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  
  
  FLonly_FL_201probes2 <- ggplot2::ggplot(InfiniumClust_FLonlySamples_WithPurity, 
                                          aes(x = factor(cluster), y = purity, fill=factor(translocation) )) +
    geom_boxplot(size = 1, outlier.size = 3) +
    ggtitle("FL samples") +
    scale_y_continuous(name = "Purity") +
    # geom_jitter(aes(color = factor(stage)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    labs(fill = "translocation") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position = "right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  
  # 27 Nov 2019
  # plotting purity_OICR_withNormal
  plot(purity_OICR_withNormal, col = c(rep(1, 10), rep(2, 155), rep(3, 5)), type = "p")
  PurityWithSTAGEandTYPE <- data.frame(purity = purity_OICR_withNormal, stage = ClinicalFile_T1$STAGE, type = ClinicalFile_T1$TYPE)
  p8 <- ggplot2::ggplot(PurityWithSTAGEandTYPE, aes(x = factor(stage), y = purity, color = factor(type))) +
    geom_boxplot(size = 2, outlier.size = 3) +
    ggtitle("Purity generated based on 201 probes and 170 patients") +
    scale_y_continuous(name = "Purity") +
    geom_jitter(aes(color = factor(type)), size = 2, alpha = 1, width = 0.2) +
    labs(x = "Stage") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Type") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  
  # Outcome analysis
  SurvivalFile_FLonly <- SurvivalFile[match(substr(colnames(BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1), 4,5) == "FL")]), 1, 9), SurvivalFile$LY_FL_ID), ]
  BetaMatrix = BetaMatrix_T1[, c(11:165)]
  ClinicalFile = ClinicalFile_T1[c(11:165), ]
  SurvivalFile = SurvivalFile_FLonly
  ClusterLabels = mclust::map(InfiniumClust_FLonlySamples$Z)
  # Then see 15_SurivivalAnalysis.R
  
  
  # running for FL only samples - repeating 10x to see stability
  InfiniumClust_FLonlySamples_10xrepeat <- InfiniumClust_FLonlySamples_10xrepeat_labels <- list()
  for (i in 1:10) {
    set.seed(i)
    InfiniumClust_FLonlySamples_10xrepeat[[i]] <- InfiniumPurify::InfiniumClust(tumor.data = beta.sel[, c(11:165)], purity = purity_OICR_withNormal[c(11:165)], K=2, maxiter = 100, tol = 0.0001)
    InfiniumClust_FLonlySamples_10xrepeat_labels[[i]] <- table(mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[i]]$Z)) 
  }
  
  #saveRDS(InfiniumClust_FLonlySamples_10xrepeat, file = "InfiniumClust_FLonlySamples_10xrepeat.rds")
  #InfiniumClust_FLonlySamples_10xrepeat_testing <- readRDS(file = "InfiniumClust_FLonlySamples_10xrepeat.rds")
  
  ClusterMemByRun <- cbind(mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[1]]$Z),
                           mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[2]]$Z),
                           mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[3]]$Z),
                           mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[4]]$Z),
                           mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[5]]$Z),
                           mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[6]]$Z),
                           mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[7]]$Z),
                           mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[8]]$Z),
                           mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[9]]$Z),
                           mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[10]]$Z))
  colnames(ClusterMemByRun) <- paste0("Run", c(1:10))
  rownames(ClusterMemByRun) <- paste0("Obs", c(1:155))
  
  gplots::heatmap.2(ClusterMemByRun, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')
  ggplot(data.frame(InfiniumClust_FLonlySamples_10xrepeat[[1]]$Z), aes(x=X1, y=X2)) + geom_point()
  
  d <- density(InfiniumClust_FLonlySamples_10xrepeat[[1]]$Z) # returns the density data
  plot(d) # plots the results
  
  gplots::heatmap.2(InfiniumClust_FLonlySamples_10xrepeat[[1]]$Z, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')
  
  
  # plotting posterior probabilities based on figure adapted from 21_ProportionVisualization.r
  BarPlot <- list()
  for (i in 1:10) {
    tableProp2 <- as.data.frame(cbind(colnames(beta.sel[, c(11:165)]),
                                      mclust::map(InfiniumClust_FLonlySamples_10xrepeat[[i]]$Z), 
                                      InfiniumClust_FLonlySamples_10xrepeat[[i]]$Z[, 1], 
                                      InfiniumClust_FLonlySamples_10xrepeat[[i]]$Z[, 2]))
    names(tableProp2) <- c("Sample", "Cluster", "Probability 1", "Probability 2")
    tableProp2$Cluster <- factor(tableProp2$Cluster, levels = c("1", "2"), labels = c("C1", "C2"))
    tableProp2Melt <- reshape::melt(tableProp2, id.vars = c("Sample", "Cluster"), 
                                    measure.vars = c("Probability 1", "Probability 2"))
    tableProp2Melt$value <- as.numeric(levels(tableProp2Melt$value))[tableProp2Melt$value]
    # Stacked + percent
    BarPlot[[i]] <- ggplot2::ggplot(tableProp2Melt, aes(fill = variable, y = value, x = Sample)) +
                                    geom_bar(position = "fill", stat = "identity") +
                                    theme(axis.text.x = element_text(angle = 90)) + 
                                    coord_cartesian(ylim = c(0, 1)) +
                                    # labs(y = "Posterior probability") +
                                    # y axis tick mark lables
                                    scale_y_continuous(name = "Posterior probability", limits = c(0: 1)) +
                                    scale_fill_discrete(name = "Cluster")
  }
  # plotting all together
  gridExtra::grid.arrange(BarPlot[[1]], BarPlot[[2]], BarPlot[[3]], BarPlot[[4]], BarPlot[[5]],
                          BarPlot[[6]], BarPlot[[7]], BarPlot[[8]], BarPlot[[9]], BarPlot[[10]],
                          ncol = 5, nrow = 2)

  
  ggplot(tableProp2Melt, aes(x = Sample)) +
    geom_density()
  
  
  # runnign 3_DensityPlot.R
  MethylationDensityPlot(ClusterLabels = mclust::map(InfiniumClust_FLonlySamples$Z),
                         PlotWithinAnnotationCategories = "No", 
                         ClinicalCategoryToVisualize = "TYPE", 
                         BetaMatrix = BetaMatrix_T1[, c(11:165)], 
                         FigureGenerate = "Yes",  
                         ImageName = "InfiniumPurify_FlOnly_2C", 
                         PNGorPDF = "png")
  
  # Runnign proportion visualization 
  ProportionVisualization(CategoryToVisualize = "Relation_to_Island", 
                          ClusterLabels = mclust::map(InfiniumClust_FLonlySamples$Z), 
                          BetaMatrix = BetaMatrix_T1[, c(11:165)], 
                          AnnotationFile = AnnotationFile, 
                          ClinicalFile = ClinicalFile_T1[c(11:165), ], 
                          ClinicalCategory = "STAGE", 
                          PlotWithinCategories = "Yes", 
                          FigureGenerate = "Yes", 
                          PNGorPDF = "png", 
                          ImageName = )
  
  # Looking at density plots of 2Cluster FL only InfiniumClust model
  Density_RelationToIsland_TYPE_InfiniumClus <- MethylationDensityPlot(AnnotationCategoryToVisualize = "Relation_to_Island", 
                                                                       ClusterLabels = mclust::map(InfiniumClust_FLonlySamples$Z),
                                                                       PlotWithinAnnotationCategories = "Yes", 
                                                                       ClinicalCategoryToVisualize = "TYPE", 
                                                                       BetaMatrix = BetaMatrix_T1[, c(11:165)], 
                                                                       AnnotationFile = AnnotationFile, 
                                                                       ClinicalFile = ClinicalFile_T1, 
                                                                       SampleSheet = sheet, 
                                                                       FigureGenerate = "Yes", 
                                                                       ImageName = "InfiniumClus_FLonly_2Cs", 
                                                                       PNGorPDF = "png") 
  # DMRCate analysis - differential methylation 
  InfiniumClust_FLonlySamples_DMRCate  <- Diff_MethylatedRegions(Method = "DMRcate", BetaMatrix = BetaMatrix_T1[, c(11:165)], 
                                                                 MvalueMatrix = MvalueMatrix_T1[, c(11:165)], 
                                                                 ContrastColumnName = "CLUSTER", ClinicalFile = ClinicalFile_T1, 
                                                                 ClusterLabels = mclust::map(InfiniumClust_FLonlySamples$Z), 
                                                                 AnnotationFile = AnnotationFile, ProduceImages = "Yes", 
                                                                 DMR = 1, PNGorPDF = "png", ExpressionFile = NA)
    

  # Function extract gene regions list from Diff_MethylatedRegions
  ExtractingElements <- function(x, FILE) {
    
    example <- unlist(strsplit(as.character(FILE$overlapping.genes[x]), "[,]"))
    unique_elements <- unique(trimws(sub("\\-.*", "", example)))
    
    GettingTable <- function(i, n, numb_unique_elements){
      element_matrix <- matrix(ncol=2, nrow=numb_unique_elements)
      element_matrix[i,1] <- unique_elements[i]
      element_matrix[i,2] <- FILE$Stouffer[n]
      return(element_matrix)
    }
    
    element_matrix <- GettingTable(i=c(1:length(unique_elements)), n=x, numb_unique_elements=length(unique_elements))
    return(element_matrix)
  }
                                 
  # Looking at all regions          
  Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast <- sapply(c(1:500), function(i) ExtractingElements(i, FILE = InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata))
  Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast, quote=FALSE)
  write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast, file="GProfiler_Input_overlapping_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_InfiniumClust1vs2.csv")
  
  # How many regions are in negative FC vs. positive FC
  length(which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff < 0 )) # negative foldchange; 4989
  length(which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff > 0 )) # positive foldchange; 16535
  
  # Looking at regions based on FC
  # positive FC
   table_Diff_MethylatedRegions_CLUSTER_AllProbes_InfiniumClust_positive <- InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata[InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff > 0, ]
   #  dim(table_Diff_MethylatedRegions_CLUSTER_AllProbes_InfiniumClust_positive) 16535     8
   # Function extract gene regions list from Diff_MethylatedRegions
   ExtractingElements <- function(x, FILE) {
     
     example <- unlist(strsplit(as.character(FILE$overlapping.genes[x]), "[,]"))
     unique_elements <- unique(trimws(sub("\\-.*", "", example)))
     
     GettingTable <- function(i, n, numb_unique_elements){
       element_matrix <- matrix(ncol=2, nrow=numb_unique_elements)
       element_matrix[i,1] <- unique_elements[i]
       element_matrix[i,2] <- FILE$Stouffer[n]
       return(element_matrix)
     }
     
     element_matrix <- GettingTable(i=c(1:length(unique_elements)), n=x, numb_unique_elements=length(unique_elements))
     return(element_matrix)
   }
   # Looking at positive FC regions          
   Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_positive <- sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_InfiniumClust_positive))
   Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_positive <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_positive, quote=FALSE)
   write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_positive, file="GProfiler_Input_overlapping_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_InfiniumClust1vs2_positive.csv")
   
   
   table_Diff_MethylatedRegions_CLUSTER_AllProbes_InfiniumClust_negative <- InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata[InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff < 0, ]
   # dim(table_Diff_MethylatedRegions_CLUSTER_AllProbes_InfiniumClust_negative) 4989    8
   # Looking at negative FC regions          
   Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_negative <- sapply(c(1:500), function(i) ExtractingElements(i, FILE = table_Diff_MethylatedRegions_CLUSTER_AllProbes_InfiniumClust_negative))
   Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_negative <- do.call(rbind, Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_negative, quote=FALSE)
   write.csv(Unique_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_1stContrast_negative, file="GProfiler_Input_overlapping_genes_Table_Diff_MethylatedRegions_CLUSTER_AllProbes_InfiniumClust1vs2_negative.csv")
   
   
  # Negative FC - getting top hits
  testing <- sort(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff[which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff < 0 )], decreasing = TRUE)
  head(testing)
  # 2.063910e-05 -2.589850e-05 -4.903444e-05 -4.908837e-05
  InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$overlapping.genes[which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff == head(testing)[1])]
  # "MGAT1, SNORD95"
  InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$overlapping.genes[which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff == head(testing)[2])]
  # "RAD23A"
  InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$overlapping.genes[which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff == head(testing)[3])]
  # "MRPS22, SNORA81, SNORD66, SNORD2, SNORD5, Metazoa_SRP, SNORA18, U4"
  InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$overlapping.genes[which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff == head(testing)[4])]
  # "THEM6, CTD-2292P10.4, CTD-2292P10.2"
  
  # Positive FC
  testing2 <- sort(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff[which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff > 0 )], decreasing = TRUE)
  head(testing2)
  # 0.1694719 0.1621857 0.1572766 0.1532893
  InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$overlapping.genes[which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff == head(testing2)[1])]
  # "TLE1"
  InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$overlapping.genes[which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff == head(testing2)[2])]
  # "SNORA30, FAM201A, ANKRD18A"
  InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$overlapping.genes[which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff == head(testing2)[3])]
  # "SNORA7, GFRA2"
  InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$overlapping.genes[which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff == head(testing2)[4])]
  InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$overlapping.genes[which(InfiniumClust_FLonlySamples_DMRCate$GRangesObject[[1]]@elementMetadata$meandiff == head(testing2)[5])]
  # "PLXNA4"
  
  ############################# #
  # running for all samples, DLBCL, FL, RLN, after removing low purity FL samples 
  # which(purity_OICR_withNormal < 0.25)
  length(which(purity_OICR_withNormal < 0.25)) # 14
  # Among these are the 5 RLN samples and 9 FL samples; ONLY removing FL samples 
  InfiniumClust_HighPurity_161samples <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_T1[, - which(purity_OICR_withNormal < 0.25)[1:9]], purity = purity_OICR_withNormal[- which(purity_OICR_withNormal < 0.25)[1:9]], K = 2, maxiter = 100, tol = 0.0001)
  table(mclust::map(InfiniumClust_HighPurity_161samples$Z))
  # 1  2 
  # 81 80
  
  InfiniumClust_HighPurity_161samples_data <- data.frame(cluster = mclust::map(InfiniumClust_HighPurity_161samples$Z), 
                                                         purity = purity_OICR_withNormal[- which(purity_OICR_withNormal < 0.25)[1:9]],
                                                         type = ClinicalFile_T1$TYPE[- which(purity_OICR_withNormal < 0.25)[1:9]],
                                                         stage = ClinicalFile_T1$STAGE[- which(purity_OICR_withNormal < 0.25)[1:9]],
                                                         translocation = ClinicalFile_T1$TRANSLOC_14_18[- which(purity_OICR_withNormal < 0.25)[1:9]])
  
  
  
  # plotting high purity samples by cluster and type
  p9 <- ggplot2::ggplot(InfiniumClust_HighPurity_161samples_data, aes(x = factor(cluster), y = purity, color=factor(type))) +
    geom_boxplot(size = 2, outlier.size = 3) +
    ggtitle("Plot - high purity samples (161) with 100 iterations, 0.0001 tolerance") +
    scale_y_continuous(name = "Purity") +
    geom_jitter(aes(color = factor(type)), size = 2, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face = "bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio = 1, legend.position = "right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  
  table(InfiniumClust_HighPurity_161samples_data$cluster, InfiniumClust_HighPurity_161samples_data$type)                                                                                                                                                   
  # DLBCL FL RLN
  # 1    10 71   0
  # 2     0 75   5
  chisq.test(table(InfiniumClust_HighPurity_161samples_data$cluster, InfiniumClust_HighPurity_161samples_data$type))
  # X-squared = 15.104, df = 2, p-value = 0.0005251
  
  table(InfiniumClust_HighPurity_161samples_data$cluster, ClinicalFile_T1$STAGE[- which(purity_OICR_withNormal < 0.25)[1:9]])                                                                                                                                                   
  # ADVANCED LIMITED
  # 1       48      23
  # 2       29      46
  chisq.test(table(InfiniumClust_HighPurity_161samples_data$cluster, ClinicalFile_T1$STAGE[- which(purity_OICR_withNormal < 0.25)[1:9]]))
  # X-squared = 11.121, df = 1, p-value = 0.0008535
  
  table(InfiniumClust_HighPurity_161samples_data$cluster, ClinicalFile_T1$TRANSLOC_14_18[- which(purity_OICR_withNormal < 0.25)[1:9]])                                                                                                                                                   
  #    0  1
  # 1 14 30
  # 2 15  8
  chisq.test(table(InfiniumClust_HighPurity_161samples_data$cluster, ClinicalFile_T1$TRANSLOC_14_18[- which(purity_OICR_withNormal < 0.25)[1:9]]))
  # X-squared = 5.5704, df = 1, p-value = 0.01827
  
  # plotting high purity samples by cluster and stage
  p10 <- ggplot2::ggplot(InfiniumClust_HighPurity_161samples_data, aes(x = factor(cluster), y = purity, color=factor(stage))) +
    geom_boxplot(size = 2, outlier.size = 3) +
    ggtitle("Plot - high purity samples (161) with 100 iterations, 0.0001 tolerance") +
    scale_y_continuous(name = "Purity") +
    geom_jitter(aes(color = factor(stage)), size = 2, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  # plotting high purity samples by cluster and type
  p11 <- ggplot2::ggplot(InfiniumClust_HighPurity_161samples_data, aes(x = factor(cluster), y = purity, color = factor(translocation))) +
    geom_boxplot(size = 2, outlier.size = 3) +
    ggtitle("Plot - high purity samples (161) with 100 iterations, 0.0001 tolerance") +
    scale_y_continuous(name = "Purity") +
    geom_jitter(aes(color = factor(translocation)), size = 2, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  
  # Do clusters differ in terms of purity for FL?
  ClusterLabels_TwoClusters_FLonly <- InfiniumClust_HighPurity_161samples_data$cluster[which(InfiniumClust_HighPurity_161samples_data$type == "FL")]                              
  Purity_TwoClusters_FLonly <- InfiniumClust_HighPurity_161samples_data$purity[which(InfiniumClust_HighPurity_161samples_data$type == "FL")] 
  stats::t.test(x = Purity_TwoClusters_FLonly[which(ClusterLabels_TwoClusters_FLonly == 1)], 
                y = Purity_TwoClusters_FLonly[which(ClusterLabels_TwoClusters_FLonly == 2)]) 
  # t = 1.6196, df = 143.92, p-value = 0.1075
  
  
  # Do clusters differ in terms of outcome?
  
  
  SurvivalFile_161samples <- SurvivalFile[match(substr(colnames(BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1[, - which(purity_OICR_withNormal < 0.25)[1:9]]), 4,5) == "FL")]), 1, 9), SurvivalFile$LY_FL_ID), ]
  
  Cluster_Outcome_p11 <- SurvivalAnalysis(BetaMatrix = BetaMatrix_T1[, which(substr(colnames(BetaMatrix_T1[, - which(purity_OICR_withNormal < 0.25)[1:9]]), 4,5) == "FL")], 
                                          ClinicalFile = ClinicalFile_T1[which(substr(colnames(BetaMatrix_T1[, - which(purity_OICR_withNormal < 0.25)[1:9]]), 4,5) == "FL"), ],
                                          SurvivalFile = SurvivalFile_161samples, 
                                          ClusterLabels = mclust::map(InfiniumClust_HighPurity_161samples$Z), 
                                          FigureGenerate = "Yes", 
                                          PNGorPDF = "png") 
  
  write.csv(InfiniumClust_HighPurity_161samples$Z, "Probabilities_161Samples_6Dec2019.csv")
  
  
  # Clustering of all samples, DLBCL, FL, RLN, after adjusting FL samples to have average purity
  mean(purity_OICR_withNormal[c(11:165)]) # 0.5576975
  MeanFLPurity_purity_OICR_withNormal <- purity_OICR_withNormal
  MeanFLPurity_purity_OICR_withNormal[c(11:165)] <- 0.56
  InfiniumClust_MeanPurity_Allsamples <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_T1, purity = MeanFLPurity_purity_OICR_withNormal, K = 2, maxiter = 100, tol = 0.0001)
  table(mclust::map(InfiniumClust_MeanPurity_Allsamples$Z))
  # 1  2 
  # 86 84
  
  # save data to a data frame
  MeanFLPurity_purity_OICR_withNormal_data <- data.frame(cluster = mclust::map(InfiniumClust_MeanPurity_Allsamples$Z), 
                                                         purity = MeanFLPurity_purity_OICR_withNormal,
                                                         type = ClinicalFile_T1$TYPE,
                                                         stage = ClinicalFile_T1$STAGE,
                                                         translocation = ClinicalFile_T1$TRANSLOC_14_18)
  
  # Compare with earlier clustering:
  table(InfiniumClust_HighPurity_161samples_data$cluster, MeanFLPurity_purity_OICR_withNormal_data$cluster[- which(purity_OICR_withNormal < 0.25)[1:9]])
  # 1  2
  # 1  7 74
  # 2 70 10
  
  table(InfiniumClust_AllSamples_201probes_100iterations_100iterations_Tol0.0001_WithPurity$cluster, MeanFLPurity_purity_OICR_withNormal_data$cluster)
  # 1  2
  # 1 64 19
  # 2 22 65
  
  # plotting high purity samples by cluster and type
  p12 <- ggplot2::ggplot(MeanFLPurity_purity_OICR_withNormal_data, aes(x = factor(cluster), y = purity, color=factor(type))) +
    geom_boxplot(size = 2, outlier.size = 3) +
    ggtitle("Plot - mean FL purity samples and others with 100 iters, 0.0001 tol") +
    scale_y_continuous(name = "Purity") +
    geom_jitter(aes(color = factor(type)), size = 2, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  table(MeanFLPurity_purity_OICR_withNormal_data$cluster, MeanFLPurity_purity_OICR_withNormal_data$type) 
  # DLBCL FL RLN
  # 1     0 81   5
  # 2    10 74   0
  chisq.test(table(MeanFLPurity_purity_OICR_withNormal_data$cluster, MeanFLPurity_purity_OICR_withNormal_data$type))
  # X-squared = 15.295, df = 2, p-value = 0.0004773
  
  table(MeanFLPurity_purity_OICR_withNormal_data$cluster, ClinicalFile_T1$STAGE)                                                                                                                                                   
  #       ADVANCED LIMITED
  # 1       25      56
  # 2       55      19
  chisq.test(table(MeanFLPurity_purity_OICR_withNormal_data$cluster, ClinicalFile_T1$STAGE))
  # X-squared = 27.533, df = 1, p-value = 1.545e-07
  
  table(MeanFLPurity_purity_OICR_withNormal_data$cluster, ClinicalFile_T1$TRANSLOC_14_18)                                                                                                                                                   
  #    0  1
  # 1 17  8
  # 2 14 31
  chisq.test(table(MeanFLPurity_purity_OICR_withNormal_data$cluster, ClinicalFile_T1$TRANSLOC_14_18))
  # X-squared = 7.4317, df = 1, p-value = 0.006409
  
  # plotting high purity samples by cluster and stage
  p10 <- ggplot2::ggplot(MeanFLPurity_purity_OICR_withNormal_data, aes(x = factor(cluster), y = purity, color=factor(stage))) +
    geom_boxplot(size = 2, outlier.size = 3) +
    ggtitle("Plot - mean FL purity samples and others with 100 iters, 0.0001 tol") +
    scale_y_continuous(name = "Purity") +
    geom_jitter(aes(color = factor(stage)), size = 2, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  # plotting high purity samples by cluster and type
  p11 <- ggplot2::ggplot(MeanFLPurity_purity_OICR_withNormal_data, aes(x = factor(cluster), y = purity, color = factor(translocation))) +
    geom_boxplot(size = 2, outlier.size = 3) +
    ggtitle("Plot - mean FL purity samples and others with 100 iters, 0.0001 tol") +
    scale_y_continuous(name = "Purity") +
    geom_jitter(aes(color = factor(translocation)), size = 2, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  
  ############################# #
  # https://bioconductor.org/packages/release/data/experiment/vignettes/FlowSorted.Blood.EPIC/inst/doc/FlowSorted.Blood.EPIC.html
  # OLOptimizedCpGs the IDOL L-DMR library for EPIC arrays
  head (IDOLOptimizedCpGs)  
  library(FlowSorted.Blood.EPIC)
  countsEPIC_170samples <- FlowSorted.Blood.EPIC::estimateCellCounts2(rgSet = Output_Data_2$RGChannelSet,
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
  names(countsEPIC_170samples)
  # "counts"         "compTable"      "normalizedData"
  
  dim(countsEPIC_170samples$counts) # 170   6
  # saveRDS(countsEPIC_170samples$counts, "EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020.rds")
  head(countsEPIC_170samples$counts)
  dim(countsEPIC_170samples$compTable) #865859     11
  length(countsEPIC_170samples$normalizedData) # 865859
  countsEPIC_170samples <- readRDS(file = "EstimatedCellCOunts_FlowSortedBloodEPIC_13Jan2020.rds")
  dim(countsEPIC_170samples) # 170   6
  
  
  Purity_170Samples <- readRDS(file = "Purity_281probes_15Nov2019.rds")
  dim(Purity_170Samples) # 170   1
  
  # See if sample order is the same
  setequal(rownames(countsEPIC_170samples$counts), rownames(Purity_170Samples)) # TRUE
  
  
  ######################################## #  
  # Perform PCA: methylation vs. immune cell subsets
  ######################################## # 
  # https://rpkgs.datanovia.com/factoextra/index.html
  # Starting principal component analysis of cell types by estimateCellCounts2
  cell.count.pca <- prcomp(countsEPIC_170samples, scale = TRUE)
  factoextra::fviz_eig(cell.count.pca) 
  # first dimension account for 51% of variability
  # second dimension account for ~20% of variability
  # third dimension account for ~14% of variability
  cell.count.pca.ind <- factoextra::get_pca_ind(cell.count.pca)
  cell.count.pca.ind <- cell.count.pca.ind$coord
  cell.count.pca.ind <- data.frame(cell.count.pca.ind)
  cell.count.pca.ind$SAMPLE_ID <- row.names(cell.count.pca.ind)
  # Retain only the first two dimensions
  cell.count.pca.ind <- cell.count.pca.ind %>%
    select(SAMPLE_ID, cell.count.Dim.1 = Dim.1, cell.count.Dim.2 = Dim.2)
  
  
  BetaMatrix.pca <- prcomp(t(BetaMatrix_T1), scale = TRUE)
  factoextra::fviz_eig(BetaMatrix.pca) 
  # first dimension account for 15% of variability
  # second dimension account for ~12.5% of variability
  # third dimension account for ~4% of variability
  # fourth dimension account for ~3% of variability
  
  # Plots based on Beta values
  BetaMatrix.pca.ind <- factoextra::get_pca_ind(BetaMatrix.pca)
  BetaMatrix.pca.ind <- BetaMatrix.pca.ind$coord
  BetaMatrix.pca.ind <- data.frame(BetaMatrix.pca.ind)
  BetaMatrix.pca.ind$SAMPLE_ID <- row.names(BetaMatrix.pca.ind)
  BetaMatrix.pca.ind <- BetaMatrix.pca.ind %>%
    select(SAMPLE_ID, BetaMatrix.Dim.1 = Dim.1, BetaMatrix.Dim.2 = Dim.2)
  
  BetaMatrixMerged.pca <- cell.count.pca.ind %>%
    left_join(BetaMatrix.pca.ind) %>%
    left_join(ClinicalFile_T1[,c("SAMPLE_ID", "TYPE")])
  
  Beta_plots_to_save <- list()
  
  Beta_plots_to_save[[1]] <- BetaMatrixMerged.pca %>%
    ggplot(aes(cell.count.Dim.1, BetaMatrix.Dim.1, col = TYPE),
           main = "PCA plot") +
    labs(x = "1st PC of cell fractions", y = "1st PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    ggtitle("Unadjusted Beta Values")
  
  Beta_plots_to_save[[2]] <- BetaMatrixMerged.pca %>%
    ggplot(aes(cell.count.Dim.1, BetaMatrix.Dim.2, col = TYPE),
           main = "PCA plot") + 
    labs(x = "1st PC of cell fractions", y = "2nd PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    ggtitle("Unadjusted Beta Values")
  
  
  gridExtra::grid.arrange(Beta_plots_to_save[[1]], Beta_plots_to_save[[2]], ncol=2)
  
  ggsave(filename = paste0("img/", date, " cor_pca_M.matrix_cell.count.png"),
         marrangeGrob(grobs = plots_to_save, ncol = 2, nrow = 1),
         device = "png", width = 25, height = 10, units = "cm")
  
  ######################### #
  # Use the 201 probes used for purity analysis
  ######################### #
  beta.sel.pca <- prcomp(t(beta.sel), scale = TRUE)
  factoextra::fviz_eig(beta.sel.pca) 
  # first dimension account for 55% of variability
  # second dimension account for ~5% of variability
  # third dimension account for ~2.5% of variability
  
  # Plots based on Beta values
  beta.sel.pca.ind <- factoextra::get_pca_ind(beta.sel.pca)
  beta.sel.pca.ind <- beta.sel.pca.ind$coord
  beta.sel.pca.ind <- data.frame(beta.sel.pca.ind)
  beta.sel.pca.ind$SAMPLE_ID <- row.names(beta.sel.pca.ind)
  beta.sel.pca.ind <- beta.sel.pca.ind %>%
    select(SAMPLE_ID, beta.sel.Dim.1 = Dim.1, beta.sel.Dim.2 = Dim.2)
  
  beta.sel.Merged.pca <- cell.count.pca.ind %>%
    left_join(beta.sel.pca.ind) %>%
    left_join(ClinicalFile_T1[,c("SAMPLE_ID", "TYPE")])
  
  
  Beta_plots_201probes_to_save <- list()
  
  Beta_plots_201probes_to_save[[1]] <- beta.sel.Merged.pca %>%
    ggplot2::ggplot(aes(cell.count.Dim.1, beta.sel.Dim.1, col = TYPE),
                    main = "PCA plot") +
    labs(x = "1st PC of cell fractions", y = "1st PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    ggtitle("Unadjusted Beta Values, 201 probes")
  
  Beta_plots_201probes_to_save[[2]] <- beta.sel.Merged.pca %>%
    ggplot2::ggplot(aes(cell.count.Dim.1, beta.sel.Dim.2, col = TYPE),
                    main = "PCA plot") + 
    labs(x = "1st PC of cell fractions", y = "2nd PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    ggtitle("Unadjusted Beta Values, 201 probes")
  
  
  gridExtra::grid.arrange(Beta_plots_201probes_to_save[[1]], Beta_plots_201probes_to_save[[2]], ncol=2)
  
  ######################### #
  # Using M value matrix with all probes
  ######################### #
  M.matrix.pca <- prcomp(t(MvalueMatrix_T1), scale = TRUE)
  factoextra::fviz_eig(M.matrix.pca)
  # first dimension account for 18.5% of variability
  # second dimension account for ~12.5% of variability
  # third dimension account for ~4% of variability
  # fourth dimension account for ~3% of variability
  # Plots based on M values
  M.matrix.pca.ind <- factoextra::get_pca_ind(M.matrix.pca)
  M.matrix.pca.ind <- M.matrix.pca.ind$coord
  M.matrix.pca.ind <- data.frame(M.matrix.pca.ind)
  M.matrix.pca.ind$SAMPLE_ID <- row.names(M.matrix.pca.ind)
  M.matrix.pca.ind <- M.matrix.pca.ind %>%
    select(SAMPLE_ID, M.matrix.Dim.1 = Dim.1, M.matrix.Dim.2 = Dim.2)
  
  merged.pca <- cell.count.pca.ind %>%
    left_join(M.matrix.pca.ind) %>%
    left_join(ClinicalFile_T1[,c("SAMPLE_ID", "TYPE")])
  
  plots_to_save <- list()
  
  plots_to_save[[1]] <- merged.pca %>%
    ggplot(aes(cell.count.Dim.1, M.matrix.Dim.1, col = TYPE),
           main = "PCA plot") +
    labs(x = "1st PC of DNA methylation", y = "1st PC of cell fractions") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() + 
    ggtitle("Unadjusted M Values")
  
  plots_to_save[[2]] <- merged.pca %>%
    ggplot(aes(cell.count.Dim.1, M.matrix.Dim.2, col = TYPE),
           main = "PCA plot") + 
    labs(x = "1st PC of DNA methylation", y = "2nd PC of cell fractions") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    ggtitle("Unadjusted M Values")
  
  par(mfrow=c(2,2))
  
  ggsave(filename = paste0("img/", date, " cor_pca_M.matrix_cell.count.png"),
         marrangeGrob(grobs = plots_to_save, ncol = 2, nrow = 1),
         device = "png", width = 25, height = 10, units = "cm")
  
  
  
  
  ######################################## #  
  # Perform PCA: methylation vs. immune cell subsets - using adjusted matrix
  ######################################## #  
  # Using Beta Matrix
  # adjust B matrix using function from https://github.com/brentp/celltypes450:
  
  # Starting principal component analysis of cell types by celltypes450, adjust.beta
  adjustedBetaMatrix_T1 <- celltypes450::adjust.beta(BetaMatrix_T1) 
  # using 361 dmrs
  # adjusting beta (this will take a while)...
  # 0.0924% of values < epsilon:
  #   0.0295% of values > 1-epsilon:
  #   these will be changed to epsilon, 1-epsilon respectively
  
  # This doesn't work!
  # cell.types = celltypes450::adjust.beta(BetaMatrix_T1, est.only=TRUE)
  # saveRDS(adjustedBetaMatrix_T1, file = "adjustedBetaMatrix_T1_celltypes450_10Jan2020.rds")
  # adjustedBetaMatrix_T1 <- readRDS(file = "adjustedBetaMatrix_T1_celltypes450_10Jan2020.rds")
  dim(adjustedBetaMatrix_T1) # 595564    170
  
  #--
  # Plot sd vs mean methylation/probe, with and without adjustment
  #--
  meta_data = ClinicalFile_T1
  B.matrix.sd.mean <- t(BetaMatrix_T1[sample(nrow(BetaMatrix_T1), 10000), ]) %>%
    data.frame() %>%
    mutate(SAMPLE_ID = row.names(.)) %>%
    left_join(meta_data[,c("SAMPLE_ID", "TYPE")], by = "SAMPLE_ID") %>%
    tidyr::pivot_longer(cols = -c(SAMPLE_ID, TYPE), names_to = "cg") %>%
    group_by(cg, TYPE) %>%
    summarize(mean = mean(value), sd = sd(value)) %>%
    mutate(matrix = "non.adjusted")
  
  B.adjusted.matrix.sd.mean <- t(adjustedBetaMatrix_T1[sample(nrow(adjustedBetaMatrix_T1), 10000), ]) %>%
    data.frame() %>%
    mutate(SAMPLE_ID = row.names(.)) %>%
    left_join(meta_data[,c("SAMPLE_ID", "TYPE")], by = "SAMPLE_ID") %>%
    tidyr::pivot_longer(cols = -c(SAMPLE_ID, TYPE), names_to = "cg") %>%
    group_by(cg, TYPE) %>%
    summarize(mean = mean(value), sd = sd(value)) %>%
    mutate(matrix = "adjusted")
  
  B.matrix.sd.mean %>%
    rbind(B.adjusted.matrix.sd.mean) %>%
    mutate(matrix = factor(matrix, levels = c("non.adjusted", "adjusted"))) %>%
    ggplot(aes(x = mean, y = sd)) +
    geom_point(alpha = 0.1) +
    facet_grid(cols = vars(TYPE), rows = vars(matrix)) +
    guides(fill = FALSE) +
    theme_bw()
  
  

  adjustedBetaMatrix_T1.pca <- stats::prcomp(t(adjustedBetaMatrix_T1))
  factoextra::fviz_eig(adjustedBetaMatrix_T1.pca) 
  # first dimension account for  ~15 of variability
  # second dimension account for ~8% of variability
  # third dimension account for ~4.5% of variability
  # fourth dimension account for ~3% of variability
  
  # Plots based on Beta values
  adjustedBetaMatrix_T1.pca.ind <- factoextra::get_pca_ind(adjustedBetaMatrix_T1.pca)
  adjustedBetaMatrix_T1.pca.ind <- adjustedBetaMatrix_T1.pca.ind$coord
  adjustedBetaMatrix_T1.pca.ind <- data.frame(adjustedBetaMatrix_T1.pca.ind)
  adjustedBetaMatrix_T1.pca.ind$SAMPLE_ID <- row.names(adjustedBetaMatrix_T1.pca.ind)
  adjustedBetaMatrix_T1.pca.ind <- adjustedBetaMatrix_T1.pca.ind %>%
    select(SAMPLE_ID, BetaMatrix.Dim.1 = Dim.1, BetaMatrix.Dim.2 = Dim.2)
  
  
  adjustedBetaMatrix_T1_merged.adj.pca <- cell.count.pca.ind %>%
    left_join(adjustedBetaMatrix_T1.pca.ind) %>%
    left_join(ClinicalFile_T1[,c("SAMPLE_ID", "TYPE")])
  
  adjustedBetaMatrix_T1_plots_to_save <- list()
  
  adjustedBetaMatrix_T1_plots_to_save[[1]] <- adjustedBetaMatrix_T1_merged.adj.pca %>%
    ggplot(aes(cell.count.Dim.1, BetaMatrix.Dim.1, col = TYPE),
           main = "PCA plot") +
    labs(x = "1st PC of cell fractions 1st", y = "PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    ggtitle("celltypes450 Adjusted Beta Values")
  
  adjustedBetaMatrix_T1_plots_to_save[[2]] <- adjustedBetaMatrix_T1_merged.adj.pca %>%
    ggplot(aes(cell.count.Dim.1, BetaMatrix.Dim.2, col = TYPE),
           main = "PCA plot") +
    labs(x = "1st PC of cell fractions", y = "2nd PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    ggtitle("celltypes450 Adjusted Beta Values")
  
  # ggsave(filename = paste0("img/", date, " cor_pca_M.matrix_cell_adj.count.png"),
  #        marrangeGrob(grobs = plots_to_save, ncol = 2, nrow = 1),
  #        device = "png", width = 25, height = 10, units = "cm")
  
  
  gridExtra::grid.arrange(Beta_plots_to_save[[1]], 
                          adjustedBetaMatrix_T1_plots_to_save[[1]],
                          Beta_plots_to_save[[2]], 
                          adjustedBetaMatrix_T1_plots_to_save[[2]], 
                          ncol=2, nrow=2)
  
  
  ######################### #
  # Use the 201 probes used for purity analysis
  ######################### #
  # adjustedbeta.sel <- celltypes450::adjust.beta(beta.sel)
  # using 1 dmrs
  # Error in B[dmrs, ] : subscript out of bounds
  # In addition: Warning message:
  # In estimate.cell.proportions(B = B, top_n = top_n, mc.cores = mc.cores,  :
  #                               We suggest you use at least 2x10^5 probes
  adjustedbeta.sel.pca <- stats::prcomp(t(adjustedbeta.sel))
  factoextra::fviz_eig(adjustedbeta.sel.pca) 
  
  
  
  
  
  
  ######################################## #  
  # using all 177 samples
  ######################################## #  
  BetaMatrix_177samples <- Output_Data_2_177samples$BetaMatrix
  probes <- iDMC[[tumor.type]]
  probes.true <- names(probes[probes == T])
  beta.sel_177 <- BetaMatrix_177samples[row.names(BetaMatrix_177samples) %in% probes.true,]
  # data(beta.emp) # An example data set for InfiniumClust and InfiniumPurify.
  purity_OICR_withNormal_177 <- getPurity(tumor.data = beta.sel_177, tumor.type = tumor.type)
  length(purity_OICR_withNormal_177)
  plot(cbind(c(1:15), purity_OICR_withNormal_177[BadSamples]), ylab = "purity", xaxt="n", xlab="", main= "Purity vs 15 minfi 'bad' samples")
  axis(1, at = 1:15, labels = names(purity_OICR_withNormal_177[BadSamples]),las=2 ,cex.axis=0.7)
  # "LY_FL_020_T1" "LY_FL_046_T1" "LY_FL_065_T1" "LY_FL_159_T1" "LY_FL_164_T1"
  #  "LY_FL_165_T1" "LY_FL_183_T1" "LY_FL_218_T1" "LY_FL_465_T1" "LY_RLN_002"  
  #  "LY_FL_066_T1" "LY_FL_149_T1" "LY_FL_168_T1" "LY_FL_287_T1" "LY_FL_297_T1"
  
  Purityof7SexSamples <- match(ClinicalFile$SAMPLE_ID[SamplesDifferSex], names(purity_OICR_withNormal))
  plot(cbind(c(1:7), purity_OICR_withNormal[Purityof7SexSamples]), ylab = "purity", xaxt="n", xlab="", main = "Purity vs 7 mismatched sex samples")
  axis(1, at = 1:7, labels = names(purity_OICR_withNormal[Purityof7SexSamples]),las=2 ,cex.axis=0.7)
  # "LY_FL_165_T1" "LY_FL_080_T1" "LY_FL_149_T1" "LY_FL_038_T1" "LY_FL_163_T1"
  #  "LY_FL_168_T1" "LY_FL_159_T1"
  
  
  
  
  ######################### #
  # Clustering with X probes
  ######################### #
  # running for all samples, DLBCL, FL, RLN
  InfiniumClust_AllSamples_605411probes_100iterations <- InfiniumPurify::InfiniumClust(tumor.data = BetaMatrix_T1, purity = purity_OICR_withNormal, K = 2, maxiter = 100, tol = 0.0001)
  table(mclust::map(InfiniumClust_AllSamples_605411probes_100iterations$Z))
  # 1  2 
  # 85 85
  InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity <- data.frame(cluster = mclust::map(InfiniumClust_AllSamples_605411probes_100iterations$Z), 
                                                                                                       purity = purity_OICR_withNormal,
                                                                                                       stage = ClinicalFile_T1$STAGE,
                                                                                                       translocation = ClinicalFile_T1$TRANSLOC_14_18)
  
  table(InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity$cluster, InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity$stage)                                                                                                                                                   
  # ADVANCED LIMITED
  # 1       51      24
  # 2       29      51
  chisq.test(table(InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity$cluster, InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity$stage))
  # X-squared = 14.38, df = 1, p-value = 0.0001494
  
  p1_605411probes <- ggplot2::ggplot(InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity, 
                                     aes(x = factor(cluster), y = purity, fill=stage)) +
    geom_boxplot(size = 1, outlier.size = 3) +
    ggtitle("Plot - all samples with 100 iterations, 0.0001 tolerance") +
    scale_y_continuous(name = "Purity (170 samples and 605,411 probes)") +
    # geom_jitter(aes(color = factor(stage)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  my_comparisons <- list(c("1", "2"))
  p2_605411probes <- ggpubr::ggviolin(InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity,
                                      x = "cluster", y = "purity", fill = "cluster", add = "boxplot", ylab=" Purity") +
    ggtitle("Plot - all samples with 100 iterations, 0.0001 tolerance") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  
  
  p3_605411probes <- ggplot2::ggplot(InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity, 
                                     aes(x = factor(cluster), y = purity, fill=factor(translocation) )) +
    geom_boxplot(size = 1, outlier.size = 3) +
    ggtitle("Plot - all samples with 100 iterations, 0.0001 tolerance") +
    scale_y_continuous(name = "Purity (170 samples and 605,411 probes)") +
    # geom_jitter(aes(color = factor(stage)), size = 3, alpha = 1, width = 0.2) +
    labs(x = "Cluster") +
    theme_bw() + theme(text = element_text(size=20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Cluster") +
    theme(aspect.ratio=1, legend.position="right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  table(InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity$cluster, InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity$translocation)                                                                                                                                                   
  #    0  1
  # 1 15 32
  # 2 16  7
  chisq.test(table(InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity$cluster, InfiniumClust_AllSamples_605411probes_100iterations_100iterations_Tol0.0001_WithPurity$translocation))
  # X-squared = 7.4119, df = 1, p-value = 0.006479
  
  
  # InfiniumDMC = Differentially Methylation Calling accounting for tumor purity
  # DM calling with normal controls
  DMC = InfiniumPurify::InfiniumDMC(tumor.data = tumor.data, normal.data = normal.data, purity = purity)
  
  
  
  ######################################## #  
  # Perform PCA: methylation vs. immune cell subsets - Before differentially Methylation Calling, on 1:165 samples
  ######################################## # 
  
  NoInfiniumDMC_BetaMatrix.pca <- prcomp(t(BetaMatrix_T1[, 1:165]), scale = TRUE)
  factoextra::fviz_eig(NoInfiniumDMC_BetaMatrix.pca) 
  # first dimension account for 14% of variability
  # second dimension account for ~12.5% of variability
  # third dimension account for ~4% of variability
  # fourth dimension account for ~3% of variability
  
  # Plots based on Beta values
  NoInfiniumDMC_BetaMatrix.pca.ind <- factoextra::get_pca_ind(NoInfiniumDMC_BetaMatrix.pca)
  NoInfiniumDMC_BetaMatrix.pca.ind <- NoInfiniumDMC_BetaMatrix.pca.ind$coord
  NoInfiniumDMC_BetaMatrix.pca.ind <- data.frame(NoInfiniumDMC_BetaMatrix.pca.ind)
  NoInfiniumDMC_BetaMatrix.pca.ind$SAMPLE_ID <- row.names(NoInfiniumDMC_BetaMatrix.pca.ind)
  NoInfiniumDMC_BetaMatrix.pca.ind <- NoInfiniumDMC_BetaMatrix.pca.ind %>%
    select(SAMPLE_ID, NoInfiniumDMC_BetaMatrix.Dim.1 = Dim.1, NoInfiniumDMC_BetaMatrix.Dim.2 = Dim.2)
  
  NoInfiniumDMC_BetaMatrixMerged.pca <- cell.count.pca.ind[c(1:165),] %>%
    left_join(NoInfiniumDMC_BetaMatrix.pca.ind) %>%
    left_join(ClinicalFile_T1[c(1:165), c("SAMPLE_ID", "TYPE")])
  
  NoInfiniumDMC_plots_to_save <- list()
  
  NoInfiniumDMC_plots_to_save[[1]] <- NoInfiniumDMC_BetaMatrixMerged.pca %>%
    ggplot(aes(cell.count.Dim.1, NoInfiniumDMC_BetaMatrix.Dim.1, col = TYPE),
           main = "PCA plot") +
    labs(x = "1st PC of cell fractions", y = "1st PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    ggtitle("Unadjusted Beta Values")
  
  NoInfiniumDMC_plots_to_save[[2]] <- NoInfiniumDMC_BetaMatrixMerged.pca %>%
    ggplot(aes(cell.count.Dim.1, NoInfiniumDMC_BetaMatrix.Dim.2, col = TYPE),
           main = "PCA plot") + 
    labs(x = "1st PC of cell fractions", y = "2nd PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    ggtitle("Unadjusted Beta Values")
  
  
  ######################################## # 
  # InfiniumPurify = Purify tumor methylomes caused by normal cell contamination.
  ## correct tumor methylome by tumor purity
  ## estimate tumor purity using 5 RLN samples x 4
  purity_InfiniumPurify_5times20 <- getPurity(tumor.data = BetaMatrix_T1[, 1:165],normal.data = BetaMatrix_T1[, rep(c(166:170), 4)], tumor.type= "DLBC")
  
  purity_InfiniumPurify_1RLN <- getPurity(tumor.data = BetaMatrix_T1[, 1:165],normal.data = BetaMatrix_T1[, rep(169, 20)], tumor.type= "DLBC")
  
  
  tumor.purified_170samples_605411probes = InfiniumPurify::InfiniumPurify(tumor.data = BetaMatrix_T1[, 1:165],
                                                                          normal.data = BetaMatrix_T1[, rep(c(166:170),4)],
                                                                          purity = purity_InfiniumPurify_5times20)
  
  ## estimate tumor purity using 1 RLN sample with low purity x 20
  tumor.purified_170samples_605411probes = InfiniumPurify::InfiniumPurify(tumor.data = BetaMatrix_T1[, 1:165],
                                                                          normal.data = BetaMatrix_T1[, rep(c(169),20)],
                                                                          purity = purity_InfiniumPurify_1RLN)
  
  
  dim(tumor.purified_170samples_605411probes) # 605411    165
  par(mfrow = c(1, 2))
  minfi::densityPlot(dat = as.matrix(tumor.purified_170samples_605411probes), main=paste0("Density of tumor.purified_170samples_605411probes"))
  minfi::densityPlot(dat = as.matrix(BetaMatrix_T1[, c(1:165)]), main = paste0("Density of BetaMatrix"))
  
  # InfiniumDMC = Differentially Methylation Calling accounting for tumor purity
  InfiniumDMC_BetaMatrix.pca <- prcomp(t(tumor.purified_170samples_605411probes), scale = TRUE)
  factoextra::fviz_eig(InfiniumDMC_BetaMatrix.pca) 
  # first dimension account for 16% of variability
  # second dimension account for ~12.% of variability
  # third dimension account for ~4% of variability
  # fourth dimension account for ~3% of variability
  
  # Plots based on Beta values
  # Extract the results for individuals
  InfiniumDMC_BetaMatrix.pca.ind <- factoextra::get_pca_ind(InfiniumDMC_BetaMatrix.pca)
  # fviz_pca_ind(InfiniumDMC_BetaMatrix.pca, col.ind = "cos2", 
  #             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  #             repel = TRUE # Avoid text overlapping (slow if many points)
  #)
  
  
  InfiniumDMC_BetaMatrix.pca.ind <- InfiniumDMC_BetaMatrix.pca.ind$coord
  InfiniumDMC_BetaMatrix.pca.ind <- data.frame(InfiniumDMC_BetaMatrix.pca.ind)
  InfiniumDMC_BetaMatrix.pca.ind$SAMPLE_ID <- row.names(InfiniumDMC_BetaMatrix.pca.ind)
  InfiniumDMC_BetaMatrix.pca.ind <- InfiniumDMC_BetaMatrix.pca.ind %>%
    select(SAMPLE_ID, InfiniumDMC_BetaMatrix.Dim.1 = Dim.1, InfiniumDMC_BetaMatrix.Dim.2 = Dim.2)
  
  InfiniumDMC_BetaMatrixMerged.pca <- cell.count.pca.ind[c(1:165),] %>%
    left_join(InfiniumDMC_BetaMatrix.pca.ind) %>%
    left_join(ClinicalFile_T1[c(1:165), c("SAMPLE_ID", "TYPE")])
  
  InfiniumDMC_Beta_plots_to_save <- list()
  
  InfiniumDMC_Beta_plots_to_save[[1]] <- InfiniumDMC_BetaMatrixMerged.pca %>%
    ggplot(aes(cell.count.Dim.1, InfiniumDMC_BetaMatrix.Dim.1, col = TYPE),
           main = "PCA plot") +
    labs(x = "1st PC of cell fractions", y = "1st PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    ggtitle("Adjusted Beta Values")
  
  InfiniumDMC_Beta_plots_to_save[[2]] <- InfiniumDMC_BetaMatrixMerged.pca %>%
    ggplot(aes(cell.count.Dim.1, InfiniumDMC_BetaMatrix.Dim.2, col = TYPE),
           main = "PCA plot") + 
    labs(x = "1st PC of cell fractions", y = "2nd PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    ggtitle("Adjusted Beta Values")
  
  
  gridExtra::grid.arrange(InfiniumDMC_Beta_plots_to_save[[1]], InfiniumDMC_Beta_plots_to_save[[2]], ncol=2)
  
  gridExtra::grid.arrange(NoInfiniumDMC_plots_to_save[[1]], 
                          InfiniumDMC_Beta_plots_to_save[[1]],
                          NoInfiniumDMC_plots_to_save[[2]], 
                          InfiniumDMC_Beta_plots_to_save[[2]], 
                          ncol=2, nrow=2)
  
  
  ######################################## # 
  # InfiniumPurify = Purify tumor methylomes caused by normal cell contamination.
  ## correct tumor methylome by tumor purity
  # Try with the RLN sample with least purity, i.e, LY_RLN_003 = 0.1201521
  tumor.purified_170samples_595564probes_normalselected = InfiniumPurify::InfiniumPurify(tumor.data = BetaMatrix_T1[, 1:165],
                                                                                         normal.data = BetaMatrix_T1[, rep(c(166:170), 4)],
                                                                                         purity = purity_OICR_withNormal)
  
  # Starting principal component analysis of InfiniumClust results
  adjustedBetaMatrix_T1.pca_InfiniumPurify <- stats::prcomp(t(tumor.purified_170samples_595564probes_normalselected))
  factoextra::fviz_eig(adjustedBetaMatrix_T1.pca_InfiniumPurify) 
  # first dimension account for  ~15 of variability
  # second dimension account for ~8% of variability
  # third dimension account for ~4.5% of variability
  # fourth dimension account for ~3% of variability
  
  # Plots based on Beta values
  adjustedBetaMatrix_T1.pca.ind_InfiniumPurify <- factoextra::get_pca_ind(adjustedBetaMatrix_T1.pca_InfiniumPurify)
  adjustedBetaMatrix_T1.pca.ind_InfiniumPurify <- adjustedBetaMatrix_T1.pca.ind_InfiniumPurify$coord
  adjustedBetaMatrix_T1.pca.ind_InfiniumPurify <- data.frame(adjustedBetaMatrix_T1.pca.ind_InfiniumPurify)
  adjustedBetaMatrix_T1.pca.ind_InfiniumPurify$SAMPLE_ID <- row.names(adjustedBetaMatrix_T1.pca.ind_InfiniumPurify)
  adjustedBetaMatrix_T1.pca.ind_InfiniumPurify <- adjustedBetaMatrix_T1.pca.ind_InfiniumPurify %>%
    select(SAMPLE_ID, BetaMatrix.Dim.1 = Dim.1, BetaMatrix.Dim.2 = Dim.2)
  
  
  adjustedBetaMatrix_T1_merged.adj.pca_InfiniumPurify <- cell.count.pca.ind[c(1:165), ] %>%
    left_join(adjustedBetaMatrix_T1.pca.ind_InfiniumPurify) %>%
    left_join(ClinicalFile_T1[c(1:165), c("SAMPLE_ID", "TYPE")])
  
  adjustedBetaMatrix_T1_plots_to_save_InfiniumPurify <- list()
  
  # library(scales)
  # show_col(hue_pal()(4))
  
  adjustedBetaMatrix_T1_plots_to_save_InfiniumPurify[[1]] <- adjustedBetaMatrix_T1_merged.adj.pca_InfiniumPurify %>%
    ggplot(aes(cell.count.Dim.1, BetaMatrix.Dim.1, col = TYPE),
           main = "InfiniumPurify PCA plot") +
    labs(x = "1st PC of cell fractions 1st", y = "PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    scale_color_manual(values = scales::hue_pal()(10)[c(1,5)]) +
    ggtitle("InfiniumPurify Adjusted Beta Values")
  
  adjustedBetaMatrix_T1_plots_to_save_InfiniumPurify[[2]] <- adjustedBetaMatrix_T1_merged.adj.pca_InfiniumPurify %>%
    ggplot(aes(cell.count.Dim.1, BetaMatrix.Dim.2, col = TYPE),
           main = "InfiniumPurify PCA plot") +
    labs(x = "1st PC of cell fractions", y = "2nd PC of DNA methylation") +
    geom_point() +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    stat_cor() +
    scale_color_manual(values = scales::hue_pal()(10)[c(1,5)]) +
    ggtitle("InfiniumPurify Adjusted Beta Values")
  
  # ggsave(filename = paste0("img/", date, " cor_pca_M.matrix_cell_adj.count.png"),
  #        marrangeGrob(grobs = plots_to_save, ncol = 2, nrow = 1),
  #        device = "png", width = 25, height = 10, units = "cm")
  
  # comparing deconvolution of InfiniumPurify (InfiniumClust) and celltype450
  # Plotting principal component analysis (PCA)
  gridExtra::grid.arrange(Beta_plots_to_save[[1]], 
                          adjustedBetaMatrix_T1_plots_to_save[[1]],
                          adjustedBetaMatrix_T1_plots_to_save_InfiniumPurify[[1]],
                          Beta_plots_to_save[[2]], 
                          adjustedBetaMatrix_T1_plots_to_save[[2]], 
                          adjustedBetaMatrix_T1_plots_to_save_InfiniumPurify[[2]],
                          ncol=3, nrow=2)
  
  
  
  
  ######################################## # 
  # 28 Nov 2019
  # RFpurify 
  # From paper RF_Purify: a novel tool for comprehensive analysis of tumor-purity in 
  # methylation array data based on random forest regression, 2019
  # R-package to predict ABSOLUTE and ESTIMATE tumor purity using a Random Forest regression model 
  
  # install.packages("devtools") 
  library(devtools)
  devtools::install_github("mwsill/RFpurify") # not working
  # prevent warnings from beeing converted to errors when calling install_github
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
  install_github('mwsill/RFpurify')
  library(RFpurify)
  library(minfiData)
  
  # load example data
  data(MsetEx)
  MsetEx
  # predict purity
  # absolute <- predict_purity(MsetEx,method="ABSOLUTE")
  # estimate <- predict_purity(MsetEx,method="ESTIMATE")
  
  
  # try by manually downloading 
  #  setwd("/Volumes/GoogleDrive/My Drive/UHN/FLOMICS/FLOMICS/Methylation/Pipeline/RFpurify-master/R")
  # source("predict_purity.R")
  # source("predict_purity_betas.R")
  
  # predict purity
  absolute <- RFpurify::predict_purity(Mset = Output_Data_2$MethylSetRaw, method = "ABSOLUTE")
  estimate <- RFpurify::predict_purity(Mset = Output_Data_2$MethylSetRaw, method = "ESTIMATE")
  
  # plot predicted ABSOLUTE purity against ESTIMATE purity
  plot(x = estimate, y = absolute, pch = 19) # ylim = c(0, 1), xlim = c(0, 1)
  
  fit <- lm(absolute~estimate)
  plot(estimate, absolute,
       xlab = "Predicted ESTIMATE purity", ylab = "Predicted ABSOLUTE purity",
       col = c(rep("red",10), rep("green", 155), rep("blue",5)),
       pch = 16, 
       main = "Predicted ESTIMATE purity against ABSOLUTE purity")
  abline(fit)
  
  fit2 <- lm(estimate~purity_OICR_withNormal)
  
  plot(estimate, purity_OICR_withNormal,
       xlab = "Predicted ESTIMATE purity", ylab = "Non-deconvoluted Purity",
       col = c(rep("red",10), rep("green", 155), rep("blue",5)),
       pch = 16, 
       main = "InfiniumPurify non-deconvoluted vs RFPurify ESTIMATE purity")
  abline(fit2)
  
  fit3 <- lm(absolute~purity_OICR_withNormal)
  
  plot(absolute, purity_OICR_withNormal,
       xlab = "Predicted ABSOLUTE purity", ylab = "Non-deconvoluted Purity",
       col = c(rep("red",10), rep("green", 155), rep("blue",5)),
       pch = 16, 
       main = "InfiniumPurify non-deconvoluted vs RFPurify ABSOLUTE purity")
  abline(fit3)
  
  # for RF purify data 
  PurityRFPurify <- data.frame(estimate = estimate, absolute = absolute, Type = ClinicalFile_T1$TYPE, Stage =  ClinicalFile_T1$STAGE)
  my_comparisons <- list(c("FL","RLN"), c("RLN","DLBCL"), c("FL","DLBCL"))
  p_type_estimate <- ggpubr::ggviolin(PurityRFPurify, x = "Type", y = "estimate", fill = "Type",
                                      add = "boxplot", ylab=" ESTIMATE") +
    ggtitle("ESTIMATE purity  by type") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  p_type_absolute <- ggpubr::ggviolin(PurityRFPurify, x = "Type", y = "absolute", fill = "Type",
                                      add = "boxplot", ylab=" ABSOLUTE") +
    ggtitle("ABSOLUTE purity  by type") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  my_comparisons <- list(c("LIMITED","ADVANCED"))
  p_stage_estimate  <- ggpubr::ggviolin(PurityRFPurify[which(PurityRFPurify$Type == "FL"),], x = "Stage", y = "estimate", fill = "Stage",
                                        add = "boxplot", ylab=" ESTIMATE") +
    ggtitle("ESTIMATE purity by stage") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  p_stage_absolute  <- ggpubr::ggviolin(PurityRFPurify[which(PurityRFPurify$Type == "FL"),], x = "Stage", y = "absolute", fill = "Stage",
                                        add = "boxplot", ylab=" ABSOLUTE") +
    ggtitle("ABSOLUTE purity by stage") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  ######################################## # 
  # 29 Nov 2019
  # trying Purity Assessment from clonal MEthylation Sites (PAMES)
  # From paper Tumor purity quantification by clonal DNA methylation signatures, 2018
  # install.packages("devtools")
  devtools::install_github("cgplab/PAMES", build_vignettes = T)  
  library(PAMES)
  browseVignettes("PAMES")
  
  lsf.str("package:PAMES")
  length(lsf.str("package:PAMES")) # 5
  
  devtools::install_github("cgplab/PAMESdata")
  library(PAMESdata)
  lsf.str("package:PAMESdata") # nothing is printed
  
  AUC_FL <- compute_AUC(tumor_table = BetaMatrix_T1[c(1:100000), c(1:165)]*100, 
                        control_table = BetaMatrix_T1[c(1:100000), c(166:170)]*100, 
                        ncores = 4, na_threshold = 0)
  info_sites_FL <- select_informative_sites(tumor_table = BetaMatrix_T1[c(1:100000), c(1:165)]*100,
                                            auc = AUC_FL, 
                                            platform = "450k")
  # [2019-11-29 15:01:08] # Select informative sites #
  # Error: Number of rows of tumor_table is not equal to the number of rows of 
  # platform_data. Be sure to use correct platform and genome version and to remove any non-'cg' probe from tumor_table.
  
  purity <- compute_purity(tumor_table = BetaMatrix_T1[, c(1:165)], 
                           list_of_sites = list(hyper=c(1, 10, 20), hypo=c(15,30,45)))
  
  
  
  
  select_informative_sites(tumor_table = BetaMatrix_T1, )
  
  ######################################## # 
  # MethylPurify
  # From paper MethylPurify: tumor purity deconvolution and differential methylation detection 
  # from single tumor DNA methylomes, 2014
  # https://www.ncbi.nlm.nih.gov/pubmed/25103624/ - tool in python
  
  
  # LUMP
  # From paper Systematic pan-cancer analysis of tumour purity, 2015
  # https://www.ncbi.nlm.nih.gov/pubmed/26634437/
  
  
  
  ######################################## # 
  # InfiniumPurify - Try InfiniumDMC: differential methylation analysis 
  # accounting for tumor purity
  
  ##  ## DM calling without normal control
  DMC = InfiniumDMC(tumor.data = BetaMatrix_T1[,1:165], purity = purity_OICR_withNormal[1:165], threshold = 0.1)
  # probability is returned
  # the function computes posterior probability to rank CpG sites.
  DMC_Sites <- which(DMC$prob < 0.05)
  length(DMC_Sites) # 101469
  rownames(DMC[DMC_Sites,])
  rownames(DMC)[DMC_Sites[1]] # "cg14725810"
  rownames(DMC)[DMC_Sites[101469]] # "cg11219407"
  
  
  ######
  # InfiniumPurify: deconvolute pure tumor methylomes
  tumor.purified = InfiniumPurify::InfiniumPurify(tumor.data = BetaMatrix_T1[, 1:165],
                                                  normal.data = BetaMatrix_T1[, rep(c(166:170),4)],
                                                  purity = purity_OICR_withNormal[1:165])
  dim(tumor.purified) # 595564    165
  
  
  # Use deconvoluted values to get purity again 
  purity_withNormal_Purified <- data.frame(purity = getPurity(tumor.data = tumor.purified, tumor.type = "DLBC"))
  fit <- lm(purity_withNormal_Purified$purity~purity_OICR_withNormal[1:165])
  
  plot(purity_OICR_withNormal[1:165], purity_withNormal_Purified$purity,
       xlab = "Non-deconvoluted Purity", ylab = "Deconvoluted Purity",
       col = c(rep("red",10), rep("green", 155)),
       ylim = c(0.4, 0.9),
       xlim = c(0.1, 0.9),
       pch = 16, 
       main = "InfiniumPurify Deconvoluted vs Non-deconvoluted Purity ")
  abline(fit)
  
  PurityWithType_Purified <- data.frame(Purity = c(purity_withNormal_Purified$purity, purity_OICR_withNormal[166:170]), Type = ClinicalFile_T1$TYPE, Stage =  ClinicalFile_T1$STAGE)
  
  
  
  my_comparisons <- list(c("FL","RLN"), c("RLN","DLBCL"), c("FL","DLBCL"))
  p_type <- ggpubr::ggviolin(PurityWithType_Purified, x = "Type", y = "Purity", fill = "Type",
                             add = "boxplot", ylab=" Purity") +
    ggtitle("Purity by type") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  my_comparisons <- list(c("LIMITED","ADVANCED"))
  p_stage <- ggpubr::ggviolin(PurityWithType_Purified[which(PurityWithType_Purified$Type == "FL"),], x = "Stage", y = "Purity", fill = "Stage",
                              add = "boxplot", ylab=" Purity") +
    ggtitle("Purity by stage") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.5, paired = FALSE)
  
  
  # plotting for "cg14725810"
  number <- 1
  DMC_cg14725810 <- data.frame(RawValues = BetaMatrix_T1[DMC_Sites[number], ], 
                               PurifiedValues = c(tumor.purified[DMC_Sites[number], ], rep(NA, 5)),
                               Type = factor(ClinicalFile_T1$TYPE, levels = c("RLN", "FL", "DLBCL")))
  
  DMC_cg14725810_melt <- melt(DMC_cg14725810, id.vars = "Type")
  
  plot1 <- ggplot(DMC_cg14725810_melt, aes(x = Type, y =  value, color = variable)) +
    geom_boxplot(size=2, outlier.size = 3) +
    scale_y_continuous(name = "Methylation Beta Value") +
    ggtitle(rownames(DMC)[DMC_Sites[number]]) +
    #geom_jitter(aes(color = variable), size = 2, alpha = 0.95, width = 0.3, na.rm = FALSE) +
    theme_bw() +  scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Type of Data") +
    theme(text = element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    theme(aspect.ratio = 1, legend.position = "right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  number <- 101469 # Tried 8, 101469
  rownames(DMC)[DMC_Sites[number]] 
  DMC_cg16629967 <- data.frame(RawValues = BetaMatrix_T1[DMC_Sites[number], ], 
                               PurifiedValues = c(tumor.purified[DMC_Sites[number], ], rep(NA, 5)),
                               Type = factor(ClinicalFile_T1$TYPE, levels = c("RLN", "FL", "DLBCL")))
  
  DMC_cg16629967_melt <- melt(DMC_cg16629967, id.vars = "Type")
  plot2 <- ggplot(DMC_cg16629967_melt, aes(x = Type, y =  value, color = variable)) +
    geom_boxplot(size=2, outlier.size = 3) +
    scale_y_continuous(name = "Methylation Beta Value") +
    ggtitle(rownames(DMC)[DMC_Sites[number]]) +
    #geom_jitter(aes(color = variable), size = 2, alpha = 0.95, width = 0.3, na.rm = FALSE) +
    theme_bw() +  scale_color_manual(values = c("#fddbc7","#b2182b","red"), name = "Type of Data") +
    theme(text = element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    theme(aspect.ratio = 1, legend.position = "right", panel.background = element_rect(colour = "black", size=1.5),  axis.title =  element_text(face = "bold"))
  
  
}








