# Updated 3 Aug 2021
# Updated 10 June 2019
# Function: Survival analysis
# Author: Anjali Silva

# Input:
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns.
# ClinicalFile: File with patient sample names, and categories. The patient order will be made to be that 
#               of BetaMatrix.  
# SurvivalFile: File with patient sample names, and categories (AGE_AT_DIAGNOSIS, ANN_ARBOR_STAGE, etc). 
#               The patient order will be made to be that of BetaMatrix.  
# ClusterLabels: Vector of integers providing cluster labels.
# FigureGenerate: Produce images or not, options = "Yes" or "No"; default "Yes". 
# PNGorPDF: Output format of the image, options = "png" or "pdf"; default "png".

# Output:
# TTPSurvivalAdvanced: Outcome from using survfit() function on "Advance" patients for time to progression (TTP).
# TTPSurvivalLimited: Outcome from using survfit() function on "Limited" patients for time to progression (TTP).

# Visuals saved to img folder
# 15_SurvivalAnalysis_survival_TTP_advanced
# 15_SurvivalAnalysis_survival_TTP_limited

SurvivalAnalysis15 <- function(BetaMatrix, 
                               ClinicalFile,
                               SurvivalFile, 
                               ClusterLabels = "NA", 
                               FigureGenerate = "No", 
                               PNGorPDF = "png",
                               ImageName = date()) {
  
  # Loading needed packages
  library(ggplot2)
  library(RColorBrewer)
  library(survival)
  library(survminer)
  
  pathNow <- getwd()
  
  # Checking user input
  if (length(which(is.na(ClusterLabels) == TRUE)) >= 1) {
    stop("No cluster labels provided");
  }
  
  # Identify entries that correspond to FL from BetaMatrix, then find those entries in Surivival File
  Cluster <- ClusterLabels[which(substr(colnames(BetaMatrix), 4, 5) == "FL")]
  
  # Ordering Survival File by the same order as patients in Beta Matrix
  matchedIDs <- match(substr(colnames(BetaMatrix[, which(substr(colnames(BetaMatrix), 4, 5) == "FL")]), 1, 9), 
                      SurvivalFile$LY_FL_ID)
  if (sum(is.na(matchedIDs) ) > 0) {
    # if NAs are present
    SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile[matchedIDs[- which(is.na(matchedIDs) == TRUE)], ]
    Cluster <- ClusterLabels[-which(is.na(matchedIDs) == TRUE)]
  } else {
    SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile[matchedIDs, ]
  }
  
  # Ordering Clinical File by the same order as patients in Beta Matrix
  matchedIDs <- match(colnames(BetaMatrix[, which(substr(colnames(BetaMatrix), 4, 5) == "FL")]), 
                      ClinicalFile$SAMPLE_ID)
  if (sum(is.na(matchedIDs)) > 0) {
    ClinicalFile_OrderedbyBetaMatrixPatients <- ClinicalFile[matchedIDs[- which(is.na(matchedIDs) == TRUE)], ]
    Cluster <- ClusterLabels[- which(is.na(matchedIDs) == TRUE)]
  } else {
    ClinicalFile_OrderedbyBetaMatrixPatients <- ClinicalFile[matchedIDs, ]
  }

  # Check if all samples are in order
  if (sum(is.na(match(ClinicalFile_OrderedbyBetaMatrixPatients$SAMPLE_ID, 
                      colnames(BetaMatrix)))) > 0) {
    stop("Order of samples in ClinicalFile and BetaMatrix differ");}
  
  
  
  
  # Extracting values
  OS <- SurvivalFile_OrderedbyBetaMatrixPatients$OS
  CodeOS <- SurvivalFile_OrderedbyBetaMatrixPatients$CODE_OS
  TTP <- as.vector(SurvivalFile_OrderedbyBetaMatrixPatients$TTP)
  CodeTTP <- SurvivalFile_OrderedbyBetaMatrixPatients$CODE_TTP
  TTT <- as.vector(SurvivalFile_OrderedbyBetaMatrixPatients$TTT)
  CodeTRANSF <- SurvivalFile_OrderedbyBetaMatrixPatients$CODE_TRANSF # in place of CODE_TTT
  CodeDDS <- as.vector(SurvivalFile_OrderedbyBetaMatrixPatients$CODE_DSS) # disease specific survival 
  CodePFS <- as.vector(SurvivalFile_OrderedbyBetaMatrixPatients$CODE_PFS)
  Stage <- ClinicalFile_OrderedbyBetaMatrixPatients$STAGE
  ANNARBORstage <- SurvivalFile_OrderedbyBetaMatrixPatients$ANN_ARBOR_STAGE
  FLIPI <- SurvivalFile_OrderedbyBetaMatrixPatients$FLIPI_BINARY
  
  
  # if an entry is empty in OS, eliminate that from all categories
  # OSemptyentry <- which(is.na(OS) == TRUE)
  # if (length(OSemptyentry) >= 1) {
  #   OS <- OS[- OSemptyentry]
  #   CodeOS <- CodeOS[- OSemptyentry]
  #   TTP <- TTP[-OSemptyentry]
  #   CodeTTP <- CodeTTP[- OSemptyentry]
  #   Cluster <- Cluster[- OSemptyentry]
  #    SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile_OrderedbyBetaMatrixPatients[- OSemptyentry, ]
  #   Stage <- Stage[- OSemptyentry]
  #   ANNARBORstage <- ANNARBORstage[- OSemptyentry]
  #    FLIPI <- FLIPI[- OSemptyentry]
  #  }
  
  # if an entry is empty in CodeOS, eliminate that from all categories
  # CodeOSemptyentry <- which(is.na(CodeOS) == TRUE)
  # if (length(CodeOSemptyentry) >= 1) {
  #   CodeOS <- CodeOS[- CodeOSemptyentry]
  #   OS <- OS[-CodeOSemptyentry]
  #    TTP <- TTP[- CodeOSemptyentry]
  #   CodeTTP <- CodeTTP[- CodeOSemptyentry]
  #   Cluster <- Cluster[- CodeOSemptyentry]
  #   SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile_OrderedbyBetaMatrixPatients[- CodeOSemptyentry, ]
  #   Stage <- Stage[- CodeOSemptyentry]
  #   ANNARBORstage <- ANNARBORstage[- CodeOSemptyentry]
  #   FLIPI <- FLIPI[- CodeOSemptyentry]
  # }
  
  # if an entry is empty in TTP, eliminate that from all categories
  # TTPemptyentry <- which(is.na(TTP) == TRUE)
  # if (length(TTPemptyentry) >= 1) {
  #   CodeOS <- CodeOS[- TTPemptyentry]
  #   OS <- OS[-TTPemptyentry]
  #   TTP <- TTP[- TTPemptyentry]
  #   CodeTTP <- CodeTTP[- TTPemptyentry]
  #   Cluster <- Cluster[- TTPemptyentry]
  #   SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile_OrderedbyBetaMatrixPatients[- TTPemptyentry, ]
  #   Stage <- Stage[- TTPemptyentry]
  #   ANNARBORstage <- ANNARBORstage[- TTPemptyentry]
  #   FLIPI <- FLIPI[- TTPemptyentry]
  # }
  
  # if an entry is empty in CodeTTP, eliminate that from all categories
  # CodeTTPemptyentry <- which(is.na(CodeTTP) == TRUE)
  # if (length(CodeTTPemptyentry) >= 1) {
  #   CodeOS <- CodeOS[- CodeTTPemptyentry]
  #   OS <- OS[- CodeTTPemptyentry]
  #   TTP <- TTP[- CodeTTPemptyentry]
  #   CodeTTP <- CodeTTP[- CodeTTPemptyentry]
  #   Cluster <- Cluster[- CodeTTPemptyentry]    
  #   SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile_OrderedbyBetaMatrixPatients[- CodeTTPemptyentry, ]
  #  Stage <- Stage[- CodeTTPemptyentry]
  #  ANNARBORstage <- ANNARBORstage[- CodeTTPemptyentry]
  #   FLIPI <- FLIPI[- CodeTTPemptyentry]
  # }
  
  # if an entry is empty in Stage, eliminate that from all categories
  # Stageemptyentry <- which(is.na(Stage) == TRUE)
  # if (length(Stageemptyentry) >= 1) {
  #   CodeOS <- CodeOS[- Stageemptyentry]
  #   OS <- OS[- Stageemptyentry]
  #   TTP <- TTP[- Stageemptyentry]
  #   CodeTTP <- CodeTTP[- Stageemptyentry]
  #   Cluster <- Cluster[- Stageemptyentry]    
  #    SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile_OrderedbyBetaMatrixPatients[- Stageemptyentry, ]
  #   Stage <- Stage[- Stageemptyentry]
  #   ANNARBORstage <- ANNARBORstage[- Stageemptyentry]
  #   FLIPI <- FLIPI[- Stageemptyentry]
  # }
  
  # if an entry is empty in ANNARBORstage, eliminate that from all categories
  # ANN_ARBOR_Stageemptyentry <- which(is.na(ANNARBORstage) == TRUE)
  # if (length(ANN_ARBOR_Stageemptyentry) >= 1) {
  #   CodeOS <- CodeOS[- ANN_ARBOR_Stageemptyentry]
  #   OS <- OS[- ANN_ARBOR_Stageemptyentry]
  #   TTP <- TTP[- ANN_ARBOR_Stageemptyentry]
  #   CodeTTP <- CodeTTP[- ANN_ARBOR_Stageemptyentry]
  #    Cluster <- Cluster[- ANN_ARBOR_Stageemptyentry]    
  #   SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile_OrderedbyBetaMatrixPatients[- ANN_ARBOR_Stageemptyentry, ]
  #   Stage <- Stage[- ANN_ARBOR_Stageemptyentry]
  #   ANNARBORstage <- ANNARBORstage[- ANN_ARBOR_Stageemptyentry]
  #   FLIPI <- FLIPI[- ANN_ARBOR_Stageemptyentry]
  # }
  
  # if an entry is empty in FLIPI, eliminate that from all categories
  # FLIPIemptyentry <- which(is.na(FLIPI) == TRUE)
  # if (length(FLIPIemptyentry) >= 1) {
  #   CodeOS <- CodeOS[- FLIPIemptyentry]
  #   OS <- OS[- FLIPIemptyentry]
  #   TTP <- TTP[- FLIPIemptyentry]
  #   CodeTTP <- CodeTTP[- FLIPIemptyentry]
  #   Cluster <- Cluster[- FLIPIemptyentry]    
  #   SurvivalFile_OrderedbyBetaMatrixPatients <- SurvivalFile_OrderedbyBetaMatrixPatients[- FLIPIemptyentry, ]
  #   Stage <- Stage[- FLIPIemptyentry]
  #   ANNARBORstage <- ANNARBORstage[- FLIPIemptyentry]
  #   FLIPI <- FLIPI[- FLIPIemptyentry]
  # }
  
  # FLIPI entries as "LOW/INTERMEDIATE, "LOW_INTERMEDIATE", "LOW-INTERMEDIATE"
  FLIPI[which(FLIPI == "LOW-INTERMEDIATE")] <- "LOW/INTERMEDIATE"
  FLIPI[which(FLIPI == "LOW_INTERMEDIATE")] <- "LOW/INTERMEDIATE"
  FLIPI <- factor(FLIPI)
  
  # making a table with information, given cluster labels
  # replace 1 with "one", 2 with "two" and 3 with "three"
  # Cluster[which(Cluster == 1)] <- "one"
  # Cluster[which(Cluster == 2)] <- "two"
    
  # Cluster <- factor(Cluster, levels = c("one", "two"))
  # Cluster <- relevel(Cluster, ref = "two")
  # Cluster <- factor(Cluster)
  
  # Define colours:
  coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                      '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                      '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                      '#000075', '#808080')
  
  if (length(levels(Cluster)) == 2 || max(Cluster) == 2) {
    # Cluster <- relevel(Cluster, ref = "2")
    # setting colour palette to be two colours
    colourPalette = coloursBarPlot[1:2]
    # provide reference group as cluster three
  } else if (length(levels(Cluster)) == 3 || max(Cluster) == 3) {
      # Cluster[which(Cluster == 3)] <- "three"
      # setting colour palette to be three colours
      colourPalette = coloursBarPlot[1:3]
      # Cluster <- factor(Cluster, levels = c("one", "three", "two"))
      # Cluster <- relevel(Cluster, ref = "three")
      #Cluster <- relevel(Cluster, ref = "3")
  } else if (length(levels(Cluster)) == 4 || max(Cluster) == 4) {
      # Cluster[which(Cluster == 4)] <- "four"
      # setting colour palette to be four colours
      colourPalette = coloursBarPlot[1:4]
      # Cluster <- factor(Cluster, levels = c("one", "two", "three", "four"))
      # Cluster <- relevel(Cluster, ref = "four")
      # Cluster <- relevel(Cluster, ref = "4")
  } else if (length(levels(Cluster)) == 5 || max(Cluster) == 5) {
    colourPalette = coloursBarPlot[1:5]
  } else {
    stop("\n Only upto 6 clusters are supported.")
  }
    
    
    # tablePatientsAll <- as.data.frame(cbind(as.numeric(CodeOS), as.numeric(OS), as.numeric(CodeTTP), as.numeric(TTP), Cluster, as.numeric(ANNARBORstage), as.numeric(FLIPI)))
    #  rownames(tablePatientsAll) <- SurvivalFile_OrderedbyBetaMatrixPatients$LY_FL_ID
    # colnames(tablePatientsAll) <- c("CodeOS", "OS", "CodeTTP", "TTP", "Cluster", "ANNARBORstage", "FLIPI")
    
    # For all
    tablePatientsAll <- data.frame(event = as.numeric(CodeDDS),  #Code TTT is Code DDS
                                   time = as.numeric(TTT))
    rownames(tablePatientsAll) <- SurvivalFile_OrderedbyBetaMatrixPatients$LY_FL_ID
    colnames(tablePatientsAll) <- c("CodeDDS", "TTT")
    # based on TTP create survival curves for all patients
    survTTPall <- survival::survfit(Surv(time = TTT, 
                                         event = CodeDDS) ~ 1, 
                                    data = tablePatientsAll)
    # print(summary(survTTPall, times = c(5, 10)))
    cat("\n", colnames(tablePatientsAll)[1], "is generated \n")
    print(survTTPall)
    summary(survTTPall, times = c(5, 10))
    
  
  
    # By Cluster; for TTT see below 
    tablePatientsAll <- data.frame(event = as.numeric(CodeTTP), 
                                   time = as.numeric(TTP), 
                                   Cluster = Cluster)
    rownames(tablePatientsAll) <- SurvivalFile_OrderedbyBetaMatrixPatients$LY_FL_ID
    colnames(tablePatientsAll) <- c("event", "time", "Cluster")
    # cat("\n tablePatientsAll is generated \n")
    # print(head(tablePatientsAll))
    
    # based on TTP create survival curves for all patients
    survTTPall <- survival::survfit(Surv(time = time, 
                                        event = event) ~ Cluster, 
                                      data = tablePatientsAll)
    # print(summary(survTTPall, times = c(5, 10)))
    cat("\n", colnames(tablePatientsAll)[1], "is generated \n")
    print(survTTPall)
    summary(survTTPall, times = c(5, 10))
    
    
    
    # plot CodeTTP and TTP
    pathNow <- getwd()
    if (FigureGenerate == FigureGenerate) {
      
      tablePatientsAll <- data.frame(event = as.numeric(CodeTTP), 
                                     time = as.numeric(TTP), 
                                     Cluster = Cluster)
      survTTPall <- survival::survfit(Surv(time = time, 
                                           event = event) ~ Cluster, 
                                      data = tablePatientsAll)
      
      survTTPallPlot <- survminer::ggsurvplot(survTTPall, 
                                        # title = "Disease−specific Survival",
                                        # title = "Overall Survival",
                                        title = "Time to Progression",
                                        font.main = "bold", 
                                        xlim = c(0, 10),
                                        xlab = "Time (years)", 
                                        ylab = "Proportion of Patients Surviving",
                                        pval = TRUE, 
                                        pval.coord = c(0, 0.05), 
                                        conf.int = F,
                                        break.time.by = 1, 
                                        ncensor.plot = F,
                                        risk.table = TRUE, # Add risk table
                                        risk.table.pos = "out", 
                                        censor = T, 
                                        risk.table.y.text.col = TRUE,
                                        risk.table.col = "strata", # Change risk table color by groups
                                        palette = colourPalette,
                                        risk.table.height = 0.25,
                                        risk.table.y.text = FALSE,
                                        legend = "top",
                                        legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTP_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5)  
      # dev.off() 
      cat("\n Printed All image")
      
      
      
      # making a table with information, seperating stages
      ListAdvanced <- which(Stage == "ADVANCED")
      CodeTTPadvanced <- as.numeric(CodeTTP[ListAdvanced])
      TTPadvanced <- as.numeric(TTP[ListAdvanced])
      ClusterAdvanced <- as.numeric(Cluster[ListAdvanced])
      
      tableTTPadvanced <- data.frame(event = as.numeric(CodeTTPadvanced), 
                                     time = as.numeric(TTPadvanced), 
                                     Cluster = ClusterAdvanced)
      survTTPadvanced <- survival::survfit(Surv(time = time, 
                                           event = event) ~ Cluster, 
                                      data = tableTTPadvanced)
      
      survTTPadvancedPlot <- survminer::ggsurvplot(survTTPadvanced, 
                                              # title = "Disease−specific Survival",
                                              # title = "Overall Survival",
                                              title = "Time to Progression - Advanced",
                                              font.main = "bold", 
                                              xlim = c(0, 10),
                                              xlab = "Time (years)", 
                                              ylab = "Proportion of Patients Surviving",
                                              pval = TRUE, 
                                              pval.coord = c(0, 0.05), 
                                              conf.int = F,
                                              break.time.by = 1, 
                                              ncensor.plot = F,
                                              risk.table = TRUE, # Add risk table
                                              risk.table.pos = "out", 
                                              censor = T, 
                                              risk.table.y.text.col = TRUE,
                                              risk.table.col = "strata", # Change risk table color by groups
                                              palette = colourPalette,
                                              risk.table.height = 0.25,
                                              risk.table.y.text = FALSE,
                                              legend = "top",
                                              legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysisTTP_Advanced_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5)  
      

      
      ListLimited <- which(Stage == "LIMITED")
      CodeTTPLimited <- as.numeric(CodeTTP[ListLimited])
      TTPLimited <- as.numeric(TTP[ListLimited])
      ClusterLimited <- as.numeric(Cluster[ListLimited])
      
      tableTTPLimited <- data.frame(event = as.numeric(CodeTTPLimited), 
                                     time = as.numeric(TTPLimited), 
                                     Cluster = ClusterLimited)
      survTTPLimited <- survival::survfit(Surv(time = time, 
                                                event = event) ~ Cluster, 
                                           data = tableTTPLimited)
      
      survTTPLimitedPlot <- survminer::ggsurvplot(survTTPLimited, 
                                                   # title = "Disease−specific Survival",
                                                   # title = "Overall Survival",
                                                   title = "Time to Progression - Limited",
                                                   font.main = "bold", 
                                                   xlim = c(0, 10),
                                                   xlab = "Time (years)", 
                                                   ylab = "Proportion of Patients Surviving",
                                                   pval = TRUE, 
                                                   pval.coord = c(0, 0.05), 
                                                   conf.int = F,
                                                   break.time.by = 1, 
                                                   ncensor.plot = F,
                                                   risk.table = TRUE, # Add risk table
                                                   risk.table.pos = "out", 
                                                   censor = T, 
                                                   risk.table.y.text.col = TRUE,
                                                   risk.table.col = "strata", # Change risk table color by groups
                                                   palette = colourPalette,
                                                   risk.table.height = 0.25,
                                                   risk.table.y.text = FALSE,
                                                   legend = "top",
                                                   legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysisTTP_Limited_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5)  
      
    }
    
    # plot CodeOS and OS
    if (FigureGenerate == FigureGenerate) {
      
      tablePatientsOSAll <- data.frame(event = as.numeric(CodeOS), 
                                     time = as.numeric(OS), 
                                     Cluster = Cluster)
      survOSall <- survival::survfit(Surv(time = time, 
                                           event = event) ~ Cluster, 
                                      data = tablePatientsOSAll)

      kmTTPall <- survminer::ggsurvplot(survOSall, 
                                        # title = "Disease−specific Survival",
                                        title = "Overall Survival",
                                        # title = "Time to Progression",
                                        font.main = "bold", 
                                        xlim = c(0, 10),
                                        xlab = "Time (years)", 
                                        ylab = "Proportion of Patients Surviving",
                                        pval = TRUE, 
                                        pval.coord = c(0, 0.05), 
                                        conf.int = F,
                                        break.time.by = 1, 
                                        ncensor.plot = F,
                                        risk.table = TRUE, # Add risk table
                                        risk.table.pos = "out", 
                                        censor = T, 
                                        risk.table.y.text.col = TRUE,
                                        risk.table.col = "strata", # Change risk table color by groups
                                        palette = colourPalette,
                                        risk.table.height = 0.25,
                                        risk.table.y.text = FALSE,
                                        legend = "top",
                                        legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_OS_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5)  
      # dev.off() 
      
      ListAdvanced <- which(Stage == "ADVANCED")
      CodeOSadvanced <- as.numeric(CodeOS[ListAdvanced])
      OSadvanced <- as.numeric(OS[ListAdvanced])
      ClusterAdvanced <- as.numeric(Cluster[ListAdvanced])
      
      tableOSadvanced <- data.frame(event = as.numeric(CodeOSadvanced), 
                                     time = as.numeric(OSadvanced), 
                                     Cluster = ClusterAdvanced)
      survOSadvanced <- survival::survfit(Surv(time = time, 
                                                event = event) ~ Cluster, 
                                           data = tableOSadvanced)
      
      survOSadvancedPlot <- survminer::ggsurvplot(survOSadvanced, 
                                        # title = "Disease−specific Survival",
                                        title = "Overall Survival - Advanced",
                                        # title = "Time to Progression",
                                        font.main = "bold", 
                                        xlim = c(0, 10),
                                        xlab = "Time (years)", 
                                        ylab = "Proportion of Patients Surviving",
                                        pval = TRUE, 
                                        pval.coord = c(0, 0.05), 
                                        conf.int = F,
                                        break.time.by = 1, 
                                        ncensor.plot = F,
                                        risk.table = TRUE, # Add risk table
                                        risk.table.pos = "out", 
                                        censor = T, 
                                        risk.table.y.text.col = TRUE,
                                        risk.table.col = "strata", # Change risk table color by groups
                                        palette = colourPalette,
                                        risk.table.height = 0.25,
                                        risk.table.y.text = FALSE,
                                        legend = "top",
                                        legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_OS_Advanced_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5) 
      
      
      ListLimited <- which(Stage == "LIMITED")
      CodeOSlimited <- as.numeric(CodeOS[ListLimited])
      OSlimited <- as.numeric(OS[ListLimited])
      ClusterLimited <- as.numeric(Cluster[ListLimited])
      
      tableOSlimited <- data.frame(event = as.numeric(CodeOSlimited), 
                                    time = as.numeric(OSlimited), 
                                    Cluster = ClusterLimited)
      survOSlimited <- survival::survfit(Surv(time = time, 
                                               event = event) ~ Cluster, 
                                          data = tableOSlimited)
      
      survOSlimitedPlot <- survminer::ggsurvplot(survOSlimited, 
                                       # title = "Disease−specific Survival",
                                       title = "Overall Survival - Limited",
                                       # title = "Time to Progression",
                                       font.main = "bold", 
                                       xlim = c(0, 10),
                                       xlab = "Time (years)", 
                                       ylab = "Proportion of Patients Surviving",
                                       pval = TRUE, 
                                       pval.coord = c(0, 0.05), 
                                       conf.int = F,
                                       break.time.by = 1, 
                                       ncensor.plot = F,
                                       risk.table = TRUE, # Add risk table
                                       risk.table.pos = "out", 
                                       censor = T, 
                                       risk.table.y.text.col = TRUE,
                                       risk.table.col = "strata", # Change risk table color by groups
                                       palette = colourPalette,
                                       risk.table.height = 0.25,
                                       risk.table.y.text = FALSE,
                                       legend = "top",
                                       legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_OS_Limited_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5) 
      
    }
    
    # plot CODE_DSS and OS    
    if (FigureGenerate == FigureGenerate) {
      
      tablePatientsDDSAll <- data.frame(event = as.numeric(CodeDDS), 
                                     time = as.numeric(OS), 
                                     Cluster = Cluster)
      survDDSall <- survival::survfit(Surv(time = time, 
                                           event = event) ~ Cluster, 
                                      data = tablePatientsDDSAll)

      
      kmTTPall <- survminer::ggsurvplot(survDDSall, 
                                        title = "Disease−specific Survival",
                                        # title = "Overall Survival",
                                        # title = "Time to Progression",
                                        font.main = "bold", 
                                        xlim = c(0, 10),
                                        xlab = "Time (years)", 
                                        ylab = "Proportion of Patients Surviving",
                                        pval = TRUE, 
                                        pval.coord = c(0, 0.05), 
                                        conf.int = F,
                                        break.time.by = 1, 
                                        ncensor.plot = F,
                                        risk.table = TRUE, # Add risk table
                                        risk.table.pos = "out", 
                                        censor = T, 
                                        risk.table.y.text.col = TRUE,
                                        risk.table.col = "strata", # Change risk table color by groups
                                        palette = colourPalette,
                                        risk.table.height = 0.25,
                                        risk.table.y.text = FALSE,
                                        legend = "top",
                                        legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_DSS_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5)      
      # dev.off() 
      cat("\n Printed All image")
      
      
      ListAdvanced <- which(Stage == "ADVANCED")
      CodeDDSadvanced <- as.numeric(CodeDDS[ListAdvanced])
      OSadvanced <- as.numeric(OS[ListAdvanced])
      ClusterAdvanced <- as.numeric(Cluster[ListAdvanced])
      
      tableDDSadvanced <- data.frame(event = as.numeric(CodeDDSadvanced), 
                                    time = as.numeric(OSadvanced), 
                                    Cluster = ClusterAdvanced)
      survDDSadvanced <- survival::survfit(Surv(time = time, 
                                               event = event) ~ Cluster, 
                                          data = tableDDSadvanced)
      
      survDDSadvancedPlot <- survminer::ggsurvplot(survDDSadvanced, 
                                        title = "Disease−specific Survival - Advanced",
                                        # title = "Overall Survival",
                                        # title = "Time to Progression",
                                        font.main = "bold", 
                                        xlim = c(0, 10),
                                        xlab = "Time (years)", 
                                        ylab = "Proportion of Patients Surviving",
                                        pval = TRUE, 
                                        pval.coord = c(0, 0.05), 
                                        conf.int = F,
                                        break.time.by = 1, 
                                        ncensor.plot = F,
                                        risk.table = TRUE, # Add risk table
                                        risk.table.pos = "out", 
                                        censor = T, 
                                        risk.table.y.text.col = TRUE,
                                        risk.table.col = "strata", # Change risk table color by groups
                                        palette = colourPalette,
                                        risk.table.height = 0.25,
                                        risk.table.y.text = FALSE,
                                        legend = "top",
                                        legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_DSS_Advanced_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5) 
      
      
      ListLimited <- which(Stage == "LIMITED")
      CodeDDSlimited <- as.numeric(CodeDDS[ListLimited])
      OSlimited <- as.numeric(OS[ListLimited])
      ClusterLimited <- as.numeric(Cluster[ListLimited])
      
      tableDDSlimited <- data.frame(event = as.numeric(CodeDDSlimited), 
                                     time = as.numeric(OSlimited), 
                                     Cluster = ClusterLimited)
      survDDSlimited <- survival::survfit(Surv(time = time, 
                                                event = event) ~ Cluster, 
                                           data = tableDDSlimited)
      
      tableDDSlimitedPlot <- survminer::ggsurvplot(survDDSlimited, 
                                                   title = "Disease−specific Survival - Limited",
                                                   # title = "Overall Survival",
                                                   # title = "Time to Progression",
                                                   font.main = "bold", 
                                                   xlim = c(0, 10),
                                                   xlab = "Time (years)", 
                                                   ylab = "Proportion of Patients Surviving",
                                                   pval = TRUE, 
                                                   pval.coord = c(0, 0.05), 
                                                   conf.int = F,
                                                   break.time.by = 1, 
                                                   ncensor.plot = F,
                                                   risk.table = TRUE, # Add risk table
                                                   risk.table.pos = "out", 
                                                   censor = T, 
                                                   risk.table.y.text.col = TRUE,
                                                   risk.table.col = "strata", # Change risk table color by groups
                                                   palette = colourPalette,
                                                   risk.table.height = 0.25,
                                                   risk.table.y.text = FALSE,
                                                   legend = "top",
                                                   legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_DSS_Limited_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5) 
    }
    
    # plot CodeTRANSF and TTT       
    if (FigureGenerate == FigureGenerate) { #  plot time to transformation as cumulative incidence or similar metric.
      tablePatientsAll <- data.frame(event = as.numeric(CodeTRANSF), 
                                     time = as.numeric(TTT), 
                                     Cluster = Cluster)
      survTTTall <- survival::survfit(Surv(time = time, 
                                           event = event) ~ Cluster, 
                                      data = tablePatientsAll)

      # grDevices::dev.off()
      # if (PNGorPDF == "png") {
      #  grDevices::png(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTT_", ImageName, ".", PNGorPDF))
      # }
      # if (PNGorPDF == "pdf") {
      #   grDevices::pdf(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTT_", ImageName, ".", PNGorPDF))
      # }
      
      kmTTPall <- survminer::ggsurvplot(survTTTall, 
                                        title = "Cumulative Incidence of Transformation",
                                        data = tablePatientsAll, 
                                        font.main = "bold", 
                                        xlim = c(0, 10),
                                        xlab = "Time (years)", 
                                        ylab = "Cumulative Hazard",
                                        pval = TRUE, 
                                        fun = "cumhaz", 
                                        ylim = c(0, 0.4),
                                        pval.coord = c(0, 0.4), 
                                        conf.int = F,
                                        break.time.by = 1, 
                                        ncensor.plot = F,
                                        risk.table = TRUE, # Add risk table
                                        risk.table.pos = "out", 
                                        censor = T, 
                                        risk.table.y.text.col = TRUE,
                                        risk.table.col = "strata", # Change risk table color by groups
                                        palette = colourPalette,
                                        risk.table.height = 0.25,
                                        risk.table.y.text = FALSE,
                                        legend = "top",
                                        legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTT_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5)
      
      
      ListAdvanced <- which(Stage == "ADVANCED")
      CodeTRANSFadvanced <- as.numeric(CodeTRANSF[ListAdvanced])
      TTTadvanced <- as.numeric(TTT[ListAdvanced])
      ClusterAdvanced <- as.numeric(Cluster[ListAdvanced])
      
      tableTTTadvanced <- data.frame(event = as.numeric(CodeTRANSFadvanced), 
                                     time = as.numeric(TTTadvanced), 
                                     Cluster = ClusterAdvanced)
      survTTTadvanced <- survival::survfit(Surv(time = time, 
                                                event = event) ~ Cluster, 
                                           data = tableTTTadvanced)
      
      survTTTadvancedPlot <- survminer::ggsurvplot(survTTTadvanced, 
                                        title = "Cumulative Incidence of Transformation - Advanced",
                                        data = tablePatientsAll, 
                                        font.main = "bold", 
                                        xlim = c(0, 10),
                                        xlab = "Time (years)", 
                                        ylab = "Cumulative Hazard",
                                        pval = TRUE, 
                                        fun = "cumhaz", 
                                        ylim = c(0, 0.4),
                                        pval.coord = c(0, 0.4), 
                                        conf.int = F,
                                        break.time.by = 1, 
                                        ncensor.plot = F,
                                        risk.table = TRUE, # Add risk table
                                        risk.table.pos = "out", 
                                        censor = T, 
                                        risk.table.y.text.col = TRUE,
                                        risk.table.col = "strata", # Change risk table color by groups
                                        palette = colourPalette,
                                        risk.table.height = 0.25,
                                        risk.table.y.text = FALSE,
                                        legend = "top",
                                        legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTT_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5)
      
      
      ListLimited <- which(Stage == "LIMITED")
      CodeTRANSFlimited <- as.numeric(CodeTRANSF[ListLimited])
      TTTlimited <- as.numeric(TTT[ListLimited])
      ClusterLimited <- as.numeric(Cluster[ListLimited])
      
      tableTTTlimited <- data.frame(event = as.numeric(CodeTRANSFlimited), 
                                     time = as.numeric(TTTlimited), 
                                     Cluster = ClusterLimited)
      survTTTlimited <- survival::survfit(Surv(time = time, 
                                                event = event) ~ Cluster, 
                                           data = tableTTTlimited)
      
      survTTTlimitedPlot <- survminer::ggsurvplot(survTTTlimited, 
                                                   title = "Cumulative Incidence of Transformation - Limited",
                                                   data = tablePatientsAll, 
                                                   font.main = "bold", 
                                                   xlim = c(0, 10),
                                                   xlab = "Time (years)", 
                                                   ylab = "Cumulative Hazard",
                                                   pval = TRUE, 
                                                   fun = "cumhaz", 
                                                   ylim = c(0, 0.4),
                                                   pval.coord = c(0, 0.4), 
                                                   conf.int = F,
                                                   break.time.by = 1, 
                                                   ncensor.plot = F,
                                                   risk.table = TRUE, # Add risk table
                                                   risk.table.pos = "out", 
                                                   censor = T, 
                                                   risk.table.y.text.col = TRUE,
                                                   risk.table.col = "strata", # Change risk table color by groups
                                                   palette = colourPalette,
                                                   risk.table.height = 0.25,
                                                   risk.table.y.text = FALSE,
                                                   legend = "top",
                                                   legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTT_", ImageName, ".", PNGorPDF), 
                      width = 8, height = 5.5)
      
      
      
    }

    
    return(invisible(NULL))
  }


 StageWiseSurvivalAnalysis15 <- function() {
    
    # Added pairwise difference calculation on 3 March 2020
    # Citation: http://www.sthda.com/english/wiki/survminer-0-3-0#pairwise-comparisons-for-survival-curves
    # Pairwise survdiff
    TTPallPairwiseDiff <- survminer::pairwise_survdiff(Surv(time = as.numeric(TTT), 
                                                             event = as.numeric(CodeTRANSF)) ~ Cluster, 
                                                       data = tablePatientsAll)
    cat("\n Pariwise Differences for all patients:\n")
    print(TTPallPairwiseDiff)
    # Symbolic number coding
    # symnum(TTPallPairwiseDiff$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    #        symbols = c("****", "***", "**", "*", "+", " "),
    #        abbr.colnames = FALSE, na = "")
    

    
    
    # making a table with information, seperating stages
    ListAdvanced <- which(Stage == "ADVANCED")
    ListLimited <- which(Stage == "LIMITED")
    
    
    
    CodeDDSadvanced <- as.numeric(CodeDDS[ListAdvanced])
    CodePFSadvanced <- as.numeric(CodePFS[ListAdvanced])
    OSadvanced <- as.numeric(OS[ListAdvanced])
    CodeOSadvanced <- as.numeric(CodeOS[ListAdvanced])
    TTPadvanced <- as.numeric(TTP[ListAdvanced])
    CodeTTPadvanced <- as.numeric(CodeTTP[ListAdvanced])
    TTTadvanced <- as.numeric(TTT[ListAdvanced])
    CodeTTTadvanced <- as.numeric(CodeTRANSF[ListAdvanced])
    ClusterAdvanced  <- Cluster[ListAdvanced]
    ANNARBORStageadvanced <- ANNARBORstage[ListAdvanced]
    FLIPIadvanced  <- FLIPI[ListAdvanced]
    TablePatientsAdvanced <- as.data.frame(cbind(as.numeric(CodeDDS[ListAdvanced]),
                                                 as.numeric(CodePFS[ListAdvanced]),
                                                 as.numeric(OS[ListAdvanced]), 
                                                 as.numeric(CodeOS[ListAdvanced]),
                                                 as.numeric(TTP[ListAdvanced]), 
                                                 as.numeric(CodeTTP[ListAdvanced]), 
                                                 as.numeric(TTT[ListAdvanced]), 
                                                 as.numeric(CodeTRANSF[ListAdvanced]),                                                 
                                                 Cluster[ListAdvanced], 
                                                 as.numeric(ANNARBORstage[ListAdvanced]), 
                                                 as.numeric(FLIPI[ListAdvanced])))
    rownames(TablePatientsAdvanced) <- ClinicalFile_OrderedbyBetaMatrixPatients$SAMPLE_ID[ListAdvanced]
    colnames(TablePatientsAdvanced) <- c("CodeDDSadvanced", "CodePFSadvanced", "OSadvanced", 
                                         "CodeOSadvanced", "TTPadvanced", "CodeTTPadvanced", 
                                         "TTTadvanced", "CodeTTTadvanced", "ClusterAdvanced", 
                                         "ANNARBORStageadvanced", "FLIPIadvanced")
    # cat("\n Finished TablePatientsAdvanced \n")
    
    
    CodeDDSlimited <- as.numeric(CodeDDS[ListLimited])
    CodePFSlimited <- as.numeric(CodePFS[ListLimited])
    OSlimited <- as.numeric(OS[ListLimited])
    CodeOSlimited <- as.numeric(CodeOS[ListLimited])
    TTPlimited <- as.numeric(TTP[ListLimited])
    CodeTTPlimited <- as.numeric(CodeTTP[ListLimited])
    TTTlimited <- as.numeric(TTT[ListLimited])
    CodeTTTlimited <- as.numeric(CodeTRANSF[ListLimited])
    Clusterlimited  <- Cluster[ListLimited]
    ANNARBORStagelimited <- ANNARBORstage[ListLimited]
    FLIPIlimited  <- FLIPI[ListLimited]
    TablePatientsLimited <- as.data.frame(cbind(as.numeric(CodeDDS[ListLimited]),
                                                 as.numeric(CodePFS[ListLimited]),
                                                 as.numeric(OS[ListLimited]), 
                                                 as.numeric(CodeOS[ListLimited]),
                                                 as.numeric(TTP[ListLimited]), 
                                                 as.numeric(CodeTTP[ListLimited]), 
                                                 as.numeric(TTT[ListLimited]), 
                                                 as.numeric(CodeTRANSF[ListLimited]),                                                 
                                                 Cluster[ListLimited], 
                                                 as.numeric(ANNARBORstage[ListLimited]), 
                                                 as.numeric(FLIPI[ListLimited])))
    rownames(TablePatientsLimited) <- ClinicalFile_OrderedbyBetaMatrixPatients$SAMPLE_ID[ListLimited]
    colnames(TablePatientsLimited) <- c("CodeDDSlimited", "CodePFSlimited", "OSlimited", 
                                         "CodeOSlimited", "TTPlimited", "CodeTTPlimited", 
                                         "TTTlimited", "CodeTTTlimited", "Clusterlimited", 
                                         "ANNARBORStagelimited", "FLIPIlimited")
    # cat("\n Finished tablePatientsLimited \n")
    
  
    # based on OS create survival curves for advanced patients
    # surv_os_advanced <- survival::survfit(Surv(time = OSadvanced, 
    #                                            event = CodeOSadvanced)~ClusterAdvanced, 
    #                                       data = TablePatientsAdvanced)
    # print(summary(surv_os_advanced))
    # log rank test
    # survdiff_os_advanced <- survival::survdiff(Surv(time = OSadvanced, 
    #                                                event = CodeOSadvanced)~ClusterAdvanced, 
    #                                           data = TablePatientsAdvanced) 
    # pval1Advanced <- 1 - stats::pchisq(survdiff_os_advanced$chisq, length(survdiff_os_advanced$n) - 1)
    # pval1Advanced <- round(pval1Advanced, 3) # 0.44
    
    # based on OS create survival curves for limited patients
    # surv_os_limited <- survival::survfit(Surv(time = OSlimited, 
    #                                           event = CodeOSlimited)~ClusterLimited, 
    #                                      data = tablePatientsLimited)
    # print(summary(surv_os_limited))
    # log rank test
    # survdiff_os_limited <- survival::survdiff(Surv(time = OSlimited, 
    #                                                event = CodeOSlimited)~ClusterLimited, 
    #                                           data = tablePatientsLimited) 
    # pval1Limited <- 1 - stats::pchisq(survdiff_os_limited$chisq,length(survdiff_os_limited$n) - 1)
    # pval1Limited <- round(pval1Limited, 3) # 0.44
    
    # based on TTP create survival curves for advanced patients
    survTTPadvanced <- survival::survfit(Surv(time = as.numeric(TTTadvanced), 
                                              event = as.numeric(CodeTTTadvanced)) ~ ClusterAdvanced, 
                                         data = TablePatientsAdvanced)
    # print(summary(survTTPadvanced))
    # log rank test
    survdiffTTPadvanced <- survival::survdiff(Surv(time = as.numeric(TTTadvanced), 
                                                   event = as.numeric(CodeTTTadvanced)) ~ ClusterAdvanced, 
                                              data = TablePatientsAdvanced) 
    pval1Advanced <- 1 - stats::pchisq(survdiffTTPadvanced$chisq, 
                                        length(survdiffTTPadvanced$n) - 1)
    pval1Advanced <- round(pval1Advanced, 3) # 0.478
    
    
    
  
    if (FigureGenerate == "Yes") {
      cat("\n Printing advanced stage image")
      if (PNGorPDF == "png") {
        grDevices::png(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTP_advanced.", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTP_advanced.", PNGorPDF))
      }
      
      kmTTPadvanced <- survminer::ggsurvplot(survTTPadvanced, 
                                               title = "59 Samples - TTT by SNF advanced stage: CodeTRANSF with TTT",
                                               font.main = "bold", 
                                               xlab = "TTT (years)", 
                                               ylab = "Proportion of patients surviving",
                                               pval = TRUE, 
                                               pval.coord = c(0, 0.05), 
                                               conf.int = F,
                                               break.time.by = 1, 
                                               ncensor.plot = F,
                                               risk.table = TRUE, # Add risk table
                                               risk.table.pos = "out", 
                                               censor = T, 
                                               risk.table.y.text.col = TRUE,
                                               risk.table.col = "strata", # Change risk table color by groups
                                               #linetype = "strata", # Change line type by groups
                                               #surv.median.line = "hv", # Specify median survival
                                               #ggtheme = theme_bw(), # Change ggplot2 theme
                                               palette = colourPalette,
                                               risk.table.height = 0.25,
                                               risk.table.y.text = FALSE,
                                               legend = "top",
                                               legend.title = "Cluster")
                                               #,legend.labs = c("C=1","C=2"))
      grDevices::dev.off() 
    }
    
    if (FigureGenerate == "Yes") { #  plot time to transformation as cumulative incidence or similar metric.
      survTTPadvanced <- survival::survfit(Surv(time = as.numeric(TTTadvanced), 
                                           event = as.numeric(CodeTTTadvanced)) ~ ClusterAdvanced, 
                                           data = TablePatientsAdvanced)
      
      kmTTPadvanced <- survminer::ggsurvplot(survTTPadvanced, 
                                             title = "59 Samples - TTT by SNF advanced stage: CodeTRANSF with TTT",
                                             data = tablePatientsAll, 
                                             font.main = "bold", 
                                             xlab = "TTT (years)", 
                                             ylab = "Cumulative Hazard",
                                             pval = TRUE, 
                                             fun = "cumhaz", 
                                             ylim = c(0, 1),
                                             pval.coord = c(0, 1), 
                                             conf.int = F,
                                             break.time.by = 1, 
                                             ncensor.plot = F,
                                             risk.table = TRUE, # Add risk table
                                             risk.table.pos = "out", 
                                             censor = T, 
                                             risk.table.y.text.col = TRUE,
                                             risk.table.col = "strata", # Change risk table color by groups
                                             palette = colourPalette,
                                             risk.table.height = 0.25,
                                             risk.table.y.text = FALSE,
                                             legend = "top",
                                             legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTT_Advanced_", ImageName, ".png"), 
                      width = 8, height = 5.5)
      
    }
    
    
    # Added pairwise difference calculation on 3 March 2020
    # Citation: http://www.sthda.com/english/wiki/survminer-0-3-0#pairwise-comparisons-for-survival-curves
    # Pairwise survdiff
    TTPadvancedPairwiseDiff <- survminer::pairwise_survdiff(Surv(time = as.numeric(TTTadvanced), 
                                                                 event = as.numeric(CodeTTTadvanced)) ~ ClusterAdvanced, 
                                                            data = TablePatientsAdvanced)
    cat("\n Pariwise differences for advanced-stage patients:\n")
    print(TTPadvancedPairwiseDiff)
    # Symbolic number coding
    # symnum(TTPadvancedPairwiseDiff$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    #        symbols = c("****", "***", "**", "*", "+", " "),
    #        abbr.colnames = FALSE, na = "")
    

  
  
  
    # based on TTP create survival curves for limited patients
    survTTPlimited <- survfit(Surv(time = as.numeric(TTTlimited), 
                                   event = as.numeric(CodeTTTlimited)) ~ ClusterLimited, 
                              data = tablePatientsLimited)
    # summary(survTTPlimited)  
    # log rank test
    survdiffTTPlimited <- survdiff(Surv(time = as.numeric(TTTlimited), 
                                          time2 = as.numeric(CodeTTTlimited)) ~ ClusterLimited, 
                                     data = tablePatientsLimited) 
    pval1Limited <- 1 - stats::pchisq(q = survdiffTTPlimited$chisq, 
                                       df = length(survdiffTTPlimited$n) - 1)
    pval1Limited <- round(pval1Limited, 3) # 0.478
    
    if (FigureGenerate == "Yes") {
      cat("\n Printing limited stage image")
      if (PNGorPDF == "png") {
        grDevices::png(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTP_limited.", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTP_limited.", PNGorPDF))
      }
      
      kmTTPlimited <- survminer::ggsurvplot(survTTPlimited, 
                                           title = "42 Samples - TTT by SNF limited stage: CodeTRANSF with TTT",
                                           font.main = "bold", 
                                           xlab = "TTT (years)", 
                                           ylab = "Proportion of patients surviving",
                                           pval = TRUE, 
                                           pval.coord = c(0, 0.05), 
                                           conf.int = F,
                                           break.time.by = 1, 
                                           ncensor.plot = F,
                                           risk.table = TRUE, # Add risk table
                                           risk.table.pos = "out", 
                                           censor = T, 
                                           risk.table.y.text.col = TRUE,
                                           risk.table.col = "strata", # Change risk table color by groups
                                           #linetype = "strata", # Change line type by groups
                                           #surv.median.line = "hv", # Specify median survival
                                           #ggtheme = theme_bw(), # Change ggplot2 theme
                                           palette = colourPalette,
                                           risk.table.height = 0.25,
                                           risk.table.y.text = FALSE,
                                           legend = "top",
                                           legend.title = "Cluster")
                                   #,legend.labs = c("C=1","C=2"))
      grDevices::dev.off() 
    }
    
    if (FigureGenerate == "Yes") { #  plot time to transformation as cumulative incidence or similar metric.
      survTTPlimited <- survival::survfit(Surv(time = as.numeric(TTTlimited), 
                                                event = as.numeric(CodeTTTlimited)) ~ Clusterlimited, 
                                           data = TablePatientsLimited)
      
      kmTTPadvanced <- survminer::ggsurvplot(survTTPlimited, 
                                             title = "42 Samples - TTT by SNF limited stage: CodeTRANSF with TTT",
                                             data = tablePatientsAll, 
                                             font.main = "bold", 
                                             xlab = "TTT (years)", 
                                             ylab = "Cumulative Hazard",
                                             pval = TRUE, 
                                             fun = "cumhaz", 
                                             ylim = c(0, 1),
                                             pval.coord = c(0, 1), 
                                             conf.int = F,
                                             break.time.by = 1, 
                                             ncensor.plot = F,
                                             risk.table = TRUE, # Add risk table
                                             risk.table.pos = "out", 
                                             censor = T, 
                                             risk.table.y.text.col = TRUE,
                                             risk.table.col = "strata", # Change risk table color by groups
                                             palette = colourPalette,
                                             risk.table.height = 0.25,
                                             risk.table.y.text = FALSE,
                                             legend = "top",
                                             legend.title = "Cluster")
      ggplot2::ggsave(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTT_Limited_", ImageName, ".png"), 
                      width = 8, height = 5.5)
    
    }
    
    
    # Added pairwise difference calculation on 3 March 2020
    # Citation: http://www.sthda.com/english/wiki/survminer-0-3-0#pairwise-comparisons-for-survival-curves
    # Pairwise survdiff
    TTPlimitedPairwiseDiff <- survminer::pairwise_survdiff(Surv(time = as.numeric(TTPlimited), 
                                                                event = as.numeric(CodeTTPlimited)) ~ ClusterLimited, 
                                                           data = tablePatientsLimited)
    cat("\n Pariwise differences for limited-stage patients:\n")
    print(TTPlimitedPairwiseDiff)
    # Symbolic number coding
    # symnum(TTPlimitedPairwiseDiff$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    #        symbols = c("****", "***", "**", "*", "+", " "),
    #        abbr.colnames = FALSE, na = "")
    
    
    # Coxp model based on TTP adjusted for cluster labels only for all patients 
    NotRun <- function() {
      coxph_TTP_all <- coxph(Surv(time = as.numeric(TTP), 
                                  time2 = as.numeric(CodeTTP)) ~ Cluster, 
                             data = tablePatientsAll)
      # cat("\n Cox model based on TTP adjusted for cluster labels only for all patients \n")
      # summary(coxph_TTP_all)
      # cat("\n Test the proportional hazards assumption for a Cox regression model fit \n")
      # cox.zph(coxph_TTP_all)
      
      # Coxp model based on TTP adjusted for cluster labels and STAGE for all patients 
      coxph_TTP_all_stage <- coxph(Surv(time = as.numeric(TTP), 
                                        time2 = as.numeric(CodeTTP)) ~ Cluster + Stage, 
                                   data = tablePatientsAll)
      # cat("\n Cox model based on TTP adjusted for cluster labels and stage for all patients \n")
      # summary(coxph_TTP_all_stage)
      # cat("\n Test the proportional hazards assumption for a Cox regression model fit \n")
      # cox.zph(coxph_TTP_all_stage)
      
      # Coxp model based on TTP adjusted for cluster labels and FLIPI for all patients 
      coxph_TTP_all_flipi <- coxph(Surv(time = as.numeric(TTP), 
                                        time2 = as.numeric(CodeTTP)) ~ Cluster + FLIPI, 
                                   data = tablePatientsAll)
      # cat("\n Cox model based on TTP adjusted for cluster labels and flipi for all patients \n")
      # summary(coxph_TTP_all_flipi)
      # cat("\n Test the proportional hazards assumption for a Cox regression model fit \n")
      # cox.zph(coxph_TTP_all_flipi)
      
      # Coxp model based on TTP adjusted for cluster labels and FLIPI + STAGE for all patients 
      coxph_TTP_all_flipi_stage <- coxph(Surv(time = as.numeric(TTP), 
                                              time2 = as.numeric(CodeTTP)) ~ Cluster + 
                                              FLIPI + Stage, data = tablePatientsAll)
      # cat("\n Cox model based on TTP adjusted for cluster labels, flipi and stage for all patients \n")
      # summary(coxph_TTP_all_flipi_stage)
      # cat("\n Test the proportional hazards assumption for a Cox regression model fit \n")
      # cox.zph(coxph_TTP_all_flipi_stage)
      
      # Coxp model based on TTP adjusted for stage for advanced patients
      coxph_TTP_advanced <- coxph(Surv(time = TTPadvanced, time2 = CodeTTPadvanced) ~ ClusterAdvanced + 
                                       FLIPIadvanced , data = TablePatientsAdvanced)
      # cat("\n Cox model based on TTP adjusted for stage for advanced patients \n")
      # summary(coxph_TTP_advanced)
      # cat("\n Test the proportional hazards assumption for a Cox regression model fit \n")
      # cox.zph(coxph_TTP_advanced)
    
      # Coxp model based on TTP adjusted for flipi + stage for limited patients
      coxph_TTP_limited <- coxph(Surv(time = TTPlimited, time2 = CodeTTPlimited) ~ ClusterLimited + 
                                   FLIPIlimited + ANNARBORStagelimited , data = tablePatientsLimited)
      # cat("\n Cox model based on TTP adjusted for stage for limited patients \n")
      # print(summary(coxph_TTP_limited))
      # cat("\n Test the proportional hazards assumption for a Cox regression model fit \n")
      # print(cox.zph(coxph_TTP_limited))
      
      # Coxp model based on TTP adjusted for flipi for limited patients
      coxph_TTP_flipi_limited <- coxph(Surv(time = TTPlimited, time2 = CodeTTPlimited) ~ ClusterLimited + 
                                         FLIPIlimited, data = tablePatientsLimited)
      # cat("\n Cox model based on TTP adjusted for flipi only for limited patients \n")
      # print(summary(coxph_TTP_flipi_limited))
      # cat("\n Test the proportional hazards assumption for a Cox regression model fit \n")
      # print(cox.zph(coxph_TTP_flipi_limited))
    }
    
    
    # An alternate - drawing by cluster
    
    tablePatientsAll <- as.data.frame(cbind(as.numeric(CodeTTP), 
                                            as.numeric(TTP), 
                                            Stage,
                                            Cluster))
    rownames(tablePatientsAll) <- SurvivalFile_OrderedbyBetaMatrixPatients$LY_FL_ID
    colnames(tablePatientsAll) <- c("CodeTTP", "TTP", "Stage", "Cluster")
    # cat("\n tablePatientsAll is generated \n")
    # print(head(tablePatientsAll))
    
    # based on TTP create survival curves for all patients
    survTTPallStage <- survival::survfit(Surv(time = as.numeric(TTP), 
                                         event = as.numeric(CodeTTP)) ~ Stage, 
                                    data = tablePatientsAll)
    # print(summary(survTTPall, times = c(5, 10)))
    cat("\n survTTPall is generated \n")
    print(survTTPall)
    
    # making a table with information, seperating stages
    ListOne <- which(Cluster == 1)
    ListTwo <- which(Cluster == 2)
    
    OSC1 <- as.numeric(OS[ListOne])
    OSC2 <- as.numeric(OS[ListTwo])
    CodeOSC1 <- as.numeric(CodeOS[ListOne])
    CodeOSC2 <- as.numeric(CodeOS[ListTwo])
    TTPC1 <- as.numeric(TTP[ListOne])
    TTPC2 <- as.numeric(TTP[ListTwo])
    CodeTTPC1 <- as.numeric(CodeTTP[ListOne])
    CodeTTPC2 <- as.numeric(CodeTTP[ListTwo])
    StageC1 <- Stage[ListOne]
    StageC2 <- Stage[ListTwo]
    
    ANNARBORStageC1 <- ANNARBORstage[ListOne]
    ANNARBORStageC2 <- ANNARBORstage[ListTwo]
    FLIPIC1  <- FLIPI[ListOne]
    FLIPIC2  <- FLIPI[ListTwo]
    TablePatientsC1 <- as.data.frame(cbind(as.numeric(CodeOS[ListOne]), 
                                           as.numeric(OS[ListOne]), 
                                           as.numeric(CodeTTP[ListOne]), 
                                           as.numeric(TTP[ListOne]), 
                                           Stage[ListOne], 
                                           as.numeric(ANNARBORstage[ListOne]), 
                                           as.numeric(FLIPI[ListOne])))
    rownames(TablePatientsC1) <- ClinicalFile_OrderedbyBetaMatrixPatients$SAMPLE_ID[ListOne]
    colnames(TablePatientsC1) <- c("CodeOSC1", "OSC1", "CodeTTPC1", 
                                   "TTPC1", "StageC1", 
                                   "ANNARBORStageC1", "FLIPIC1")
    
    TablePatientsC2 <- as.data.frame(cbind(as.numeric(CodeOS[ListTwo]), 
                                           as.numeric(OS[ListTwo]), 
                                           as.numeric(CodeTTP[ListTwo]), 
                                           as.numeric(TTP[ListTwo]), 
                                           Stage[ListTwo], 
                                           as.numeric(ANNARBORstage[ListTwo]), 
                                           as.numeric(FLIPI[ListTwo])))
    rownames(TablePatientsC2) <- ClinicalFile_OrderedbyBetaMatrixPatients$SAMPLE_ID[ListTwo]
    
    
    survTTPC1 <- survival::survfit(Surv(time = as.numeric(TTPC1), 
                                        event = as.numeric(CodeTTPC1)) ~ StageC1, 
                                   data = TablePatientsC1)
    survTTPC2 <- survival::survfit(Surv(time = as.numeric(TTPC2), 
                                        event = as.numeric(CodeTTPC2)) ~ StageC2, 
                                   data = TablePatientsC2)
    
    
    if (FigureGenerate == "Yes") {
      cat("\n Printing C1 image")
      if (PNGorPDF == "png") {
        grDevices::png(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTP_C1", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTP_C1.", PNGorPDF))
      }
      
      kmTTPC1 <- survminer::ggsurvplot(survTTPC2, 
                                             title = "101 Samples - TTP by SNF Cluster 2",
                                             font.main = "bold", 
                                             xlab = "TTP (years)", 
                                             ylab = "Proportion of patients surviving",
                                             pval = TRUE, 
                                             pval.coord = c(0, 0.05), 
                                             conf.int = F,
                                             break.time.by = 1, 
                                             ncensor.plot = F,
                                             risk.table = TRUE, # Add risk table
                                             risk.table.pos = "out", 
                                             censor = T, 
                                             risk.table.y.text.col = TRUE,
                                             risk.table.col = "strata", # Change risk table color by groups
                                             #linetype = "strata", # Change line type by groups
                                             #surv.median.line = "hv", # Specify median survival
                                             #ggtheme = theme_bw(), # Change ggplot2 theme
                                             palette = colourPalette,
                                             risk.table.height = 0.25,
                                             risk.table.y.text = FALSE,
                                             legend = "top",
                                             legend.title = "Stage")
      #,legend.labs = c("C=1","C=2"))
      grDevices::dev.off() 
    }
    
    
    if (FigureGenerate == "Yes") {
      cat("\n Printing C2 image")
      if (PNGorPDF == "png") {
        grDevices::png(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTP_C2", PNGorPDF))
      }
      if (PNGorPDF == "pdf") {
        grDevices::pdf(paste0(pathNow, "/img/15_SurvivalAnalysis_survival_TTP_C2.", PNGorPDF))
      }
      
      kmTTPC1 <- survminer::ggsurvplot(survTTPC2, 
                                       title = "101 Samples - TTP by SNF Cluster 2",
                                       font.main = "bold", 
                                       xlab = "TTP (years)", 
                                       ylab = "Proportion of patients surviving",
                                       pval = TRUE, 
                                       pval.coord = c(0, 0.05), 
                                       conf.int = F,
                                       break.time.by = 1, 
                                       ncensor.plot = F,
                                       risk.table = TRUE, # Add risk table
                                       risk.table.pos = "out", 
                                       censor = T, 
                                       risk.table.y.text.col = TRUE,
                                       risk.table.col = "strata", # Change risk table color by groups
                                       #linetype = "strata", # Change line type by groups
                                       #surv.median.line = "hv", # Specify median survival
                                       #ggtheme = theme_bw(), # Change ggplot2 theme
                                       palette = colourPalette,
                                       risk.table.height = 0.25,
                                       risk.table.y.text = FALSE,
                                       legend = "top",
                                       legend.title = "Stage")
      #,legend.labs = c("C=1","C=2"))
      grDevices::dev.off() 
    }
    
    
    
    RESULTS <- list(TTPSurvivalAdvanced = survTTPadvanced,
                    TTPSurvivalLimited = survTTPlimited,
                    TTPAll = survTTPall
                    #TTPCoxphAdvanced = coxph_TTP_advanced,
                    #TTPCoxphLimited = coxph_TTP_limited
                    )
    
    class(RESULTS) <- "SurvivalAnalysis_ASilva"
    return(RESULTS)
    
 }
 # [END]
