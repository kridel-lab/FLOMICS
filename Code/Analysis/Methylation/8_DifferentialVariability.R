# Updated 3 Aug 2021
# Updated 26 March 2019
# Function: Differential variability with no summarization. Rather than testing for differences 
#           in mean methylation, test for differences between group variances.
# http://bioconductor.org/packages/devel/bioc/vignettes/missMethyl/inst/doc/missMethyl.html#testing-for-differential-methylation-using
# Author: Anjali Silva

# Input:
# ClinicalFile: File with patient sample names, and categories. 
# MvalueMatrix: A matrix of probes x patients, with probes in rows and patients as columns.
# ProbeOrGene: Indicate whether the MvalueMatrix has probes or genes. ProbeOrGene options include "Gene" or "Probe".
# ContrastColumnName: Type of differential analysis to be done. ContrastColumnName options 
#                     include "STAGE", "SEX", "SITE_BIOPSY", "TYPE_BIOPSY", "COO", "TRANSLOC_14_18".
# PNGorPDF: Output format of the image, options = "png" or "pdf".


# Output:
# SummaryOfResults: Results output from Multiple Testing Across Genes and Contrasts using
#                   limma package.
# MArrayLM: Empirical Bayes Statistics from applying eBayes function.
# topP: Top 100 differentially methylated features (probes or genes as defined) based on P value.
# topLogFC: Top 100 differentially methylated features (probes or genes as defined) based on LogFC.
# topGOList: GO and KEGG results for top 50 differentially methylated regions.

# Visuals saved to img folder
# 6_VolcanoPlot_",ContrastColumnName,ProbeOrGene,"_toplogFC.p*

DifferentialVariability8 <- function(ClinicalFile, 
                                     MvalueMatrix, 
                                     ProbeOrGene, 
                                     ContrastColumnName, 
                                     PNGorPDF) {
  
  # Loading needed packages
  library(limma)
  library(missMethyl)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(ggplot2)
  library(ggrepel)
  library(stringi)
  
  if(nrow(ClinicalFile) != ncol(MvalueMatrix)) {
    cat("\n ***ClinicalFile and MvalueMatrix does not have the same number of patients.\n")
  }
  
  # Remove empty entries in the clinical file
  ClinicalFile[ClinicalFile == ""] <- NA
  
  # # # # # # # # # # # # # # # # # # # # # # # #
  # Differential methylation between ADVANCED-LIMITED
  # Select advanced and limited stage cases only from the FL cases
  if(ContrastColumnName == "STAGE") {
    stage_status <- factor(stringi::stri_trim(ClinicalFile[which(ClinicalFile[ , 
                                              which(colnames(ClinicalFile) == "TYPE")] == "FL"),
                                              which(colnames(ClinicalFile) == "STAGE")]), 
                                              levels = c("ADVANCED", "LIMITED"))
    design_epic2 <- stats::model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fitvar <- missMethyl::varFit(data = MvalueMatrix[, which(ClinicalFile[ , 
                                        which(colnames(ClinicalFile) == "TYPE")] == "FL")], 
                                        design = design_epic2, coef = c(1, 2))
    summary <- summary(decideTests(fitvar))
    topDV <- missMethyl::topVar(fitvar, coef = 2)
  } else if(ContrastColumnName == "SEX") {
    stage_status <- factor(stringi::stri_trim(ClinicalFile[which(ClinicalFile[ , 
                                              which(colnames(ClinicalFile) == "TYPE")] == "FL"),
                                              which(colnames(ClinicalFile) == "SEX")]), 
                                              levels = c("M", "F"), 
                                              labels = c("Male", "Female"))
    design_epic2 <- stats::model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fitvar <- missMethyl::varFit(data = MvalueMatrix[ , which(ClinicalFile[ , 
                                        which(colnames(ClinicalFile) == "TYPE")] == "FL")], 
                                        design = design_epic2, coef = c(1, 2))
    summary <- summary(decideTests(fitvar))
    topDV <- missMethyl::topVar(fitvar, coef = 2)
  } else if(ContrastColumnName == "SITE_BIOPSY") {
    stage_status <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[ , 
                                              which(colnames(ClinicalFile) == "SITE_BIOPSY")] == TRUE))),
                                              which(colnames(ClinicalFile) == "SITE_BIOPSY")]), 
                                              levels = c("EN", "LN"))
    design_epic2 <- stats::model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fitvar <- missMethyl::varFit(data = MvalueMatrix[ , - c(which(is.na(ClinicalFile[ , 
                                        which(colnames(ClinicalFile) == "SITE_BIOPSY")] == TRUE)))],
                                        design = design_epic2, coef = c(1, 2))
    summary <- summary(decideTests(fitvar))
    topDV <- missMethyl::topVar(fitvar, coef = 2)
  } else if(ContrastColumnName == "TYPE_BIOPSY"){
    stage_status <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[ , 
                                              which(colnames(ClinicalFile) == "TYPE_BIOPSY")]) == "TRUE")), 
                                              which(colnames(ClinicalFile) == "TYPE_BIOPSY")]), 
                                              levels = c("CORE", "TISSUE"))
    design_epic2 <- stats::model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fitvar <- missMethyl::varFit(data = MvalueMatrix[ , - c(which(is.na(ClinicalFile[ , 
                                        which(colnames(ClinicalFile) == "TYPE_BIOPSY")]) == "TRUE"))], 
                                        design = design_epic2, coef = c(1, 2))
    summary <- summary(decideTests(fitvar))
    topDV <- missMethyl::topVar(fitvar, coef = 2)
  } else if(ContrastColumnName == "COO") {
    # https://support.bioconductor.org/p/44216/ # contrasts for multiple comparisons
    stage_status <- factor(stringi::stri_trim(ClinicalFile[which(ClinicalFile[ , 
                                              which(colnames(ClinicalFile) == "TYPE")] == "DLBCL"),
                                              which(colnames(ClinicalFile) == "COO")]), 
                                              levels = c("GCB", "ABC"))
    design_epic2 <- stats::model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fitvar <- missMethyl::varFit(data = MvalueMatrix[ , which(ClinicalFile[, 
                                        which(colnames(ClinicalFile) == "TYPE")] == "DLBCL")], 
                                        design = design_epic2, coef = c(1, 2))
    summary <- summary(decideTests(fitvar))
    topDV <- missMethyl::topVar(fitvar, coef = 2)
  } else if(ContrastColumnName == "TRANSLOC_14_18") {
    stage_status <- factor(stringi::stri_trim(ClinicalFile[which(! is.na(ClinicalFile[ , 
                                              which(colnames(ClinicalFile) == "TRANSLOC_14_18")]) == TRUE),
                                              which(colnames(ClinicalFile) == "TRANSLOC_14_18")]), 
                                              levels = c("0", "1"), 
                                              labels = c("NoTranslocation", "Translocation"))
    design_epic2 <- stats::model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fitvar <- missMethyl::varFit(data = MvalueMatrix[ , which(! is.na(ClinicalFile[ , 
                                 which(colnames(ClinicalFile) == "TRANSLOC_14_18")]) == TRUE)], 
                                 design = design_epic2, coef = c(1, 2))
    summary <- summary(decideTests(fitvar))
    topDV <- missMethyl::topVar(fitvar, coef = 2)
  #}else if(ContrastColumnName=="INSTITUTION"){
    #  # https://support.bioconductor.org/p/44216/
    # stage_status <- factor(stri_trim(ClinicalFile[,which(colnames(ClinicalFile)=="INSTITUTION")]), levels = c("BCCA", "UHN","JGH"))
    
    #stage_status_bcca_uhn <- factor(stage_status[-which(stage_status=="JGH")])
    #design_epic_bcca_uhn <- model.matrix (~0 + stage_status_bcca_uhn)
    #colnames(design_epic_bcca_uhn) <- levels(stage_status_bcca_uhn)
    #fitvar_bcca_uhn <- varFit(data=MvalueMatrix[,-which(stage_status=="JGH")], design = design_epic_bcca_uhn, coef = c(1,2))
    #summary <- summary(decideTests(fitvar_bcca_uhn))
    #topDV <- topVar(fitvar_bcca_uhn, coef=2)
    #}else if(ContrastColumnName=="EPIC_QC"){
    #stage_status <- factor(stri_trim(ClinicalFile[,which(colnames(ClinicalFile)=="EPIC_QC")]), levels = c("Fair","Good","Poor"))
    #design_epic2 <- model.matrix (~0 + stage_status)
    #colnames(design_epic2) <- levels(stage_status)
    #fit.reduced <- lmFit(MvalueMatrix,design_epic2)
    #contrasts <- makeContrasts("Good-Fair", "Poor-Good", "Poor-Fair", levels = design_epic2)
    #}else if(ContrastColumnName=="TYPE"){
    #stage_status <- factor(stri_trim(ClinicalFile[,which(colnames(ClinicalFile)=="TYPE")]), levels = c("DLBCL","FL","RLN"))
    #design_epic2 <- model.matrix (~0 + stage_status)
    #colnames(design_epic2) <- levels(stage_status)
    #fit.reduced <- lmFit(MvalueMatrix,design_epic2)
    #contrasts <- makeContrasts("FL-DLBCL", "RLN-DLBCL", "RLN-FL", levels = design_epic2)
  }else{
    cat("\n ***Warning: Not a valid ContrastColumnName.\n")
  }
  

  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()

  # gene ontology analysis
  if(ProbeOrGene == "Probe") {
    # gene ontology analysis 
    sigCpGs <- rownames(topDV) # names of top 50 CpG sites
    gstEpic <- gometh(sig.cpg = sigCpGs, 
                       all.cpg = rownames(MvalueMatrix), 
                       collection = c("GO","KEGG"), 
                       array.type = "EPIC")
    topGOList <- topGSA(gstEpic)
  
  } else if(ProbeOrGene == "Gene") {
    gstEpic <- "For GO output, input of probes (not genes) needed"
    topGOList <- "For GO output, input of probes (not genes) needed"
  }


  RESULTS <- list(AllVariabilityResults = fitvar,
                  SummaryVariabilityResults = summary,
                  topVariabilityResults = topDV,
                  GOResults = gstEpic,
                  topGOResults = topGOList)
  
  class(RESULTS) <- "DifferentialVariability_ASilva"
  return(RESULTS)
}
# [END]
