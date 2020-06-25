# Updated 12 Feb 2019
# https://bioconductor.org/packages/devel/bioc/vignettes/missMethyl/inst/doc/missMethyl.html#testing-for-differential-methylation-using
# Function: Differential methylation based on probe level (with no summarization)
# Author: Anjali Silva

# Input:
# ClinicalFile: File with patient sample names, and categories. 
# MvalueMatrix: A matrix of probes x patients, with probes in rows and patients as columns.
# BetaMatrix: 
# ProbeOrGene: Indicate whether the MvalueMatrix has probes or genes. ProbeOrGene options include "Gene" or "Probe".
# ContrastColumnName: Type of differential analysis to be done. ContrastColumnName options include "STAGE", "SEX", "SITE_BIOPSY", "TYPE_BIOPSY", "INSTITUTION", "COO", "TYPE", "TRANSLOC_14_18", "EPIC_QC".
# ProduceImages: Produce images or not, options = "Yes" or "No"
# PNGorPDF: Output format of the image, options = "png" or "pdf"

# Output:
# SummaryOfResults: results output from Multiple Testing Across Genes and Contrasts using limma package
# MArrayLM: Empirical Bayes Statistics from applying limma::eBayes function
# topP: Table of all features from Linear Model Fit (probes or genes as defined) based on P value using limma::eBayes() and limma::topTable() functions
# topPGOTesting: All results from gene ontology enrichment for the vector of CpG sites from topP using gometh() function
# topPGOList: 20 most enriched pathways for the vector of CpG sites from topP using topGO() function
# topLogFC: Table of all features from Linear Model Fit (probes or genes as defined) based on LogFC
# topLogFCTesting: All results from gene ontology enrichment for the vector of CpG sites from topLogFC using gometh() function
# topLogFCGOList: 20 most enriched pathways for the vector of CpG sites from topLogFC using topGO() function

# Visuals saved to img folder
# 7_VolcanoPlot_", ContrastColumnName, ProbeOrGene, "_toplogFC.p*



DifferentialMethylation <- function(ClinicalFile, 
                                    MvalueMatrix, 
                                    BetaMatrix, 
                                    ProbeOrGene, 
                                    ContrastColumnName, 
                                    RGChannelSet, 
                                    SampleSheet, 
                                    ProduceImages = "Yes", 
                                    PNGorPDF = "png",
                                    AnnotationFile = NA) { 
  
  # Loading needed packages
  # RegularPckgs=c("limma","missMethyl","IlluminaHumanMethylationEPICanno.ilm10b2.hg19","ggplot2","ggrepel","stringi")
  # "stringi" package used to remove extra white space, if present, e.g. "EN " to "EN"
  library(missMethyl)
  library(limma)
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
    stage_status <- factor(stri_trim( ClinicalFile[which(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE")] == "FL"), 
                                                   which(colnames(ClinicalFile) == "STAGE")]), levels = c("ADVANCED", "LIMITED"))
    design_epic2 <- model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fit.reduced <- limma::lmFit(MvalueMatrix[,which(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE")] == "FL")], design_epic2)
    contrasts <- limma::makeContrasts("ADVANCED-LIMITED", levels = design_epic2)
  } else if(ContrastColumnName == "SEX") {
    stage_status <- factor(stri_trim(ClinicalFile[which(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE")] == "FL"), 
                                                  which(colnames(ClinicalFile) == "SEX")]), levels = c("M", "F"), labels = c("Male", "Female"))
    design_epic2 <- model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fit.reduced <- limma::lmFit(MvalueMatrix[,which(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE")] == "FL")], design_epic2)
    contrasts <- limma::makeContrasts("Female-Male", levels = design_epic2)
  } else if(ContrastColumnName == "SITE_BIOPSY") { 
    stage_status <- factor(stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[ , which(colnames(ClinicalFile) == "SITE_BIOPSY")] == TRUE))), 
                                                  which(colnames(ClinicalFile) == "SITE_BIOPSY")]), levels = c("EN", "LN"))
    design_epic2 <- model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fit.reduced <- limma::lmFit(MvalueMatrix[ , - c(which(is.na(ClinicalFile[ , which(colnames(ClinicalFile) == "SITE_BIOPSY")] == TRUE)))], design_epic2)
    contrasts <- limma::makeContrasts("LN-EN", levels = design_epic2)
  } else if(ContrastColumnName == "TYPE_BIOPSY") {
    stage_status <- factor(stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE_BIOPSY")]) == "TRUE")),
                                                  which(colnames(ClinicalFile) == "TYPE_BIOPSY")]), levels = c("CORE", "TISSUE"))
    design_epic2 <- model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fit.reduced <- limma::lmFit(MvalueMatrix[ , - c(which(is.na(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE_BIOPSY")]) == "TRUE"))], design_epic2)
    contrasts <- limma::makeContrasts("TISSUE-CORE", levels = design_epic2)
  } else if(ContrastColumnName == "INSTITUTION") { 
    # https://support.bioconductor.org/p/44216/
    stage_status <- factor(stri_trim(ClinicalFile[ , which(colnames(ClinicalFile) == "INSTITUTION")]), levels = c("BCCA", "UHN","JGH"))
    design_epic2 <- model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fit.reduced <- limma::lmFit(MvalueMatrix,design_epic2)
    contrasts <- limma::makeContrasts("UHN-BCCA", "JGH-BCCA", "JGH-UHN", levels = design_epic2)
  } else if(ContrastColumnName == "EPIC_QC") { 
    stage_status <- factor(stri_trim(ClinicalFile[,which(colnames(ClinicalFile) == "EPIC_QC")]), levels = c("Fair", "Good", "Poor"))
    design_epic2 <- model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fit.reduced <- limma::lmFit(MvalueMatrix,design_epic2)
    contrasts <- limma::makeContrasts("Good-Fair", "Poor-Good", "Poor-Fair", levels = design_epic2)
  } else if(ContrastColumnName == "COO") { 
    # https://support.bioconductor.org/p/44216/ # contrasts for multiple comparisons
    stage_status <- factor(stri_trim(ClinicalFile[which(ClinicalFile[,which(colnames(ClinicalFile) == "TYPE")] == "DLBCL"), 
                                                  which(colnames(ClinicalFile) == "COO")]), levels = c("GCB", "ABC"))
    design_epic2 <- model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fit.reduced <- limma::lmFit(MvalueMatrix[,which(ClinicalFile[,which(colnames(ClinicalFile) == "TYPE")] == "DLBCL")], design_epic2)
    contrasts <- limma::makeContrasts("ABC-GCB", levels = design_epic2)
  } else if(ContrastColumnName == "TYPE") {
    stage_status <- factor(stri_trim(ClinicalFile[,which(colnames(ClinicalFile) == "TYPE")]), levels = c("DLBCL", "FL", "RLN"))
    design_epic2 <- model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fit.reduced <- limma::lmFit(MvalueMatrix,design_epic2)
    # contrasts <- limma::makeContrasts("FL-DLBCL", "RLN-DLBCL", "RLN-FL", levels = design_epic2)
    contrasts <- limma::makeContrasts("DLBCL-FL", "DLBCL-RLN", "FL-RLN", levels = design_epic2)
  } else if(ContrastColumnName == "TRANSLOC_14_18") { 
    stage_status <- factor(stri_trim(ClinicalFile[which(!is.na(ClinicalFile[ , which(colnames(ClinicalFile) == "TRANSLOC_14_18")]) == TRUE),
                                                  which(colnames(ClinicalFile) == "TRANSLOC_14_18")]), 
                           levels = c("0", "1"), labels = c("NoTranslocation", "Translocation"))
    design_epic2 <- model.matrix (~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    fit.reduced <- limma::lmFit(MvalueMatrix[, which(! is.na(ClinicalFile[, which(colnames(ClinicalFile) == "TRANSLOC_14_18")]) == TRUE)], design_epic2)
    contrasts <- limma::makeContrasts("Translocation-NoTranslocation", levels = design_epic2)
  } else if(ContrastColumnName == "CLUSTER") {
    if (length(ClinicalFile[, which(colnames(ClinicalFile) == "CLUSTER")]) != ncol(BetaMatrix)) {
      stop("ContrastColumnName is set to CLUSTER, but ClusterLabels !=  number of observations");}
    ClusterLabels <- ClinicalFile[, which(colnames(ClinicalFile) == "CLUSTER")]
    if (length(unique(ClusterLabels)) == 5) {
      stage_status <- factor(ClusterLabels, levels = c("1", "2", "3", "4", "5"), labels = c("one", "two", "three", "four", "five"))
    } else if (length(unique(ClusterLabels)) == 4) {
      stage_status <- factor(ClusterLabels, levels = c("1", "2", "3", "4"), labels = c("one", "two", "three", "four"))
    } else if(length(unique(ClusterLabels)) == 3) {
      stage_status <- factor(ClusterLabels, levels = c("1", "2", "3"), labels = c("one", "two", "three"))
    } else if(length(unique(ClusterLabels)) == 2) {
      stage_status <- factor(ClusterLabels, levels = c("1", "2"), labels = c("one", "two"))
    } else{
      stop("ClusterLabels have more than 5 clusters. Currently only upto 4 clusters supported");}
    
    design_epic2 <- stats::model.matrix(~ 0 + stage_status)
    colnames(design_epic2) <- levels(stage_status)
    if(length(which(is.na(ClinicalFile[ , which(colnames(ClinicalFile) == "CLUSTER")] == TRUE))) > 0) {
      fit.reduced <- limma::lmFit(MvalueMatrix[ , - c(which(is.na(ClinicalFile[ , which(colnames(ClinicalFile) == "CLUSTER")] == TRUE)))], design_epic2)
    } else {
      fit.reduced <- limma::lmFit(MvalueMatrix, design_epic2)
    }
  
    if (length(unique(ClusterLabels)) == 5) {
      contrasts <- limma::makeContrasts("one-two","one-three","one-four", "one-five", "two-three","two-four", "two-five", "three-four", "three-five", "four-five", levels = design_epic2)
    } else if (length(unique(ClusterLabels)) == 4) {
      contrasts <- limma::makeContrasts("one-two","one-three","one-four","two-three","two-four","three-four", levels = design_epic2)
    } else if (length(unique(ClusterLabels)) == 3) {
      contrasts <- limma::makeContrasts("one-two","one-three","two-three", levels = design_epic2)
    } else if(length(unique(ClusterLabels)) == 2) {
      contrasts <- limma::makeContrasts("one-two", levels = design_epic2)
    }
  } else { 
    cat("\n ***Warning: Not a valid ContrastColumnName.\n")
  }
  
  fit.reduced <- limma::eBayes(contrasts.fit(fit.reduced, contrasts))
  fit.reduced_MultipleTesting <- limma::decideTests(fit.reduced, adjust.method = "BH", p.value = 0.01)
  Significant_Probes <- rownames(fit.reduced_MultipleTesting)[which(fit.reduced_MultipleTesting != 0)]
  Summary_Results <- summary(limma::decideTests(fit.reduced, adjust.method = "BH", p.value = 0.01))
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  position_20genes <- c(1:10)
  
  # All methylated CpGs for the given contrast extracted using limma::topTable
  # sorted by P value
  if(length(levels(stage_status)) >= 3) { 
    # need to adjust for multiple testing  
    # order probes by P value and then select probes with adj.P.Val < 0.01
    top_P_Coef1 <- limma::topTable(fit.reduced, adjust.method="BH", coef = 1, number = nrow(MvalueMatrix), sort.by="p")
    top_P_Coef11 <- top_P_Coef1[which(top_P_Coef1[, 5] < 0.01), ]
    top_P_Coef2 <- limma::topTable(fit.reduced, adjust.method = "BH", coef = 2, number = nrow(MvalueMatrix), sort.by="p")
    top_P_Coef22 <- top_P_Coef2[which(top_P_Coef2[, 5] < 0.01), ]
    top_P_Coef3 <- limma::topTable(fit.reduced, adjust.method = "BH", coef = 3, number = nrow(MvalueMatrix), sort.by="p")
    top_P_Coef33 <- top_P_Coef3[which(top_P_Coef3[, 5] < 0.01), ]
    
    top_P <- list(FirstComparison = top_P_Coef11,
                  SecondComparison = top_P_Coef22,
                  ThirdComparison = top_P_Coef33)
    
    # Among the 3 contrasts, find common probes (those that have been selected based on adj.P.Val < 0.01)
    ProbesAmongAll <- Reduce(intersect, list(rownames(top_P_Coef11),rownames(top_P_Coef22),rownames(top_P_Coef33)))
    
    if(length(ProbesAmongAll) == 0) {
      cat("\n No common", paste0(ProbeOrGene, "s"),"were found among all comparisons.")
    } else { 
      # length(ProbesAmongAll) # 1537
      
      # Iding the locations where the probes are located in each comparison 
      Locations <- lapply(c(1:length(ProbesAmongAll)), function(i) c(which(rownames(top_P_Coef11) == ProbesAmongAll[i]), 
                                                                     which(rownames(top_P_Coef22) == ProbesAmongAll[i]), 
                                                                     which(rownames(top_P_Coef33) == ProbesAmongAll[i])))
      
      # Save the locations into a matrix with ncols = number of comparisons and nrows = number of rows
      Locations2 <- matrix(do.call(c, Locations), byrow = T, ncol = length(levels(stage_status)))
      # dim(Locations2) # 1537    3
      
      # Construct a Volcano plot for FL-DLBCL
      # Volcano plot; Note 1.3 is -log10(0.05)
      # Note 2 is -log10(0.01) 
      
      # find needed genes in top_P
      # position_20genes <- match(rownames(beta.sel), rownames(top_P_Coef1))
      # top_P_Coef1[position_20genes,]
      
      # checking hypomethylated genes
      # Locations2 <- top_P_Coef33[which((top_P_Coef33$logFC < 0) == TRUE), ]
      
      # ggplot2::ggplot(top_P_Coef1, aes(logFC, -log10(adj.P.Val) )) +
      #   geom_point(color="#e0e0e0")+
      #   ggtitle(paste0("Plot of ", colnames(Summary_Results)[1])) +
      #   geom_text_repel(data= top_P_Coef11[Locations2,], 
      #                   aes(x=logFC, y=-log10(adj.P.Val), label=rownames(top_P_Coef11)[Locations2]),
      #                   fontface=3 ) + geom_point(data=top_P_Coef11[Locations2,], aes(x=logFC,
      #                                                                                 y=-log10(adj.P.Val)), colour="#b2182b", size=3) + 
      #    geom_point(data=top_P_Coef11[Locations2,], aes(x=logFC, y=-log10(adj.P.Val)),
      #              colour="#4d4d4d", size=3) + labs(color = "GeneClass", y="-log10(adj.P.Value)")+
      #   geom_hline(yintercept= 1.3, linetype="dotted") +
      #   geom_hline(yintercept= 2, linetype="dotted") +  
      #   theme_bw() + 
      #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      #   scale_color_manual(values=c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
      #   theme(aspect.ratio=1, text = element_text(size=15))
      
      # ggplot2::ggplot(top_P_Coef3, aes(logFC, -log10(adj.P.Val) )) +
      #   geom_point(color="#e0e0e0") +
      #   ggtitle(paste0("Plot of ", colnames(Summary_Results)[3])) +
      #   geom_hline(yintercept= 2, linetype="dotted") + 
      #   labs(y="-log10(adj. P-value)") +
      #   theme_bw() + 
      #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      #   scale_color_manual(values=c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
      #   theme(aspect.ratio=1, text = element_text(size=15))
      
      if(ProduceImages == "Yes") { 
        volcanoplot <- ggplot2::ggplot(top_P_Coef1, aes(logFC, -log10(adj.P.Val))) +
                                       geom_point(color = "#e0e0e0") +
                                       ggtitle(paste0("Plot of ", colnames(Summary_Results)[1])) +
                                       geom_text_repel(data = top_P_Coef11[Locations2[c(1:10), 1], ], 
                                                       aes(x = logFC, 
                                                           y = -log10(adj.P.Val),
                                                           label = rownames(top_P_Coef11)[Locations2[c(1:10), 1]]),
                                                      fontface = 3 ) + 
                                       geom_point(data = top_P_Coef11[Locations2[c(1:5), 1], ], 
                                                  aes(x = logFC, y = -log10(adj.P.Val)), 
                                                  colour = "#b2182b", size = 3) + 
                                       geom_point(data = top_P_Coef11[Locations2[c(6:10), 1], ], 
                                                  aes(x = logFC, y = -log10(adj.P.Val)),
                                                  colour = "#4d4d4d", size = 3) + 
                                       labs(color = "GeneClass", y = "-log10(adj.P.Value)")+
                                       geom_hline(yintercept = 1.3, linetype = "dotted") +
                                       geom_hline(yintercept = 2, linetype = "dotted") +  
                                       theme_bw() +  
                                       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                                       scale_color_manual(values = c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
                                       theme(aspect.ratio = 1, text = element_text(size=15))
        ggsave(paste0(pathNow, "/img/7_VolcanoPlot_", ContrastColumnName, "_1stComparison_topPValue.", PNGorPDF))
        
        # Construct a Volcano plot for RLN-DLBCL
        volcanoplot <- ggplot2::ggplot(top_P_Coef2, aes(logFC, -log10(adj.P.Val))) +
                                       geom_point(color = "#e0e0e0") +
                                       ggtitle(paste0("Plot of ", colnames(Summary_Results)[2])) +
                                       geom_text_repel(data = top_P_Coef22[Locations2[c(1:10), 2], ], 
                                                      aes(x = logFC, 
                                                          y = -log10(adj.P.Val), 
                                                          label = rownames(top_P_Coef22)[Locations2[c(1:10), 2]]),
                                                          fontface = 3 ) + 
                                       geom_point(data = top_P_Coef22[Locations2[c(1:5),2], ], 
                                                  aes(x = logFC, y = -log10(adj.P.Val)), 
                                                  colour="#b2182b", size=3) + 
                                       geom_point(data = top_P_Coef22[Locations2[c(6:10), 2], ], 
                                                  aes(x = logFC, y = -log10(adj.P.Val)),
                                                  colour = "#4d4d4d", size = 3) + 
                                       labs(color = "GeneClass", y="-log10(adj.P.Value)") +
                                       geom_hline(yintercept = 1.3, linetype = "dotted") +
                                       geom_hline(yintercept = 2, linetype = "dotted") +  
                                       theme_bw() + 
                                       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                                       scale_color_manual(values = c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
                                       theme(aspect.ratio = 1, text = element_text(size = 15))
        ggsave(paste0(pathNow, "/img/7_VolcanoPlot_", ContrastColumnName, "_2ndComparison_topPValue.", PNGorPDF))
        
        # Construct a Volcano plot for RLN-FL
        volcanoplot <- ggplot2::ggplot(top_P_Coef3, aes(logFC, -log10(adj.P.Val))) +
                                      geom_point(color = "#e0e0e0")+
                                      ggtitle(paste0("Plot of", colnames(Summary_Results)[3])) +
                                      geom_text_repel(data = top_P_Coef33[Locations2[c(1:10), 3], ], 
                                                      aes(x = logFC, 
                                                          y = -log10(adj.P.Val), 
                                                          label = rownames(top_P_Coef33)[Locations2[c(1:10), 3]]), 
                                                          fontface = 3 ) + 
                                      geom_point(data = top_P_Coef33[Locations2[c(1:5), 3], ], 
                                                aes(x = logFC, y = -log10(adj.P.Val)), 
                                                colour = "#b2182b", size = 3) + 
                                      geom_point(data = top_P_Coef33[Locations2[c(6:10),3], ], 
                                                 aes(x = logFC, y = -log10(adj.P.Val)),
                                                 colour = "#4d4d4d", size = 3) + 
                                                 labs(color = "GeneClass", 
                                                      y = "-log10(adj. P.Value)") +
                                      geom_hline(yintercept = 1.3, linetype = "dotted") +
                                      geom_hline(yintercept = 2, linetype = "dotted") +  
                                      theme_bw() + 
                                      theme(panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank()) + 
                                      scale_color_manual(values = c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
                                      theme(aspect.ratio = 1, text = element_text(size=15))
        ggsave(paste0(pathNow,"/img/7_VolcanoPlot_",ContrastColumnName,"_3rdComparison_topPValue.",PNGorPDF))
      }
    }
  } else if(length(levels(stage_status)) == 2) { 
    top_P <- limma::topTable(fit.reduced, coef = 1, number = nrow(MvalueMatrix), sort.by="p")
    top_P_Coef11 <- top_P[which(top_P$adj.P.Val < 0.01), ]
    
    # Determine hyper and hypomethylation
    # length(which((top_P_Coef22$logFC > 0) == TRUE))
    # length(which((top_P_Coef22$logFC < 0) == TRUE))
    if (all(is.na(AnnotationFile) != TRUE)) {
      testing_up <- match(rownames(top_P_Coef11[which((top_P_Coef11$logFC > 0) == TRUE), ]), 
                          AnnotationFile$V1)
      (table(AnnotationFile$Relation_to_Island[testing_up]) / sum(table(AnnotationFile$Relation_to_Island[testing_up]))) * 100
      geneRegion <- sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Group)
      geneRegion[which(geneRegion == "")] <- "NA"
      table(geneRegion[testing_up])
      
      
      geneRegion_breakdown_up <- data.frame(table(AnnotationFile$Relation_to_Island[testing_up]))
      geneRegion_breakdown_up <- data.frame(table(geneRegion[testing_up]))
      
      # Using plot_ly to create pie-chart for gene region
      geneRegionType <- plotly::plot_ly(geneRegion_breakdown_up, 
                                        labels = ~Var1, values = ~Freq, 
                                        marker = list(colors = c("blue","orange","green","red", "purple", "brown")),
                                        type = 'pie',
                                        textposition = 'outside',textinfo = 'label+percent') %>%
                                        layout(title = '', showlegend = FALSE, 
                                               xaxis = list(showgrid = FALSE, 
                                                            zeroline = FALSE, 
                                                            showticklabels = FALSE),
                                               yaxis = list(showgrid = FALSE, 
                                                            zeroline = FALSE, 
                                                            showticklabels = FALSE))
                                      
      # Create heatmap
      # Match significant probes with Beta Matrix
      matchBetaSigProb_up <- match(rownames(top_P_Coef11[which((top_P_Coef11$logFC > 0) == TRUE), ]), 
                                   rownames(BetaMatrix_T1))
      write.csv(x = BetaMatrix_T1[matchBetaSigProb_up, ], 
                file = "matchBetaSigProb_up.csv")
      heatmap(x = BetaMatrix_T1[matchBetaSigProb_up, ])
      
      
      
      testing_down <- match(rownames(top_P_Coef11[which((top_P_Coef11$logFC < 0) == TRUE), ]), AnnotationFile$V1)
      (table(AnnotationFile$Relation_to_Island[testing_down]) / sum(table(AnnotationFile$Relation_to_Island[testing_down]))) * 100
      #table(AnnotationFile$chr[testing_down])
      #table(AnnotationFile$Phantom4_Enhancers[testing_down])
      #table(AnnotationFile$Phantom5_Enhancers[testing_down])
      #table(AnnotationFile$X450k_Enhancer[testing_down])
      #table(AnnotationFile$DNase_Hypersensitivity_NAME[testing_down])
      #table(AnnotationFile$DNase_Hypersensitivity_Evidence_Count[testing_down])
      #table(AnnotationFile$UCSC_RefGene_Group[testing_down])
      
      geneRegion <- sub("\\;.*", "", AnnotationFile$UCSC_RefGene_Group)
      geneRegion[which(geneRegion == "")] <- "NA"
      table(geneRegion[testing_down])
      
      geneRegion_breakdown_down <- data.frame(table(AnnotationFile$Relation_to_Island[testing_down]))
      geneRegion_breakdown_down <- data.frame(table(geneRegion[testing_down]))
      # Using plot_ly to create pie-chart for gene region
      geneRegionType <- plotly::plot_ly(geneRegion_breakdown_down, labels = ~Var1, values = ~Freq, 
                                        marker = list(colors = c("blue","orange","green","red", "purple", "brown")),
                                        type = 'pie',
                                        textposition = 'outside',textinfo = 'label+percent') %>%
        layout(title = '', showlegend = FALSE, 
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
      
      
      # }
      
      
      # rank by p value
      top_P <- limma::topTable(fit.reduced, coef=1, number = nrow(MvalueMatrix), sort.by="p")
      # length(which(top_P$adj.P.Val < 0.01)) # of significant probes
      
      # find needed genes in top_P
      # position_20genes <- match(rownames(beta.sel), rownames(top_P))
      # length(which(is.na(position_20genes) == TRUE)) # 201 -length(which(is.na(position_20genes) == TRUE))
      # top_P[position_20genes,]
      
      
      # Volcano plot; Note 1.3 is -log10(0.05)
      volcanoplot <- ggplot2::ggplot(top_P, aes(logFC, - log10(adj.P.Val) )) +
                                     geom_point(color = "#e0e0e0") +
                                     ggtitle(paste0("Plot of ", colnames(Summary_Results))) +
                                     geom_text_repel(data = top_P[position_20genes, ], 
                                                    aes(x = logFC, 
                                                        y = -log10(adj.P.Val), 
                                                        label = rownames(top_P)[position_20genes]),
                                                        fontface = 3 ) + 
                                     geom_point(data = top_P[position_20genes[1:5], ], 
                                                aes(x = logFC, 
                                                    y = -log10(adj.P.Val)), 
                                                colour = "#b2182b", 
                                                size = 2) + 
                                     geom_point(data = top_P[position_20genes[6:10], ], 
                                                aes(x = logFC, y = -log10(adj.P.Val)),
                                                colour = "#4d4d4d", size = 2) + 
                                     labs(color = "GeneClass", y = "-log10(adj.P.Value)")+
                                     geom_hline(yintercept =  1.3, linetype = "dotted") +  
                                     theme_bw() + 
                                     theme(panel.grid.major = element_blank(), 
                                           panel.grid.minor = element_blank()) + 
                                     scale_color_manual(values=c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
                                     theme(aspect.ratio = 1, text = element_text(size = 15))
      
      # plotting based on 201 probes
      #   ggplot2::ggplot(top_P, aes(logFC, -log10(adj.P.Val) )) +
      #   geom_point(color="#e0e0e0")+
      #   ggtitle(paste0("Plot of ", colnames(Summary_Results), " S_Shore")) +
      #   geom_point(data=top_P[position_20genes,], aes(x=logFC, y=-log10(adj.P.Val)),
      #                colour="#4d4d4d", size=2) + labs(color = "GeneClass", y="-log10(adj.P.Value)")+
      #     geom_hline(yintercept= 1.3, linetype="dotted") +  
      #     geom_hline(yintercept= 2, linetype="dotted") +  
      #    theme_bw() + 
      #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      #     scale_color_manual(values=c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
      #     theme(aspect.ratio=1, text = element_text(size=15))
      if(ProduceImages == "Yes") {
        ggsave(paste0(pathNow, "/img/7_VolcanoPlot_", ContrastColumnName, "_topPValue.", PNGorPDF))
      }
    }
    
    # sorted by log fold change
    if(length(levels(stage_status)) >= 3) { 
      
      top_logFC_Coef1 <- limma::topTable(fit.reduced, 
                                         adjust.method  = "BH", 
                                         coef = 1, 
                                         number = nrow(MvalueMatrix), 
                                         sort.by = "logFC")
      top_logFC_Coef11 <- top_logFC_Coef1[which(top_logFC_Coef1[ , 5] < 0.01), ]
      top_logFC_Coef2 <- limma::topTable(fit.reduced, 
                                         adjust.method = "BH", 
                                         coef = 2, 
                                         number = nrow(MvalueMatrix), 
                                         sort.by = "logFC")
      top_logFC_Coef22 <- top_logFC_Coef2[which(top_logFC_Coef2[ , 5] < 0.01), ]
      top_logFC_Coef3 <- limma::topTable(fit.reduced, 
                                         adjust.method = "BH", 
                                         coef = 3, 
                                         number = nrow(MvalueMatrix), 
                                         sort.by = "logFC")
      top_logFC_Coef33 <- top_logFC_Coef3[which(top_logFC_Coef3[ , 5] < 0.01), ]
      
      top_logFC <- list(FirstComparison = top_logFC_Coef11,
                        SecondComparison = top_logFC_Coef22,
                        ThirdComparison = top_logFC_Coef33)
      
      # Among the 3 contrasts, find common probes (those that have been selected based on adj.P.Val < 0.01)
      ProbesAmongAll_logFC <- Reduce(intersect, 
                                     list(rownames(top_logFC_Coef11),
                                          rownames(top_logFC_Coef22),
                                          rownames(top_logFC_Coef33)))
      # length(ProbesAmongAll) # 1537
      
      if (length(ProbesAmongAll_logFC) == 0) {
        cat("\n No common", paste0(ProbeOrGene, "s"), "were found among all comparisons.")
      } else {
        
        # Iding the locations where the probes are located in each comparison 
        Locations_logFC <- lapply(c(1:length(ProbesAmongAll_logFC)), 
                                  function(i) c(which(rownames(top_logFC_Coef11) == ProbesAmongAll_logFC[i]), 
                                                which(rownames(top_logFC_Coef22) == ProbesAmongAll_logFC[i]), 
                                                which(rownames(top_logFC_Coef33) == ProbesAmongAll_logFC[i])))
        
        # Save the locations into a matrix with ncols = number of comparisons and nrows = number of rows
        Locations2_logFC <- matrix(do.call(c, Locations_logFC), byrow = T, ncol = length(levels(stage_status)))
        # dim(Locations2) # 1537    3
        
        # Volcano plot; Note 1.3 is -log10(0.05)
        if (ProduceImages == "Yes") {
          volcanoplot <- ggplot2::ggplot(top_logFC_Coef1, aes(logFC, -log10(adj.P.Val) )) +
                                         geom_point(color = "#e0e0e0") +
                                         ggtitle(paste0("Plot of ", colnames(Summary_Results)[1])) +  
                                         geom_text_repel(data = top_logFC_Coef1[Locations2_logFC[c(1:10), 1], ], 
                                                        aes(x = logFC, 
                                                            y = -log10(adj.P.Val), 
                                                            label = rownames(top_logFC_Coef11)[Locations2_logFC[c(1:10), 1]]),
                                                            fontface = 3 ) + 
                                         geom_point(data = top_logFC_Coef11[Locations2_logFC[c(1:5),1], ], 
                                                   aes(x = logFC, 
                                                       y = - log10(adj.P.Val)), 
                                                   colour = "#b2182b", 
                                                   size = 2) + 
                                         geom_point(data=top_logFC_Coef11[Locations2_logFC[c(6:10),1], ], 
                                                   aes(x = logFC, y = - log10(adj.P.Val)),
                                                   colour = "#4d4d4d", size=2) + 
                                         labs(color = "GeneClass", y="-log10(adj.P.Value)")+
                                         geom_hline(yintercept = 1.3, linetype = "dotted") +  
                                         theme_bw() + 
                                         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                                         scale_color_manual(values = c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
                                         theme(aspect.ratio = 1, text = element_text(size = 15))
          ggsave(paste0(pathNow, "/img/7_VolcanoPlot_", ContrastColumnName, "_1stComparison_FC.", PNGorPDF))
          
          
          volcanoplot <- ggplot2::ggplot(top_logFC_Coef2, aes(logFC, - log10(adj.P.Val) )) +
                                         geom_point(color = "#e0e0e0") +
                                         ggtitle(paste0("Plot of ", colnames(Summary_Results)[2])) +
                                         geom_text_repel(data = top_logFC_Coef2[Locations2_logFC[c(1:10), 2], ], 
                                                        aes(x = logFC, 
                                                            y = - log10(adj.P.Val), 
                                                            label = rownames(top_logFC_Coef22)[Locations2_logFC[c(1:10), 2]]),
                                                        fontface = 3 ) + 
                                         geom_point(data = top_logFC_Coef22[Locations2_logFC[c(1:5),2], ], 
                                                    aes(x = logFC, 
                                                        y = -log10(adj.P.Val)), 
                                                    colour = "#b2182b", size = 2) + 
                                         geom_point(data = top_logFC_Coef22[Locations2_logFC[c(6:10), 2], ], 
                                                    aes(x = logFC, y = -log10(adj.P.Val)),
                                                    colour = "#4d4d4d", size = 2) + 
                                         labs(color = "GeneClass", y = "-log10(adj.P.Value)")+
                                         geom_hline(yintercept = 1.3, linetype = "dotted") +  
                                         theme_bw() + 
                                         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                                         scale_color_manual(values = c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
                                         theme(aspect.ratio = 1, text = element_text(size = 15))
          ggsave(paste0(pathNow, "/img/7_VolcanoPlot_", ContrastColumnName, "_2ndComparison_FC.", PNGorPDF))
          
          
          volcanoplot <- ggplot2::ggplot(top_logFC_Coef3, aes(logFC, -log10(adj.P.Val) )) +
                                         geom_point(color = "#e0e0e0")+
                                         ggtitle(paste0("Plot of ", colnames(Summary_Results)[3])) +
                                         geom_text_repel(data = top_logFC_Coef3[Locations2_logFC[c(1:10), 3], ], 
                                                        aes(x = logFC, 
                                                            y = - log10(adj.P.Val), 
                                                            label = rownames(top_logFC_Coef33)[Locations2_logFC[c(1:10), 3]]),
                                                        fontface = 3 ) + 
                                         geom_point(data = top_logFC_Coef33[Locations2_logFC[c(1:5), 3], ], 
                                                    aes(x = logFC, 
                                                        y = - log10(adj.P.Val)), colour = "#b2182b", size = 2) + 
                                         geom_point(data = top_logFC_Coef33[Locations2_logFC[c(6:10),3], ], 
                                                    aes(x = logFC, y = -log10(adj.P.Val)),
                                                    colour = "#4d4d4d", size = 2) + 
                                         labs(color = "GeneClass", y = "-log10(adj.P.Value)")+
                                         geom_hline(yintercept = 1.3, linetype="dotted") +  
                                         theme_bw() + 
                                         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                                         scale_color_manual(values = c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
                                         theme(aspect.ratio = 1, text = element_text(size = 15))
          ggsave(paste0(pathNow, "/img/7_VolcanoPlot_", ContrastColumnName, "_3rdComparison_FC.", PNGorPDF))
        }
      }
    } else if (length(levels(stage_status)) == 2) {
      top_logFC <- limma::topTable(fit.reduced, coef = 1, number = nrow(MvalueMatrix), sort.by = "logFC")
      
      if (ProduceImages == "Yes") {
        volcanoplot <- ggplot2::ggplot(top_logFC, aes(logFC, -log10(adj.P.Val) )) +
                                       geom_point(color = "#e0e0e0")+
                                       ggtitle(paste0("Plot of ", colnames(Summary_Results))) +
                                       geom_text_repel(data = top_logFC[position_20genes, ], 
                                                      aes(x = logFC, 
                                                          y = -log10(adj.P.Val), 
                                                          label = rownames(top_logFC)[position_20genes]),
                                                          fontface = 3 ) + 
                                       geom_point(data = top_logFC[position_20genes[1:5], ], aes(x = logFC,
                                                                                                 y = -log10(adj.P.Val)), 
                                                  colour = "#b2182b", size = 2) + 
                                       geom_point(data = top_logFC[position_20genes[6:10], ], 
                                                  aes(x=logFC, y=-log10(adj.P.Val)),
                                                  colour = "#4d4d4d", size = 2) + 
                                       labs(color = "GeneClass", y = "-log10(adj.P.Value)") +
                                       geom_hline(yintercept = 1.3, linetype = "dotted") +  
                                       theme_bw() + 
                                       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                                       scale_color_manual(values = c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
                                       theme(aspect.ratio = 1, text = element_text(size=15))
        ggsave(paste0(pathNow,"/img/7_VolcanoPlot_",ContrastColumnName,ProbeOrGene,"_toplogFC.",PNGorPDF))
        
        # ggplot2::ggplot(top_logFC, aes(logFC, -log10(adj.P.Val) )) +
        #   geom_point(color="#e0e0e0") +
        #   ggtitle(paste0("Plot of ", colnames(Summary_Results)[1])) +
        #   geom_hline(yintercept= 2, linetype="dotted") + 
        #   labs(y="-log10(adj. P-value)") +
        #   theme_bw() + 
        #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        #   scale_color_manual(values=c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
        #   theme(aspect.ratio=1, text = element_text(size=15))
        
      }
    }
  }
    
    
    if(ProbeOrGene == "Probe" && length(levels(stage_status)) == 2) { 
      # gene ontology analysis based on topP
      sigCpGs <- rownames(top_P) # names of top 50 CpG sites
      gst_epic <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(MvalueMatrix), 
                         collection = c("GO","KEGG"), array.type = "EPIC")
      topPGOList <- topGSA(gst_epic)
      
      # gene ontology analysis based on logFC
      sigCpGs <- rownames(top_logFC) # names of top 50 CpG sites
      gst_epic_FC <- gometh(sig.cpg = sigCpGs, all.cpg = rownames(MvalueMatrix), 
                            collection = c("GO","KEGG"), array.type = "EPIC")
      topGO_list_FC <- topGSA(gst_epic)
    } else if(ProbeOrGene == "Gene") {
      gst_epic <- "For GO output, input of probes (not genes) needed"
      gst_epic_FC <- "For GO output, input of probes (not genes) needed"
      topPGOList <- "For GO output, input of probes (not genes) needed"
      topGO_list_FC <- "For GO output, input of probes (not genes) needed"
    } else if (length(levels(stage_status)) == 3) {
      gst_epic <- "Code needs to be adjusted for multiple comparisons"
      gst_epic_FC <- "Code needs to be adjusted for multiple comparisons"
      topPGOList <- "Code needs to be adjusted for multiple comparisons"
      topGO_list_FC <- "Code needs to be adjusted for multiple comparisons"
    }
    
    # Plotting
    if (ProduceImages == "Yes" && length(levels(stage_status)) == 2) {
      if (PNGorPDF == "pdf") {
        pdf(paste0(pathNow, "/img/7_Histogram_", ContrastColumnName, "_Pvalue.pdf"), 
            width = 50, height = 50, pointsize = 50)
      } else {
        png(paste0(pathNow, "/img/7_Histogram_", ContrastColumnName, "_Pvalue.png"))
      }
      hist(top_P$P.Value, breaks = 50)
      dev.off()
      
      if (PNGorPDF == "pdf") {
        pdf(paste0(pathNow, "/img/7_Histogram_", ContrastColumnName, "_logFC.pdf"), 
            width = 50, height = 50, pointsize = 50)
      } else {
        png(paste0(pathNow, "/img/7_Histogram_", ContrastColumnName, "_logFC.png"))
      }
      hist(top_logFC$logFC, breaks = 50)
      dev.off()
    }
    
    # making 1-D Scatter Plots
    cpgs <- rownames(top_logFC)
    # par(mfrow=c(2,2))
    # pdf(paste0(pathNow,"/img/3_ScatterPlots_logFC_probes.pdf"), width = 1600, height = 1600, pointsize = 25)
    # for(i in 1:4){
    # stripchart(beta[rownames(beta)==cpgs[i],]~design_epic2[,1],method="jitter",
    # group.names=tolower(levels(stage_status)),pch=16,cex=1.5,col=c(4,2),ylab="Beta values",
    # vertical=TRUE,cex.axis=1.5,cex.lab=1.5)
    # title(cpgs[i],cex.main=1.5)
    # }
    # dev.off()
    
    # Assembling a list of probes shared among comparisons, if multiple comparisons were present
    if (length(levels(stage_status)) == 3) {
      MultipleComparisonCommonProbes = list(ListPValues = ProbesAmongAll, 
                                            ListFCValues = ProbesAmongAll_logFC)
    } else if (length(levels(stage_status)) == 2) {
      MultipleComparisonCommonProbes = NA
    }
    
    # Added 15 Oct 2019
    # Removal of unwanted variation when performing a differential methylation
    # RUVfit uses the RUV-inverse method by default
    # k = number of components of unwanted variation to remove, where 0 ≤ k < no.samples
    # 2 cases: if negative controls are known: 
    #          if negative controls are not known apriori: identified via RUVm
    # http://bioconductor.org/packages/devel/bioc/vignettes/missMethyl/inst/doc/missMethyl.html#testing-for-differential-methylation-using
    
    
    removing_unwanted_variation <- function() {
      # Performing RUV - Remove Unwanted Variation from missMethyl
      # This function has not been generalized. 
      # Ref: http://bioconductor.org/packages/devel/bioc/vignettes/missMethyl/inst/doc/missMethyl.html#testing-for-differential-methylation-using
      
      # Stage 1 involves performing a differential methylation analysis
      # using RUV-inverse (by default) and the 613 Illumina negative controls
      # (INCs) as negative control features. This will produce a list of CpGs
      # ranked by p-value according to their level of association with the
      # factor of interest. This list can then be used to identify a set of
      # empirical control probes (ECPs), which will capture more of the 
      # unwanted variation than using the INCs alone.
      
      # Negative control features: probes/genes/etc. that are known a priori
      # to not truly be associated with the biological factor of interest, 
      # but are affected by unwanted variation. 
      
      # setup the factor of interest
      stage_status <- factor(stri_trim(ClinicalFile[which(ClinicalFile[ , which(colnames(ClinicalFile) == "TYPE")] == "FL"),
                                                     which(colnames(ClinicalFile) == "STAGE")]), 
                             levels = c("ADVANCED", "LIMITED"))
      grp <- factor(stri_trim(ClinicalFile[which(ClinicalFile[,which(colnames(ClinicalFile) == "TYPE")] == "FL"),
                                            which(colnames(ClinicalFile) == "STAGE")]), labels = c(1, 0))
      
      # extract Illumina negative control data
      INCs <- missMethyl::getINCs(RGChannelSet)
      #head(INCs)
      
      # Added step: removed samples from INCs
      colnames(INCs) <- Output_Data_1$SampleSheet$Sample_Name
      
      # Remove KRI prefix
      # Doesn't work, why?
      #colnames(INCs[, which(substr(colnames(INCs), 1, 3) == "KRI")])   <- as.character(substr(colnames(INCs[, which(substr(colnames(INCs), 1, 3) == "KRI")]), 5, 16))
      
      # Try another method - remove KRI prefix
      colnames <- Output_Data_1$SampleSheet$Sample_Name
      # Remove KRI prefix
      colnames[which(substr(colnames, 1, 3) == "KRI")] <- substr(colnames[which(substr(colnames, 1, 3) == "KRI")], 5, 16)
      colnames[which(substr(colnames, 10, 13) == "")] <- paste0(colnames[which(substr(colnames, 10, 13) == "")], "_T1")
      colnames(INCs) <- colnames
      # Manually remove "LY_FL_159_T1"
      INCs <- INCs[, - which(colnames(INCs) == "LY_FL_159_T1")]
      # Manually rename "LY_FL_159_T1_rep" to "LY_FL_159_T1"
      colnames(INCs)[which(colnames(INCs) == "LY_FL_159_T1_rep")] <- "LY_FL_159_T1"
      # Update INCs so that it has the same samples as Output_Remove_2_MvalueMatrix_T1
      INCs <- INCs[, match(colnames(MvalueMatrix), colnames(INCs))]
      
      # Because the contrast is STAGE, need to account for that in INCs **NEED GENERATLIZATION
      INCs <- INCs[, ! is.na(ClinicalFile$STAGE)]
      # dim(INCs) # 411 155
      MvalueMatrix <- MvalueMatrix[, ! is.na(ClinicalFile$STAGE)]
      # dim(MvalueMatrix) # 559110    155
      
      # Add negative control data to M-values
      Mc <- rbind(MvalueMatrix, INCs)
      
      # Create vector marking negative controls in data matrix
      ctl1 <- rownames(Mc) %in% rownames(INCs)
      cat("FALSE and  TRUE: ", table(ctl1))
      #  FALSE   TRUE 
      # 559110    411 
      # Samples in  experiment is not greater than the number of Illumina negative controls (170 < 411)
      
      # Stage 1 analysis
      # Remove unwanted variation when testing for differential methylation
      # ctl	= logical vector, length == nrow(Y). Features that are to be used as 
      # negative control variables are indicated as TRUE, all other features are FALSE.
      rfit1 <- missMethyl::RUVfit(Y = Mc, X = grp, ctl = ctl1) 
      # Post-process and summarize the results of call to RUVfit.
      rfit2 <- missMethyl::RUVadj(Y = Mc, fit = rfit1)
      
      # Designate the CpGs that are least associated with the factor
      # of interest based on FDR-adjusted p-value as ECPs
      # Table of top-ranked differentially methylated CpGs obatained from a differential methylation analysis using RUV
      # p.BH = numeric, cutoff value for Benjamini-Hochberg adjusted p-values.
      top1 <- missMethyl::topRUV(rfit2, num = Inf, p.BH = 1)
      # head(top1)
      # dim(top1) # 559521      9
      
      # Retain non signficant CpGs
      # F.p.BH = Benjamini-Hochberg adjusted p-values for testing all of the factors of interest simultaneously.
      # p.BH_X1	= Benjamini-Hochberg adjusted p-values for the factor of interest.
      ctl2 <- rownames(MvalueMatrix) %in% rownames(top1[top1$p.BH_X1.0 > 0.5, ])
      cat("FALSE and  TRUE: ", table(ctl2))
      #   FALSE   TRUE 
      #284598 274512
      
      # Stage 2 involves performing a second differential methylation 
      # analysis on the original data using RUV-inverse (by default) and the ECPs.
      
      # Perform RUV adjustment and fit
      rfit3 <- missMethyl::RUVfit(Y = MvalueMatrix, X = grp, ctl = ctl2) # Stage 2 analysis
      
      rfit4 <- missMethyl::RUVadj(Y = MvalueMatrix, fit = rfit3)
      
      # Look at table of top results
      top2 <- missMethyl::topRUV(rfit4) # Retain only probes (CpGs) with a p values < 0.05 
      # dim(top2)
      
      length(which(rfit4$C$F.p.BH < 0.05)) # 108593
      MvalueMatrix_rfit4 <- MvalueMatrix[rfit4$C$F.p.BH < 0.05,]
      BetaMatrix_rfit4 <- BetaMatrix[rfit4$C$F.p.BH < 0.05, ! is.na(ClinicalFile$STAGE)]
      # saveRDS(BetaMatrix_rfit4, 
      #        file = paste0(Path,"/7_BetaMatrix_rfit4_108593probes.rds"))
      
      BetaMatrix_rfit4 <- readRDS(file = '7_BetaMatrix_rfit4_108593probes.rds')
      
      
      # Visualising the effect of RUVm adjustment
      # getAdj function can be used to extract the adjusted values from the RUVm fit object produced by RUVfit
      # NOTE: The adjusted values should only be used for visualisations - it is NOT recommended that they are used in any downstream analysis.
      Madj <- missMethyl::getAdj(Y = MvalueMatrix, fit = rfit3) 
      
      par(mfrow = c(1, 2))
      limma::plotMDS(MvalueMatrix, labels = grp, col = as.integer(factor(grp)),
                     main = "Unadjusted", gene.selection = "common")
      legend("topleft",legend = c("Adv","Lim"), pch = 16, cex = 1, col = 1:2)
      limma::plotMDS(Madj, labels = grp, col = as.integer(factor(grp)),
                     main = "Adjusted: RUV-inverse", gene.selection = "common")
      legend("topleft", legend = c("Adv", "Lim"), pch = 16, cex = 1, col = 1:2)
      
      
      # Perform RPMM
      # Performs beta latent class modeling using recursively-partitioned mixture model
      rpmm <- RPMM::blcTree(t(BetaMatrix_rfit4), verbose=0)
      # Get weight matrix and show first few rows
      rpmmWeightMatrix <- RPMM::blcTreeLeafMatrix(rpmm)
      # Get class assignments and compare with tissue
      rpmmClass <- RPMM::blcTreeLeafClasses(rpmm) # classes are assigned based on terminal nodes
      rpmmClass_numbers <- c(rpmmClass) 
      
      # compare with previous results based on standard deviation
      point2SD_Probes_4ClusterModel <- read.csv(paste0(pathNow,
                                                       "/DataFiles/5394_point2SD_Probes_4ClusterModel_RPMM_FLonly_Labels.csv"), 
                                                row.names = 1)
      
      
      OrderNew <- match(colnames(BetaMatrix_rfit4), rownames(point2SD_Probes_4ClusterModel))
      table(rpmmClass_numbers, point2SD_Probes_4ClusterModel$V2[OrderNew])
      table(rpmmClass_numbers, ClinicalFile$STAGE[which(ClinicalFile$TYPE == "FL")])
      table(rpmmClass_numbers, ClinicalFile$TRANSLOC_14_18[which(ClinicalFile$TYPE == "FL")])
      
      # Analyze pathways 
      #
      # install.packages("BioMethyl_1.1.tar.gz", repos = NULL)  # Not working 
      # source("BioMethyl-masterR_Oct2019/calDEG.R")
      # source("BioMethyl-masterR_Oct2019/calExpr.R")
      # source("BioMethyl-masterR_Oct2019/calGSEA.R")
      # source("BioMethyl-masterR_Oct2019/filterMethyData.R")
      # source("BioMethyl-masterR_Oct2019/referCancerType.R")
      # Example - data not working 
      # BiocManager::install("ENmix")
      # data(MethData)
      # dat <- filterMethyData(MethData)
      # myExpr <- calExpr(dat, "BRCA", Example = T, SaveOut = F, "")
      # samp_1 <- colnames(myExpr)[1:10]
      # samp_2 <- colnames(myExpr)[11:20]
      # mydf <- calDEG(myExpr, Sample_1 = samp_1, Sample_2 = samp_2, SaveOut = F, "")
      # xx <- calGSEA(myExpr, mydf, DEGthr = c(0, 0.01), Sample_1 = samp_1, Sample_2 = samp_2, OutFile = "", GeneSet = "C2")    
      
      
      BiocManager::install("methylGSA")
      library(methylGSA)
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      data(cpgtoy)
      head(cpg.pval, 20)
      methylGSA::methylglm
      
      library(BiocParallel)
      pvalues_rfit4 <- rfit4$C$F.p.BH
      names(pvalues_rfit4) <- rownames(rfit4$C)
      
      # Implement logistic regression adjusting for number of probes in 
      # enrichment analysis
      res_p <- methylGSA::methylglm(cpg.pval = pvalues_rfit4, 
                                   array.type = "EPIC",
                                   GS.type = "KEGG", parallel = TRUE)
      head(res_p, 5)
      
      # Enrichment analysis after adjusting multiple p-values of each gene by
      # Robust Rank Aggregation
      res2 = methylGSA::methylRRA(cpg.pval = pvalues_rfit4, 
                                  array.type = "EPIC",
                                  method = "GSEA")
      head(res2, 5)
    }
    
    
    
    RESULTS <- list(SummaryOfResults = Summary_Results,
                    SignificantProbes = Significant_Probes,
                    MArrayLM = fit.reduced,
                    topP = top_P,
                    topPGOTesting=gst_epic,
                    topPGOList = topPGOList,
                    topLogFC = top_logFC,
                    topLogFCTesting= gst_epic_FC,
                    topLogFCGOList = topGO_list_FC,
                    MultipleComparisonCommonProbes = MultipleComparisonCommonProbes)
    
    class(RESULTS) <- "DifferentialMethylation_ASilva"
    return(RESULTS)
    
  }
  