#----------------------------------------------------------------------------------
#Gene Expression - Differential Gene Expression RNA-Seq
#----------------------------------------------------------------------------------

#Sarah Russell
#Date: Aug 24, 2020

#----------------------------------------------------------------------------------
#PACKAGES
#----------------------------------------------------------------------------------

# Loading needed packages
library(dplyr)
library(data.table)
library(edgeR)
library(EnhancedVolcano)
library(ggpubr)
library(grDevices)
library(KEGGprofile)
library(limma)
library(RColorBrewer)
library(stringi)
library(pheatmap)
library(readxl)

#----------------------------------------------------------------------------------
#DATA
#----------------------------------------------------------------------------------
date = Sys.Date()
setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_gene/")
pathNow <- getwd()
RNAseqCountMatrix=data.frame(fread("STAR_quantmode_counts_matrix_FL_136_patients.txt"))
AnnotationFileEnsemblToGenes=fread("gene_IDs.csv")
ClinicalFile=as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.xlsx"))

samples <- data.frame(names = c("gene",ClinicalFile$rna_seq_file_sample_ID)) %>%
  mutate(names = as.character(names))
RNAseqCountMatrix <- RNAseqCountMatrix %>% select(one_of(samples$names))
rownames(RNAseqCountMatrix)=RNAseqCountMatrix$gene
RNAseqCountMatrix$gene=NULL

#----------------------------------------------------------------------------------
#ANALYSIS
#----------------------------------------------------------------------------------

DifferentialExpressionRNAseq <- function(RNAseqCountMatrix,
                                         ContrastColumnName = "CLUSTER",
                                         AnnotationFileEnsemblToGenes,
                                         ClinicalFile,
                                         LogFCcutoff = 2,
                                         FDRcutoff = 0.05,
                                         ProduceImages = "Yes",
                                         ProduceHeatMap = "No",
                                         PNGorPDF = "pdf") {

#### Keeping all ensembl genes for now, filter for protein coding later

#  create a DGEList object:
exprsRNAseq <- edgeR::DGEList(counts = RNAseqCountMatrix)
exprsRNAseqNormFacs <- edgeR::calcNormFactors(exprsRNAseq)


# Set contrasts
###This part differs from original code based off different input/Clinical file formats. Code below is similar 
    ### to contrats set during Telescope edgeR analysis
  if(ContrastColumnName == "STAGE") {
    sample_info = as.data.table(filter(ClinicalFile, STAGE %in% c("ADVANCED", "LIMITED")))
    z = which(colnames(exprsRNAseqNormFacs) %in% sample_info$rna_seq_file_sample_ID)
    exprsRNAseqNormFacs = exprsRNAseqNormFacs[,z]
    con_sample_info = sample_info[order(match(rna_seq_file_sample_ID, colnames(exprsRNAseqNormFacs)))]
    print(con_sample_info$rna_seq_file_sample_ID == colnames(exprsRNAseqNormFacs))
    group=con_sample_info$STAGE
    statusContrastName=factor(group,levels = c("ADVANCED", "LIMITED"))

    designRNAseq <- stats::model.matrix (~0 + statusContrastName)
    colnames(designRNAseq) <- levels(statusContrastName)
    contrasts <- limma::makeContrasts("ADVANCED-LIMITED", levels = designRNAseq)

} else if(ContrastColumnName == "CLUSTER") {
    con_sample_info = ClinicalFile[order(match(rna_seq_file_sample_ID, colnames(exprsRNAseqNormFacs)))]
    print(con_sample_info$rna_seq_file_sample_ID == colnames(exprsRNAseqNormFacs))
    group=as.character(con_sample_info$Cluster)
    statusContrastName=factor(group,levels = c("1", "2"))

    designRNAseq <- stats::model.matrix (~0 + statusContrastName)
    colnames(designRNAseq) <- c("Cluster1","Cluster2")
    contrasts <- limma::makeContrasts("Cluster1-Cluster2", levels = designRNAseq)
}

  # Estimate dispersion.
  # Estimate the genewise dispersion estimates over all genes, allowing for a possible
  # abundance trend.
  # Allowing gene-specific dispersion is necessary in order that differential
  # expression is not driven by outliers. Therefore the tagwise dispersions are
  # strongly recommended in model fitting and testing for differential expression.
  # This method estimates common dispersion, trended dispersions and tagwise dispersions
  # in one run and is recommended.

  dispEstimates <- edgeR::estimateDisp(y = exprsRNAseqNormFacs,
                                       design = designRNAseq,
                                       robust = TRUE)

  # The dispersion estimates can be viewed in a BCV plot
  # Plot the genewise biological coefficient of variation (BCV) against gene abundance
  # (in log2 counts per million).
  # edgeR::plotBCV(dispEstimates)

  # quasi-likelihood (QL) dispersions can be estimated and visualized
  dispEstimateQL <- edgeR::glmQLFit(y = dispEstimates,
                                    design = designRNAseq,
                                    robust = TRUE)
  # plotQLDisp(dispEstimateQL)

  # Fit generalized linear model
  # Such a model is an extension of classical linear models to non-normally distributed
  # response data.
  fitGenLModel <- edgeR::glmQLFit(y = dispEstimates,
                                  design = designRNAseq)

  # Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests



# Differential gene expression analysis
 if(ProduceImages == "Yes") {
   if (PNGorPDF == "pdf") {
      pdf(paste("/cluster/home/srussell/logFCsVSaverageCountSize",ContrastColumnName, ".pdf",sep="_"), width=60, height=60, pointsize=50)
   }
   
   if (ncol(contrasts) == 1) {
     par(mfrow = c(1, 1))
   } else if (ncol(contrasts) == 2) {
     par(mfrow = c(1, 2))
   } else if (ncol(contrasts) <= 4) {
     par(mfrow = c(2, 2))
   } else if (ncol(contrasts) <= 6) {
     par(mfrow = c(2, 3))
   }
 }


 
 # Fit model and plot log CPM against logFC
 topTagsPvalue <- filterFCandFDR <- IDsavelogFCandFDRup <- IDsavelogFCandFDRdown <-
   IDsavelogFCandFDRupGeneList <- IDsavelogFCandFDRdownGeneList <- list()
 for (i in 1:ncol(contrasts)) {
   # Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
   # Conduct genewise statistical tests for a given coefficient or contrast.
   fitGenLModQLTest <- edgeR::glmQLFTest(glmfit = fitGenLModel, 
                                         contrast = contrasts[, i])
   
   if(ProduceImages == "Yes") {
     # Plot all the logFCs against average count size, highlighting the DE genes
     limma::plotMD(fitGenLModQLTest, main = paste0("Contrast: ", colnames(contrasts)[i]))
     graphics::abline(h = c(-1, 1), col = "blue")
   }
   
   # Extracts the most differentially expressed genes ranked by p-value 
   # Issue present with duplicated issues. Not sure why.
   topTagsPvalue[[i]] <- data.frame(edgeR::topTags(fitGenLModQLTest, 
                                                   n = Inf, 
                                                   adjust.method = "BH", 
                                                   sort.by = "PValue")) %>% 
                                                   mutate(ensemblGeneId = row.names(.))
   
   
   # AllDifferentilalENSEMBLids[67]
   # Filter based on prespecified logFC and FDR cutoff values. 
   # Note for logFC, the absolute value is taken.
   filterFCandFDR[[i]] <- topTagsPvalue[[i]] %>% dplyr::filter(abs(logFC) > LogFCcutoff & FDR < FDRcutoff) 
   cat("\n Contrast:", colnames(contrasts)[i], "has total of", nrow(filterFCandFDR[[i]]), 
       "ENSEMBL IDs that passed logFC > ", LogFCcutoff, "and FDR < ", FDRcutoff, ".")
   
   # Seperate filtered values based on + or - logFC
   IDsavelogFCandFDRup[[i]] <- filterFCandFDR[[i]] %>% dplyr::filter(logFC > 0) 
   IDsavelogFCandFDRupGeneList[[i]] <- IDsavelogFCandFDRup[[i]]$ensemblGeneId
   cat("\n Contrast:", colnames(contrasts)[i], "has", nrow(IDsavelogFCandFDRup[[i]]), 
       "ensembl IDs (+ve)logFC and significant.")
   IDsavelogFCandFDRdown[[i]] <- filterFCandFDR[[i]] %>% dplyr::filter(logFC < 0) 
   IDsavelogFCandFDRdownGeneList[[i]] <- IDsavelogFCandFDRdown[[i]]$ensemblGeneId
   cat("\n Contrast:", colnames(contrasts)[i], "has", nrow(IDsavelogFCandFDRdown[[i]]), 
       "ensembl IDs (-ve)logFC and significant.")
 }
 if(ProduceImages == "Yes") {
  grDevices::dev.off()
}


  # name list with contrast names
  names(topTagsPvalue) <- names(filterFCandFDR) <- names(IDsavelogFCandFDRup) <- 
    names(IDsavelogFCandFDRdown) <- colnames(contrasts)
  
  if(ProduceImages == "Yes") {
    # Enhanced volcano plot
    for (i in 1:ncol(contrasts)) {
      par(mfrow = c(1, 1))
      
      # https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
      EnhancedVolcano::EnhancedVolcano(toptable = topTagsPvalue[[i]],
                                       lab = topTagsPvalue[[i]]$ensemblGeneId,
                                       x = 'logFC',
                                       y = 'PValue',
                                       pCutoff = c(0.01, 0.05),
                                       FCcutoff = 2.0,
                                       title = paste0("VolcanoPlot Contrast:", colnames(contrasts)[i]),
                                       xlim = c(-15, 15),
                                       gridlines.major = FALSE,
                                       gridlines.minor = FALSE)
      
      if (PNGorPDF == "pdf") {
        ggplot2::ggsave(paste0("/cluster/home/srussell/EnhancedVolcano_Contrast_", colnames(contrasts)[i], ".pdf"))
      } 
    } 
    
    # GGplot2 volcano plot
    for (i in 1:ncol(contrasts)) {
      grDevices::dev.off()
      par(mfrow = c(1, 1))
      VolcanoPlot <- ggplot2::ggplot(topTagsPvalue[[i]]) +
                                    geom_point(aes(x = logFC, y= -log10(PValue))) +
                                    ggtitle(paste0("Contrast: ", colnames(contrasts)[i])) +
                                    labs(y = "-log10(adjusted P.Value)", x = "logFC") +
                                    geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
                                               col = "#4d9221", size = 1) +
                                    geom_hline(yintercept = -log10(0.001), linetype = "dashed", 
                                               col = "#762a83", size = 1) +
                                    theme_bw()  +
                                    theme(text = element_text(size=15), 
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(), 
                                          panel.background = element_rect(colour = "black", size = 1.5), 
                                          axis.title =  element_text(face = "bold"),
                                          aspect.ratio = 1, legend.position = "none", 
                                          axis.text.x = element_text(face = "bold"), 
                                          axis.text.y = element_text(face = "bold")) 
      
      ggplot2::ggsave(paste0("/cluster/home/srussell/VolcanoPlot_GGplot_Contrast_", colnames(contrasts)[i], ".png"))
    }
  }
  
  RESULTS <- list(TopTagsEachContrast = topTagsPvalue,
                    AllDEIDsTable = filterFCandFDR)
    class(RESULTS) <- "DifferentialExpressionRNAseq"
    return(RESULTS)      
  }
