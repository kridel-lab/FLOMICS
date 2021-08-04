# Updated 3 Aug 2021
# Updated 16 January 2020
# Function: Visualize pie chart of probes from beta matrix based on: probe type entries,
#           CpG context entries,  CpG context type, Fantom5 entries, transcription binding
#           site type, open chromatin regions, DNAse I hypersensitive region. 
# Author: Anjali Silva

# Input:
# BetaMatrix: Matrix of beta values for probes x patients, with probes in rows and patients as columns.
# AnnotationFile: Matrix of annotations for all the probes found in BetaMatrix. 
#                 It is of size probes x annotations.
# ProduceImages: Produce images or not, options = "Yes" or "No"; default "Yes" 
# PNGorPDF: Output format of the image, options = "png" or "pdf"; default "png"
# ImageName: Character string indicating the name for the image.

PieCharts30 <- function(BetaMatrix, 
                      AnnotationFile = NA, 
                      ProduceImages = "Yes", 
                      PNGorPDF = "png", 
                      ImageName = NA) {
  library("dplyr")  
  library("plotly")
  library("stringr")
  library("grDevices")
  
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  # Matching rownames in methylation file with corresponding location in annotation file
  match_ids_methylation_beta <- match(rownames(BetaMatrix), AnnotationFile$V1)
  
  # Bringing annotation file in the same order as BetaMatrix
  AnnotationFile2 <- data.frame(AnnotationFile[match_ids_methylation_beta, ])
  
  # Looking at chromosomes - creating a table
  chromosomeEntries <- data.frame(table(AnnotationFile2$chr))
  # Using ggplots to create pie-chart for chromosomes - not used as labels overlap
  # pie = ggplot(Chromosome_breakdown, aes(x = "", y = proportion, fill = Var1)) + 
  #       geom_bar(stat = "identity", width = 1) + 
  #       # Convert to pie (polar coordinates) and add labels
  #       coord_polar("y", start = 0) + 
  #       geom_text(aes(label = paste0(Var1, "(", round(proportion), "%)"), x = 1.4), position = position_stack(vjust = 0.5)) +
         # Remove labels and add title
  #       labs(x = NULL, y = NULL, fill = "Chromosomes", title = "Chromosomes") +
  #       theme_classic() + theme(axis.line = element_blank(),
  #                                     axis.text = element_blank(),
  #                                     axis.ticks = element_blank(),
  #                                     plot.title = element_text(hjust = 0.5, color = "#666666"))
  
  # Using plot_ly to create pie-chart for chromosomes
  chromosomeType <- plotly::plot_ly(chromosomeEntries, 
                                    labels = ~Var1, 
                                    values = ~Freq, 
                                    type = 'pie',
                                    textposition = 'outside', 
                                    textinfo = 'label+percent') %>%
                                    layout(title = '', showlegend = FALSE, 
                                    xaxis = list(showgrid = FALSE, 
                                                 zeroline = FALSE, showticklabels = FALSE),
                                    yaxis = list(showgrid = FALSE, 
                                                 zeroline = FALSE, showticklabels = FALSE))
  # Produce image 
  if(ProduceImages == "Yes") {
    chromosomeType
     if (PNGorPDF == "png") {
      grDevices::png(paste0(pathNow,"/img/30_PieChart_",
                            ImageName, "_chromosomeType.", PNGorPDF))
    }
    if (PNGorPDF == "pdf") { 
      pdf(paste0(pathNow,"/img/30_PieChart_", ImageName, 
                 "_chromosomeType.", PNGorPDF), width = 10, height = 20)
    }
  } 
 
   
  # Looking at probe type - creating a table
  probeTypeEntries <- data.frame(table(AnnotationFile2$Type))
  # Using plot_ly to create pie-chart for probe type 
  probeType <- plotly::plot_ly(probeTypeEntries, labels = ~Var1, values = ~Freq, type = 'pie',
                       textposition = 'outside',textinfo = 'label+percent') %>%
                       layout(title = '', showlegend = FALSE, 
                       xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                       yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  # Produce image 
  if(ProduceImages == "Yes") {
    probeType
    if (PNGorPDF == "png") {
      grDevices::png(paste0(pathNow,"/img/30_PieChart_", 
                            ImageName, "_probeType.", PNGorPDF))
    }
    if (PNGorPDF == "pdf") { 
      pdf(paste0(pathNow,"/img/30_PieChart_", ImageName, 
                 "_probeType.", PNGorPDF), width = 10, height = 20)
    }
  } 
  
  # Looking at probe CpGcontext
  CpGcontextEntries <- data.frame(table(AnnotationFile2$Relation_to_Island))
  # Using plot_ly to create pie-chart for chromosomes
  CpGcontextType <- plotly::plot_ly(CpGcontextEntries, 
                                    labels = ~Var1, values = ~Freq, 
                                    marker = list(colors = c("#a6cee3", "#1f78b4", 
                                                             "#fb9a99", "#fdbf6f", 
                                                             "#ff7f00", "#cab2d6")),
                                    type = 'pie',
                                    textposition = 'outside', textinfo = 'label+percent') %>%
                                    layout(title = '', showlegend = FALSE, 
                                    xaxis = list(showgrid = FALSE, 
                                                 zeroline = FALSE, showticklabels = FALSE),
                                    yaxis = list(showgrid = FALSE, 
                                                 zeroline = FALSE, showticklabels = FALSE))
  # Produce image 
  
  if(ProduceImages == "Yes") {
    if (PNGorPDF == "png") {
        grDevices::png(paste0(pathNow,"/img/30_PieChart_", 
                              ImageName, "_CpGcontextType.", PNGorPDF))
    }
    if (PNGorPDF == "pdf") { 
      pdf(paste0(pathNow,"/img/30_PieChart_", ImageName, 
                 "_CpGcontextType.", PNGorPDF), width = 10, height = 20)
    }
  }
  
  
  # Looking at probe Fantom5
  Fantom5entries <- AnnotationFile2$Phantom5_Enhancers
  Fantom5entries[which(AnnotationFile2$Phantom5_Enhancers != "")] <- 
    "Phantom 5 enhancer"
  Fantom5entries[which(AnnotationFile2$Phantom5_Enhancers == "")] <- 
    "No Phantom 5 enhancer"
  Fantom5Breakdown <- data.frame(table(Fantom5entries))
  # Using plot_ly to create pie-chart for Fantom5
  Fantom5Type <- plotly::plot_ly(Fantom5Breakdown, labels = ~Fantom5entries, values = ~Freq, type = 'pie',
                         textposition = 'outside',textinfo = 'label+percent') %>%
                         layout(title = '', showlegend = FALSE, 
                         xaxis = list(showgrid = FALSE, 
                                      zeroline = FALSE, showticklabels = FALSE),
                         yaxis = list(showgrid = FALSE, 
                                      zeroline = FALSE, showticklabels = FALSE))
  # Produce image 
  if(ProduceImages == "Yes") {
    if (PNGorPDF == "png") {
      grDevices::png(paste0(pathNow,"/img/30_PieChart_", 
                            ImageName, "_Fantom5Type.", PNGorPDF))
    }
    if (PNGorPDF == "pdf") { 
      pdf(paste0(pathNow,"/img/30_PieChart_", ImageName, 
                 "_Fantom5Type.", PNGorPDF), width = 10, height = 20)
    }
  }
  
  
  # Looking at probe transcription factor binding sites (TFBS)
  TFBSentries <- AnnotationFile2$TFBS_NAME
  TFBSentries[which(AnnotationFile2$Phantom5_Enhancers != "")] <- 
    "Transcription binding site"
  TFBSentries[which(AnnotationFile2$Phantom5_Enhancers == "")] <- 
    "Outside transcription binding site"
  TFBSBreakdown <- data.frame(table(TFBSentries))
  # Using plot_ly to create pie-chart for TFBS
  TFBSBType <- plotly::plot_ly(TFBSBreakdown, labels = ~TFBSentries, values = ~Freq, type = 'pie',
                       textposition = 'outside',textinfo = 'label+percent') %>%
                       layout(title = '', showlegend = FALSE, 
                       xaxis = list(showgrid = FALSE, 
                                    zeroline = FALSE, showticklabels = FALSE),
                       yaxis = list(showgrid = FALSE, 
                                    zeroline = FALSE, showticklabels = FALSE))
  # Produce image 
  if(ProduceImages == "Yes") {
    if (PNGorPDF == "png") {
      grDevices::png(paste0(pathNow,"/img/30_PieChart_", 
                            ImageName, "_TFBSBType.", PNGorPDF))
    }
    if (PNGorPDF == "pdf") { 
      pdf(paste0(pathNow,"/img/30_PieChart_", ImageName, 
                 "_TFBSBType.", PNGorPDF), width = 10, height = 20)
    }
  }
  
  
  # Looking at probe open chromatin region
  ChromatinEntries <- AnnotationFile2$OpenChromatin_NAME
  ChromatinEntries[which(AnnotationFile2$OpenChromatin_NAME != "")] <- "Open chromatin regions"
  ChromatinEntries[which(AnnotationFile2$OpenChromatin_NAME == "")] <- "Non-open chromatin regions"
  
  ChromatinBreakdown <- data.frame(table(ChromatinEntries))
  # Using plot_ly to create pie-chart for TFBS
  ChromatinType <- plotly::plot_ly(ChromatinBreakdown, 
                                   labels = ~ChromatinEntries, values = ~Freq, type = 'pie',
                                   textposition = 'outside',textinfo = 'label+percent') %>%
                                   layout(title = '', showlegend = FALSE, 
                                   xaxis = list(showgrid = FALSE, 
                                                zeroline = FALSE, showticklabels = FALSE),
                                   yaxis = list(showgrid = FALSE, 
                                                zeroline = FALSE, showticklabels = FALSE))
  # Produce image 
  if(ProduceImages == "Yes") {
    if (PNGorPDF == "png") {
      grDevices::png(paste0(pathNow,"/img/30_PieChart_", 
                            ImageName, "_ChromatinType.", PNGorPDF))
    }
    if (PNGorPDF == "pdf") { 
      pdf(paste0(pathNow,"/img/30_PieChart_", ImageName, 
                 "_ChromatinType.", PNGorPDF), width = 10, height = 20)
    }
  }
  
  
  # Looking at probe DNAse I hypersensitive region
  DNAseIHypEntries <- AnnotationFile2$DNase_Hypersensitivity_NAME
  DNAseIHypEntries[which(AnnotationFile2$OpenChromatin_NAME != "")] <- 
    "DNAse I hypersensitive region"
  DNAseIHypEntries[which(AnnotationFile2$OpenChromatin_NAME == "")] <- 
    "Outside a DNAse I hypersensitive region"
  DNAseIHypBreakdown <- data.frame(table(DNAseIHypEntries))
  # Using plot_ly to create pie-chart for TFBS
  DNAseIHypType <- plotly::plot_ly(DNAseIHypBreakdown, 
                                   labels = ~DNAseIHypEntries, values = ~Freq, type = 'pie',
                                   textposition = 'outside',textinfo = 'label+percent') %>%
                                   layout(title = '', showlegend = FALSE, 
                                   xaxis = list(showgrid = FALSE, 
                                                zeroline = FALSE, showticklabels = FALSE),
                                   yaxis = list(showgrid = FALSE, 
                                                zeroline = FALSE, showticklabels = FALSE))
  
  # Produce image 
  if(ProduceImages == "Yes") {
    if (PNGorPDF == "png") {
      grDevices::png(paste0(pathNow,"/img/30_PieChart_", 
                            ImageName, "_DNAseIHypType.", PNGorPDF))
    }
    if (PNGorPDF == "pdf") { 
      pdf(paste0(pathNow,"/img/30_PieChart_", ImageName, 
                 "_DNAseIHypType.", PNGorPDF), width = 10, height = 20)
    }
  }
  
  
  # Looking at probe gene region
  geneRegion <- sub("\\;.*", "", AnnotationFile2$UCSC_RefGene_Group)
  geneRegion[which(geneRegion == "")] <- "NA"
  geneRegion_breakdown <- data.frame(table(geneRegion))
  # Using plot_ly to create pie-chart for gene region
  geneRegionType <- plotly::plot_ly(geneRegion_breakdown, 
                                    labels = ~geneRegion, values = ~Freq, 
                                    marker = list(colors = c("#8dd3c7", "#ffffb3", 
                                                             "#bebada","#fb8072", 
                                                             "#80b1d3", "#fdb462", 
                                                             "#fccde5", "#d9d9d9")),
                                    type = 'pie',
                                    textposition = 'outside',textinfo = 'label+percent') %>%
                                    layout(title = '', showlegend = FALSE, 
                                    xaxis = list(showgrid = FALSE, 
                                                 zeroline = FALSE, showticklabels = FALSE),
                                    yaxis = list(showgrid = FALSE, 
                                                 zeroline = FALSE, showticklabels = FALSE))
  # Produce image 
  if(ProduceImages == "Yes") {
    if (PNGorPDF == "png") {
      grDevices::png(paste0(pathNow,"/img/30_PieChart_", 
                            ImageName, "_geneRegionType.", PNGorPDF))
    }
    if (PNGorPDF == "pdf") { 
      pdf(paste0(pathNow,"/img/30_PieChart_", ImageName, 
                 "_geneRegionType.", PNGorPDF), width = 10, height = 20)
    }
  }
  

  
  return(NULL)
}
# [END]
