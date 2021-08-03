# Updated 3 Aug 2021
# Updated 3 June 2019
# Function: bumphunter package looks for genomic regions that are differentially methylated 
#          between two conditions; DMRcate extracts the most differentially methylated 
#          regions (DMRs) and variably methylated regions (VMRs) from IlluminaR Infinium
#          BeadChip Array samples via kernel smoothing.
# Ref: https://htmlpreview.github.io/?https://github.com/hansenlab/tutorial.450k/blob/master/inst/doc/methylation450k.html#bumphunter-to-find-differentially-methylated-regions-dmrs
# Ref: http://bioconductor.org/packages/release/bioc/vignettes/DMRcate/inst/doc/DMRcate.pdf
#      https://rpubs.com/anavoj/133334

# Input
# Method: Method to be used, options = "bumphunter", "DMRcate".
# BetaMatrix: A matrix of beta values for probes (or genes) x patients, with probes (or genes) 
#             in rows and patients as columns.
# ContrastColumnName: Type of differential analysis to be done. ContrastColumnName options 
#                     include "STAGE", "SEX", "SITE_BIOPSY", "TYPE_BIOPSY", "INSTITUTION", "COO", 
#                     "TYPE", "TRANSLOC_14_18", "EPIC_QC", "CLUSTER".
# ClinicalFile: File with patient sample names, and categories. 
# ClusterLabels: If "ContrastColumnName" == "CLUSTER", vector of integers should be provided.
# AnnotationFile: Matrix of annotations for all the probes found in BetaMatrix. It is of size
#                 probes x annotations.
# ProduceImages: Produce images or not, options = "Yes" or "No"
# DMR: One integer (only) indicating which DMR to be plotted for DMRcate. This only affects plotting. 
# PNGorPDF: Output format of the image, options = "png" or "pdf".
# ExpressionFile: A string indicating the name of expression file in text format, if present.
#                 Otherwise, set to NA. 

# Output
# Bumphunter_Output
# ****** Under construction

# DMRcate_output
# AnnotatedCpGs: A list length of the number of contrasts in ContrastColumnName, that contains 
#                annotated CpGs with their chromosome position and test statistic.
# DifferentiallyMethylatedRegions: A list length of the number of contrasts in ContrastColumnName, 
#                                 that contains significantly differentially (or variable) methylated
#                                 regions.
# GRangesObject: A list length of the number of contrasts in ContrastColumnName, that contains 
#                corresponding GRanges object. 


# Visuals saved to img folder
# 9_DifferentiallyMethylatedRegions_DMR_*


# Note to self 
#************** FIX line 379

DiffMethylatedRegions9 <- function(Method = "DMRcate", 
                                  BetaMatrix,
                                  MvalueMatrix, 
                                  ContrastColumnName = "TYPE", 
                                  ClinicalFile, 
                                  ClusterLabels = NA, 
                                  AnnotationFile, 
                                  ProduceImages = "Yes", 
                                  DMR = 1, 
                                  PNGorPDF = "png", 
                                  ExpressionFile = "differential_gene_expression_tum_env_rcd1Aug2019.txt") {
  
  # Loading needed packages
  library(stringi)
  library(limma)
  
  # Remove empty entries in the clinical file
  ClinicalFile[ClinicalFile == ""] <- NA
  
  if(! is.na(ExpressionFile)) {
    TumEnvGeneExp <- read.delim(file = ExpressionFile) 
  }
  
  # defining functions for each method
  bumphunter_function <- function(BetaMatrix, 
                                  ContrastColumnName, 
                                  ClinicalFile, 
                                  AnnotationFile, 
                                  ProduceImages, 
                                  PNGorPDF) {
  # ********** Under construction in terms of defining  model matrix and methylation data track ************
  cat("\n Running bumphunter... Note: bumps are identified but plotting under construction")  
  library(bumphunter)  
    
  # Getting the path of the file, which should contain a folder called "img"
  pathNow <- getwd()
  
  # Define  phenotype of interest
  # One way to define model matrix
  if(ContrastColumnName == "TYPE") {
    stage_status <- factor(stringi::stri_trim(ClinicalFile[ , 
                           which(colnames(ClinicalFile) == "TYPE")]), 
                           levels = c("DLBCL","FL","RLN"))
    design_epic2 <- stats::model.matrix (~stage_status)
    
    # obtaining only the rows in Beta matrix from annotation file 
    matchingrows <- match(rownames(BetaMatrix),  AnnotationFile$X)
    
    #Found 58 bumps.
    DLBCLvsFL <- bumphunter(object = BetaMatrix, 
                            design = design_epic2, 
                            cutoff = 0.2, 
                            B = 100, 
                            nullMethod = "bootstrap",
                            coef = 2,
                            type = "Beta",
                            chr = as.vector(AnnotationFile$chr[matchingrows]),
                            pos = as.vector(AnnotationFile$pos[matchingrows]))
    
    # Found 268 bumps.
    DLBCLvsRLN <- bumphunter(object = BetaMatrix, 
                             design = design_epic2, 
                             cutoff = 0.2, 
                             coef = 3,
                             B = 100, 
                             nullMethod = "bootstrap",
                             type = "Beta",
                             chr = as.vector(AnnotationFile$chr[matchingrows]),
                             pos = as.vector(AnnotationFile$pos[matchingrows]))

    
    type2 <- relevel(stage_status, "FL")
    designEpic3 <- stats::model.matrix(~ type2)
    
    # Found 190 bumps.
    FLvsRLN <- bumphunter(object = BetaMatrix, 
                          design = designEpic3, 
                          cutoff = 0.2, 
                          coef=3,
                          B = 100,
                          nullMethod = "bootstrap",
                          type = "Beta",
                          chr = as.vector(AnnotationFile$chr[matchingrows]),
                          pos = as.vector(AnnotationFile$pos[matchingrows]))
    
    if(ProduceImages == "Yes") {
      if (PNGorPDF == "png") {
        png(paste0(pathNow,"/img/9_DifferentiallyMethylatedRegions_Bumphunter.",PNGorPDF))
      }
      if (PNGorPDF == "pdf") { 
        pdf(paste0(pathNow,"/img/9_DifferentiallyMethylatedRegions_Bumphunter.",PNGorPDF))
      }
    
    genome <- "hg19"
    # NOTE: Using the non-permuted 
    dmr <- bh_dmrs2$table[1, ]
    chrom <- dmr$chr
    start <- dmr$start
    end <- dmr$end
    minbase <- start - 0.25 * (end - start)
    maxbase <- end + 0.25 * (end - start)
    pal <- c("#E41A1C", "#377EB8")
    
    # Start building the tracks
    iTrack <- IdeogramTrack(genome = genome, 
                            chromosome = dmr$chr, 
                            name = "")
    
    gTrack <- GenomeAxisTrack(col = "black", 
                              cex = 1, 
                              name = "", 
                              fontcolor = "black")
    
    # NOTE: This track takes a little while to create
    rTrack <- UcscTrack(genome = genome, 
                        chromosome = chrom, 
                        track = "ccdsGene",
                        from = minbase,
                        to = maxbase, 
                        trackType = "GeneRegionTrack",
                        rstarts = "exonStarts", 
                        rends = "exonEnds", 
                        gene = "name",
                        symbol = "name2", 
                        transcript = "name",
                        strand = "strand",
                        fill = "darkblue",
                        stacking = "squish", 
                        name = "RefSeq",
                        showId = TRUE, 
                        geneSymbol = TRUE)
    
    ##   'arg' should be one of "1000G Ph1 Accsbl", "1000G Ph1 Vars", "1000G Ph3 Accsbl", "1000G Ph3 Vars", "5% Lowest S", "acembly", "AceView Genes", "Affy Exon Array", "Affy GNF1H", "Affy RNA Loc", "Affy U133", "Affy U133Plus2", "Affy U95", "affyExonArray", "affyGnf1h", "affyU133", "affyU133Plus2", "affyU95", "All SNPs(138)", "All SNPs(141)", "All SNPs(142)", "All SNPs(144)", "All SNPs(146)", "All SNPs(147)", "All SNPs(150)", "All SNPs(151)", "Allen Brain", "allenBrainAli", "allHg19RS_BW", "Alt Haplotypes", "altLocations", "Assembly", "AUGUSTUS", "augustusGene", "BAC End Pairs", "bacEndPairs", "Broad ChromHMM", "Broad Histone", "BU ORChID", "Burge RNA-seq", "burgeRnaSeqGemMapperAlign", "Caltech RNA-seq", "Cand. Gene Flow", "CCDS", "ccdsGene", "CD34 DnaseI", "CGAP SAGE", "cgapSage", "chainSelf", "Chromosome Band", "Chromosome Band (Ideogram)", "ClinGen CNVs", "clinvar", "ClinVar Variants", "cnvDevDelay", "Common SNPs(138)", "Common SNPs(141)", "Common SNPs(142)", "Common SNPs(144)", 
    
    
    ##### ?????????????????????????????????????
    # methylation data track
    gr <- granges(GRset.funnorm)
    gr$beta <- getBeta(GRset.funnorm)
    
    methTrack <- DataTrack(range = gr,
                           groups = targets$status,
                           genome = genome,
                           chromosome = chrom, 
                           ylim = c(-0.05, 1.05), 
                           col = pal,
                           type = c("a","p"), 
                           name = "DNA Meth.\n(beta value)",
                           background.panel = "white", 
                           legend = TRUE, 
                           cex.title = 0.8,
                           cex.axis = 0.8, 
                           cex.legend = 0.8)
    
    # DMR position data track
    dmrTrack <- AnnotationTrack(start = start, 
                                end = end, 
                                genome = genome, 
                                name = "DMR",
                                chromosom = chrom)
    
    # Finally, plot the tracks
    tracks <- list(iTrack, gTrack, methTrack, dmrTrack, rTrack)
    sizes <- c(2, 2, 5, 2, 3) # set up the relative sizes of the tracks
    plotTracks(tracks, 
               from = minbase, 
               to = maxbase, 
               showTitle = TRUE, 
               add53 = TRUE,
               add35 = TRUE, 
               grid = TRUE, 
               lty.grid = 3, 
               sizes = sizes, 
               length(tracks))
    
    dev.off()
    }
    
    RESULTS <- list(DLBCLvsFL = DLBCLvsFL,
                    DLBCLvsRLN = DLBCLvsRLN,
                    FLvsRLN = FLvsRLN)
    
    class(RESULTS) <- "DifferentiallyMethylatedRegions_Bumphunter"
    return(RESULTS)
  }

  if(ContrastColumnName == "STAGE") {
    stage_status <- factor(stringi::stri_trim(ClinicalFile[which(ClinicalFile[ , 
                           which(colnames(ClinicalFile) == "TYPE")] == "FL"), 
                           which(colnames(ClinicalFile) == "STAGE")]), 
                           levels = c("ADVANCED", "LIMITED"))
    design_epic2 <- stats::model.matrix (~stage_status)
    
    # obtaining only the rows in Beta matrix from annotation file 
    matchingrows <- match(rownames(BetaMatrix[which(ClinicalFile[ , 
                          which(colnames(ClinicalFile) == "TYPE")] == "FL"), ]), AnnotationFile$X)
    
    AdvVsLimited <- bumphunter(object = BetaMatrix[which(ClinicalFile[, 
                               which(colnames(ClinicalFile) == "TYPE")] == "FL"), ], 
                               design = design_epic2, 
                               cutoff = 0.2, 
                               B = 100, 
                               coef = 2,
                               type = "Beta",
                               chr = as.vector(AnnotationFile$chr[matchingrows]),
                               pos = as.vector(AnnotationFile$pos[matchingrows]))
  }
  
  RESULTS <- list(Bumphunter_Output = bh_dmrs2)
  
  class(RESULTS) <- "DifferentiallyMethylatedRegions_Bumphunter"
  return(RESULTS)
  
  }
  
  DMRcate_function <- function(BetaMatrix, 
                               ContrastColumnName, 
                               ClinicalFile, 
                               AnnotationFile, 
                               ProduceImages, 
                               PNGorPDF) {
    
    cat("\n Running DMRcate... ")
    library(DMRcate)
    library(Gviz)
    
    # ********************** Example **********************
    #
    # data(dmrcatedata)
    # myMs <- logit2(myBetas)
    # nrow(snpsall)
    # nrow(myMs)
    # myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05) # 9382 probes
    # nrow(myMs.noSNPs)
    # patient <- factor(sub("-.*", "", colnames(myMs)))
    # type2 <- factor(sub(".*-", "", colnames(myMs)))
    # Releveling 
    # type3<-relevel(type2, "Tumour")
    # design3 <- stats::model.matrix(~~patient+ type3)
    # design <- stats::model.matrix(~patient + type2)
    # myannotation2 <- cpg.annotate("array", myMs.noSNPs, what="M", arraytype = "450K",
    #                              analysis.type="differential", design=design, coef=39)
    # dmrcoutput2 <- dmrcate(myannotation2, lambda=1000, C=2)
    # results.ranges2 <- extractRanges(dmrcoutput2, genome = "hg19")
    # results.ranges2
    # groups2 <- c(Tumour="magenta", Normal="forestgreen")
    # cols2 <- groups2[as.character(type2)]
    # samps2 <- c(1:6, 38+(1:6))
    # DMR.plot(ranges=results.ranges2, dmr=1, CpGs=myBetas, what="Beta", arraytype = "450K",
    #          phen.col=cols2, genome="hg19", samps=samps2)
    #
    # ********************** Example ends **********************
    
    
    # Set contrasts 
    if(ContrastColumnName == "STAGE") {
      # check if NA values are present, if so remove those patients 
      if(length(which(is.na(ClinicalFile[, which(colnames(ClinicalFile) == "STAGE")]) == TRUE)) > 0) {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[ , 
                                                  which(colnames(ClinicalFile) == "STAGE")]) == TRUE)), 
                                                  which(colnames(ClinicalFile) == "STAGE")]), 
                                                  levels = c("ADVANCED", "LIMITED"))
        MvalueMatrix <- MvalueMatrix[, which(is.na(ClinicalFile[, 
                                       which(colnames(ClinicalFile) == ContrastColumnName)]) == FALSE)]
      } else {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[ , 
                               which(colnames(ClinicalFile) == "STAGE")]), 
                               levels = c("ADVANCED", "LIMITED"))
      }
      design_epic2 <- stats::model.matrix (~0 + stage_status)
      colnames(design_epic2) <- levels(stage_status)
      contrasts <- limma::makeContrasts("ADVANCED-LIMITED", levels = design_epic2)
    } else if(ContrastColumnName == "SEX") {
      if(length(which(is.na(ClinicalFile[, which(colnames(ClinicalFile) == "SEX")]) == TRUE)) > 0) {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[-c(which(is.na(ClinicalFile[, 
                                                  which(colnames(ClinicalFile) == "SEX")]) == TRUE)), 
                                                  which(colnames(ClinicalFile) == "SEX")]),
                                                  levels = c("M", "F"), 
                                                  labels = c("Male", "Female"))
        MvalueMatrix <- MvalueMatrix[, which(is.na(ClinicalFile[, 
                        which(colnames(ClinicalFile) == ContrastColumnName)]) == FALSE)]
      } else {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[, 
                               which(colnames(ClinicalFile) == "SEX")]),  
                               levels = c("M", "F"), labels = c("Male","Female"))
      }
      design_epic2 <- stats::model.matrix(~0 + stage_status)
      colnames(design_epic2) <- levels(stage_status)
      contrasts <- limma::makeContrasts("Female-Male", levels = design_epic2)
    } else if(ContrastColumnName == "SITE_BIOPSY") {
      if(length(which(is.na(ClinicalFile[,which(colnames(ClinicalFile) == "SITE_BIOPSY")]) == TRUE)) > 0) {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[- 
                               c(which(is.na(ClinicalFile[, which(colnames(ClinicalFile) == "SITE_BIOPSY")]) == TRUE)), 
                               which(colnames(ClinicalFile) == "SITE_BIOPSY")]),  levels = c("EN", "LN"))
        MvalueMatrix <- MvalueMatrix[, which(is.na(ClinicalFile[, 
                                    which(colnames(ClinicalFile) == ContrastColumnName)]) == FALSE)]
      }else {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[, 
                                                  which(colnames(ClinicalFile) == "SITE_BIOPSY")] == TRUE))), 
                                                  which(colnames(ClinicalFile) == "SITE_BIOPSY")]), 
                                                  levels = c("EN", "LN"))
      }
      design_epic2 <- stats::model.matrix(~0 + stage_status)
      colnames(design_epic2) <- levels(stage_status)
      contrasts <- limma::makeContrasts("LN-EN", levels = design_epic2)
    } else if(ContrastColumnName == "TYPE_BIOPSY") {
      if(length(which(is.na(ClinicalFile[,which(colnames(ClinicalFile) == "TYPE_BIOPSY")]) == TRUE)) > 0) {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[-c(which(is.na(ClinicalFile[, 
                                                  which(colnames(ClinicalFile) == "TYPE_BIOPSY")]) == TRUE)), 
                                                  which(colnames(ClinicalFile) == "TYPE_BIOPSY")]), 
                                                  levels = c("CORE", "TISSUE"))
        MvalueMatrix <- MvalueMatrix[, which(is.na(ClinicalFile[, 
                                                   which(colnames(ClinicalFile) == ContrastColumnName)]) == FALSE)]
      } else {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[-c(which(is.na(ClinicalFile[, 
                                                  which(colnames(ClinicalFile) == "TYPE_BIOPSY")] == TRUE))), 
                                                  which(colnames(ClinicalFile) == "TYPE_BIOPSY")]), 
                                                  levels = c("CORE", "TISSUE"))
      }
      design_epic2 <- stats::model.matrix(~0 + stage_status)
      colnames(design_epic2) <- levels(stage_status)
      contrasts <- limma::makeContrasts("TISSUE-CORE", levels = design_epic2)
    } else if(ContrastColumnName == "INSTITUTION") {
      # https://support.bioconductor.org/p/44216/
      if(length(which(is.na(ClinicalFile[,which(colnames(ClinicalFile) == "INSTITUTION")]) == TRUE)) > 0) {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[,
                               which(colnames(ClinicalFile) == "INSTITUTION")]) == TRUE)), 
                               which(colnames(ClinicalFile) == "INSTITUTION")]), 
                               levels = c("BCCA", "UHN","JGH"))
        MvalueMatrix <- MvalueMatrix[, which(is.na(ClinicalFile[, 
                        which(colnames(ClinicalFile) == ContrastColumnName)]) == FALSE)]
      } else {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[, which(colnames(ClinicalFile) == "INSTITUTION")]), 
                               levels = c("BCCA", "UHN","JGH"))
      }
      design_epic2 <- stats::model.matrix(~0 + stage_status)
      colnames(design_epic2) <- levels(stage_status)
      contrasts <- limma::makeContrasts("UHN-BCCA", "JGH-BCCA", "JGH-UHN", levels = design_epic2)
    } else if(ContrastColumnName == "EPIC_QC") {
      if(length(which(is.na(ClinicalFile[,which(colnames(ClinicalFile) == "EPIC_QC")]) == TRUE)) > 0) {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[, 
                                                  which(colnames(ClinicalFile) == "EPIC_QC")]) == TRUE)), 
                                                  which(colnames(ClinicalFile) == "EPIC_QC")]), 
                                                  levels = c("Fair", "Good", "Poor"))
        MvalueMatrix<- MvalueMatrix[, which(is.na(ClinicalFile[, 
                       which(colnames(ClinicalFile) == ContrastColumnName)]) == FALSE)]
      } else {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[, 
                        which(colnames(ClinicalFile) == "EPIC_QC")]), 
                        levels = c("Fair","Good","Poor"))
      }
      design_epic2 <- stats::model.matrix(~0 + stage_status)
      colnames(design_epic2) <- levels(stage_status)
      contrasts <- limma::makeContrasts("Good-Fair", "Poor-Good", "Poor-Fair", 
                                        levels = design_epic2)
    } else if(ContrastColumnName == "COO") {
      # https://support.bioconductor.org/p/44216/ # contrasts for multiple comparisons
      if(length(which(is.na(ClinicalFile[, which(colnames(ClinicalFile) == "COO")]) == TRUE)) > 0) {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[-c(which(is.na(ClinicalFile[, 
                                                  which(colnames(ClinicalFile) == "COO")]) == TRUE)), 
                                                  which(colnames(ClinicalFile) == "COO")]), 
                                                  levels = c("GCB", "ABC"))
        MvalueMatrix <- MvalueMatrix[, which(is.na(ClinicalFile[, 
                        which(colnames(ClinicalFile) == ContrastColumnName)]) == FALSE)]
      } else {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[which(ClinicalFile[, 
                                                  which(colnames(ClinicalFile) == "TYPE")] == "DLBCL"), 
                                                  which(colnames(ClinicalFile) == "COO")]), 
                                                  levels = c("GCB", "ABC"))
      }
      design_epic2 <- stats::model.matrix(~0 + stage_status)
      colnames(design_epic2) <- levels(stage_status)
      contrasts <- limma::makeContrasts("ABC-GCB", levels = design_epic2)
    } else if(ContrastColumnName == "TYPE") {
      if(length(which(is.na(ClinicalFile[,which(colnames(ClinicalFile) == "TYPE")]) == TRUE)) > 0) {
        if (length(unique(ClinicalFile[,which(colnames(ClinicalFile) == "TYPE")])) == 2) {
          stage_status <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[ ,
                                                    which(colnames(ClinicalFile) == "TYPE")]) == TRUE)), 
                                                    which(colnames(ClinicalFile) == "TYPE")]), 
                                                    levels = c("DLBCL", "FL", "RLN"))
        } else {
          stage_status <- factor(stringi::stri_trim(ClinicalFile[- c(which(is.na(ClinicalFile[, 
                                                    which(colnames(ClinicalFile) == "TYPE")]) == TRUE)), 
                                                    which(colnames(ClinicalFile) =="TYPE")]), 
                                                    levels = c("DLBCL", "FL", "RLN"))}
        MvalueMatrix<- MvalueMatrix[, which(is.na(ClinicalFile[,
                       which(colnames(ClinicalFile) == ContrastColumnName)]) == FALSE)]
      }else{
        if (length(unique(ClinicalFile[, which(colnames(ClinicalFile) == "TYPE")])) == 2) {
          stage_status <- factor(stringi::stri_trim(ClinicalFile[, which(colnames(ClinicalFile) == "TYPE")]))
        } else {
            stage_status <- factor(stringi::stri_trim(ClinicalFile[,which(colnames(ClinicalFile) == "TYPE")]), 
                                   levels = c("DLBCL","FL","RLN"))
        }
      }
      design_epic2 <- stats::model.matrix(~ 0 + stage_status)
      colnames(design_epic2) <- levels(stage_status)
      if (length(unique(ClinicalFile[,which(colnames(ClinicalFile) == "TYPE")])) == 2) {
        contrasts <- limma::makeContrasts(paste0(unique(ClinicalFile[, which(colnames(ClinicalFile) == "TYPE")])[1], "-", 
                                                 unique(ClinicalFile[, which(colnames(ClinicalFile) == "TYPE")])[2]), 
                                          levels = design_epic2)
      } else { 
        contrasts <- limma::makeContrasts("DLBCL-FL", "DLBCL-RLN", "FL-RLN", levels = design_epic2)
      }
    } else if (ContrastColumnName == "TRANSLOC_14_18") {
      if (length(which(is.na(ClinicalFile[,which(colnames(ClinicalFile) == "TRANSLOC_14_18")]) == TRUE)) > 0) {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[-c(which(is.na(ClinicalFile[, 
                                                  which(colnames(ClinicalFile) == "TRANSLOC_14_18")]) == TRUE)),
                                                  which(colnames(ClinicalFile) == "TRANSLOC_14_18")]),  
                                                  levels = c("0","1"), 
                                                  labels = c("NoTranslocation", "Translocation"))
        MvalueMatrix <- MvalueMatrix[ , which(is.na(ClinicalFile[,
                                    which(colnames(ClinicalFile) == ContrastColumnName)]) == FALSE)]
      } else {
        stage_status <- factor(stringi::stri_trim(ClinicalFile[which(! is.na(ClinicalFile[, 
                                                  which(colnames(ClinicalFile) == "TRANSLOC_14_18")]) == TRUE),
                                                  which(colnames(ClinicalFile) == "TRANSLOC_14_18")]), 
                                                  levels = c("0","1"), 
                                                  labels = c("NoTranslocation", "Translocation"))
      }
      design_epic2 <- stats::model.matrix(~ 0 + stage_status)
      colnames(design_epic2) <- levels(stage_status)
      contrasts <- limma::makeContrasts("Translocation-NoTranslocation", 
                                        levels = design_epic2)
    } else if(ContrastColumnName == "CLUSTER") {
      if (sum(is.na(ClusterLabels)) > 0) {
        stop("ContrastColumnName is set to CLUSTER, but \n
             ClusterLabels not provided. Provide the ClusterLabels");}
      if (length(unique(ClusterLabels)) == 5) {
        stage_status <- factor(ClusterLabels, levels = c("1", "2", "3", "4", "5"), 
                               labels = c("one", "two", "three", "four", "five"))
      } else if (length(unique(ClusterLabels)) == 4) {
        stage_status <- factor(ClusterLabels, levels = c("1", "2", "3", "4"), 
                               labels = c("one", "two", "three", "four"))
      } else if(length(unique(ClusterLabels)) == 3) {
        stage_status <- factor(ClusterLabels, levels = c("1", "2", "3"), 
                               labels = c("one", "two", "three"))
      } else if(length(unique(ClusterLabels)) == 2) {
        stage_status <- factor(ClusterLabels, levels = c("1", "2"), 
                               labels = c("one", "two"))
      } else{
        stop("ClusterLabels have more than 5 clusters. \n
             Currently only upto 4 clusters supported");}

      design_epic2 <- stats::model.matrix(~0 + stage_status)
      colnames(design_epic2) <- levels(stage_status)
      if (length(unique(ClusterLabels)) == 5) {
        contrasts <- limma::makeContrasts("one-two", "one-three",
                                          "one-four", "one-five", 
                                          "two-three", "two-four", 
                                          "two-five", "three-four", 
                                          "three-five", "four-five", 
                                          levels = design_epic2)
      } else if (length(unique(ClusterLabels)) == 4) {
        contrasts <- limma::makeContrasts("one-two", "one-three",
                                          "one-four", "two-three",
                                          "two-four", "three-four", 
                                          levels = design_epic2)
      } else if (length(unique(ClusterLabels)) == 3) {
        contrasts <- limma::makeContrasts("one-two", "one-three", "two-three",
                                          levels = design_epic2)
      } else if(length(unique(ClusterLabels)) == 2) {
        contrasts <- limma::makeContrasts("one-two", 
                                          levels = design_epic2)
      }
      
    } else {
      cat("\n ***Warning: Not a valid ContrastColumnName.\n")
    }
    cat("\n Contrast has been set")
    
    cat("\n Filter out probes 2 nucleotides or closer \n 
        to a SNP that have a minor allele frequency greater than 0.05")
    # Measurements on the array may be confounded by proximity to SNPs, and cross-hybridisation to other areas of the genome
    # So filter out probes 2 nucleotides or closer to a SNP that have a minor allele frequency greater than 0.05
    myMs.noSNPs <- DMRcate::rmSNPandCH(MvalueMatrix, dist = 2, mafcut = 0.05) # 557951 probes
    
    cat("\n After filtering,", nrow(MvalueMatrix) - nrow(myMs.noSNPs), 
        "probes were removed from the original", 
        nrow(MvalueMatrix), "probes. \n Now there are", 
        nrow(myMs.noSNPs), "probes.")
    cat("\n Setting contrast for ContrastColumnName:", ContrastColumnName)
    
    cpgannotation <- function(i) {
      cat("\n Contrast considered is:", colnames(contrasts)[i], "\n")
      AnnotateCpG <- DMRcate::cpg.annotate(datatype = "array",
                                           object = myMs.noSNPs,
                                           what = "M",
                                           arraytype = "EPIC",
                                           analysis.type = "differential",
                                           design = design_epic2,
                                           contrasts = TRUE,
                                           cont.matrix = contrasts,
                                           coef = colnames(contrasts)[i], # column name in cont.matrix
                                           fdr = 0.05) 
      return(AnnotateCpG)
    }
    
    cat("\n Setting function for cpgannotation")
    cat("\n Annotate CpGs with their chromosome position and test statistic")
    # Annotate CpGs with their chromosome position and test statistic
    myannotations <- lapply(c(1:length(colnames(contrasts))), 
                            function(x) cpgannotation(x))
    names(myannotations) <- c(colnames(contrasts))
    print(myannotations)
    
    cat("\n Computes a kernel estimate against a null comparison to identify 
        significantly differentially (or variable) methylated regions.")
    # Computes a kernel estimate against a null comparison to identify significantly differentially (or variable) methylated regions.
    # dmrcoutput <- lapply(c(1:length(colnames(contrasts))), function(j) DMRcate::dmrcate(object = myannotations[[j]], lambda = 1000, C = 2))
    # 12 March 2020 - removed lapply to add tryCatch, so if error, the loop will still go on with out haltering.
    dmrcoutput <- list()
    for (j in 1:length(colnames(contrasts))) {
      dmrcoutput[[j]] <- NA
      tryCatch({
        cat("\n Running dmrcate for contrast: ", colnames(contrasts)[j])
        dmrcoutput[[j]] <- DMRcate::dmrcate(object = myannotations[[j]], lambda = 1000, C = 2)
        if(typeof(dmrcoutput[[j]]) == "logical") {stop(paste0("Contrast: ", 
                                  colnames(contrasts)[j], "didn't work."))}
      }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) 
    }
    
    names(dmrcoutput) <- c(colnames(contrasts))
    cat("\n DMRCate output:\n")
    print(dmrcoutput)
    
    cat("\n Convert the DMR list to a GRanges object")
    # Convert the DMR list to a GRanges object, which uses the genome argument to annotate overlapping promoter regions (+/- 2000 bp from TSS)
    # results.ranges <- lapply(c(1:length(colnames(contrasts))), function(k) DMRcate::extractRanges(dmrcoutput[[k]], genome = "hg19"))
    # 12 March 2020 - removed lapply to add tryCatch, so if error, the loop will still go on with out haltering.
    results.ranges <- list()
    for (k in 1:length(colnames(contrasts))) {
      results.ranges[[k]] <- NA
      tryCatch({
        cat("\n Running DMRcate::extractRanges for contrast: ", colnames(contrasts)[k])
        results.ranges[[k]] <- DMRcate::extractRanges(dmrcoutput[[k]], genome = "hg19")
        if(typeof(results.ranges[[k]]) == "logical") {stop(paste0("Contrast: ", 
                                      colnames(contrasts)[k], "didn't work."))}
      }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) 
      
      # when ranked based on adjusted p values, the top DMRs overlapped with the promoters of genes 
      # orderedbyStouffer <- results.ranges[[k]][order(results.ranges[[k]]@elementMetadata@listData$Stouffer)]
      
    }
    
    names(results.ranges) <- c(colnames(contrasts))
    cat("\n Done converting the DMR list to a GRanges object")
    
    
    # Obtain genes corresponding to each DMRcate region for given contrast
    geneVector <- chromosomeVector <- list()
    for (k in 1:length(colnames(contrasts))) {
      # Getting genes falling in DMRs
      # This is given by: results.ranges[[1]][, 8]$overlapping.genes 
      geneVector[[k]] <- unique(unlist(strsplit(results.ranges[[k]]@elementMetadata$overlapping.genes, 
                                split=", ")))
      chromosomeVector[[k]] <- unlist(lapply(c(1:length(results.ranges[[k]])), 
                              function(x) results.ranges[[k]]@seqnames[x]@values))
    }
    names(geneVector) <- c(colnames(contrasts))
    
    # Produce image 
    if(ProduceImages == "Yes") {
      cat("\n Produce image ")
      
      # Obtaining path to save images
      pathNow <- getwd()
      
      figuregenerator_DMR <- function(r) {
        if(DMR > length(results.ranges[[r]]$no.cpgs)) {
          stop("DMR argument cannot be greater than ", length(results.ranges[[r]]$no.cpgs)); }

          groups_FL_RLN <- c("1" = "blue", "2" = "green", "3"="maroon")
          types_FL_RLN <- stage_status
          cols_FL_RLN <- groups_FL_RLN[types_FL_RLN]
          DMRcate::DMR.plot(ranges = results.ranges[[r]], 
                            dmr = 1, 
                            CpGs = BetaMatrix[,which(stage_status %in% unlist(strsplit(colnames(contrasts)[r], "-")) == TRUE)],
                            what = "Beta", arraytype = "EPIC",
                            phen.col = cols_FL_RLN, genome="hg19",
                            samps = c(1:10, 40:50))
                   #samps=c(1:length(which(stage_status %in% unlist(strsplit(colnames(contrasts)[r], "-")) == TRUE))))

        if (PNGorPDF == "png") {
          png(paste0(pathNow,"/img/9_DifferentiallyMethylatedRegions_DMR_", 
                     DMR, "_comparison_", colnames(contrasts)[r], ".", PNGorPDF))
        }
        if (PNGorPDF == "pdf") { 
          pdf(paste0(pathNow,"/img/9_DifferentiallyMethylatedRegions_DMR_", 
                     DMR, "_comparison_", colnames(contrasts)[r], ".", PNGorPDF), width = 10, height = 20)
        }
        
        # old code 
        # groups_FL_RLN <- c("1" = "blue", "2" = "green", "3"="maroon")
        # types_FL_RLN <- c(ClinicalFile[which(ClinicalFile[,which(colnames(ClinicalFile) ==ContrastColumnName)] %in% unlist(strsplit(colnames(contrasts)[r], "-")) == TRUE),which(colnames(ClinicalFile) == ContrastColumnName)])
        # cols_FL_RLN <- groups_FL_RLN[types_FL_RLN]
        
        # DMR.plot(ranges = results.ranges[[r]],
        #         dmr = DMR,
        #          CpGs = BetaMatrix[,which(ClinicalFile[,which(colnames(ClinicalFile) == ContrastColumnName)] %in% unlist(strsplit(colnames(contrasts)[r], "-")) == TRUE)],
        #          what = "Beta",
        #          arraytype = "EPIC",
        #          phen.col = cols_FL_RLN,
        #          genome = "hg19",
        #          samps = c(1:length())
        # )
      }
      
      figuregeneratorVolcano <- function(p) {
        
        new_work <- function(TumEnvGeneExp = TumEnvGeneExp) {
          # Annotating points with gene expression data provided on 1 Aug 2019
          TumEnvGeneExp2$gene
          
          colvector <- rep("#e0e0e0", length = length(results.ranges[[p]]$overlapping.promoters))
          for (x in 1:length(results.ranges[[p]]$overlapping.promoters)) {
            example <- unlist(strsplit(as.character(results.ranges[[p]]$overlapping.promoters[x]), "[,]"))
            unique_elements <- unique(trimws(sub("\\-.*", "", example)))
            
            matchingelements <- match(unique_elements, TumEnvGeneExp2$gene)
            # Look in to each element that was matched
            for(n in 1:length(matchingelements)) {
              # Look in to only matched elements that are not NA values
              if(! is.na(matchingelements[n])) {
                # If the matched element is microenvironment, give it one colour
                # Otherwise, give another colour
                if (TumEnvGeneExp2$purity[matchingelements[n]] == "microenvironment") {
                  colvector[x] = 1
                } else if(TumEnvGeneExp2$purity[matchingelements[n]] == "tumour") {
                  colvector[x] = 2
                }
              }
            }
          }
          
          # to see how many genes were not found in gene expression dataset
          length(colvector) - length(which(colvector == "#e0e0e0"))
          # out of total genes 31837, 14911 were not found
    
          # P-value (Stouffer) cutoff on differentially methylated regions (default is 0.05)
          table <- data.frame(FoldChange = results.ranges[[p]]$meanbetafc, Stouffer = -log10(results.ranges[[p]]$Stouffer))
          
          ggplot(table, aes(FoldChange, Stouffer, color = colvector)) +
            geom_point(size = 1) + 
            scale_color_manual(labels = c("none", "high in microenvironment", "high in tumour"), values = c("#e0e0e0", "3", "2")) +
            theme_bw() +
            guides(color = guide_legend("Legend"))+
            ggtitle(paste0("Plot of one-two")) +
            # geom_hline(yintercept = 0, linetype="dotted") +
            # geom_hline(yintercept= -log10(0.01), linetype="dotted") +  
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
        
          RobertsCode <- function() {
            library(dplyr)
            library(ggplot2)
            library(ggpubr)
            
            diff_exprs <- read.table(ExpressionFile, sep = "\t", header = T)
            
            df <- read.csv("Table_Diff_MethylatedRegions_CLUSTER_AllProbes_C1vsC2_2RPMM.csv") %>%
              mutate(GENE = gsub("\\-.*", "", overlapping.promoters)) %>%
              left_join(diff_exprs, by = "GENE") %>%
              mutate(direction = ifelse(adj.P.Val < 0.05 & logFC < 0 & Stouffer < 1e-50, "microenvironment",
                                        ifelse(adj.P.Val < 0.05 & logFC > 0 & Stouffer < 1e-50, "tumor", "not significant"))) %>%
              filter(! is.na(direction))
            
            ##
            
            ggplot(df, aes(x = meanbetafc, y = -log10(Stouffer), color = direction, alpha = 0.4)) +
              geom_point() +
              theme_classic() +
              scale_color_manual(values=c("#E69F00", "gray90", "#56B4E9"))
          }
          
        }
        
        # Print figures only if results.ranges have results stored in this. Help avoid errors. 
        if (! is.na(results.ranges[[p]][1])) {
          # P-value (Stouffer) cutoff on differentially methylated regions (default is 0.05)
          table <- data.frame(FoldChange = results.ranges[[p]]$meandiff, 
                              Stouffer = -log10(results.ranges[[p]]$Stouffer))
          
          # volcano plot
          volcanoplot <- ggplot2::ggplot(table, aes(FoldChange, Stouffer)) +
            geom_point(color = "#e0e0e0") + 
            #coord_cartesian(xlim = c(, 20)) +
            ggtitle(paste0("Plot of ", colnames(contrasts)[p])) +
            geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
            #geom_hline(yintercept= -log10(0.01), linetype="dotted") +  
            theme_bw() + 
            labs(x = "Mean difference", y = "-log10(Stouffer)") + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            scale_color_manual(values = c("#e0e0e0", "#b2182b", "#4d4d4d")) + 
            theme(aspect.ratio = 1, text = element_text(size = 15))
          ggplot2::ggsave(paste0(pathNow, "/img/9_DifferentiallyMethylatedRegions_comparison_", colnames(contrasts)[p], "_VolcanoPlot.", PNGorPDF))
          
          # old plot
          # if (PNGorPDF=="png"){
          #   png(paste0(pathNow,"/img/9_DifferentiallyMethylatedRegions_DMR_",DMR,"_comparison_",colnames(contrasts)[p],"VolcanoPlot.",PNGorPDF))
          # }
          # if (PNGorPDF=="pdf"){ 
          #   pdf(paste0(pathNow,"/img/9_DifferentiallyMethylatedRegions_DMR_",DMR,"_comparison_",colnames(contrasts)[p],"VolcanoPlot.",PNGorPDF), width = 10, height = 20)
          # }
          # plot(results.ranges[[p]]$meanbetafc, -log10(results.ranges[[p]]$Stouffer))
          # dev.off()
        }
      }
      
      if (length(levels(stage_status)) == 2) {
        numbComparisons = 1
      } else if (length(levels(stage_status)) == 3) {
        numbComparisons = 3
      } else if (length(levels(stage_status)) == 4) {
        numbComparisons = 6
      } else if (length(levels(stage_status)) == 5) {
        numbComparisons = 10
      }
      
      for (u in rev(1: numbComparisons)) {
        figuregeneratorVolcano(p = u)
        #figuregenerator_DMR(r = u) # only plot first 30 samples
      }
    }
    
    # GO analysis
    notrun <- function() {
      enrichmentGO <- list()
      for (k in 1:length(colnames(contrasts))) {
        enrichmentGO[[k]] <- missMethyl::goregion(regions = results.ranges[[1]][1:length(results.ranges[[k]])], 
                                                  all.cpg = rownames(myMs.noSNPs),
                                                  collection = "GO", 
                                                  array.type = "EPIC")
        enrichmentGO[[k]] <- enrichmentGO[[k]][order(enrichmentGO[[k]]$P.DE),]
        dim(enrichmentGO[[k]]) # 22750     6
        head(enrichmentGO[[k]], 10)
        

        
        
        table(enrichmentGO[[k]]$TERM[enrichmentGO[[k]]$P.DE <0.01])
        
        
        library(GOplot)
        circ <- GOplot::circle_dat(EC$david, EC$genelist)
        
        
        library(enrichplot)
        library(DOSE)
        data(geneList)
        de <- names(geneList)[abs(geneList) > 2]
        
        edo <- enrichDGN(de)
        barplot(edo, showCategory=20)
        
        
        library(DOSE)
        data(geneList)
        de <- names(geneList)[1:100]
        x <- enrichDO(de)
        barplot(x)
        
      }
    }
    
    RESULTS <- list(AnnotatedCpGs = myannotations,
                    DifferentiallyMethylatedRegions = dmrcoutput,
                    GRangesObject = results.ranges,
                    GenesFallingInRegions = geneVector)
    class(RESULTS) <- "DifferentiallyMethylatedRegions_DMRcate_ASilva"
    return(RESULTS)
  } 
  
  
  if (Method == "bumphunter") {
    RESULTSFinal <- bumphunter_function(BetaMatrix = BetaMatrix, 
                                        ContrastColumnName = ContrastColumnName, 
                                        ClinicalFile = ClinicalFile, 
                                        AnnotationFile = AnnotationFile, 
                                        ProduceImages = ProduceImages, 
                                        PNGorPDF = PNGorPDF)
  } else if (Method == "DMRcate") {
    RESULTSFinal <- DMRcate_function(BetaMatrix = BetaMatrix, 
                                     ContrastColumnName = ContrastColumnName, 
                                     ClinicalFile = ClinicalFile, 
                                     AnnotationFile = AnnotationFile, 
                                     ProduceImages = ProduceImages, 
                                     PNGorPDF = PNGorPDF)
  }

  class(RESULTSFinal) <- "DifferentiallyMethylatedRegions_ASilva"
  return(RESULTSFinal)
}
# [END]
