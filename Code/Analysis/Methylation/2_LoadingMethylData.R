# Updated 5 Feb 2020
# Updated 8 Nov 2019
# Updated 7 Feb 2019
# Function: Load methylation IDAT files using sheet from 1_QCRemoveSamples.R to carry out 
#           normalization, then generate Mvalues and Beta values using mainly functions from minfi
#           and missMethyl packages.
# based on https://bioconductor.org/packages/release/bioc/vignettes/missMethyl/inst/doc/missMethyl.html
# and https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html
# other links: https://support.bioconductor.org/p/77176/
# https://bioconductor.org/packages/release/data/experiment/vignettes/BeadArrayUseCases/inst/doc/BeadArrayUseCases.pdf
# Author: Anjali Silva and Robert Kridel

# Input:
# Sheet: A data frame generated as the output from 1st step: 1_QCRemoveSamples(). 
#        It should be of size observations x 8 colnames. Observations may include patient IDs.
#        Colnames may include "Sample_Name" Sample_Well" "Sample_Plate" "Sample_Group" 
#        "Sample_Type" "Array" "Slide" "Basename". E.g., Sheet_updSamples_170_Ordered_T1. 
# ClinicalFile: A data frame with observations (e.g., patient IDs) along rows and categories 
#               (e.g., SEX SITE_BIOPSY, etc.) on columns. The number of observations must match 
#               the number of observations in 'Sheet' argument. 
# AnnotationFile: A data frame (data table) of annotations for all the probes found in BetaMatrix.
#                 It is of size probes x annotations.
# FigureGenerate: Should figures be generated, default = "No". Options: "Yes", "No"
# TableGenerate: Should tables be generated, default = "No". Options: "Yes", "No"
# PNGorPDF: Should figures be generated in PNG or PDF format. Options: "png", "pdf"; 
#           Default is "png"

# Output:
# PhenoData: Phenotypic data found within. 
# RGChannelSet: RGChannelSet object of raw data from the IDAT files.
# TotalProbes: Number of total probes from original reading of IDAT files.
# FilteredProbes: Number of total probes resulting after filtering out poor quality 
#                 probes using detectionP() function.
# MethylSetRaw: MethylSet without normalization.
# MethylSetWithNorm: MethylSet with Subset-quantile Within Array Normalisation.
# BetaMatrix: Matrix of Beta values of size probes x updated samples (refer 
#             MissingSamplesClinicalFile).
# MvalueMatrix: Matrix of M values of size probes x updated samples (refer 
#               MissingSamplesClinicalFile).
# PredictedSex: A character vector of "M" and "F" values indicating predicted sex 
#               for samples from minfi::minfiQC.


LoadingMethylData <- function(Sheet, 
                              ClinicalFile, 
                              AnnotationFile,
                              FigureGenerate = "No", 
                              TableGenerate = "No", 
                              PNGorPDF = "png") {
  
  
  ########################################## #
  # Loading needed packages
  ########################################## #
  library(minfi)
  library(missMethyl)
  library(Biobase)
  library(shinyMethyl)
  library(data.table)
  library(dplyr)
  library(IlluminaHumanMethylationEPICmanifest)
  
  # Getting the path to save files, which should contain a folder called "img"
  pathNow <- getwd()
  
  
  # Defining the path to reading the an Illumina methylation sample sheet
  # sheet <- read.metharray.sheet(base =file.path(SampleSheet_Path, "extdata"))
  RGset <- minfi::read.metharray.exp(targets = as.data.frame(Sheet), force = TRUE)
  # class(RGset) # "RGChannelSet"
  # RGChannelSet: raw data from the IDAT files; this data is organized at the probe 
  # (not CpG locus) level. This data has two channels: Red and Green.
  # dim(RGset) # 1051539     170
  colnames(RGset) <- Sheet$Sample_Name
  
  # Associated pheno data, i.e., information contained in the targets object from RGset
  pd <- Biobase::pData(RGset)
  
  # Creation of a MethylSet without normalization
  # Converts the Red/Green channel for an Illumina methylation array into methylation 
  # signal, without using any normalization.
  mSet <- minfi::preprocessRaw(RGset)
  # class(mSet) # MethylSet
  # mSet is MethylSet objects contains only the methylated and unmethylated signals. 
  # MethylSet : data organized by the CpG locus level, but not mapped to a genome.
  # This data has two channels: Meth (methylated) and Unmeth (unmethylated).
  # dim(mSet) # 865859    170
  
  GMSet <- minfi::mapToGenome(mSet)
  # class(GMSet) # "GenomicMethylSet"
  # GenomicMethylSet : like a MethylSet, but mapped to a genome.
  # Could in principle be smaller in case the annotation 
  # you use says that some probes do not map to the genome
  # dim(GMSet) # 865859    170
  
  # To see what type of annotation packages is used by array
  AnnotationPackage <- BiocGenerics::annotation(RGset)
  #                          array                     annotation 
  # "IlluminaHumanMethylationEPIC"                 "ilm10b4.hg19" 
  
  
  
  ########################################## #
  # Quality control
  ########################################## #
  QCanalysis <- minfi::minfiQC(object = GMSet, fixOutliers = TRUE, verbose = FALSE)
  # names(QCanalysis) # "object" "qc"  
  # For getQC, a DataFrame with two columns: mMed and uMed which are the chipwide medians 
  # of the Meth and Unmeth channels.
  # QCanalysis$qc$predictedSex
  
  estSex <- minfi::getSex(GMSet)
  GMSetNew <- minfi::addSex(GMSet, sex = estSex)
  # densityPlot(mSet, sampGroups = ClinicalFile$TYPE)
  
  # Plotting QC plots
  if (FigureGenerate == "Yes") {
    if (PNGorPDF == "pdf") {
      pdf(paste0(pathNow, "/img/1_QCplots.pdf"), width = 60, height = 60, pointsize = 50)
    } else {
      png(paste0(pathNow,"/img/1_QCplots.png"))
    }
    par(mfrow = c(1,2))
    # which uses the log median intensity in both the methylated (M) and unmethylated (U) 
    # channels
    plotQC(QCanalysis$qc, badSampleCutoff = 10.5)
    
    performing_internal_analysis <- function() {
      # Performing internal QC analysis
      
      # First with 170 samples (then with 177 samples)
      # plot(QCanalysis$qc$mMed,QCanalysis$qc$uMed, xlim = c(8,13), ylim=c(8,13))
      # QCvalues <- cbind(QCanalysis$qc$mMed, QCanalysis$qc$uMed)
      # BadSamples <- rowSums(QCvalues < 10.5)  == 2
      # rownames(QCanalysis$qc)[BadSamples]
      # bad samples when 170 were used = "LY_FL_020_T1" "LY_FL_046_T1" "LY_FL_065_T1" "LY_FL_149_T1" "LY_FL_159_T1" "LY_FL_164_T1"
      # [7] "LY_FL_165_T1" "LY_FL_168_T1" "LY_FL_183_T1" "LY_FL_218_T1" "LY_FL_287_T1" "LY_FL_297_T1"
      # [13] "LY_FL_465_T1" "LY_RLN_002" 
      # ClinicalFile_T1$SAMPLE_ID[BadSamples] # double check if IDs are the same, yes
      # ClinicalFile_T1$INSTITUTION[BadSamples]
      # ClinicalFile_T1$TYPE[BadSamples]
      # ClinicalFile_T1$EPIC_QC[BadSamples]
      # ClinicalFile_T1$STAGE[BadSamples]
      # ClinicalFile_T1$TRANSLOC_14_18[BadSamples]
      
      # bad samples when 177 were used = 
      # [1] "LY_FL_020_T1" "LY_FL_046_T1" "LY_FL_065_T1" "LY_FL_159_T1" "LY_FL_164_T1" "LY_FL_165_T1"
      # [7] "LY_FL_183_T1" "LY_FL_218_T1" "LY_FL_465_T1" "LY_RLN_002"   "LY_FL_066_T1" "LY_FL_149_T1"
      # [13] "LY_FL_168_T1" "LY_FL_287_T1" "LY_FL_297_T1"
      Match_Bad177 <- match(rownames(QCanalysis$qc)[BadSamples], ClinicalFile$SAMPLE_ID)
      ClinicalFile$SAMPLE_ID[Match_Bad177] # samples match
      
      table(ClinicalFile$SEX[Match_Bad177] )
      table(ClinicalFile$TYPE[Match_Bad177])
      table(ClinicalFile$INSTITUTION[Match_Bad177])
      table(ClinicalFile$STAGE[Match_Bad177])
      table(ClinicalFile$TRANSLOC_14_18[Match_Bad177])
      which(ClinicalFile$EPIC_QC == "Bad") #108
      which(ClinicalFile$EPIC_QC == "Poor") #108
      # sex analysis
      Match_clinicalSamples177 <- match(ClinicalFile$SAMPLE_ID, colnames(GMSetNew))
      table(GMSetNew[, Match_clinicalSamples177]$predictedSex, ClinicalFile$SEX)
      # comparing which samples differ
      SamplesDifferSex <- which(GMSetNew[, Match_clinicalSamples177]$predictedSex != ClinicalFile$SEX)
      ClinicalFile$SEX[SamplesDifferSex]
      ClinicalFile$SAMPLE_ID[SamplesDifferSex]
      table(ClinicalFile$EPIC_QC[SamplesDifferSex])
      # "LY_FL_165_T1" "LY_FL_080_T1" "LY_FL_149_T1" "LY_FL_038_T1" "LY_FL_163_T1"
      #  "LY_FL_168_T1" "LY_FL_159_T1"
    
      # see samples shared between minfi 'bad' samples and mismatched sex samples
      Reduce(intersect, list(ClinicalFile$SAMPLE_ID[Match_Bad177], 
                             ClinicalFile$SAMPLE_ID[SamplesDifferSex]))
      # "LY_FL_159_T1" "LY_FL_165_T1" "LY_FL_149_T1" "LY_FL_168_T1"
      
      SamplesDifferSex_mSet <- match(ClinicalFile$SAMPLE_ID[SamplesDifferSex], colnames(mSet))
      
      par(mfrow = c(2,2))
      densityPlot(mSet[, BadSamples], main = "Raw, 15 bad samples")
      densityPlot(mSet[, SamplesDifferSex_mSet], main = "Raw, 7 samples with sex difference")
      densityPlot(mSet, main = "Raw, all samples")
    }
    # 
    GMSetNew <- minfi::addSex(GMSet, sex = estSex)
    plotSex(GMSetNew, id = NULL)
    dev.off()
  }
  
  # mdsPlot Multi-dimensional scaling plots giving an overview of similarities and
  # differences between samples.
  # Multi-dimensional scaling (MDS) plots showing a 2-d projection of distances between samples
  # mdsPlot(dat = RGset, 
  #        numPositions = 2000,
  #        main = "2D projection of distances between samples using 2000 most variable positions")
  
  # Plot the overall density distribution of beta values and the density distributions of the Infinium I and II probe types.
  # plotBetasByType(data = mSet[,1])
  
  # Plot control probe signals
  # controlStripPlot(rgSet = RGset)
  
  # Density bean plots of methylation Beta values.
  # densityBeanPlot(dat = mSet, sampGroups = NULL, sampNames = NULL, main = NULL,
  #                 pal = brewer.pal(8, "Dark2"), numPositions = 10000)
  
  # QC report for Illumina Infinium Human Methylation 450k arrays
  # qcReport(rgSet = RGset, sampNames = NULL, sampGroups = NULL, pdf = "qcReport.pdf",
  #         maxSamplesPerPage = 24, controls = c("BISULFITE CONVERSION I",
  #                                              "BISULFITE CONVERSION II", "EXTENSION", "HYBRIDIZATION",
  #                                              "NON-POLYMORPHIC", "SPECIFICITY I", "SPECIFICITY II", "TARGET REMOVAL"))
  
  ########################################## #
  # Runninc shinyMethyl
  ########################################## #
  # myShinyMethylSet <- shinyMethyl::shinySummarize(RGset)
  # runShinyMethyl(myShinyMethylSet)
  
  
  
  ########################################## #
  # Normalization
  ########################################## #
  cat("\n Subset-quantile within array normalization \n")
  # Subset-quantile within array normalization based on missMethyl
  # SWAN corrects for the technical differences between the Infinium I and II assay designs and produces a smoother overall β value distribution
  mSetSw <- missMethyl::SWAN(mSet, verbose = TRUE)
  # class(mSetSw) # MethylSet object
  # MethylSet : data organized by the CpG locus level, but not mapped to a genome.
  # This data has two channels: Meth (methylated) and Unmeth (unmethylated).
  
  # Do QC analysis based on SWAN normalization based on missMethyl
  # GMSet2 <- minfi::mapToGenome(mSetSw)
  # QCanalysis2 <- minfi::minfiQC(object = GMSet2, fixOutliers = FALSE, verbose = FALSE)
  # plotQC(QCanalysis2$qc)
  # Error in kmeans(dd, centers = c(min(dd), max(dd))) : 
  # initial centers are not distinct
  
  # Subset-quantile within array normalization based on minfi
  # datSwan <- minfi::preprocessSWAN(RGset, mSet = mSet)
  # class(datSwan) # MethylSet object
  # dim(datSwan)
  
  # Plotting density distributions
  if (FigureGenerate == "Yes") {
    if (PNGorPDF == "pdf") {
      pdf(paste0(pathNow, "/img/1_DensityPlots_RawVsSWAN_LY_DLC_001.pdf"), 
          width = 60, height = 60, pointsize = 50)
    } else {
      png(paste0(pathNow,"/img/1_DensityPlots_RawVsSWAN_LY_DLC_001.png"))
    }
    par(mfrow = c(1,2))
    # Sample LY_DLC_001
    densityByProbeType(mSet[, 1], main = "Raw, Sample LY_DLC_001")
    densityByProbeType(mSetSw[, 1], main = "missMethyl SWAN")
    #densityByProbeType(mSet[, 1], main = "Raw, Sample LY_DLC_001")
    #densityByProbeType(datSwan[, 1], main = "minfi SWAN")
    dev.off()
    
    if (PNGorPDF == "pdf") {
      pdf(paste0(pathNow, "/img/1_DensityPlots_RawVsSWAN_LY_FL_008_T1.pdf"), 
          width = 60, height = 60, pointsize = 50)
    } else {
      png(paste0(pathNow,"/img/1_DensityPlots_RawVsSWAN_LY_FL_008_T1.png"))
    }
    par(mfrow = c(1,2))
    # Sample LY_FL_008_T1
    densityByProbeType(mSet[, 13], main = "Raw, Sample LY_FL_008_T1")
    densityByProbeType(mSetSw[, 13], main = "missMethyl SWAN")
    #densityByProbeType(mSet[, 13], main = "Raw, Sample LY_FL_008_T1")
    #densityByProbeType(datSwan[, 13], main = "minfi SWAN")
    dev.off()
    
    if (PNGorPDF == "pdf") {
      pdf(paste0(pathNow, "/img/1_DensityPlots_RawVsSWAN_LY_RLN_001.pdf"), 
          width = 60, height = 60, pointsize = 50)
    } else {
      png(paste0(pathNow,"/img/1_DensityPlots_RawVsSWAN_LY_RLN_001.png"))
    }
    par(mfrow = c(1,2))
    # Sample LY_RLN_001
    densityByProbeType(mSet[, 166], main = "Raw, Sample LY_RLN_001")
    densityByProbeType(mSetSw[, 166], main = "missMethyl SWAN")
    #densityByProbeType(mSet[, 166], main = "Raw, Sample LY_RLN_001")
    #densityByProbeType(datSwan[, 166], main = "minfi SWAN")
    dev.off()
  }
  
  
  ########################################## #
  # Filter out poor quality probes 
  ########################################## #
  
  # Poor quality probes can be filtered out based on the detection p-value less than 0.01
  detP <- minfi::detectionP(RGset)
  keep <- base::rowSums(detP < 0.01) == ncol(RGset)
  mSetSw_new <- mSetSw[keep, ]
  # dim(mSetSw_new) # 628221    177
  
  ########################################## #
  # Filtering probes affected by SNPs
  ########################################## #
  
  # https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#methylset-and-ratioset
  # Probe, CpG and SBE correspond the SNPs present inside the probe body, at
  # the CpG interrogation and at the single nucleotide extension respectively.
  # The columns with rs give the names of the SNPs while the columns with maf 
  # gives the minor allele frequency of the SNPs based on the dbSnp database. 
  
  # GRset_SNP <- addSnpInfo(RGset)
  # Error in .isGenomicOrStop(object) : 
  # object is of class 'RGChannelSet', but needs to be of class 'GenomicMethylSet' or 
  # 'GenomicRatioSet'
  # GRset_SNP <- dropLociWithSnps(GRset_SNP, snps=c("SBE","CpG"), maf=0)
  
  # We strongly recommend to drop the probes that contain either a SNP at the CpG 
  # interrogation (CpG_rs) or at the single nucleotide extension (SBE_rs),for any minor 
  # allele frequency
  probes.SNPs <- minfi::getSnpInfo(mSetSw_new) 
  CpG_rs <- which(! is.na(probes.SNPs$CpG_rs) == "TRUE")
  length(CpG_rs) # 19,671
  SBE_rs <- which(! is.na(probes.SNPs$SBE_rs) == "TRUE")
  length(SBE_rs) # 10,554
  # Probe_rs <- which(! is.na(probes.SNPs$Probe_rs) == "TRUE")
  #length(Probe_rs) # 109,785
  # Testing <- dropLociWithSnps(GMSet, snps=c("SBE","CpG"), maf=0) # 30,435 probes
  
  # Find common SNPs between CpG_rs, SBE_rs
  # CommonSNPs <- Reduce(intersect, list(CpG_rs, SBE_rs))
  AllSNPs <- union(CpG_rs, SBE_rs) # 20,073
  # length(CommonSNPs) # 10152
  mSetSw_new <- mSetSw_new[- AllSNPs, ]
  # dim(mSetSw_new) # 618069    170
  
  # method shown in minfi tutorial 
  # https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#qc-plot
  # TestSNPs <- minfi::getSnpInfo(GMSet) 
  # Only provides SNPs present inside the probe body, nothing for 
  # SNPs at the CpG interrogation and at the single nucleotide extension
  # GRset <-  minfi::addSnpInfo(GMSet)
  # GRset <-  minfi::dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)
  
  
  
  ########################################## #
  # Filtering probes corresponding sex chromosomes (X and Y)
  ########################################## #
  
  # Obtaining same probes as that remain in BetaMatrix after removal of probes based on other criteria
  obtain_probes_annotationFile_order <- match(rownames(mSetSw_new), AnnotationFile$V1)
  # AnnotationFile$V1[obtain_probes_annotationFile_order] # checked order is correct
  # See which of these probes are from sex chromosomes
  remove_sex_csome_X <- which(AnnotationFile$chr[obtain_probes_annotationFile_order] == "chrX")
  # length(remove_sex_csome_X) # 12541
  remove_sex_csome_Y <- which(AnnotationFile$chr[obtain_probes_annotationFile_order] == "chrY")
  # length(remove_sex_csome_Y) #  117
  mSetSw_new <- mSetSw_new[- c(remove_sex_csome_X, remove_sex_csome_Y), ]
  # dim(mSetSw_new) # 605411    170
  
  
  ########################################## #
  # Extracting Beta and M-values
  ########################################## #
  meth <- minfi::getMeth(mSetSw_new)
  unmeth <- minfi::getUnmeth(mSetSw_new)
  Mval <- base::log2((meth + 100) / (unmeth + 100))
  # dim(Mval) # 605411    170
  
  # M-values; M = logit(Beta) = log( Meth / Unmeth )
  # Mval <- minfi::getM(mSetSw_new)
  
  # Beta values; Beta = Meth / (Meth + Unmeth + offset)
  beta <- minfi::getBeta(mSetSw_new)
  # densityPlot(mSetSw_new, sampGroups = ClinicalFile$TYPE)
  # dim(beta)  # 605411    170
  
  ########################################## #
  # Extracting copy number values
  ########################################## #
  # get copy number values which are defined as the sum of the methylation and 
  # unmethylation channel.
  # copynumber <- minfi::getCN(mSetSw_new)
  # dim(copynumber) #  605411    170
  
  
  ########################################## #
  # Remove probes that show value of 0 across all samples
  ########################################## #
  ZeroRows_beta <- which(rowSums(beta) == 0) # no probes were present with rowSums 0
  if(length(ZeroRows_beta) >0 ) {
    # Remove these probes from Beta matrix
    beta <- beta[- ZeroRows_beta, ]
  }
  # dim(beta) # 605411    170
  
  # cat("\n Writing Results")
  RESULTS <- list(PhenoData = pd,
                  ArrayAnnotation = AnnotationPackage[2],
                  RGChannelSet = RGset,
                  GenomicMethylSet = GMSetNew, 
                  TotalProbes = nrow(detP),
                  FilteredProbes = nrow(detP) - nrow(mSetSw_new),
                  MethylSetRaw = mSet,
                  MethylSetWithNorm = mSetSw_new,
                  BetaMatrix = beta,
                  MvalueMatrix = Mval,
                  PredictedSex = QCanalysis$qc$predictedSex)
  
  class(RESULTS) <- "LoadingMethylData_ASilva"
  return(RESULTS)
}


