# Updated 3 Aug 2021
# Created 18 May 2021
# Description:
# The scriptscRNAseqExpressionset46 is a script to generate single cell data expression set.
# The scriptBulkRNAseqExpressionset46 is a script 
# The scBulkRNAseqBisqueDecomposition46 is a function for reference-based decomposition 
# The plotBisqueData is a function for plotting bisque output

# Adapted from: https://cran.r-project.org/web/packages/BisqueRNA/vignettes/bisque.html


# Extract count matrix from seurat and generate scRNAseq expression set - not a function, a script
scriptscRNAseqExpressionset46 <- function() {
  
  library("Biobase")
  library("BisqueRNA")
  # install.packages("Seurat") 
  library("Seurat")
  
  # single-cell data expression set
  # Karin's output from seurat is located here:
  # output <- "/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April82021/"
  # With respect to Anjali's computer
  output <- "/Users/anjalisilva/UHN/kridel-lab - FLOMICS/Analysis-Files/Seurat/April82021/"
  #sc.eset <- readRDS(paste(output, "sc_eset_for_bisque.rds", sep=""))

# Adapted from Karin's code
# https://github.com/kridel-lab/FLOMICS/blob/master/Code/BioinformaticsProcessing/snRNAseq/004_prepare_seurat_object_for_bisque.R
  scCountMatrix <- readRDS(paste0(output,"pc_genes_only_no_seurat_integrated_dim_20_2000_2021-04-08_samples_clusters.rds"))
  get.cell.names <- function(obj) base::colnames(obj)
  get.ident <- function(obj) Seurat::Idents(object = obj)
  get.raw.data <- function(obj) Seurat::GetAssayData(object = obj[["RNA"]], slot = "counts")
  delimiter <- "_"
  position <- 2
  individual.ids <- base::sapply(base::strsplit(get.cell.names(scCountMatrix),
                                                delimiter), `[[`, position)
  
  base::names(individual.ids) <- get.cell.names(scCountMatrix)
  individual.ids <- base::factor(individual.ids)
  n.individuals <- base::length(base::levels(individual.ids))
  base::message(base::sprintf("Split sample names by \"%s\"", delimiter),
                base::sprintf(" and checked position %i.", position),
                base::sprintf(" Found %i individuals.", n.individuals))
  base::message(base::sprintf("Example: \"%s\" corresponds to individual \"%s\".",
                              get.cell.names(scCountMatrix)[1], individual.ids[1]))
  sample.ids <- base::names(get.ident(scCountMatrix))
  sc.pheno <- base::data.frame(check.names = F, 
                               check.rows = F,
                               stringsAsFactors = F,
                               row.names = sample.ids,
                               SubjectName = individual.ids,
                               cellType = get.ident(scCountMatrix))
  sc.meta <- base::data.frame(labelDescription = base::c("SubjectName",
                                                         "cellType"),
                              row.names = base::c("SubjectName",
                                                  "cellType"))
  sc.pdata <- methods::new("AnnotatedDataFrame",
                           data = sc.pheno,
                           varMetadata = sc.meta)
  
  sc.data <- base::as.matrix(get.raw.data(scCountMatrix)[, sample.ids, drop = F])
  sc.eset <- Biobase::ExpressionSet(assayData = sc.data,
                                    phenoData = sc.pdata)
  # sc.eset$SubjectName # 1,2,3 is respectively "LY_FL_062_T1", "LY_FL_064_T1", "LY_FL_076_T1"
  # length(sc.eset$SubjectName) # 18676
  # assayData: 29354 features, 18676 samples 
  
  RESULTS <- list(sc.eset = sc.eset)
  class(RESULTS) <- "scriptscRNAseqExpressionset46"
  return(RESULTS)
}
scOutput <- scriptscRNAseqExpressionset46()



# Extract bulk RNAseq expression set - not a function, a script
scriptBulkRNAseqExpressionset46 <- function(RNAseqProtCode) {
# Bulk RNA-seq data can be converted from a matrix (columns are samples, rows are genes)
# to an ExpressionSet as follows:
# Karin's code, if needed: https://github.com/kridel-lab/FLOMICS/blob/master/Code/Analysis/snRNAseq/bisque_001_sarah_final_code.R

  scRNAseqSamples <- c("LY_FL_062_T1", "LY_FL_064_T1", "LY_FL_076_T1")
  matchBulkRNAseq <- match(scRNAseqSamples, colnames(RNAseqProtCode))
  RNAseqProtCode3Samples <- RNAseqProtCode
  colnames(RNAseqProtCode3Samples)[matchBulkRNAseq] <- c("1", "2", "3")
  bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(RNAseqProtCode3Samples))
  
  RESULTS <- list(bulk.eset = bulk.eset)
  class(RESULTS) <- "scriptBulkRNAseqExpressionset46"
  return(RESULTS)
}
rnaSeqOutput <- scriptBulkRNAseqExpressionset46(RNAseqProtCode = RNAseqProtCode)




# Reference-based decomposition - function
# Input:
# bulk.eset: Bulk RNAseq output from function Biobase::ExpressionSet
# sc.eset: scRNAseq output from function Biobase::ExpressionSet
scBulkRNAseqBisqueDecomposition46 <- function(bulk.eset,
                                              sc.eset) {
  library("Biobase")
  library("BisqueRNA")
  library("Seurat")
  library("grDevices")

  
  # Reference-based decomposition
  refDecomp <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = bulk.eset, 
                                                      sc.eset = sc.eset, 
                                                      markers = NULL, 
                                                      use.overlap = TRUE)
  # Decomposing into 16 cell types.
  # Found 3 samples with bulk and single-cell expression.
  # Remaining 98 bulk samples will be decomposed.
  # Using 15149 genes in both bulk and single-cell expression. (OLD)
  # Using 25947 genes in both bulk and single-cell expression.
  
  # Converting single-cell counts to CPM and filtering zero variance genes.
  # Filtered 1 zero variance genes.
  # Converting bulk counts to CPM and filtering unexpressed genes.
  # Filtered 0 unexpressed genes.
  # Generating single-cell based reference from 18676 cells.
  
  # Learning bulk transformation from overlapping samples.
  # Applying transformation to bulk samples and decomposing.
  # Dropped an additional 468 genes for which a transformation could not be learned.
  
  # names(refDecomp)
  # "bulk.props"       "sc.props"         "rnorm"            "genes.used"      
  # [5] "transformed.bulk"
  
  # A list. Slot bulk.props contains a matrix of cell type proportion estimates with
  # cell types as rows and individuals as columns. Slot sc.props contains a matrix
  # of cell type proportions estimated directly from counting single-cell data. Slot 
  # rnorm contains Euclidean norm of the residuals for each individual's proportion 
  # estimates. Slot genes.used contains vector of genes used in decomposition. Slot 
  # transformed.bulk contains the transformed bulk expression used for decomposition. 
  
  #scRNAseqSamples <- c("LY_FL_062_T1", "LY_FL_064_T1", "LY_FL_076_T1")
  #matchBulkRNAseq <- match(scRNAseqSamples, colnames(refDecomp$bulk.props)) # not present
  
  #refDecomp$bulk.props
  #summary(refDecomp$bulk.props)
  
  #refDecomp$sc.props
  #summary(refDecomp$sc.props)
  
  # A list is returned with decomposition estimates in slot bulk.props.
  # ref.based.estimates <- refDecomp$bulk.props
  # knitr::kable(ref.based.estimates, digits = 2)
  

  RESULTS <- list(EntireOutput = refDecomp,
                  RefEstimatesBulkProp = refDecomp$bulk.props,
                  RefEstimatesBulkscRNAseqProp = cbind(refDecomp$bulk.props, refDecomp$sc.props))
  class(RESULTS) <- "scBulkRNAseqBisqueDecomposition46"
  return(RESULTS)
}

BisqueDecompOutput <- scBulkRNAseqBisqueDecomposition46(bulk.eset = rnaSeqOutput$bulk.eset,
                                  sc.eset = scOutput$sc.eset)





# plotBisqueData is a script
# Input:
# RefEstimates
# ClinicalFile: File with patient sample names, and categories. The patient order will be made to be that 
#               of RefEstimates. Should include column "CLUSTER".

RefEstimates = BisqueDecompOutput$RefEstimatesBulkscRNAseqProp
plotBisqueData <- function(RefEstimates,
                           ClinicalFile,
                           PNGorPDF) {
  
  library(reshape2)
  
  allCellTypes <- rownames(RefEstimates) # Total number of cells
  # Assign cells 
  Bcells <- c(0, 1, 11, 12, 2, 4)
  Tcells <- c(15, 5, 6, 7, 8, 9)
  Othercells <- as.vector(setdiff(allCellTypes, c(Bcells, Tcells)), mode = "numeric")
  
  # determine cell
  cellOrganizer <- function(ColumnToFill, Variable) {
    ColumnToFill[which(Variable %in% Bcells)] <- "Bcells"
    ColumnToFill[which(Variable %in% Tcells)] <- "Tcells"
    ColumnToFill[which(Variable %in% Othercells)] <- "Other"
    ColumnToFill <- factor(ColumnToFill, levels = c("Bcells", "Tcells", "Other"))
    return(ColumnToFill)
  }
  
  
  # Determine color given category to visualize
  colorFunction <- function(CategoryToVisualize) {
    # Change violin plot colors by groups
    if(CategoryToVisualize == "TYPE") {
      ColourChoice <- c("#d6604d", "#5aae61", "#2166ac")
    } else if(CategoryToVisualize == "STAGE") {
      ColourChoice <- c("#762a83","#c2a5cf")
    } else if(CategoryToVisualize == "SEX") {
      ColourChoice <- c("#b35806","#fdb863")
    } else if(CategoryToVisualize == "SITE_BIOPSY") {
      ColourChoice <- c("#f1b6da","#c51b7d")
    } else if(CategoryToVisualize == "TYPE_BIOPSY") {
      ColourChoice <- c("#a6dba0", "#878787")
    } else if(CategoryToVisualize == "TRANSLOC_14_18") {
      ColourChoice <- c("#e0e0e0","#878787")
    } else {
      ColourChoice <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                        '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080')
    }
    return(ColourChoice)
  }
  
  # Determine number of pairwise comparisons
  comparisonFunction <- function(numberComparisons) {  
    # Define comparisons
    if(numberComparisons == 2) {
      comparisonOptions <- list(as.character(utils::combn(1:2, 2)))
    } else if(numberComparisons >= 3) {
        comparisonOptions <- list() # empty list to save data
        comparisonsPresent <- utils::combn(0:numberComparisons, 2)
        for(i in 1:ncol(comparisonsPresent)) {            
          comparisonOptions[[i]] <- as.character(comparisonsPresent[ , i])
        }
    }
    return(comparisonOptions)
  }
  
  
  # Convert RefEstimates to data frame
  tRefEstimates <- data.frame(SAMPLE_ID = rownames(t(RefEstimates)),
                              t(RefEstimates))
  
  #tRefEstimates$SAMPLE_ID[99:101] <- c("LY_FL_062_T1", "LY_FL_064_T1", "LY_FL_076_T1")
  
  # Combine with Clinical file data 
  tRefEstimatesClinicalFile <- merge(tRefEstimates, 
                                     ClinicalFile, by = "SAMPLE_ID")
  

  
  # STAGE
  dataSTAGE <- data.frame(reshape2::melt(tRefEstimatesClinicalFile[, 
                          c("SAMPLE_ID", paste0("X", allCellTypes), "STAGE")],
                          id = c("STAGE", "SAMPLE_ID", "SAMPLE_ID")))
  dataSTAGE$STAGE <- factor(dataSTAGE$STAGE)
  dataSTAGE$variable <- factor(substr(dataSTAGE$variable, 2, 3))
  dataSTAGE$variable <- factor(dataSTAGE$variable, levels = c(0, 1:(length(allCellTypes) - 1)))
  dataSTAGE$CellType <- factor(cellOrganizer(ColumnToFill = dataSTAGE$SAMPLE_ID, 
                                               Variable = dataSTAGE$variable))
  dataSTAGE$value <- as.numeric(dataSTAGE$value)
  
  
  pValues <- dataSTAGE %>% group_by(variable) %>%
    dplyr::summarize(pval = wilcox.test(value ~ STAGE)$p.value) %>%
    data.frame() %>% dplyr::mutate(p.adjst = stats::p.adjust(as.numeric(pval), method = "bonferroni")) 
  dataSTAGE <- merge(dataSTAGE, pValues, by = "variable")

  STAGEPlot <- ggpubr::ggboxplot(dataSTAGE, 
                                   x = "variable", 
                                   y = "value", 
                                   width = 0.8,
                                   size = 1,
                                   color = "STAGE",
                                   ylab = "Fraction of cells", 
                                   xlab = "Cluster",
                                   font.label = list(size = 20, color = "black"), 
                                   add = "jitter",
                                   palette = colorFunction("STAGE")[1:2]) +
                                   ggtitle(paste0("CLUSTER - 10 Feb 2021 labels (98 deconvoluted; no sex X and Y)")) +
                                   facet_wrap(~ CellType, scales = "free") +
                                   geom_text(aes(label = paste0("P=", round(p.adjst, 3)), y = 0.65), size = 2)
  pathNow <- getwd()
  ggplot2::ggsave(paste0(pathNow, "/img/46_scRNAseqBisque_STAGE_", ImageName, ".", PNGorPDF),
                  width = 3.5, height = 3.5)
  
  
  
  
  
  # Cluster
  dataCLUSTER <- data.frame(reshape2::melt(tRefEstimatesClinicalFile[, 
                            c("SAMPLE_ID", paste0("X", allCellTypes), "CLUSTER")],
                            id = c("CLUSTER", "SAMPLE_ID")))
  dataCLUSTER$CLUSTER <- factor(dataCLUSTER$CLUSTER)
  dataCLUSTER$variable <- factor(substr(dataCLUSTER$variable, 2, 3))
  dataCLUSTER$variable <- factor(dataCLUSTER$variable, levels = c(0, 1:length(allCellTypes) - 1))
  dataCLUSTER$CellType <- factor(cellOrganizer(ColumnToFill = dataCLUSTER$SAMPLE_ID, 
                                        Variable = dataCLUSTER$variable))
  dataCLUSTER$value <- as.numeric(dataCLUSTER$value)
  
  # identical(unique(snf_analysis$SAMPLE_ID), unique(dataCLUSTER$SAMPLE_ID))
  
  # Comparing input values of Karin vs Anjali
  # par(mfrow = c(1, 2))
  # AnjaliVals <- dataCLUSTER %>% group_by(variable) %>% filter(variable == 1) %>% 
  # data.frame() %>% dplyr::select(value) 
  # plot(x = c(1:98), y = AnjaliVals$value,  
  #      ylim = c(0, 0.6),
  #      xlab= "Patient", ylab = "Cell Proportion", 
  #      main = "Anjali for Cell Population 1")
  
  # KarinVals <- snf_analysis %>% group_by(cell_type) %>% filter(cell_type == 1) %>% 
  # data.frame() %>% dplyr::select(value) 
  # plot(x = c(1:98), y = KarinVals$value,
  #      ylim = c(0, 0.6),
  #      xlab = "Patient", ylab = "Cell Proportion", 
  #      main = "Karin for Cell Population 1")

  pValues <- dataCLUSTER %>% group_by(variable) %>%
    dplyr::summarize(pval = wilcox.test(value ~ CLUSTER)$p.value) %>%
    data.frame() %>% dplyr::mutate(p.adjst = stats::p.adjust(as.numeric(pval), method = "bonferroni")) 
  dataCLUSTER <- merge(dataCLUSTER, pValues, by = "variable")

  
  clusterPlot <- ggpubr::ggboxplot(dataCLUSTER, 
                                 x = "variable", 
                                 y = "value", 
                                 width = 0.8,
                                 size = 1,
                                 color = "CLUSTER",
                                 ylab = "Fraction of cells", 
                                 xlab = "Cluster",
                                 font.label = list(size = 20, color = "black"), 
                                 add = "jitter",
                                 palette = colorFunction("CLUSTER")[1:2]) +
    ggtitle(paste0("CLUSTER - 10 Feb 2021 labels (98 deconvoluted; no sex X and Y)")) +
    facet_wrap(~ CellType, scales = "free") +
    geom_text(aes(label = paste0("P=", round(p.adjst, 3)), y = 0.65), size = 2)
  
  
  pathNow <- getwd()
  ggplot2::ggsave(paste0(pathNow, "/img/46_SNFCluster_", ImageName, ".", PNGorPDF),
                  width = 3.5, height = 3.5)
  
  return(invisible(NULL))
}
# [END]
