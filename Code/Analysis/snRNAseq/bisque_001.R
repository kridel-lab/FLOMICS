#-------------------------------------------------------------------------------
#Run BISQUE on FL samples
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

#in terminal go to UHN FLOMICS folder or set this as working directory in
#Rstudio
#/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS for example

source("/Users/kisaev/github/FLOMICS/Code/Analysis/load_scripts_data_KI.R")

#-------------------------------------------------------------------------------
#Bisque analysis
#-------------------------------------------------------------------------------

#1. load bulk seq data into object
#Bulk RNA-seq data can be converted from a matrix
#(columns are samples, rows are genes) to an ExpressionSet as follows:
colnames(exp)[1] = "ensgene"
exp = merge(exp, genes_class, by = "ensgene")
exp = as.data.frame(exp)
rownames(exp) = exp$symbol
exp$ensgene = NULL
exp$symbol = NULL
exp = as.matrix(exp)
#change sample names in bulk rnaseq data to match them here (1,2,3,4)
z = which(colnames(exp) %in% c("LY_FL_062_T1", "LY_FL_064_T1", "LY_FL_076_T1", "LY_FL_227_T1"))
colnames(exp)[z] = c(1, 2, 3) #LY_FL_227_T1 not actually included in bulk rnaseq

bulk.eset <- Biobase::ExpressionSet(assayData = exp)

#2. Single-cell data requires additional information in the ExpressionSet,
#specificially cell type labels and individual labels. Individual labels
#indicate which individual each cell originated from.

#To add this information, Biobase requires it to be stored in a data frame format.
#Assuming we have character vectors of cell type labels (cell.type.labels) and individual labels (individual.labels), we can convert
#scRNA-seq data (with counts also in matrix format) as follows:

#can automatically convert from seurat object
sc.eset = readRDS("Analysis-Files/Seurat/sc_eset_for_bisque.rds")

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=TRUE)
ref.based.estimates <- res$bulk.props
saveRDS(ref.based.estimates, file="Analysis-Files/Seurat/bisque_decomposed_samples.rds")
