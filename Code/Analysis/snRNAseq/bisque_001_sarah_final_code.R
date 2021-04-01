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

#prep bulk data
rownames(exp) = exp$symbol
exp$ensgene = NULL
exp$symbol = NULL
exp$biotype = NULL
exp = as.matrix(exp)

#2. Single-cell data requires additional information in the ExpressionSet,
#specificially cell type labels and individual labels. Individual labels
#indicate which individual each cell originated from.

#To add this information, Biobase requires it to be stored in a data frame format.
#Assuming we have character vectors of cell type labels (cell.type.labels) and individual labels (individual.labels), we can convert
#scRNA-seq data (with counts also in matrix format) as follows:

#can automatically convert from seurat object
#sc.eset = readRDS("Analysis-Files/Seurat/seurat_objects/sc_eset_for_bisque_rmFL277dim20_prelim.rds")
output="/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/April2021/"
sc.eset=readRDS(paste(output, "sc_eset_for_bisque.rds", sep=""))

get_bisque_res = function(tier){

  #RUN for TIER1
  if(tier=="tier1"){

    patients = t1_patients
    print(tier)
    print(dim(patients))
    z = which(colnames(exp) %in% patients$rna_seq_file_sample_ID)
    exp_dat=exp[,z]

    z = which(colnames(exp_dat) %in% c("LY_FL_062_T1", "LY_FL_064_T1", "LY_FL_076_T1")) #, "LY_FL_227_T1"))
    colnames(exp_dat)[z] = c(1, 2, 3) #LY_FL_227_T1 not actually included in bulk rnaseq

    bulk.eset <- Biobase::ExpressionSet(assayData = exp_dat)

    res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=TRUE)
    ref.based.estimates <- res$bulk.props
    saveRDS(ref.based.estimates, file=paste("Analysis-Files/Seurat/", "tier1", "_bisque_decomposed_samples.rds", sep=""))
  }

  if(tier=="tier2"){
    patients = t2_patients
    print(tier)
    print(dim(patients))
    z = which(colnames(exp) %in% patients$rna_seq_file_sample_ID)
    exp_dat=exp[,z]

    z = which(colnames(exp_dat) %in% c("LY_FL_062_T1", "LY_FL_064_T1", "LY_FL_076_T1")) #, "LY_FL_227_T1"))
    colnames(exp_dat)[z] = c(1, 2, 3) #LY_FL_227_T1 not actually included in bulk rnaseq
    bulk.eset <- Biobase::ExpressionSet(assayData = exp_dat)

    res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=TRUE)
    ref.based.estimates <- res$bulk.props
    saveRDS(ref.based.estimates, file=paste("Analysis-Files/Seurat/", "tier2", "_bisque_decomposed_samples.rds", sep=""))
  }

  if(tier=="tier3"){
    patients = t3_patients
    print(tier)
    print(dim(patients))
    z = which(colnames(exp) %in% patients$rna_seq_file_sample_ID)
    exp_dat=exp[,z]

    z = which(colnames(exp_dat) %in% c("LY_FL_062_T1", "LY_FL_064_T1", "LY_FL_076_T1")) #, "LY_FL_227_T1"))
    colnames(exp_dat)[z] = c(1, 2, 3) #LY_FL_227_T1 not actually included in bulk rnaseq
    bulk.eset <- Biobase::ExpressionSet(assayData = exp_dat)

    res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=TRUE)
    ref.based.estimates <- res$bulk.props
    saveRDS(ref.based.estimates, file=paste("Analysis-Files/Seurat/", "tier3", "_bisque_decomposed_samples.rds", sep=""))
  }

  print("done")
}

tiers = c("tier1", "tier2", "tier3")
get_bisque_res("tier1")
get_bisque_res("tier2")
get_bisque_res("tier3")

#3. using full data as before and set up bulk rnaseq data for each tier based on quality of rna-seq

#change sample names in bulk rnaseq data to match them here (1,2,3,4)
z = which(colnames(exp) %in% c("LY_FL_062_T1", "LY_FL_064_T1", "LY_FL_076_T1")) #, "LY_FL_227_T1"))
colnames(exp)[z] = c(1, 2, 3) #LY_FL_227_T1 not actually included in bulk rnaseq

bulk.eset <- Biobase::ExpressionSet(assayData = exp)

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=TRUE)
ref.based.estimates <- res$bulk.props
saveRDS(ref.based.estimates, file="Analysis-Files/Seurat/bisque_decomposed_samples.rds")
