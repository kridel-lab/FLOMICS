#-------------------------------------------------------------------------------
#generate TMM matrix
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

y <- DGEList(counts=exp)
y$samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

dge <- calcNormFactors(y, method = "TMM")
dge$samples
tmm <- cpm(dge)

write.table(tmm, file="/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/RNAseq/TMM/tmm_matrix_via_star_counts.txt",
quote=F, row.names=F)
