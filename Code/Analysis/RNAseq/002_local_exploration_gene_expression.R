#local exploration of gene expression across clusters and stages
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors = F)
#avoid scientific notation
options(scipen=999)
#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr",
"mclust", "data.table", "plyr",
"ggrepel", "stringr", "maftools", "ggpubr", "readxl", "skimr",
 "edgeR", "annotables", "EnvStats")
library(gridExtra)
lapply(packages, require, character.only = TRUE)

#date
date=Sys.Date()

#directory with FLOMICS related matrices
setwd("/Users/kisaev/github/FLOMICS/Data")

#----------------------------------------------------------------------
#data filtering
#----------------------------------------------------------------------

#gene annotations
#UCSC gene classes - only protein coding genes
genes_class = as.data.table(grch37)
genes_class = as.data.table(filter(genes_class, biotype == "protein_coding"))
genes_class = as.data.table(filter(genes_class, !(is.na(entrez))))
genes_class = unique(genes_class[,c("ensgene", "symbol")])

exp <- fread("2020-06-18STAR_quantmode_counts_matrix_FL_136_patients.txt")

#these are raw counts
#normalize gene counts using EdgeR

#set up edgeR?
exp = as.data.frame(exp)
rownames(exp) = exp[,1]
exp = exp[,-1]

z = which(pats$Tumor_Sample_Barcode %in% colnames(exp))
pats = pats[z,]

y <- DGEList(counts=exp, group=pats$CLUSTER_InfiniumClust)

#Usually a gene is required to have a count of 5-10 in a library to be
#considered expressed in that library. Users should also filter with
#count-per-million (CPM) rather than filtering on the
#counts directly, as the latter does not account for differences in
#library sizes between samples.

#We filter out lowly expressed genes using the following commands:
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

#get TMM
y <- calcNormFactors(y)
y <- estimateDisp(y)

#get matrix of TMM values
norm_counts = as.data.frame(t(cpm(y)))

get_gene_name = function(gene){
  name = genes_class$symbol[which(genes_class$ensgene == gene)]
  return(name)
}

colnames(norm_counts) = sapply(colnames(norm_counts), get_gene_name)
norm_counts$Tumor_Sample_Barcode = rownames(norm_counts)
z = which(duplicated(colnames(norm_counts)))
norm_counts = norm_counts[,-z]
norm_counts = as.data.table(norm_counts)
norm_counts = merge(norm_counts, pats, by="Tumor_Sample_Barcode")
norm_counts = as.data.table(norm_counts)

norm_counts$STAGE[which(str_detect(norm_counts$Tumor_Sample_Barcode, "RLN"))] = "RLN"
norm_counts$STAGE[which(str_detect(norm_counts$Tumor_Sample_Barcode, "DLC"))] = "DLBCL"
norm_counts$CLUSTER_InfiniumClust = factor(norm_counts$CLUSTER_InfiniumClust)

get_gene_boxplot = function(gene){
  z = which(colnames(norm_counts) %in% c(gene, "Tumor_Sample_Barcode", "CLUSTER_InfiniumClust", "STAGE"))
  gene_dat = norm_counts[,..z]
  colnames(gene_dat)[2] = "gene_counts"
  p = ggboxplot(x="CLUSTER_InfiniumClust", y="gene_counts", add="jitter",
  fill="STAGE", gene_dat)+ylab(gene)+stat_n_text()+ theme_minimal()
  print(p)
}
