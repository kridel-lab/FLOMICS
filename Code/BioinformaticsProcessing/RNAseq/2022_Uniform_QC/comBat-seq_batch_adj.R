##Ting Liu-Sep,2022
##ComBat-seq is a batch effect adjustment tool for bulk RNA-seq count data. It is an improved model based on the popular ComBat
#https://github.com/zhangyuqing/ComBat-seq

#module load R/4.1.0

packages <- c("edgeR", "plyr", "dplyr", "AnnotationDbi", "org.Hs.eg.db", "GeneOverlap", "ggpubr", "gplots", "ComplexHeatmap", "ggplot2", "gProfileR", "gprofiler2", "forcats","stringr")
lapply(packages, require, character.only = TRUE)
library("devtools")
library("sva")

rm(list=ls())

setwd("/Users/tingliu/Desktop/UHN_projects/FLOMICS/RNA_seq_2022_rerun/total_365_samples/batch_QC/passed_290/passed_htseq/")
#setwd("/cluster/projects/kridelgroup/FLOMICS/2022_combined_RNA-Seq_analysis/batch_QC/Sep2022_passed_290/passed_htseq/")
output="/Users/tingliu/Desktop/UHN_projects/FLOMICS/RNA_seq_2022_rerun/total_365_samples/batch_QC/passed_290/"

date <- Sys.Date()
# Read in file names
#files <- list.files(pattern ="_gene_counts.txt$")
files <- list.files()
##define a function to get the prefix
get_prefix=function(file_name){
  file_prefix=unlist(strsplit(file_name, "\\_LY"))[1]
  return(file_prefix)
}

#ini the group arry
group=c()
for(i in files){
  cat (i,"\n")
  prefix=get_prefix(i)
  group=c(group,prefix)
}

# Generate DGElist object
df <- readDGE(files, header = FALSE, group = group)
df$samples # shows counts per library

#plot the PCA before adj
raw_df <- calcNormFactors(df)
pdf(paste(output, "diff_exp_MDS_plot_before_adj_290.pdf", sep=""))
plotMDS(raw_df, labels=group, col=c(rep("green",208),rep("blue",14), rep("red",68)))
dev.off()

#check num of 0 count samples
table(rowSums(df$counts==0)==290)   #summary the genes that have 0 count in 290

#ComBat-seq takes untransformed, raw count matrix as input. Same as ComBat, it requires a known batch variable.
#so we don't use filterByExpr here
count_matrix<-as.matrix(df$counts)
batch <-group
adjusted <- ComBat_seq(count_matrix, batch=batch, group=NULL)
df$counts<-adjusted

##filterByExpr function automaticaly filter low exps genes
#df_2=df #save the df_2 to compare normalised lcmp
keep <- filterByExpr(df)
table(keep)  #auto filtering
df <- df[keep, keep.lib.sizes=FALSE]

# Normalize.
# Normalization may not actually be required:
# Normalization issues arise only to the extent that technical factors
# have sample-specific effects.
df <- calcNormFactors(df)

#design <- model.matrix(~group)
#df <- estimateDisp(df,design)

#pdf(paste(output, "diff_exp_MDS_plot.pdf", sep=""))
#plotMDS(df, labels=group, col=c(rep("#999999",208),rep("#E69F00",14), rep("#56B4E9",68)))
#dev.off()

#plot the PCA for adjusted samples
pdf(paste(output, "diff_exp_MDS_plot_adj_290.pdf", sep=""))
plotMDS(df, labels=group, col=c(rep("green",208),rep("blue",14), rep("red",68)))
dev.off()

