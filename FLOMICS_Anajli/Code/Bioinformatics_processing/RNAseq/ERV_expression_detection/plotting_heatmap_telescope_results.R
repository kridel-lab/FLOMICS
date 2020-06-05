##Author - Karin Isaev 
##Figure XX

#----SET UP----------------------------------------------------------------------------------

options(stringsAsFactors=F)
setwd("/Users/kisaev/FLOMICS_Anjali/FLOMICS_Anajli/Data")

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
              "ggrepel", "stringr", "maftools", "ggpubr")
lapply(packages, require, character.only = TRUE)

library(gProfileR)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(pheatmap)
library(EnvStats)

date = Sys.Date()

#----DATA-----------------------------------------------------------------------------------

cpm_matrix = fread("CPM_exprs_matrix_norm_filt_ERVs_FL_109_cases.csv", data.table=FALSE)
rownames(cpm_matrix) = cpm_matrix$ERV_ids
cpm_matrix$ERV_ids = NULL

sample_info = fread("2020-01-16_TELESCOPE_OUTPUT_WITH_SAMPLE_ANNOTATION.csv")
sample_info = as.data.table(filter(sample_info, SAMPLE_ID %in% colnames(cpm_matrix)))

diff_exp = fread("2020-01-16_Telescope_ERVs_differentially_ADVANCED_vs_LIMITED.csv")
z = which(rownames(cpm_matrix) %in% diff_exp$transcript)
cpm_matrix = cpm_matrix[z,]

#----ANALYSIS-------------------------------------------------------------------------------

cpm_matrix = as.matrix(cpm_matrix)
sample_info = unique(sample_info[,c("SAMPLE_ID", "STAGE")])
sample_info = sample_info[order(match(SAMPLE_ID, colnames(cpm_matrix)))]
sample_info$STAGE = as.factor(sample_info$STAGE)

sample_info$SAMPLE_ID == colnames(cpm_matrix)

my_sample_col <- data.frame(sample = sample_info$STAGE)
rownames(my_sample_col) = colnames(cpm_matrix)

#pdf("Deeptools_bam_files_ChIPseq_correlation.pdf", width=9)
cpm_matrix = log1p(cpm_matrix)
pheatmap(cpm_matrix, annotation_col = my_sample_col, scale = "row", main="Limited vs Advanced ERVs")
#dev.off()

#test diff expression using boxplot
gene="ERV316A3_2q21.2b"

get_plot = function(gene){
  get_exp_values = as.data.frame((cpm_matrix[which(rownames(cpm_matrix)==gene),]))
  get_exp_values$SAMPLE_ID = rownames(get_exp_values)
  get_exp_values = join(get_exp_values, sample_info)
  colnames(get_exp_values)[1] = "gene_exp"
  get_exp_values$gene_exp = log1p(get_exp_values$gene_exp)
  g = ggboxplot(get_exp_values, x="STAGE", y="gene_exp") + stat_n_text() + ylab("log1p(CPM)") + stat_compare_means() + ggtitle(gene)
  print(g)
}

pdf("testing_diff_exp_genes_boxplots.pdf")
llply(unique(rownames(cpm_matrix)), get_plot)
dev.off()
      
