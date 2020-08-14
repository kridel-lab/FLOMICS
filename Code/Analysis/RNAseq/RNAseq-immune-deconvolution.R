#----------------------------------------------------------------------
#Karin Isaev
#Use RNA-seq count matrix (first normalize to TPM)
#Run through immune deconvolution tools to estimate fraction of immune
#cells present in each sample
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors = F)
#avoid scientific notation
options(scipen=999)
#load libraries
packages <- c("dplyr", "readr", "ggplot2", "tidyr",
"data.table", "plyr",
"ggrepel", "stringr", "maftools", "ggpubr", "readxl", "skimr",
 "edgeR", "annotables", "EnvStats")
library(gridExtra)
lapply(packages, require, character.only = TRUE)
library(limSolve)
library(immunedeconv) #<- main package with tools for immune deconvolution

#date
date=Sys.Date()

#getwd() --> FLOMICS teams folder
#cd /Users/kisaev/UHN/kridel-lab - Documents/FLOMICS

#----------------------------------------------------------------------
#data filtering
#----------------------------------------------------------------------

#gene annotations
#UCSC gene classes - only protein coding genes
genes_class = as.data.table(grch37)
genes_class = as.data.table(filter(genes_class, biotype == "protein_coding"))
genes_class = as.data.table(filter(genes_class, !(is.na(entrez))))
genes_class = unique(genes_class[,c("ensgene", "symbol")])
#keep only one ens id per gene name
z = which(duplicated(genes_class$symbol))
genes_class = genes_class[-z,]

#load in count matrix produced via STAR
exp <- fread("RNAseq/counts/STAR_quantmode_counts_matrix_FL_136_patients.txt")
#change ensembl ids to gene names as required for package to work
colnames(exp)[1] = "ensgene"
exp = merge(exp, genes_class, by = "ensgene")
exp = as.data.frame(exp)
rownames(exp) = exp$symbol
exp$symbol = NULL
exp$ensgene = NULL

# All FLOMICS samples included - load sample information
all.samples.DNAseq.FLOMICS <- fread("metadata/sample_annotations_rcd6Nov2019.csv")

#sample info with rna-seq qc
rnaseq_qc = fread("metadata/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. prepare normalize count matrix - TMM
y <- DGEList(counts=exp)
keep <- filterByExpr(y, min.count=10)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
tmm = cpm(y)

#run analysis - save plots and output from analysis

run_immdeco = function(exp_matrix, tool_used, qc_data){
  #2. try running a tool
  res = deconvolute(tmm, tool_used, tumor=TRUE)
  immune_cells = as.data.frame(res)

  #tool specific plotting
  if(tool_used == "quantiseq"){
  tool_plot = res %>%
    gather(sample, fraction, -cell_type) %>%
    # plot as stacked bar chart
    ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
      geom_bar(stat='identity') +
      coord_flip() +
      scale_fill_brewer(palette="Paired") +
      scale_x_discrete(limits = rev(levels(res)))+
              theme(axis.text.y = element_text(color = "grey20", size = 4))
  }

  #3. only keep 132 patients used in RNA-seq in the end
  z = which(colnames(immune_cells) %in% rnaseq_qc$rna_seq_file_sample_ID)
  immune_cells = immune_cells[,c(1,z)]

  #4. add tag based on whether sample is limited advanced and FL vs DLBCL
  immune_cells = melt(as.data.table(immune_cells))
  colnames(immune_cells)[2] = "rna_seq_file_sample_ID"
  immune_cells = merge(immune_cells, rnaseq_qc, by = "rna_seq_file_sample_ID")
  immune_cells$Cluster = factor(immune_cells$Cluster)

  #5. plot distribution of each cell type frequency across disease and stages

  pdf(paste("Analysis-Files/Immune-Deconvolution/", date, "_", tool_used, "_results.pdf", sep=""), width=15)
  g1 = ggboxplot(filter(immune_cells, !(is.na(STAGE))), x="cell_type", y="value", fill="STAGE") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste("Tool used=", tool_used, " via TMM values calcualted by EdgeR", sep=""))+
  stat_compare_means(aes(group = STAGE), label = "p.signif")

  g2 = ggboxplot(immune_cells, x="cell_type", y="value", fill="TYPE") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste("Tool used=", tool_used, " via TMM values calcualted by EdgeR", sep=""))+
  stat_compare_means(aes(group = TYPE), label = "p.signif")

  g3 = ggboxplot(immune_cells, x="cell_type", y="value", fill="Cluster",palette = "jco") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste("Tool used=", tool_used, " via TMM values calcualted by EdgeR", sep=""))+
  stat_compare_means(aes(group = Cluster), label = "p.signif")

  print(g1)
  print(g2)
  print(g3)
  print(tool_used)

  dev.off()

  #6. save file
  file_name=paste("Analysis-Files/Immune-Deconvolution/", date, "_", tool_used, "_results.txt", sep="")
  write.table(immune_cells, file_name, quote=F, row.names=F, sep=";")

  print("done analysis")

}

run_immdeco(tmm, "quantiseq", rnaseq_qc)
run_immdeco(tmm, "mcp_counter", rnaseq_qc)
