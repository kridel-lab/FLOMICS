#----------------------------------------------------------------------
#Karin Isaev
#Use RNA-seq TPM matrix from Kallisto
#Plot gene expression versus cell type fraction
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors = F)
options(scipen=999) #avoid scientific notation

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "tidyr",
"data.table", "plyr",
"ggrepel", "stringr", "maftools", "ggpubr", "readxl", "skimr",
 "edgeR", "annotables", "EnvStats", "gridExtra")
lapply(packages, require, character.only = TRUE)

#date
date=Sys.Date()

args = commandArgs(trailingOnly = TRUE) #patient ID
geneA = args[1]
print(geneA) #this should be gene name provided as input
geneB = args[2]
print(geneB)

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

#load in results from kallisto
tpm = fread("RNAseq/counts/2020-09-01_kallisto_gene_based_counts.txt", data.table=F)
colnames(tpm)[1] = "ensgene"
tpm = merge(tpm, genes_class, by = "ensgene")
tpm = as.data.frame(tpm)
rownames(tpm) = tpm$symbol
tpm$symbol = NULL
tpm$ensgene = NULL

# All FLOMICS samples included - load sample information
all.samples.DNAseq.FLOMICS <- fread("metadata/sample_annotations_rcd6Nov2019.csv")

#sample info with rna-seq qc
rnaseq_qc = fread("metadata/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")

#mutation calls
muts = fread("DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl.csv")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

get_mutation_gene_exp_correlation = function(gene_mut, gene_exp){

  print(gene_mut)
  print(gene_exp)

  gene_dat = t(muts[muts$Hugo_Symbol == gene_mut,])
  gene_dat = as.data.frame(gene_dat)
  colnames(gene_dat) = gene_dat[1,]
  gene_dat$rna_seq_file_sample_ID = rownames(gene_dat)
  gene_dat = gene_dat[-1,]
  colnames(gene_dat)[1] = "gene_mut"

  exp_dat = t(tpm[which(rownames(tpm) == gene_exp),])
  exp_dat = as.data.frame(exp_dat)
  if(!(dim(exp_dat)[2] == 0)){

  exp_dat$rna_seq_file_sample_ID = rownames(exp_dat)
  colnames(exp_dat)[1] = "gene"

  z = which(!(str_detect(exp_dat$rna_seq_file_sample_ID, "T1")))
  exp_dat$rna_seq_file_sample_ID[z] = paste(exp_dat$rna_seq_file_sample_ID[z], "_T1", sep="")

  all_dat = merge(gene_dat, exp_dat, by = "rna_seq_file_sample_ID")
  z = which(!(str_detect(rnaseq_qc$rna_seq_file_sample_ID, "T1")))
  rnaseq_qc$rna_seq_file_sample_ID[z] = paste(rnaseq_qc$rna_seq_file_sample_ID[z], "_T1", sep="")
  all_dat = merge(all_dat, rnaseq_qc, by="rna_seq_file_sample_ID")

  #create dataset for plotting
  plot = all_dat %>%
    select(rna_seq_file_sample_ID, gene, gene_mut, STAGE, TYPE)

  #z = which(str_detect(plot$rna_seq_file_sample_ID, "DLC"))
  #plot$STAGE[z] = "DLBCL"
  #z = which(str_detect(plot$rna_seq_file_sample_ID, "RLN"))
  #plot = plot[-z,]

  plot$STAGE = factor(plot$STAGE, levels=c("LIMITED", "ADVANCED"))
  plot$gene_mut = factor(plot$gene_mut, levels=c("0", "1"))

  #binary gene expression vs continuous cell fraction
  g2=ggboxplot(plot, x = "gene_mut", y = "gene", palette = c("#00AFBB", "#E7B800"),
  fill = "gene_mut", facet.by = "STAGE") +
  ylab(paste(gene_exp, "TPM")) +
  xlab(paste(gene_mut, "mutation"))
  g2=ggpar(g2, font.tickslab = c(5, "plain", "black"))+
   stat_compare_means() + stat_n_text()
  print(g2)

  #do not split by stages
  g3=ggboxplot(plot, x = "gene_mut", y = "gene", palette = c("#00AFBB", "#E7B800"),
  fill = "gene_mut") +
  ylab(paste(gene_exp, "TPM")) +
  xlab(paste(gene_mut, "mutation"))
  g3=ggpar(g3, font.tickslab = c(5, "plain", "black"))+
   stat_compare_means() + stat_n_text()
  print(g3)

}}

pdf(paste(geneA, "mut", geneB, "exp.pdf", sep="_"))
get_mutation_gene_exp_correlation(geneA, geneB)
dev.off()
