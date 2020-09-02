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

#results from CIBERSORT abs
cibersort=fread("Analysis-Files/Immune-Deconvolution/2020-09-01_cibersort_abs_results.txt")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

get_gene_immune_correlation = function(gene){

  gene_dat = t(tpm[which(rownames(tpm) == gene),])
  gene_dat = as.data.frame(gene_dat)
  gene_dat$rna_seq_file_sample_ID = rownames(gene_dat)

  cibersort_temp = cibersort
  cibersort_temp = merge(cibersort_temp, gene_dat, by = "rna_seq_file_sample_ID")
  colnames(cibersort_temp)[ncol(cibersort_temp)] = "gene"

  #create dataset for plotting
  plot = cibersort_temp %>%
    select(rna_seq_file_sample_ID, cell_type, value, gene, STAGE, TYPE) %>%
      melt(measure = "value")

  z = which(str_detect(plot$rna_seq_file_sample_ID, "DLC"))
  plot$STAGE[z] = "DLBCL"
  z = which(str_detect(plot$rna_seq_file_sample_ID, "RLN"))
  plot = plot[-z,]

  #correlation between cell type value and FOXP3 measure
  g=ggscatter(plot, x = "value", y = "gene",size=1,add = "reg.line",
  color = "black", shape = 21) +
            ylab(paste(gene, "TPM")) +
  stat_cor(aes(color=STAGE), size=2) + xlab("Cell type fraction")
  g=ggpar(g, font.tickslab = c(5, "plain", "black"))
  g=facet(g, facet.by="cell_type", ,
         panel.labs.font = list(color = "black", size=6))
  print(g)

}

get_gene_immune_correlation("FOXP3")
get_gene_immune_correlation("EZH2")
