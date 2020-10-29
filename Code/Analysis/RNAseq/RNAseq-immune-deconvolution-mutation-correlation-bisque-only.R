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

#args = commandArgs(trailingOnly = TRUE) #patient ID
#index = args[1]
#print(index) #this should be gene name provided as input
#gene_name=index

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
full_cells=fread("Analysis-Files/Immune-Deconvolution/2020-09-16_full_cell_types_results.txt")
some_cells=fread("Analysis-Files/Immune-Deconvolution/2020-09-16_subset_cell_types_results.txt")

#mutation calls
muts = fread("DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl.csv")

#results from bisque
bisque = readRDS("Analysis-Files/Seurat/bisque_decomposed_samples.rds")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

get_mutation_immune_correlation = function(gene, im_data){

  print(gene)
  gene_dat = t(muts[muts$Hugo_Symbol == gene,])
  gene_dat = as.data.frame(gene_dat)
  colnames(gene_dat) = gene_dat[1,]
  gene_dat$rna_seq_file_sample_ID = rownames(gene_dat)
  gene_dat = gene_dat[-1,]
  colnames(gene_dat)[1] = "gene_mut"

  temp_dat = im_data
  z = which(str_detect(colnames(im_data), "T2"))
  im_data = im_data[,-z]
  z = which(str_detect(colnames(im_data), "RLN"))
  im_data = im_data[,-z]

  z = which(!(str_detect(colnames(im_data), "T1")))
  colnames(im_data)[z] = paste(colnames(im_data)[z], "_T1", sep="")

  im_data = melt(im_data)
  colnames(im_data)[2] = "rna_seq_file_sample_ID"

  im_data = merge(im_data, gene_dat, by = "rna_seq_file_sample_ID")
  colnames(im_data)[2] = "cell_type"
  #create dataset for plotting
  plot = im_data %>%
    select(rna_seq_file_sample_ID, cell_type, value, gene_mut) %>%
      melt(measure = "value")

  #plot$STAGE = factor(plot$STAGE, levels=c("LIMITED", "ADVANCED"))
  plot$gene_mut = factor(plot$gene_mut, levels=c("0", "1"))

  #binary gene expression vs continuous cell fraction
  g2=ggboxplot(plot, x = "gene_mut", y = "value",
  fill = "gene_mut", palette = c("#00AFBB", "#E7B800")) +
  ylab("Cell type fraction") +
  xlab(gene)
  g2=ggpar(g2, font.tickslab = c(5, "plain", "black"))
  g2=facet(g2, facet.by="cell_type", ,
         panel.labs.font = list(color = "black", size=6))+
         stat_compare_means(aes(group = gene_mut), label = "p.signif")
  print(g2)

  #if(!( 0 %in% c(table(plot$gene_mut, plot$STAGE)))){
  if((!(table(plot$gene_mut)[1] == 0))& (!(table(plot$gene_mut)[2] == 0))){
  res <- as.data.table(plot %>% group_by(cell_type) %>%
         do(w = wilcox.test(value~gene_mut, data=., paired=FALSE)) %>%
         dplyr::summarise(cell_type, Wilcox_Pval = w$p.value))

  meds = as.data.table(plot %>% group_by(cell_type) %>%
          dplyr::summarize(median=median(value)))

  res = merge(res, meds)
  res$gene = gene
  head(res)
  return(res)
#}
}
}

genes = as.list(unique(muts$Hugo_Symbol))
methods = "bisque"

#obtain summary of associations between gene mutations and immune cell values across methods
print(methods)
file_name = paste("Analysis-Files/Immune-Deconvolution/", "BISQUE", "_","full_immune_cells_vs_mutations", ".pdf", sep="")
pdf(file_name, width=12, height=12)
all_res = as.data.table(ldply(llply(genes, get_mutation_immune_correlation, bisque, .progress="text")))
dev.off()
all_res = all_res[order(Wilcox_Pval)]
all_res$fdr = p.adjust(all_res$Wilcox_Pval, method="fdr")

all_results = all_res
file_name = paste("Analysis-Files/Immune-Deconvolution/", "all_methods_", "full_immune_cells_vs_mutations", ".txt", sep="")
write.table(all_results, file_name, quote=F, row.names=F, sep="\t")
