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
"data.table", "plyr", "openxlsx",
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

#cell types vs immune cells summary
muts_immune = fread("Analysis-Files/Immune-Deconvolution/all_methods_full_immune_cells_vs_mutations.txt")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. remove cells types that we don't need
#muts_immune = as.data.table(filter(muts_immune, !(cell_type %in% c("Endothelial cell",
#"Myeloid dendritic cell",  "Common myeloid progenitor" ,
#"Mast cell", "Myeloid dendritic cell activated", "microenvironment score", "Eosinophil", "Neutrophil",
#"immune score", "Plasmacytoid dendritic cell", "Macrophage M1", "Macrophage M2",
#"Monocyte", "Macrophage", "Granulocyte-monocyte progenitor", "Mast cell resting",
#"Macrophage M0", "Mast cell activated", "uncharacterized cell", "Myeloid dendritic cell resting"))))

#2. filter those comparisons where the median was 0 or around 0
muts_immune = as.data.table(filter(muts_immune, median > 0.05))
muts_immune$fdr = NULL
muts_immune = muts_immune[order(Wilcox_Pval)]

#3. summarize which genes are associated with which cell types
muts_immune$cell_gene = paste(muts_immune$gene, muts_immune$cell_type, sep="_")
write.xlsx(muts_immune, file="Analysis-Files/Immune-Deconvolution/Cell_Type_Mutation_Associations_CleanedUp.xlsx")

#4. keep only those with sig association in at least one method
muts_immune = as.data.table(filter(muts_immune, Wilcox_Pval < 0.1))
gene_cell = as.data.table(table(muts_immune$cell_gene))
gene_cell = gene_cell[order(-N)]

#5. summarize in tile plot
t = as.data.table(table(muts_immune$cell_type))
t=t[order(-N)]
muts_immune$cell_type = factor(muts_immune$cell_type, levels=t$V1)

tt = as.data.table(table(muts_immune$gene))
tt=tt[order(-N)]
muts_immune$gene = factor(muts_immune$gene, levels=tt$V1)

pdf("Analysis-Files/Immune-Deconvolution/Immune_Cells_Mutations_Associations_Summary_sig.pdf")
ggplot(muts_immune, aes(x = cell_type, y = gene, fill = -log10(Wilcox_Pval))) +
      geom_tile(aes(fill = -log10(Wilcox_Pval)), colour = "grey50") +
      labs(x = "Cell Type", y = "Gene with mutation", fill = "-log10(Wilcox P value)") +
      theme_classic() + rotate_x_text()+
      scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
dev.off()
