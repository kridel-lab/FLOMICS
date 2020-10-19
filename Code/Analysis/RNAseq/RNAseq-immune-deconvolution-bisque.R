#----------------------------------------------------------------------
#Karin Isaev
#Use RNA-seq count matrix (first normalize to TPM)
#Run through immune deconvolution tools to estimate fraction of immune
#cells present in each sample
#----------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

#in terminal go to UHN FLOMICS folder or set this as working directory in
#Rstudio
#/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS for example

source("/Users/kisaev/github/FLOMICS/Code/Analysis/load_scripts_data_KI.R")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

bisque = readRDS("Analysis-Files/Seurat/bisque_decomposed_samples.rds")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#visualize results from bisque
immune_cells = as.data.frame(bisque)
immune_cells$cell_type = rownames(immune_cells)

#add tag based on whether sample is limited advanced and FL vs DLBCL
immune_cells = as.data.table(immune_cells)
immune_cells = melt((immune_cells))
colnames(immune_cells)[2] = "rna_seq_file_sample_ID"
immune_cells = merge(immune_cells, rnaseq_qc, by = "rna_seq_file_sample_ID")
immune_cells$Cluster = factor(immune_cells$Cluster)

head(immune_cells)

#plot distribution of each cell type frequency across disease and stages
#set up file for plotting
pdf(paste("Analysis-Files/Immune-Deconvolution/", date, "BISQUE_results.pdf", sep=""), width=15)
g1 = ggboxplot(filter(immune_cells, !(is.na(STAGE))), x="cell_type", y="value", fill="STAGE") +
theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
stat_compare_means(aes(group = STAGE), label = "p.signif")

g2 = ggboxplot(immune_cells, x="cell_type", y="value", fill="TYPE") +
theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
stat_compare_means(aes(group = TYPE), label = "p.signif")

g3 = ggboxplot(immune_cells, x="cell_type", y="value", fill="Cluster",palette = "jco") +
theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
stat_compare_means(aes(group = Cluster), label = "p.signif")

print(g1)
print(g2)
print(g3)

dev.off()
