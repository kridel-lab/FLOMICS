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

labels = fread("/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/Cluster Labels/InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv")
colnames(labels)[2] = "SAMPLE_ID"
rnaseq_qc = merge(rnaseq_qc, labels, by="SAMPLE_ID")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

get_bisque_summ = function(dat){
  #visualize results from bisque
  immune_cells = as.data.frame(dat)
  immune_cells$cell_type = rownames(immune_cells)

  #add tag based on whether sample is limited advanced and FL vs DLBCL
  immune_cells = as.data.table(immune_cells)
  immune_cells = melt((immune_cells))
  colnames(immune_cells)[2] = "rna_seq_file_sample_ID"
  immune_cells = merge(immune_cells, rnaseq_qc, by = "rna_seq_file_sample_ID")
  immune_cells$InfinumClust = factor(immune_cells$InfinumClust)
  immune_cells$SNFClust = factor(immune_cells$SNFClust)
  immune_cells$tSeqClust = factor(immune_cells$tSeqClust)
  patients_dat = unique(immune_cells[,c("SAMPLE_ID", "STAGE", "TYPE", "SNFClust", "InfinumClust", "tSeqClust")])
  head(immune_cells)

  #plot distribution of each cell type frequency across disease and stages
  #set up file for plotting
  pdf(paste("Analysis-Files/Immune-Deconvolution/", date, "_BISQUE_results.pdf", sep=""), width=15)
  g1 = ggboxplot(filter(immune_cells, !(is.na(STAGE))), x="cell_type", y="value", fill="STAGE") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_compare_means(aes(group = STAGE), label = "p.signif")+
  ggtitle(paste(names(table(patients_dat$STAGE)[1]), table(patients_dat$STAGE)[1],
  names(table(patients_dat$STAGE)[2]), table(patients_dat$STAGE)[2]))

  g2 = ggboxplot(immune_cells, x="cell_type", y="value", fill="TYPE") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_compare_means(aes(group = TYPE), label = "p.signif")+
  ggtitle(paste(names(table(patients_dat$TYPE)[1]), table(patients_dat$TYPE)[1],
    names(table(patients_dat$TYPE)[2]), table(patients_dat$TYPE)[2],
    names(table(patients_dat$TYPE)[3]), table(patients_dat$TYPE)[3]))

  g3 = ggboxplot(immune_cells, x="cell_type", y="value", fill="InfinumClust",palette = "jco") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_compare_means(aes(group = Cluster), label = "p.signif")+
  ggtitle(paste(names(table(patients_dat$InfinumClust)[1]), table(patients_dat$InfinumClust)[1],
 names(table(patients_dat$InfinumClust)[2]), table(patients_dat$InfinumClust)[2]))

  g4 = ggboxplot(immune_cells, x="cell_type", y="value", fill="SNFClust",palette = "jco") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_compare_means(aes(group = Cluster), label = "p.signif")+
  ggtitle(paste(names(table(patients_dat$SNFClust)[1]), table(patients_dat$SNFClust)[1],
 names(table(patients_dat$SNFClust)[2]), table(patients_dat$SNFClust)[2]))

  g5 = ggboxplot(filter(immune_cells, !(is.na(SNFClust))), x="cell_type", y="value", fill="SNFClust",palette = "jco") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_compare_means(aes(group = Cluster), label = "p.signif")+
  ggtitle(paste(names(table(patients_dat$SNFClust)[1]), table(patients_dat$SNFClust)[1],
 names(table(patients_dat$SNFClust)[2]), table(patients_dat$SNFClust)[2]))

  g6 = ggboxplot(immune_cells, x="cell_type", y="value", fill="tSeqClust",palette = "jco") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_compare_means(aes(group = Cluster), label = "p.signif")+
  ggtitle(paste(names(table(patients_dat$tSeqClust)[1]), table(patients_dat$tSeqClust)[1],
 names(table(patients_dat$tSeqClust)[2]), table(patients_dat$tSeqClust)[2]))

  g7 = ggboxplot(filter(immune_cells, !(is.na(tSeqClust))), x="cell_type", y="value", fill="tSeqClust",palette = "jco") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_compare_means(aes(group = Cluster), label = "p.signif")+
  ggtitle(paste(names(table(patients_dat$tSeqClust)[1]), table(patients_dat$tSeqClust)[1],
 names(table(patients_dat$tSeqClust)[2]), table(patients_dat$tSeqClust)[2]))

  print(g1)
  print(g2)
  print(g3)
  #print(g4)
  print(g5)
  #print(g6)
  print(g7)

  dev.off()
}

get_bisque_summ(bisque)
