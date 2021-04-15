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

bisque_T1 = readRDS("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/Bisque/tier1_bisque_decomposed_samples.rds")
bisque_T2 = readRDS("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/Bisque/tier2_bisque_decomposed_samples.rds")
bisque_T3 = readRDS("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/Bisque/tier3_bisque_decomposed_samples.rds")

old_labels = fread("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Cluster Labels/InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv")
colnames(old_labels)[2] = "SAMPLE_ID"
print(table(old_labels$SNFClust))

#labels updated after mutations for n=31 plos medicine patients were reanalzyed
labels = fread("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Cluster Labels/InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv")
colnames(labels)[2] = "SAMPLE_ID"
labels$SNFClust = labels$SNFClust10Feb2021
print(table(labels$SNFClust))

#all RNA-seq info
rnaseq_qc = merge(rnaseq_qc, labels, by="SAMPLE_ID")

#mutation data
muts = fread("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl.csv")

dat = bisque_T3
z = which(str_detect(colnames(dat), "FL"))
dat = dat[,z]

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
immune_cells$cell_facet=""

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

  temp_dat = immune_cells

  z = which(!(str_detect(temp_dat$rna_seq_file_sample_ID, "T1")))
  temp_dat$rna_seq_file_sample_ID[z] = paste(temp_dat$rna_seq_file_sample_ID[z], "_T1", sep="")

  #create dataset for plotting
  plot = temp_dat %>%
    select(rna_seq_file_sample_ID, SNFClust10Feb2021, cell_type, value) %>%
      melt(measure = "value")

  #add gene mutation info
  plot = merge(plot, gene_dat)
  plot$gene_mut = factor(plot$gene_mut, levels=c("0", "1"))

  #binary mutation status vs continuous cell fraction
  g2=ggboxplot(plot, x = "gene_mut", y = "value",
  fill = "gene_mut", palette = c("#00AFBB", "#E7B800")) +
  ylab("Cell type fraction") +
  xlab(gene)
  g2=ggpar(g2, font.tickslab = c(5, "plain", "black"))
  g2=facet(g2, facet.by="cell_type", ,
         panel.labs.font = list(color = "black", size=6))+
         stat_compare_means(aes(group = gene_mut), label = "p.signif")+ylim(c(0,1))
  print(g2)

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
file_name = paste("Analysis-Files/Immune-Deconvolution/", "all_methods_", "full_immune_cells_vs_mutations", ".csv", sep="")
write.csv(all_results, file_name, quote=F, row.names=F)
