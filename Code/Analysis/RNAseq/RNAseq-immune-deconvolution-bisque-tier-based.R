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

bisque_T1 = readRDS("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/tier1_bisque_decomposed_samples.rds")
bisque_T2 = readRDS("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/tier2_bisque_decomposed_samples.rds")
bisque_T3 = readRDS("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Analysis-Files/Seurat/tier3_bisque_decomposed_samples.rds")

old_labels = fread("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Cluster Labels/InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv")
colnames(old_labels)[2] = "SAMPLE_ID"
print(table(old_labels$SNFClust))

#labels updated after mutations for n=31 plos medicine patients were reanalzyed
labels = fread("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Cluster Labels/InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv")
colnames(labels)[2] = "SAMPLE_ID"
labels$SNFClust = labels$SNFClust10Feb2021
print(table(labels$SNFClust))

rnaseq_qc = merge(rnaseq_qc, labels, by="SAMPLE_ID")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

get_bisque_summ = function(dat, tier){
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
  immune_cells$cell_facet=""

#  z=which(immune_cells$cell_type %in% c("B cells_0", "B cells_1", "B cells_2", "naive B or malignant B_9",
#"proliferating B cell_11", "memory B cell_12", "Cluster 13"))
#  immune_cells$cell_facet[z] = "B cells"

  z=which(immune_cells$cell_type %in% c(0, 1, 2, 12, 11, 16))
  immune_cells$cell_facet[z] = "B cells"

#  z=which(immune_cells$cell_type %in% c("Tfh cells_3", "CD8 T cells_4", "CD4 Treg cells_5",
#"naive T cells_7", "memory T cells_8", "proliferating T cell_17"))
#  immune_cells$cell_facet[z] = "T cells"

  z=which(immune_cells$cell_type %in% c(4, 5, 6, 7, 8, 17))
  immune_cells$cell_facet[z] = "T cells"

  immune_cells$cell_facet[immune_cells$cell_facet==""] = "The others"

  immune_cells = immune_cells[order(cell_facet, cell_type)]
  immune_cells$cell_type = factor(immune_cells$cell_type, levels=unique(immune_cells$cell_type))

#  immune_cells$cell_type = factor(immune_cells$cell_type, levels=c("B cells_0", "B cells_1",
#"B cells_2", "naive B or malignant B_9", "proliferating B cell_11", "memory B cell_12",
#"Cluster 13", "Tfh cells_3", "CD8 T cells_4", "CD4 Treg cells_5", "naive T cells_7", "memory T cells_8",
#"proliferating T cell_17", "Cluster 6", "macrophage or monocyte_10", "stromal cells_14",
#"endothelial cells_15", "Cluster 16", "macrophage or monocyte_18" ,"Cluster 19"))

  immune_cells$cell_facet = factor(immune_cells$cell_facet, levels=unique(immune_cells$cell_facet))

  #plot distribution of each cell type frequency across disease and stages
  #set up file for plotting

  #calculate p-value for each comparison
  #looking at STAGE
  stage_analysis = as.data.table(filter(immune_cells, !(is.na(STAGE))))
  res_stage_analysis = as.data.table(stage_analysis %>% group_by(cell_type) %>%
        dplyr::summarize(pval = wilcox.test(value ~ STAGE)$p.value))
  res_stage_analysis$fdr = p.adjust(res_stage_analysis$pval)
  res_stage_analysis$fdr=round(res_stage_analysis$fdr, digits=4)
  stage_analysis=merge(stage_analysis, res_stage_analysis, by="cell_type")
  z = which(duplicated(stage_analysis[,c("cell_type", "fdr")]))
  stage_analysis$fdr[z] = ""

  #looking at TYPE
  type_analysis = as.data.table(filter(immune_cells, !(is.na(TYPE))))
  res_type_analysis = as.data.table(type_analysis %>% group_by(cell_type) %>%
        dplyr::summarize(pval = kruskal.test(value ~ TYPE)$p.value))
  res_type_analysis$fdr = p.adjust(res_type_analysis$pval)
  res_type_analysis$fdr=round(res_type_analysis$fdr, digits=4)
  type_analysis=merge(type_analysis, res_type_analysis, by="cell_type")
  z = which(duplicated(type_analysis[,c("cell_type", "fdr")]))
  type_analysis$fdr[z] = ""

  #looking at InfinumClust
  infi_analysis = as.data.table(filter(immune_cells, !(is.na(InfinumClust))))
  res_infi_analysis = as.data.table(infi_analysis %>% group_by(cell_type) %>%
        dplyr::summarize(pval = wilcox.test(value ~ InfinumClust)$p.value))
  res_infi_analysis$fdr = p.adjust(res_infi_analysis$pval)
  res_infi_analysis$fdr=round(res_infi_analysis$fdr, digits=4)
  infi_analysis=merge(infi_analysis, res_infi_analysis, by="cell_type")
  z = which(duplicated(infi_analysis[,c("cell_type", "fdr")]))
  infi_analysis$fdr[z] = ""

  #looking at SNFClust
  snf_analysis = as.data.table(filter(immune_cells, !(is.na(SNFClust))))
  res_snf_analysis = as.data.table(snf_analysis %>% group_by(cell_type) %>%
        dplyr::summarize(pval = wilcox.test(value ~ SNFClust)$p.value))
  res_snf_analysis$fdr = p.adjust(res_snf_analysis$pval)
  res_snf_analysis$fdr=round(res_snf_analysis$fdr, digits=4)
  snf_analysis=merge(snf_analysis, res_snf_analysis, by="cell_type")
  z = which(duplicated(snf_analysis[,c("cell_type", "fdr")]))
  snf_analysis$fdr[z] = ""

  #looking at SNFClust
  tseq_analysis = as.data.table(filter(immune_cells, !(is.na(tSeqClust))))
  res_tseq_analysis = as.data.table(tseq_analysis %>% group_by(cell_type) %>%
        dplyr::summarize(pval = wilcox.test(value ~ tSeqClust)$p.value))
  res_tseq_analysis$fdr = p.adjust(res_tseq_analysis$pval)
  res_tseq_analysis$fdr=round(res_tseq_analysis$fdr, digits=4)
  tseq_analysis=merge(tseq_analysis, res_tseq_analysis, by="cell_type")
  z = which(duplicated(tseq_analysis[,c("cell_type", "fdr")]))
  tseq_analysis$fdr[z] = ""

  #plotting - automate for all the groups

  pdf(paste("Analysis-Files/Immune-Deconvolution/", date, "_", tier, "_BISQUE_results.pdf", sep=""), width=9,height=6)

  #STAGE
  g1 = ggboxplot(stage_analysis, x="cell_type", y="value", fill="STAGE",group="STAGE", outlier.shape=20) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste(names(table(patients_dat$STAGE)[1]), table(patients_dat$STAGE)[1],
  names(table(patients_dat$STAGE)[2]), table(patients_dat$STAGE)[2]))#+
  g1 = facet(g1, facet.by="cell_facet", scales = "free")+
     geom_text(aes(label = fdr, y=0.55), size=2)+xlab("Cell Type")+ylab("Fraction of cells")
  g1 = ggpar(g1, legend="bottom")

  #TYPE
  g2 = ggboxplot(type_analysis, x="cell_type", y="value", fill="TYPE", outlier.shape=20) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste(names(table(patients_dat$TYPE)[1]), table(patients_dat$TYPE)[1],
    names(table(patients_dat$TYPE)[2]), table(patients_dat$TYPE)[2],
    names(table(patients_dat$TYPE)[3]), table(patients_dat$TYPE)[3]))
  g2 = facet(g2, facet.by="cell_facet", scales = "free")+
       geom_text(aes(label = fdr, y=0.55), size=2)+xlab("Cell Type")+ylab("Fraction of cells")
   g2 = ggpar(g2, legend="bottom")

  #InfinumClust
  g3 = ggboxplot(infi_analysis, x="cell_type", y="value", fill="InfinumClust",palette = "jco", outlier.shape=20) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste(names(table(patients_dat$InfinumClust)[1]), table(patients_dat$InfinumClust)[1],
 names(table(patients_dat$InfinumClust)[2]), table(patients_dat$InfinumClust)[2]))
  g3 = facet(g3, facet.by="cell_facet", scales = "free")+
      geom_text(aes(label = fdr, y=0.55), size=2)+xlab("Cell Type")+ylab("Fraction of cells")
  g3 = ggpar(g3, legend="bottom")

  #SNF "1" = "#4363D8", "2" = "#F58231"
  snf_analysis$SNFClust = factor(snf_analysis$SNFClust, levels=c(1,2))
  g4 = ggboxplot(snf_analysis, x="cell_type", y="value", fill="SNFClust",palette = c("#4363D8", "#F58231"), outlier.shape=20) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste(names(table(patients_dat$SNFClust)[1]), table(patients_dat$SNFClust)[1],
 names(table(patients_dat$SNFClust)[2]), table(patients_dat$SNFClust)[2]))
  g4 = facet(g4, facet.by="cell_facet", scales = "free")+
     geom_text(aes(label = fdr, y=0.55), size=2)+xlab("Cell Type")+ylab("Fraction of cells")
  g4 = ggpar(g4, legend="bottom")

  #tSeqClust
  g5 = ggboxplot(tseq_analysis, x="cell_type", y="value", fill="tSeqClust",palette = "jco", outlier.shape=20) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste(names(table(patients_dat$tSeqClust)[1]), table(patients_dat$tSeqClust)[1],
 names(table(patients_dat$tSeqClust)[2]), table(patients_dat$tSeqClust)[2]))
  g5 = facet(g5, facet.by="cell_facet", scales = "free")+
    geom_text(aes(label = fdr, y=0.55), size=2)+xlab("Cell Type")+ylab("Fraction of cells")
  g5 = ggpar(g5, legend="bottom")

  print(g1)
  print(g2)
  print(g3)
  print(g4)
  print(g5)

  dev.off()
}

get_bisque_summ(bisque_T1, "tier1")
get_bisque_summ(bisque_T2, "tier2")
get_bisque_summ(bisque_T3, "tier3")
