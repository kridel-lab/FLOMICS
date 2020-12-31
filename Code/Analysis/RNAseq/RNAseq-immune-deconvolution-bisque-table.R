#----------------------------------------------------------------------
#Sarah Russell
#Generate descriptive statistics for for bisque results comparing stage,
#type and SNF Cluster
#----------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

#in terminal go to UHN FLOMICS folder or set this as working directory in
#Rstudio
#/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS for example

setwd("/Users/sarahrussell/UHN/kridel-lab - FLOMICS")
source("/Users/srussell/github/FLOMICS/Code/Analysis/load_scripts_data_KI.R")
library(doBy)
output("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Seurat/clustering/")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

setwd(output)
#load in dim20 and dim30 seurat objects
bisque_20 = readRDS("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Seurat/seurat_objects/bisque_decomposed_samples_rmFL277dim20_prelim.rds")
bisque_30 = readRDS("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Seurat/seurat_objects/bisque_decomposed_samples_rmFL277dim30_prelim.rds")

labels = fread("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Cluster Labels/InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv")
colnames(labels)[2] = "SAMPLE_ID"
rnaseq_qc = merge(rnaseq_qc, labels, by="SAMPLE_ID")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

get_bisque_summstat = function(dat){
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

#generate table comparing FL stage
  dat1 = data.table(compare_means(value ~ STAGE, data = filter(immune_cells, !(is.na(STAGE))), group.by = "cell_type", method = "wilcox.test", p.adjust.method = "fdr") %>%
   mutate(p.adj = format.pval(p.adj, digits = 5)) %>%
   mutate(p.adj.signif = ifelse(p.adj < .05, "sig", "NA"), comparison="STAGE"))
   dat1$group1=NULL
   dat1$group2=NULL
#calculate adv + lim mean,median
 adv=summaryBy(value ~ cell_type, data = filter(immune_cells, STAGE=="ADVANCED"),
             FUN = list(mean, median))
 lim=summaryBy(value ~ cell_type, data = filter(immune_cells, STAGE=="LIMITED"),
              FUN = list(mean, median))
 stat1=merge(adv,lim,by="cell_type")
#group1=ADV, group2=LIM
 colnames(stat1)=c("cell_type","group1.mean","group1.median","group2.mean","group2.median")
 dat1=merge(dat1,stat1,by="cell_type")
 dat1$group3.mean="NA"
 dat1$group3.median="NA"

 #generate table comparing sample type
 dat2 = data.table(compare_means(value ~ TYPE, data = immune_cells, group.by = "cell_type", method = "kruskal.test", p.adjust.method = "fdr") %>%
   mutate(p.adj = format.pval(p.adj, digits = 5)) %>%
   mutate(p.adj.signif = ifelse(p.adj < .05, "sig", "NA"),comparison="TYPE"))
 #calculate DLBCL, FL, RLN mean,median
  DLBCL=summaryBy(value ~ cell_type, data = filter(immune_cells, TYPE=="DLBCL"),
              FUN = list(mean, median))
  FL=summaryBy(value ~ cell_type, data = filter(immune_cells, TYPE=="FL"),
              FUN = list(mean, median))
  RLN=summaryBy(value ~ cell_type, data = filter(immune_cells, TYPE=="RLN"),
              FUN = list(mean, median))
  stat2=merge(DLBCL,FL,by="cell_type")
  stat2=merge(stat2,RLN,by="cell_type")
#group1=DLBCL, group2=FL, group3=RLN
  colnames(stat2)=c("cell_type","group1.mean","group1.median","group2.mean","group2.median","group3.mean","group3.median")
  dat2=merge(dat2,stat2,by="cell_type")

#generate table comparing SNFClust
  dat3 = data.table(compare_means(value ~ SNFClust, data = filter(immune_cells, !(is.na(SNFClust))), group.by = "cell_type", method = "wilcox.test", p.adjust.method = "fdr") %>%
   mutate(p.adj = format.pval(p.adj, digits = 5)) %>%
   mutate(p.adj.signif = ifelse(p.adj < .05, "sig", "NA"),comparison="SNFClust"))
   dat3$group1=NULL
   dat3$group2=NULL
#calculate DLBCL, FL, RLN mean,median
  snf1=summaryBy(value ~ cell_type, data = filter(immune_cells, SNFClust=="1"),
            FUN = list(mean, median))
  snf2=summaryBy(value ~ cell_type, data = filter(immune_cells, SNFClust=="2"),
            FUN = list(mean, median))
  stat2=merge(snf1,snf2,by="cell_type")
#group1=1, group2=2
  colnames(stat2)=c("cell_type","group1.mean","group1.median","group2.mean","group2.median")
  dat3=merge(dat3,stat2,by="cell_type")
  dat3$group3.mean="NA"
  dat3$group3.median="NA"

#combine all into one table
  dat=rbind(dat1,dat2,dat3)
  dat$.y.=NULL

  return(dat)

}
bisque_stat_20=get_bisque_summstat(bisque_20)
write.csv(bisque_stat_20, file="/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Seurat/clustering/dim20_3samples/summary_stats_20.csv",
quote=F, row.names=F)


bisque_stat_30=get_bisque_summstat(bisque_30)
write.csv(bisque_stat_30, file="/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Seurat/clustering/dim30_3samples/summary_stats_30.csv",
quote=F, row.names=F)
