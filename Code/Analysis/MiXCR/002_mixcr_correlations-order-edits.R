#----------------------------------------------------------------------
#001_mixcr_correlations.R
#Sarah Rusell
#Use MiXCR metrics and bisque results
#to plot correlation values
#R/3.6.3

#done locally
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors = F)
options(scipen=999) #avoid scientific notation
date = Sys.Date()

setwd("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Mixcr")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "tidyr",
              "data.table", "plyr","ggrepel", "stringr",
              "ggpubr", "readxl", "tidyverse","reshape2",
              "ggcorrplot","naniar","corrr","purrr","Hmisc")
lapply(packages, require, character.only = TRUE)

source("/Users/sarahrussell/github/FLOMICS/Code/Analysis/MiXCR/correlation_source.R")

#----------------------------------------------------------------------
#load data
#----------------------------------------------------------------------

# All FLOMICS samples included - load sample information
all.samples.RNAseq.FLOMICS <- fread("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/metadata/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")

#cluster labels
all.clus.labels=fread("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Cluster Labels/InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv")

#results from bisque
bisque = as.data.frame(readRDS("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Seurat/SARAH/bisque_decomposed_samples.rds"))

#mixcr diversity
diversity=as.data.frame(fread("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Mixcr/_2020-10-08_repetoire_diversity.csv"))
colnames(diversity)[1]="sample"

#tcr & bcr mixcr metrics
tcr_metrics=as.data.frame(fread("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Mixcr/_2020-10-02_tcr_metrics.csv"))
tcr_metrics$tcr_clones=tcr_metrics$tra_clones+tcr_metrics$trb_clones

bcr_metrics=as.data.frame(fread("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Mixcr/_2020-10-02_bcr_metrics.csv"))
bcr_metrics$bcr_clones=bcr_metrics$IGH_clones+bcr_metrics$IGK_clones+bcr_metrics$IGL_clones

#----------------------------------------------------------------------
#data filtering
#----------------------------------------------------------------------
metrics=merge(
  bcr_metrics %>% select(sample,total_clones,BCR_richness=bcr_clones,BCR_unique_CDR3=unique_CDR3),
  tcr_metrics %>% select(sample,TCR_richness=tcr_clones,TCR_unique_CDR3=unique_CDR3),
  by="sample", all = T
)

#keep in mind that sc samples not included in bisque: LY_FL_062_T1,LY_FL_064_T1,LY_FL_076_T1
cells=select(bisque,contains(metrics$sample))
#ncol(cells) should be 127

cells <- data.frame(t(cells))
colnames(cells)=sub("X","cluster.",colnames(cells))
celltype=colnames(cells)
cells$sample=rownames(cells)

metrics=filter(metrics,sample %in% cells$sample)
diversity=filter(diversity,sample %in% cells$sample)

labels=merge(
  all.samples.RNAseq.FLOMICS %>% select(ID=SAMPLE_ID,sample=rna_seq_file_sample_ID,TYPE,STAGE),
  all.clus.labels %>% select(ID,SNFClust),by="ID") %>%
  filter(., sample %in% cells$sample)

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

formatted_cors(metrics,diversity)

formatted_cors_gr(metrics,diversity,"STAGE")
formatted_cors_gr(metrics,diversity,"SNFClust")
