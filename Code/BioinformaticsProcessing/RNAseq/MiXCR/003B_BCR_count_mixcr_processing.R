#----------------------------------------------------------------------------------
#003B_BCR_count_mixcr_processing.R
#----------------------------------------------------------------------------------

#Sarah Russell
#Date: July 29th, 2020
#This script takes contatenated MiXCR results file
#and examines the distribution of unique clonotypes found within
#samples based on unique CDR3 sequence *for BCRs

#----------------------------------------------------------------------------------
#PACKAGES
#----------------------------------------------------------------------------------
date = Sys.Date()

library(data.table)
library(dplyr)
library(plyr)
library(tidyverse)
library(readxl)
library(ggpubr)

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS")

#----------------------------------------------------------------------------------
#DATA
#----------------------------------------------------------------------------------
###information regarding each sample and which stage of disease and cluster they are part of
sample_info = as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.xlsx"))

clonotypes=fread("_2020-07-30_MIXCR_OUTPUT_WITH_SAMPLE_ANNOTATION.csv", sep="\t")
#> dim(clonotypes)
#[1] 48292      39

#----------------------------------------------------------------------------------
#ANALYSIS
#----------------------------------------------------------------------------------
cdr3=data.frame(clonotypes$cloneCount,clonotypes$cloneFraction,clonotypes$allVHitsWithScore,clonotypes$aaSeqCDR3,clonotypes$sample,clonotypes$STAGE,clonotypes$TYPE,clonotypes$CLUSTER)
colnames(cdr3)=c("count","fraction","vhits","aaSeqCDR3","sample","stage","type","cluster")

#> dim(cdr3)
#[1] 48292       8
#write.table(cdr3, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS/", date, "CDR3_CLEANED_VHITS.csv", sep="_"), sep="\t", quote=F, row.names=F)

cdr3_TCR=as.data.table(filter(cdr3, !grepl('TR', vhits)))
dim(cdr3_TCR)
#[1] 38802       8
length(unique(cdr3_TCR$sample))
#[1] 130

#remove sequences that have frameshift or stop codon
cdr3_TCR <- cdr3_TCR[!grepl("\\~",cdr3_TCR$aaSeqCDR3),]
cdr3_TCR <- cdr3_TCR[!grepl("\\*",cdr3_TCR$aaSeqCDR3),]
cdr3_TCR <- cdr3_TCR[!grepl("\\_",cdr3_TCR$aaSeqCDR3),]

#> dim(cdr3_TCR)
#[1] 36799    9

#> length(unique(cdr3_TCR$sample))
#[1] 130


ID=as.character(unique(cdr3_TCR$sample))
get_uniq_CDR=function(ID){
    count=length(unique(cdr3_TCR[cdr3_TCR$sample==ID,]$aaSeqCDR3))
    count_matrix=as.data.frame(cbind(count,ID))
    count_matrix$type <- ifelse(grepl("DLC",count_matrix$ID)==TRUE, "DLC",
      ifelse(grepl("FL",count_matrix$ID)==TRUE,"FL","RLN"))
    return(count_matrix)
}
counts = as.data.table(ldply(llply(ID, get_uniq_CDR)))

stage=as.character(cdr3_TCR$stage[match(unique(cdr3_TCR$sample),cdr3_TCR$sample)])
cluster=as.character(cdr3_TCR$cluster[match(unique(cdr3_TCR$sample),cdr3_TCR$sample)])
all_counts=as.data.frame(cbind(counts,stage,cluster))

all_counts$count=as.character(all_counts$count)
all_counts$count=as.numeric(all_counts$count)

#boxplot of counts per type
pdf("/cluster/home/srussell/clonotypes_per_type_BCR.pdf", width=7, height=4)
ggboxplot(all_counts, x = "type", y = "count",
          title = "Number of unique T cell clones (unique CDR3 sequence) per type", ylab = "TCR clones", xlab = "Sample Type",
          color = "type", palette = "jco",
          add = "jitter")
dev.off()

#boxplot of counts per stage
pdf("/cluster/home/srussell/clonotypes_per_stage_BCR.pdf", width=7, height=4)
ggboxplot(all_counts[all_counts$type=="FL",], x = "stage", y = "count",
          title = "Number of unique T cell clones (unique CDR3 sequence) per stage", ylab = "TCR clones", xlab = "Sample Stage",
          color = "stage", palette = "jco",
          add = "jitter")
dev.off()

#boxplot of counts per stage
pdf("/cluster/home/srussell/clonotypes_per_cluster_BCR.pdf", width=7, height=4)
ggboxplot(all_counts, x = "cluster", y = "count",
          title = "Number of unique T cell clones (unique CDR3 sequence) per cluster", ylab = "TCR clones", xlab = "Sample Cluster",
          color = "cluster", palette = "jco",
          add = "jitter")
dev.off()
