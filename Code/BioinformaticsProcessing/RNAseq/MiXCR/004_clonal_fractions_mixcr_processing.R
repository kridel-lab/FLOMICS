#----------------------------------------------------------------------------------
#004_clonal_fractions_mixcr_processing.R
#----------------------------------------------------------------------------------

#Sarah Russell
#Date: August 4, 2020
#This script takes contatenated MiXCR results file
#and examines the distribution of the top clonotype (highest clonal fraction)
#found within each samples

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
# Remove duplicates based sample ID & select fo highest clonal fraction hit
top_fraction=clonotypes[!duplicated(clonotypes$sample),]

top_clones=data.frame(top_fraction$cloneCount,top_fraction$cloneFraction,top_fraction$allCHitsWithScore, top_fraction$allVHitsWithScore, top_fraction$aaSeqCDR3,top_fraction$sample,top_fraction$STAGE,top_fraction$TYPE,top_fraction$CLUSTER)
colnames(top_clones)=c("count","fraction","chits", "vhits","aaSeqCDR3","sample","stage","type","cluster")

###add receptor chain
top_clones$chain <- ifelse(grepl("TRA",top_clones$vhits)==TRUE, "TRA",
  ifelse(grepl("TRB",top_clones$vhits)==TRUE,"TRB",
  ifelse(grepl("TRD",top_clones$vhits)==TRUE,"TRD",
  ifelse(grepl("TRG",top_clones$vhits)==TRUE,"TRG",
  ifelse(grepl("IGH",top_clones$vhits)==TRUE,"IGH",
  ifelse(grepl("IGK",top_clones$vhits)==TRUE,"IGK","IGL"))))))

top_type=as.data.frame(table(top_clones$chain,top_clones$type),stringsAsFactors=F)
top_stage=as.data.frame(table(top_clones$chain,top_clones$stage),stringsAsFactors=F)
top_cluster=as.data.frame(table(top_clones$chain,top_clones$cluster),stringsAsFactors=F)

###plot histogram showing distribution of chains within all three comparison groups
pdf("/cluster/home/srussell/top_clone_type.pdf", width=5, height=4)
g=ggbarplot(top_type, x="Var1", y="Freq", label = top_type$Freq , lab.size=2, fill="Var1", facet.by=c("Var2"),
palette=c("jco"))+
theme_bw()+
theme(text = element_text(size=10)) + ylab("Number of Top Clones in each sample") + xlab("Type")
ggpar(g, legend="bottom", legend.title="Chain") + geom_hline(yintercept=0)
dev.off()

pdf("/cluster/home/srussell/top_clone_stage.pdf", width=5, height=4)
g=ggbarplot(top_stage, x="Var1", y="Freq", label = top_stage$Freq , lab.size=2, fill="Var1", facet.by=c("Var2"),
palette=c("jco"))+
theme_bw()+
theme(text = element_text(size=10)) + ylab("Number of Top Clones in each sample") + xlab("Stage")
ggpar(g, legend="bottom",legend.title="Chain") + geom_hline(yintercept=0)
dev.off()

pdf("/cluster/home/srussell/top_clone_cluster.pdf", width=5, height=4)
g=ggbarplot(top_cluster, x="Var1", y="Freq", label = top_cluster$Freq , lab.size=2, fill="Var1", facet.by=c("Var2"),
palette=c("jco"))+
theme_bw()+
theme(text = element_text(size=10)) + ylab("Number of Top Clones in each sample") + xlab("Cluster")
ggpar(g, legend="bottom",legend.title="Chain") + geom_hline(yintercept=0)
dev.off()

###########################################################################
###examine distribution of chains within all samples that have an IGH top hit

top_clones$IGHchain <- ifelse(grepl("IGHA",top_clones$chits)==TRUE, "IGHA",
  ifelse(grepl("IGHG",top_clones$chits)==TRUE,"IGHG",
  ifelse(grepl("IGHD",top_clones$chits)==TRUE,"IGHD",
  ifelse(grepl("IGHE",top_clones$chits)==TRUE,"IGHE",
  ifelse(grepl("IGHM",top_clones$chits)==TRUE,"IGHM",NA)))))

top_clones_IGH=top_clones[!is.na(top_clones$IGHchain),]
#> dim(top_clones_IGH)
#[1] 42 11

IGHtop_stage=as.data.frame(table(top_clones_IGH$IGHchain,top_clones_IGH$stage),stringsAsFactors=F)

###plot histogram of distribution of IGH chains
###only looking at ADV vs LIM, IGH to hits only found within some FL samples
pdf("/cluster/home/srussell/IGHtop_clone_stage.pdf", width=5, height=4)
g=ggbarplot(IGHtop_stage, x="Var1", y="Freq", label = IGHtop_stage$Freq , lab.size=2, fill="Var1", facet.by=c("Var2"),
palette=c("jco"))+
theme_bw()+
theme(text = element_text(size=10)) + ylab("Number of Top IGH Clones in each sample") + xlab("Stage")
ggpar(g, legend="bottom",legend.title="Chain") + geom_hline(yintercept=0)
dev.off()

###########################################################################

###plot boxplot showing distribution of top clonal fraction VALUES within all three comparison groups
pdf("/cluster/home/srussell/top_fraction_per_type.pdf", width=7, height=4)
g=ggboxplot(top_clones, x = "type", y = "fraction",
          title = "Highest clonal fraction for each sample per type", ylab = "Clonal fraction", xlab = "Sample Type",
          color = "type", palette = "jco",
          add = "jitter")
g + stat_compare_means()
dev.off()


pdf("/cluster/home/srussell/top_fraction_per_stage.pdf", width=7, height=4)
g=ggboxplot(top_clones[top_clones$type=="FL",], x = "stage", y = "fraction",
          title = "Highest clonal fraction for each sample per stage", ylab = "Clonal fraction", xlab = "Sample Stage",
          color = "stage", palette = "jco",
          add = "jitter")
g + stat_compare_means()
dev.off()

pdf("/cluster/home/srussell/top_fraction_per_cluster.pdf", width=7, height=4)
g=ggboxplot(top_clones, x = "cluster", y = "fraction",
          title = "Highest clonal fraction for each sample per cluster", ylab = "Clonal fraction", xlab = "Sample Cluster",
          color = "cluster", palette = "jco",
          add = "jitter")
g + stat_compare_means()
dev.off()
###########################################################################
###could probably summarize many of above graphs in functions
###graphs located in Teams > FLOMICS > Mixcr > clonotype_fractions
