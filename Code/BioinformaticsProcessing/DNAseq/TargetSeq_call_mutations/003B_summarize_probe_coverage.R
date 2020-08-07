#----------------------------------------------------------------------
#Summarize probe coverage from picard - part two
#----------------------------------------------------------------------

#Karin Isaev
#Started August 5th 2020
#tested on R version 3.5.0

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#this script takes output from picard collect targeted pcr metrics
#merges it with sample information
#make summary csv file with coverage per gene per patient
#simple barplot summaries of coverage per gene/probe/region
#calculate differences in coverage between coding and non-coding regions

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)
date=Sys.Date()

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
"plyr",
"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)


#----------------------------------------------------------------------
#Load most up to date data matrix from previous script "all_res"
#----------------------------------------------------------------------

all_res = fread(list.files(pattern="picard_tools_coverage")[length(list.files(pattern="picard_tools_coverage"))])

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#calculate mean coverage between coding and non-coding
c_nc = as.data.table(all_res %>% group_by(type) %>% dplyr::summarize(mean_cov = mean(mean_coverage)))
write.csv(c_nc, file=paste(date, "coding_noncoding_coverage_summary.csv", sep="_"),
  quote=F, row.names=F)

#calculate mean coverage by gene
gene_cov = as.data.table(all_res %>% group_by(gene, type) %>% dplyr::summarize(mean_cov = mean(mean_coverage)))
gene_cov = gene_cov[order(-mean_cov)]
gene_cov$gene = factor(gene_cov$gene, levels=unique(gene_cov$gene))
write.csv(gene_cov, file=paste(date, "gene_based_coverage_summary.csv", sep="_"),
  quote=F, row.names=F)

#mean target based coverage
target_cov = as.data.table(all_res %>% group_by(region, type) %>% dplyr::summarize(mean_cov = mean(mean_coverage)))
target_cov = target_cov[order(-mean_cov)]
target_cov$region = factor(target_cov$region, levels=unique(target_cov$region))
write.csv(target_cov, file=paste(date, "target_based_coverage_summary.csv", sep="_"),
  quote=F, row.names=F)

#calculate mean coverage by sample
sample_cov = as.data.table(all_res %>% group_by(External_identifier, TYPE, STAGE) %>% dplyr::summarize(mean_cov = mean(mean_coverage)))
sample_cov = sample_cov[order(-mean_cov)]
sample_cov$External_identifier = factor(sample_cov$External_identifier, levels=unique(sample_cov$External_identifier))
write.csv(sample_cov, file=paste(date, "sample_based_coverage_summary.csv", sep="_"),
  quote=F, row.names=F)

#plot barplot summary per gene coverage
pdf("/cluster/projects/kridelgroup/FLOMICS/DATA/002_coverage_detailed_summaries.pdf",
width=10)

#mean cov coding and non-coding
ggplot(c_nc, aes(x=type, y=mean_cov, fill=type))+
geom_bar(stat="identity", color="black")+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#mean gene based coverage
ggplot(gene_cov, aes(x=gene, y=mean_cov, fill=type))+
geom_bar(stat="identity", color="black")+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#mean target based coverage
ggplot(target_cov, aes(x=region, y=mean_cov))+
geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x=element_blank())

#mean sample based coverage
ggplot(sample_cov, aes(x=External_identifier, y=mean_cov, fill=STAGE))+
geom_bar(stat="identity", color="black")+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
