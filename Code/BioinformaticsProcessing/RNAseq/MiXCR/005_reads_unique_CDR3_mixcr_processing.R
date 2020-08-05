#----------------------------------------------------------------------------------
#005_reads_unique_CDR3_mixcr_processing.R
#----------------------------------------------------------------------------------

#Sarah Russell
#Date: Aug 4, 2020
#This script takes contatenated MiXCR results file
#and examines the relationship between the number of successfully
#mapped reads and unique clonotypes (unique CDR3 sequences)
#present per sample

#before processing in R:
#cat *.report | grep 'Successfully aligned reads:.*' > all_aligned.txt
#ls *.report | grep -oP '.*?(?=\.)' > all_report_names.txt

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

read_counts=read.table("all_aligned.txt")
read_counts=read_counts$V4
samples=read.table("all_report_names.txt")

reads_per_sample=cbind(samples,read_counts)
colnames(reads_per_sample)=c("sample","reads")


clonotypes=fread("_2020-07-30_MIXCR_OUTPUT_WITH_SAMPLE_ANNOTATION.csv", sep="\t")

cdr3=data.frame(clonotypes$cloneCount,clonotypes$cloneFraction,clonotypes$allVHitsWithScore,clonotypes$aaSeqCDR3,clonotypes$sample,clonotypes$STAGE,clonotypes$TYPE,clonotypes$CLUSTER)
colnames(cdr3)=c("count","fraction","vhits","aaSeqCDR3","sample","stage","type","cluster")

#remove sequences that have frameshift or stop codon
cdr3_filt <- cdr3[!grepl("\\~",cdr3$aaSeqCDR3),]
cdr3_filt <- cdr3[!grepl("\\*",cdr3$aaSeqCDR3),]
cdr3_filt <- cdr3[!grepl("\\_",cdr3$aaSeqCDR3),]

#get unique CDR3 counts
ID=as.character(unique(cdr3$sample))
get_uniq_CDR=function(ID){
    count=length(unique(cdr3_filt[cdr3_filt$sample==ID,]$aaSeqCDR3))
    count_matrix=as.data.frame(cbind(count,ID))
    count_matrix$type <- ifelse(grepl("DLC",count_matrix$ID)==TRUE, "DLC",
      ifelse(grepl("FL",count_matrix$ID)==TRUE,"FL","RLN"))
    return(count_matrix)
}
counts = as.data.table(ldply(llply(ID, get_uniq_CDR)))

stage=as.character(cdr3_filt$stage[match(unique(cdr3_filt$sample),cdr3_filt$sample)])
cluster=as.character(cdr3_filt$cluster[match(unique(cdr3_filt$sample),cdr3_filt$sample)])
all_counts=as.data.frame(cbind(counts,stage,cluster))

all_counts$count=as.character(all_counts$count)
all_counts$count=as.numeric(all_counts$count)
all_counts$ID=as.character(all_counts$ID)

reads_sample=filter(reads_per_sample,sample %in% all_counts$ID)
print(all_counts$ID==reads_sample$sample)
all_counts=cbind(all_counts,reads_sample$reads)
colnames(all_counts)[6] = "reads"

#plot scatter plot of CDR3 counts vs. read counts
pdf("/cluster/home/srussell/reads_vs_uniqueclones.pdf", width=5, height=4)
ggscatter(all_counts, x = "count", y = "reads",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          xlab = "Unique CDR3 sequence per sample",
          ylab = "Successfully aligned reads per sample",
          add.params = list(color = "blue",
                            fill = "lightgray")
          )+
  stat_cor(method = "spearman", label.x = 2500, label.y = 60000)  # Add correlation coefficient
dev.off()


pdf("/cluster/home/srussell/reads_vs_uniqueclones_filt.pdf", width=5, height=4)
ggscatter(all_counts[all_counts$reads <= 20000 & all_counts$count < 2000,], x = "count", y = "reads",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          xlab = "Unique CDR3 sequence per sample",
          ylab = "Successfully aligned reads per sample",
          add.params = list(color = "blue",
                            fill = "lightgray")
          )+
  stat_cor(method = "spearman", label.x = 600, label.y = 5000)  # Add correlation coefficient
dev.off()
###########################################################################
###could probably summarize many of above graphs in functions
###graphs located in Teams > FLOMICS > Mixcr > read_coverage
