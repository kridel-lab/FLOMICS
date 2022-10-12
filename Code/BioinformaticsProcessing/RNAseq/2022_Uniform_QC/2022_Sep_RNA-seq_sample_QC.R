##explore the STAR log file

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
date <- Sys.Date()
##
##read rrna, exonic, insertsize, starlog tables
setwd("/Users/tingliu/Desktop/UHN_projects/FLOMICS/RNA_seq_2022_rerun/total_365_samples/sample_QC/")
starlog=read.table("2022-07-27_STAR_QC_total_input_short.txt",header=T,sep = "\t") 
output="/Users/tingliu/Desktop/UHN_projects/FLOMICS/RNA_seq_2022_rerun/total_365_samples/sample_QC/"
starlog$squence_core=factor(starlog$squence_core,order=T)
rrna_perct=read.table("all_rRNA_contam_summary.txt",header=T,sep = "\t") 
bam_qc=read.table("insertsize_qualimap_summary.txt",header=T,sep = "\t") 
rna_qc=read.table("RNA_qualimap_summary.txt",header=T,sep = "\t") 
exon_perct_num=as.numeric(unlist(strsplit(rna_qc$exonic_perct, "%")))  #get rid of %, make it as numeric
##create new dataframe
exonic_perct_num<-data.frame("sample_id"=rna_qc$sample_id, "exonic_perct_num"=exon_perct_num)
#add picard coding pct cal
picard_out=read.table("/Users/tingliu/Desktop/UHN_projects/FLOMICS/RNA_seq_2022_rerun/total_365_samples/sample_QC/all_picard_RnaMetrics_summary_v2.txt",header=T,sep = "\t")

#merge all tables
all_qc_matrix= merge(merge(merge(merge(bam_qc, rrna_perct),exonic_perct_num),starlog),picard_out)
#table(all_qc_matrix$median_insert_size >= 150)
#table(all_qc_matrix$picard_RnaMetrics_perct >= 5)
#table(all_qc_matrix$exonic_perct_num >= 20)
#check how many of samples both rrna_contam_perct<=35 and picard_RnaMetrics_perct>=5
#table((all_qc_matrix$rrna_contam_perct<=35)&(all_qc_matrix$picard_RnaMetrics_perct>=5))
#table((all_qc_matrix$rrna_contam_perct<=35)&(all_qc_matrix$exonic_perct_num>=20))

#add cols to indicate if pass the QC,rrna_contam_perct<=35 for all tiers
#exonic_perct_num>20 for tier1, and get rid of an outliner, 250 samples were kept
all_qc_matrix$qc_tier1<-ifelse((all_qc_matrix$rrna_contam_perct<=35)&(all_qc_matrix$exonic_perct_num>20),'Y','N')
all_qc_matrix$qc_tier1[2]="N"   #get rid of LY_DLC_002 as its median insert size is an outliner
#picard_RnaMetrics_perct >= 5 for tier2, 290 samples in total 
all_qc_matrix$qc_tier2<-ifelse((all_qc_matrix$rrna_contam_perct<=35)&(all_qc_matrix$picard_RnaMetrics_perct>=5),'Y','N')
#picard_RnaMetrics_ALIGNED_perct >= 5 for tier3, 301 samples in total 
all_qc_matrix$qc_tier3<-ifelse((all_qc_matrix$rrna_contam_perct<=35)&(all_qc_matrix$picard_RnaMetrics_ALIGNED_perct>=5),'Y','N')

#passed_qc_matrix<-all_qc_matrix[((all_qc_matrix$rrna_contam_perct<=35)&(all_qc_matrix$picard_RnaMetrics_perct >= 5)),]
#not_passed_qc_matrix<-all_qc_matrix[all_qc_matrix$qc=='N',]
#passed_qc_matrix<-all_qc_matrix[all_qc_matrix$qc=='Y',]
#write.table(passed_qc_matrix, paste(date,"_FLOMICS_passed_qc301_matrix.txt",sep = ""), quote=F, sep = "\t", row.names = F)
#passed_qc_matrix$contm_lv<-ifelse((passed_qc_matrix$rrna_contam_perct<1),'Y','N')
#passed_qc_matrix$cod_lv<-ifelse((passed_qc_matrix$picard_RnaMetrics_perct>10),'Y','N')


write.table(all_qc_matrix, paste(date,"_FLOMICS_all_qc_matrix_v2.txt",sep = ""), quote=F, sep = "\t", row.names = F)

#all_qc_matrix %>% 
#  filter(mean_insert_size > 150) %>% 
#  filter(rrna_perct< 35) %>% 
#  filter(exonic_perct_num > 5)->passed_qc_matrix
passed_qc_matrix<-all_qc_matrix[all_qc_matrix$qc=='Y',]
not_passed_qc_matrix<-all_qc_matrix[all_qc_matrix$qc=='N',]
write.table(passed_qc_matrix, paste(date,"_FLOMICS_passed_qc301_matrix.txt",sep = ""), quote=F, sep = "\t", row.names = F)

dim(passed_qc_matrix)

##plot the uniq mapp rate between passed and not passed

setwd("/Users/tingliu/Desktop/UHN_projects/FLOMICS/RNA_seq_2022_rerun/total_365_samples/plot_passed_comp/")
pdf("pct_of_Uniquely_mapped_reads_between_passed_and_nonpassed.pdf", width=8)
my.p1 <- ggplot(data = all_qc_matrix, aes(x = qc))
my.p1 +  geom_boxplot(aes(y = Pct_Uniquely_mapped_reads, colour = qc))
dev.off()

#plot median insert size between passed and not passed
pdf("median_insert_size_diff_between_passed_and_nonpassed.pdf", width=8)
my.p2 <- ggplot(data = all_qc_matrix, aes(x = qc))
my.p2 +  geom_boxplot(aes(y = median_insert_size, colour = qc))
dev.off()

pdf("mean_insert_size_diff_between_passed_and_nonpassed.pdf", width=8)
my.p2 <- ggplot(data = all_qc_matrix, aes(x = qc))
my.p2 +  geom_boxplot(aes(y = mean_insert_size, colour = qc))
dev.off()

#plot the raw data without 
my.p1 <- ggplot(data = all_qc_matrix, aes(x = squence_core)) + geom_boxplot(aes(y = Pct_Uniquely_mapped_reads, colour = squence_core)) +
  labs(title = "% of Uniquely_mapped_reads")
my.p2 <- ggplot(data = all_qc_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = median_insert_size, colour = squence_core))+
  labs(title = "median_insert_size")
my.p3 <- ggplot(data = all_qc_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = rrna_contam_perct, colour = squence_core)) +
  labs(title = "rrna contamination perct")
my.p4 <- ggplot(data = all_qc_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = exonic_perct_num, colour = squence_core))+
  labs(title = "% exonic reigin reads")
my.p5 <- ggplot(data = all_qc_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = picard_RnaMetrics_perct, colour = squence_core))+
  labs(title = "% exonic reigin bases")
my.p6 <- ggplot(data = all_qc_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = Unmapped_reads, colour = squence_core))+
  labs(title = "% pct Unmapped_reads")

##plot the merged dig
#https://www.jqhtml.com/53242.html
setwd("/Users/tingliu/Desktop/UHN_projects/FLOMICS/RNA_seq_2022_rerun/total_365_samples/plot_passed_comp/")
pdf(paste(date, "_all_sample_plot_cores.pdf", sep=""),width=22,height=10)
ggarrange(my.p1,my.p2,my.p3,my.p4,my.p5,my.p6,ncol=3,nrow=2,labels=c("A","B","C","D","E","F"))
dev.off()


#plot passed 301 samples
my.p1 <- ggplot(data = passed_qc_matrix, aes(x = squence_core)) + geom_boxplot(aes(y = Pct_Uniquely_mapped_reads, colour = squence_core)) +
  labs(title = "% of Uniquely_mapped_reads")
my.p2 <- ggplot(data = passed_qc_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = median_insert_size, colour = squence_core))+
  labs(title = "median_insert_size")
my.p3 <- ggplot(data = passed_qc_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = rrna_contam_perct, colour = squence_core)) +
  labs(title = "rrna contamination perct")
my.p4 <- ggplot(data = passed_qc_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = exonic_perct_num, colour = squence_core))+
  labs(title = "% exonic reigin reads")
my.p5 <- ggplot(data = passed_qc_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = picard_RnaMetrics_perct, colour = squence_core))+
  labs(title = "% exonic reigin bases")
my.p6 <- ggplot(data = passed_qc_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = Unmapped_reads, colour = squence_core))+
  labs(title = "% pct Unmapped_reads")

setwd("/Users/tingliu/Desktop/UHN_projects/FLOMICS/RNA_seq_2022_rerun/total_365_samples/plot_passed_comp/")
pdf(paste(date, "_passed_samples_plot_core.pdf", sep=""),width=22,height=10)
ggarrange(my.p1,my.p2,my.p3,my.p4,my.p5,my.p6,ncol=3,nrow=2,labels=c("A","B","C","D","E","F"))
dev.off()



#plot passed 251 exon 20 samples, v2 290 samples cuse use PF_BASES to cal the picard_RnaMetrics_perct
my.p1 <- ggplot(data = passed_exon20_matrix, aes(x = squence_core)) + geom_boxplot(aes(y = Pct_Uniquely_mapped_reads, colour = squence_core)) +
  labs(title = "% of Uniquely_mapped_reads")
my.p2 <- ggplot(data = passed_exon20_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = median_insert_size, colour = squence_core))+
  labs(title = "median_insert_size")
my.p3 <- ggplot(data = passed_exon20_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = rrna_contam_perct, colour = squence_core)) +
  labs(title = "rrna contamination perct")
my.p4 <- ggplot(data = passed_exon20_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = exonic_perct_num, colour = squence_core))+
  labs(title = "% exonic reigin reads")
my.p5 <- ggplot(data = passed_exon20_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = picard_RnaMetrics_perct, colour = squence_core))+
  labs(title = "% exonic reigin bases")
my.p6 <- ggplot(data = passed_exon20_matrix, aes(x = squence_core)) +  geom_boxplot(aes(y = Unmapped_reads, colour = squence_core))+
  labs(title = "% pct Unmapped_reads")

setwd("/Users/tingliu/Desktop/UHN_projects/FLOMICS/RNA_seq_2022_rerun/total_365_samples/plot_passed_comp/")
pdf(paste(date, "_passed_exon20_251_samples_plot_core.pdf", sep=""),width=22,height=10)
ggarrange(my.p1,my.p2,my.p3,my.p4,my.p5,my.p6,ncol=3,nrow=2,labels=c("A","B","C","D","E","F"))
dev.off()

##investigate batch effect and generate the report by BatchQC
#create the gene by sample matrix
#load the expression data
#output="/Users/tingliu/Desktop/UHN_projects/FLOMICS/RNA_seq_2022_rerun/total_365_samples/htseq_out/all/"

#setwd("/Users/tingliu/UHN/kridel-lab - FLOMICS/RNAseq/2022_combined_RNAseq_total365/htseq_out/rename_counts/")
setwd("/Users/tingliu/Desktop/UHN_projects/FLOMICS/RNA_seq_2022_rerun/total_365_samples/passed_364/passed_htseq/")
packages <- c("edgeR", "plyr", "dplyr", "AnnotationDbi", "org.Hs.eg.db", "GeneOverlap", "ggpubr", "gplots", "ComplexHeatmap", "ggplot2", "gProfileR", "gprofiler2", "forcats","stringr")
lapply(packages, require, character.only = TRUE)

library(BatchQC)

# Read in file names
files <- list.files(pattern ="_gene_counts.txt$")

##define a function to get the prefix
get_prefix=function(file_name){
  file_prefix=unlist(strsplit(file_name, "\\_LY"))[1]
  return(file_prefix)
}

#ini the group arry
group=c()
for(i in files){
  cat (i,"\n")
  prefix=get_prefix(i)
  group=c(group,prefix)
}

  
# Generate DGElist object
df <- readDGE(files, header = FALSE, group = group)
df$samples # shows counts per library

#check num of 0 count samples, 364 passed QC samples in total 
table(rowSums(df$counts==0)==364)  #summary the genes that have 0 count in 365 

##filterByExpr function automaticaly filter low exps genes
keep <- filterByExpr(df)
table(keep)
df_keep <- df[keep, keep.lib.sizes=FALSE]

write.table(df_keep$counts, file="passed_df_counts.txt", quote=F, row.names = F)

setwd("/Users/tingliu/UHN/kridel-lab - FLOMICS/RNAseq/2022_combined_RNAseq_total365/htseq_out/")
metadf <- read.table("metadata_passed.txt",header=T)
exprdata<-read.table("passed_df_counts_mod.txt",header=T)


core=metadf$squence_core

batchQC(dat=as.matrix(exprdata), batch=core,
        report_file="batchqc_report.html", report_dir=".",
        report_option_binary="111111111",
        view_report=TRUE, interactive=TRUE, batchqc_output=TRUE)


##plot median_insert_size
my.p1 <- ggplot(data = all_qc_matrix[all_qc_matrix$squence_core!='TGL',], aes(x = squence_core))
my.p1 +  geom_boxplot(aes(y = median_insert_size, colour = squence_core))

boxplot(all_qc_matrix[all_qc_matrix$squence_core!='OICR',]$median_insert_size)

geom_boxplot(data=all_qc_matrix[all_qc_matrix$squence_core!='OICR',], aes(y = median_insert_size, colour = squence_core))




##sum the number of samples that mapped rate<-50
starlog[(starlog$P_Uniquely_mapped_reads<=50),]$sample
##48samples
