#----------------------------------------------------------------------------------
#007_vegan_diversity_indices.R
#----------------------------------------------------------------------------------

#Sarah Russell
#Date: Aug 20, 2020
#This script takes concatenated MiXCR results and calculates different
#diversity metrics and MiXCR summary metrics using the "vegan" package.

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
library(reshape2)
library(vegan)

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

all_clones=data.frame(clonotypes$cloneCount,clonotypes$cloneFraction,clonotypes$allCHitsWithScore, clonotypes$allVHitsWithScore, clonotypes$aaSeqCDR3,clonotypes$sample,clonotypes$STAGE,clonotypes$TYPE,clonotypes$CLUSTER,stringsAsFactors=F)
colnames(all_clones)=c("count","fraction","chits", "vhits","aaSeqCDR3","sample","stage","type","cluster")

#remove sequences that have frameshift or stop codon
all_clones_filt <- all_clones[!grepl("\\~",all_clones$aaSeqCDR3),]
all_clones_filt <- all_clones_filt[!grepl("\\*",all_clones_filt$aaSeqCDR3),]
all_clones_filt <- all_clones_filt[!grepl("\\_",all_clones_filt$aaSeqCDR3),]

###add receptor chain
all_clones_filt$chain <- ifelse(grepl("TRA",all_clones_filt$vhits)==TRUE, "TRA",
  ifelse(grepl("TRB",all_clones_filt$vhits)==TRUE,"TRB",
  ifelse(grepl("TRD",all_clones_filt$vhits)==TRUE,"TRD",
  ifelse(grepl("TRG",all_clones_filt$vhits)==TRUE,"TRG",
  ifelse(grepl("IGH",all_clones_filt$vhits)==TRUE,"IGH",
  ifelse(grepl("IGK",all_clones_filt$vhits)==TRUE,"IGK","IGL"))))))

#----------------------------------------------------------------------------------
#get unique CDR3 counts
ID=as.character(unique(all_clones_filt$sample))
get_uniq_CDR=function(ID){
    count=length(unique(all_clones_filt[all_clones_filt$sample==ID,]$aaSeqCDR3))
    count_matrix=as.data.frame(cbind(count,ID))
    count_matrix$type <- ifelse(grepl("DLC",count_matrix$ID)==TRUE, "DLBLC",
      ifelse(grepl("FL",count_matrix$ID)==TRUE,"FL","RLN"))
    return(count_matrix)
}
counts = as.data.table(ldply(llply(ID, get_uniq_CDR)))

stage=as.character(all_clones_filt$stage[match(unique(all_clones_filt$sample),all_clones_filt$sample)])
cluster=as.character(all_clones_filt$cluster[match(unique(all_clones_filt$sample),all_clones_filt$sample)])
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

#plot number of reads vs unique clones
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
#----------------------------------------------------------------------------------
#filter for TCR sequences only
all_TR=filter(all_clones_filt, chain %in% c("TRA","TRB"))
#> dim(all_TR)
#[1] 7919   10


#could also filter for all BCR sequences only 
#all_BR=filter(all_clones_filt, chain %in% c("IGH","IGL","IGK))

#calculate USR (unqiue sequencing read) ; defined as a sequence read having no identity in TRV, TRJ
#and deduced amino acid sequence of CDR3 with the other sequence reads; Concatenate string of “TRV gene name”_”
#deduced amino acid sequence of CDR3 region”_” TRJ gene name” of individual USR (for example: TRBV1_CASTRVVJFG_TRBJ2-5)
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5059964/
#not necessary to do - no difference between considering unique USRs vs. unique aaSeqCDR3
#for(k in 1:nrow(all_TR)){
#  all_TR$USR[k]=paste(all_TR$chain[k],all_TR$aaSeqCDR3[k],sep="_")
#}

CDR3_counts=data.frame(acast(all_TR, sample~aaSeqCDR3, value.var="count", fun.aggregate=sum))

#calculate diversity metrics on all TCRs - no sampling
shannon=diversity(CDR3_counts, index="shannon")
invsimpson=diversity(CDR3_counts, index="invsimpson")

count_sums=rowSums(CDR3_counts)
indices=data.frame(shannon,invsimpson,count_sums)

top_fraction=all_TR[!duplicated(all_TR$sample),]
annotations=top_fraction %>% select(sample,stage,type,cluster,chain)
rownames(annotations)=annotations$sample

top_fraction=all_clones_filt[!duplicated(all_clones_filt$sample),]
annotations=top_fraction %>% select(sample,fraction,stage,type,cluster,chain)
rownames(annotations)=annotations$sample

rownames(annotations)==rownames(indices)
clones=cbind(indices, annotations[match(rownames(indices), rownames(annotations)),])

###can plot overall TCR diversity, relative fraction of top clonotype, top clonotype chain, etc.
#----------------------------------------------------------------------------------
#examine specific alpha & beta diversity

all_TR=filter(all_clones_filt, chain %in% c("TRA","TRB"))
alpha=all_TR[all_TR$chain=="TRA",]
beta=all_TR[all_TR$chain=="TRB",]

all_alpha=data.frame(acast(alpha, sample~USR, value.var="count", fun.aggregate=sum))
all_beta=data.frame(acast(beta, sample~USR, value.var="count", fun.aggregate=sum))

#----------------------------------------------------------------------------------
#look at CDR3 counts in Samples

tcrs=all_TR %>% select(count,sample,aaSeqCDR3)
tcr_counts = aggregate(tcrs[,1],tcrs[,-1],sum)

colnames(tcr_counts)[3]="count"
tcr_counts=data.frame(tcr_counts[order(tcr_counts$sample),])

#get total counts (number of clones) per unique aaSeqCDR3 per sample
s_tcr=tcr_counts$sample
get_prob = function(sample){
  x=tcr_counts
  totalc=sum(x[x$sample==sample,]$count)
  return(totalc)
}
all_reads=unlist(llply(s_tcr,get_prob))

#get clonal fraction of those clones
tcr_counts$prob=(tcr_counts$count/all_reads)
tcr_counts$treads=all_reads

#----------------------------------------------------------------------------------
#before processing in R:
##cd /cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS
#ls *.report | grep -oP '.*?(?=\.)' > all_report_names.txt
#cat *.report | grep 'Total sequencing reads:.*' > total_reads.txt
#cat *.report | grep 'Average number of reads per clonotype:.*' | awk '{print $7}' > av_reads_clones.txt
#cat *.report | grep 'Successfully aligned reads:.*' | awk '{print $4}' > suc_reads.txt
#cat *.report | grep 'Final clonotype count:.*' | awk '{print $4}' > total_clones.txt
#cat *.report | grep -Po 'LY_.*_R1_sorted\.fastq\.gz\,|TR(A|B) chains:.*' > tcr_qc.txt


#cat *.report | grep 'TRA chains:.*' | awk '{print $3}' | sed -n '1~2p' > tra_reads.txt
#cat *.report | grep -Po '(LY_[A-Z]*_[0-9]*)|(_T[1|2]|TRA chains:.*)'
#cat *.report | grep 'TRB chains:.*' | awk '{print $3}' | sed -n '1~2p' > trb_reads.txt
#cat *.report | grep 'TRA chains:.*' | awk '{print $3}' | sed -n '2~2p' > tra_clones.txt
#cat *.report | grep 'TRB chains:.*' | awk '{print $3}' | sed -n '2~2p' > trb_clones.txt
#cat *.report| grep -E 'Total sequencing reads:.*|Successfully aligned reads:.*|((TR(A|B|D|G))|(IG(H|K|L))) chains:.*|Final clonotype count:.*|Average number of reads per clonotype:.' > allqc.txt

sample=fread("all_report_names.txt",header=F)
treads=fread("total_reads.txt",header=F)
avreads=fread("av_reads_clones.txt",header=F)
tclones=fread("total_clones.txt",header=F)
tcr_qc=fread("trc_qc.csv",header=F)
colnames(tcr_qc)=c("sample","tra_reads","trv_reads","tra_clones","trb_clones")

qc=cbind(sample,treads,avreads,tclones)
colnames(qc)=c("sample","total_reads","av_reads_clone","total_clones")
qc=merge(qc,tcr_qc,by="sample")

#clonotypes=fread("_2020-07-30_MIXCR_OUTPUT_WITH_SAMPLE_ANNOTATION.csv", sep="\t")
qc_keep=filter(qc,sample %in% unique(clonotypes$sample))

#count CDR3 sequences (abundance)
sample=as.character(unique(all_clones_filt$sample))
get_uniq_CDR=function(ID){
    count=length(unique(all_clones_filt[all_clones_filt$sample==ID,]$aaSeqCDR3))
    count_matrix=as.data.frame(cbind(ID,count))
    return(count_matrix)
}
counts = as.data.table(ldply(llply(sample, get_uniq_CDR)))
colnames(counts)[1]="sample"
colnames(counts)[2]="unique_CDR3"
counts$unique_CDR3=as.character(counts$unique_CDR3)
counts$unique_CDR3=as.numeric(counts$unique_CDR3)

qc_annos=merge(qc_annos,counts,by="sample")

#get number of unique CDR3 TCRs
sample=as.character(unique(all_TR$sample))
get_uniq_CDR=function(ID){
    count=length(unique(all_TR[all_TR$sample==ID,]$aaSeqCDR3))
    count_matrix=as.data.frame(cbind(ID,count))
    return(count_matrix)
}
counts = as.data.table(ldply(llply(sample, get_uniq_CDR)))
colnames(counts)[1]="sample"
colnames(counts)[2]="TCRunique_CDR3"
counts$TCRunique_CDR3=as.character(counts$TCRunique_CDR3)
counts$TCRunique_CDR3=as.numeric(counts$TCRunique_CDR3)
qc_annos=merge(qc_annos,counts,by="sample",all=T)
qc_annos$TCRunique_CDR3[is.na(qc_annos$TCRunique_CDR3)] <- 0

annos=all_clones_filt[!duplicated(all_clones_filt$sample),] %>% select(sample,stage,type,cluster)
qc_annos=merge(qc_keep,annos,by='sample')

#matrix of different qc metrics for TCR

#----------------------------------------------------------------------------------
#sampling of TCR CDR3 sequences

tcr_all_samples = tcr_counts %>% slice(rep(1:n(), tcr_counts$count))
count_50=filter(tcr_all_samples, treads >= 50)
#94 samples
count_100=filter(tcr_all_samples, treads >= 100)
#80 samples
count_150=filter(tcr_all_samples, treads >= 150)
#60 samples

#for whatever sampling population
ID=unique(count_50$sample)
sample_cdr3 = function(ids){
  ID=unique(count_50$sample)
  x=count_50
  y = sample_n(x[x$sample==ids,], size = 50, replace = FALSE)
  return(y)
}
samp50=as.data.table(ldply(llply(ID,sample_cdr3)))

x=table(samp50$sample,samp50$aaSeqCDR3)

shannon_50=diversity(x, index="shannon")
invsimpson_50=diversity(all50, index="invsimpson")
indices_50=data.frame(shannon_50,invsimpson_50)
indices_50$sample=rownames(indices_50)

qc=as.data.frame(merge(qc,indices_50,by="sample",all=T))

#write.csv(qc, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS/", date, "qc_annos.csv",sep="_"), quote=F, row.names=F)
##################################################################
#get CDR3 richness for top 25% of clones

#all_TR=filter(all_clones_filt, chain %in% c("TRA","TRB"))
#qc=fread("_2020-08-21_qc_annos.csv")

IDs=unique(all_TR$sample)
get_top = function(sample_name){
x=all_TR[all_TR$sample==sample_name,]
  for(i in 1:nrow(x)){
  if(i==1){
    x$perc[i]=x$count[i]
  } else {
    x$perc[i]=x$count[i]+x$perc[i-1]
  }
    x$frac[i]=x$perc[i]/sum(x$count)
  }
  x=x[x$frac <= 0.25,]
  return(x)
  }
tr_samp=as.data.table(ldply(llply(IDs,get_top)))
tr_samp$perc=NULL
tr_samp$frac=NULL

k=filter(all_TR, !sample %in% tr_samp$sample)
tr_top=rbind(tr_samp,k)

IDs=unique(tr_top$sample)
get_counts = function(sample_name){
  x=sample_name
  y=nrow(tr_top[tr_top$sample==sample_name,])
  z=data.table(sample=x,count=y)
  return(z)
}
tr_samp=as.data.table(ldply(llply(IDs,get_counts)))

annos=all_TR[!duplicated(all_TR$sample),] %>% select(sample,stage,type,cluster)
tr_samp_anno=merge(tr_samp,annos,by="sample")

##################################################################
all_TR=filter(all_clones_filt, chain %in% c("TRA","TRB"))
qc=fread("_2020-08-21_qc_annos.csv")

#Pielou's eveness = shannon diversity/log(# of unique species)
#will have to perform sampling again

tcrs=all_TR %>% select(count,sample,aaSeqCDR3)
tcr_counts = aggregate(tcrs[,1],tcrs[,-1],sum)
colnames(tcr_counts)[3]="count"
tcr_counts=data.frame(tcr_counts[order(tcr_counts$sample),])

#get total counts (number of clones) per unique aaSeqCDR3 per sample
s_tcr=tcr_counts$sample
get_prob = function(sample){
  x=tcr_counts
  totalc=sum(x[x$sample==sample,]$count)
  return(totalc)
}
all_reads=unlist(llply(s_tcr,get_prob))

#get clonal fraction of those clones
tcr_counts$prob=(tcr_counts$count/all_reads)
tcr_counts$treads=all_reads

tcr_all_samples = tcr_counts %>% slice(rep(1:n(), tcr_counts$count))
count_50=filter(tcr_all_samples, treads >= 50)
#94 samples
count_100=filter(tcr_all_samples, treads >= 100)
#80 samples
count_150=filter(tcr_all_samples, treads >= 150)
#60 samples

#for whatever sampling population
ID=unique(count_150$sample)
sample_cdr3 = function(ids){
  x=count_150
  y = sample_n(x[x$sample==ids,], size = 150, replace = FALSE)
  return(y)
}
samp150=as.data.table(ldply(llply(ID,sample_cdr3)))

x1=table(samp50$sample,samp50$aaSeqCDR3)
x2=table(samp100$sample,samp100$aaSeqCDR3)
x3=table(samp150$sample,samp150$aaSeqCDR3)

shannon_50=diversity(x1, index="shannon")
invsimpson_50=diversity(x1, index="invsimpson")
pielou_50=shannon_50/log(specnumber(x1))
indices_50=data.frame(shannon_50,invsimpson_50,pielou_50)
indices_50$sample=rownames(indices_50)
indices_50=as.data.table(indices_50)

shannon_100=diversity(x2, index="shannon")
invsimpson_100=diversity(x2, index="invsimpson")
pielou_100=shannon_100/log(specnumber(x2))
indices_100=data.frame(shannon_100,invsimpson_100,pielou_100)
indices_100$sample=rownames(indices_100)
indices_100=as.data.table(indices_100)


shannon_150=diversity(x3, index="shannon")
invsimpson_150=diversity(x3, index="invsimpson")
pielou_150=shannon_150/log(specnumber(x3))
indices_150=data.frame(shannon_150,invsimpson_150,pielou_150)
indices_150$sample=rownames(indices_150)
indices_150=as.data.table(indices_150)

qc_all=as.data.frame(merge(qc,indices_50,by="sample",all=T))
qc_all=as.data.frame(merge(qc_all,indices_100,by="sample",all=T))
qc_all=as.data.frame(merge(qc_all,indices_150,by="sample",all=T))



pdf("/cluster/home/srussell/simp_stage.pdf", width=7, height=4)
g=ggboxplot(y[y$type=="FL" & grepl("invsimpson",y$metric),], x = "stage", y = "div",
          title = "Inverse Simpson TCR CDR3 Diversity", ylab = "Inverse Simpson Index", xlab = "Sample Stage",
          color = "stage", palette = "jco", facet.by = "sampling"
          )
g + stat_compare_means() + grids(linetype = "dashed")
dev.off()

pdf("/cluster/home/srussell/pie_stage.pdf", width=7, height=4)
g=ggboxplot(y[y$type=="FL" & grepl("pielou",y$metric) & y$div > 0.9,], x = "stage", y = "div",
          title = "Pielou's Eveness - TCR CDR3 Diversity", ylab = "Pielou's Eveness Index", xlab = "Sample Stage",
          color = "stage", palette = "jco", facet.by = "sampling"
          )
g + stat_compare_means() + grids(linetype = "dashed")
dev.off()


#for whatever sampling population
ID=unique(count_50$sample)
sample_cdr3 = function(ids){
  x=count_50
  y = sample_n(x[x$sample==ids,], size = 50, replace = FALSE)
  return(y)
}
samp50=as.data.table(ldply(llply(ID,sample_cdr3)))

#for whatever sampling population
ID=unique(count_100$sample)
sample_cdr3 = function(ids){
  x=count_100
  y = sample_n(x[x$sample==ids,], size = 100, replace = FALSE)
  return(y)
}
samp100=as.data.table(ldply(llply(ID,sample_cdr3)))

#for whatever sampling population
ID=unique(count_150$sample)
sample_cdr3 = function(ids){
  x=count_150
  y = sample_n(x[x$sample==ids,], size = 150, replace = FALSE)
  return(y)
}
samp150=as.data.table(ldply(llply(ID,sample_cdr3)))



##############
#pielou

all_TR
CDR3_counts
shannon=diversity(CDR3_counts, index="shannon")
pielou=shannon/log(specnumber(CDR3_counts))
indices=data.frame(shannon,pielou)

indices$sample=rownames(indices)
indices=as.data.table(indices)

tcr_qc=fread("_2020-08-21_qc_annos.csv",header=T)
qc_all=as.data.frame(merge(tcr_qc,indices,by="sample",all=T))


pdf("/cluster/home/srussell/pielou_type.pdf", width=7, height=4)
g1=ggboxplot(qc_all[qc_all$type!="RLN",], x = "type", y = "pielou",
          title=" ", ylab = "Pielou's Evenness Index", xlab = "Sample Type",
          color = "type", palette = "jco",
          )
my_comparisons1 <- list( c("DLBCL", "FL"))
(g1 + stat_compare_means(comparisons = my_comparisons1) + grids(linetype = "dashed"))
dev.off()
#######
pdf("/cluster/home/srussell/pielou_stage.pdf", width=7, height=4)
g=ggboxplot(qc_all[qc_all$type=="FL",], x = "stage", y = "pielou",
          tile=" ",  ylab = "Pielou's Eveness Index", xlab = "Sample Stage",
          color = "stage", palette = "jco",
          )
my_comparisons2 <- list( c("ADVANCED", "LIMITED"))
(g + stat_compare_means(comparisons = my_comparisons2) + grids(linetype = "dashed"))
dev.off()
#######
pdf("/cluster/home/srussell/pielou_cluster.pdf", width=7, height=4)
g3=ggboxplot(qc_all, x = "cluster", y = "pielou",
          tile = " ", ylab = "Pielou's Eveness Index", xlab = "Sample Cluster",
          color = "cluster", palette = "jco",
          )
my_comparisons3 <- list( c(1, 2))
(g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed"))
dev.off()

######
(g1 + stat_compare_means(comparisons = my_comparisons1) + grids(linetype = "dashed")) +
(g2 + stat_compare_means(comparisons = my_comparisons2) + grids(linetype = "dashed")) +


#for whatever sampling population
ID=unique(count_50$sample)
sample_cdr3 = function(ids){
  x=count_50
  y = sample_n(x[x$sample==ids,], size = 50, replace = FALSE)
  return(y)
}
samp50=as.data.table(ldply(llply(ID,sample_cdr3)))

#for whatever sampling population
ID=unique(count_100$sample)
sample_cdr3 = function(ids){
  x=count_100
  y = sample_n(x[x$sample==ids,], size = 100, replace = FALSE)
  return(y)
}
samp100=as.data.table(ldply(llply(ID,sample_cdr3)))

#for whatever sampling population
ID=unique(count_150$sample)
sample_cdr3 = function(ids){
  x=count_150
  y = sample_n(x[x$sample==ids,], size = 150, replace = FALSE)
  return(y)
}
samp150=as.data.table(ldply(llply(ID,sample_cdr3)))





#look at other metadata groupings
#calculate metrics on new meta data with vegan package, similar to before

clonotypes=fread("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Mixcr/_2020-07-30_MIXCR_OUTPUT_WITH_SAMPLE_ANNOTATION.csv", sep="\t")
all_clones=data.frame(clonotypes$cloneCount,clonotypes$cloneFraction,clonotypes$allCHitsWithScore, clonotypes$allVHitsWithScore, clonotypes$aaSeqCDR3,clonotypes$sample,clonotypes$STAGE,clonotypes$TYPE,clonotypes$CLUSTER,stringsAsFactors=F)
colnames(all_clones)=c("count","fraction","chits", "vhits","aaSeqCDR3","sample","stage","type","cluster")

#remove sequences that have frameshift or stop codon
all_clones_filt <- all_clones[!grepl("\\~",all_clones$aaSeqCDR3),]
all_clones_filt <- all_clones_filt[!grepl("\\*",all_clones_filt$aaSeqCDR3),]
all_clones_filt <- all_clones_filt[!grepl("\\_",all_clones_filt$aaSeqCDR3),]

###add receptor chain
all_clones_filt$chain <- ifelse(grepl("TRA",all_clones_filt$vhits)==TRUE, "TRA",
                                ifelse(grepl("TRB",all_clones_filt$vhits)==TRUE,"TRB",
                                       ifelse(grepl("TRD",all_clones_filt$vhits)==TRUE,"TRD",
                                              ifelse(grepl("TRG",all_clones_filt$vhits)==TRUE,"TRG",
                                                     ifelse(grepl("IGH",all_clones_filt$vhits)==TRUE,"IGH",
                                                            ifelse(grepl("IGK",all_clones_filt$vhits)==TRUE,"IGK","IGL"))))))
all_TR=filter(all_clones_filt, chain %in% c("TRA","TRB"))

#could also filter for BCR chains
#all_BR=filter(all_clones_filt, chain %in% c("IGH","IGK","IGL"))


CDR3_counts=data.frame(acast(all_TR, sample~aaSeqCDR3, value.var="count", fun.aggregate=sum))

#calculate diversity metrics on all TCRs - no sampling
shannon=diversity(CDR3_counts, index="shannon")
invsimpson=diversity(CDR3_counts, index="invsimpson")
pielou=shannon/log(specnumber(CDR3_counts))

indices=data.frame(shannon,invsimpson,pielou)
indices$sample=rownames(indices)
indices=as.data.table(indices)

qc=qc %>% select(sample,stage,type,cluster,SITE_BIOPSY,TYPE_BIOPSY,GRADE)
indices=merge(indices,qc,by="sample")

metrics=c("shannon","invsimpson","pielou")
get_diversity = function(met){
g1=ggboxplot(indices[indices$type=="FL",], x = "SITE_BIOPSY", y = met,
             title = " ", ylab = met, xlab = "Biopsy Site",
             color = "SITE_BIOPSY", palette = "jco")
my_comparisons1 <- list( c("LN", "EN"))

g2=ggboxplot(indices[indices$type %in% c("FL","DLBCL"),], x = "TYPE_BIOPSY", y = met,
             title = " ", ylab = met, xlab = "Biopsy Type",
             color = "TYPE_BIOPSY", palette = "jco", facet.by="type")
my_comparisons2 <- list( c("TISSUE", "CORE"))

g3=ggboxplot(indices[indices$type=="FL" & indices$GRADE %in% c("1","2","3A"),], x = "GRADE", y = met,
             title = " ", ylab = met, xlab = "Sample Grade",
             color = "GRADE", palette = "jco")
my_comparisons3 <- list( c("1", "2"),c("2","3A"),c("1","3A"))
######
(g1 + stat_compare_means(comparisons = my_comparisons1) + grids(linetype = "dashed")) +
(g2 + stat_compare_means(comparisons = my_comparisons2) + grids(linetype = "dashed")) +
(g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed"))

ggsave(filename = paste("/Users/sarahrussell/R/plots/",met, "diversity.pdf",sep="_"), width=10, height=5)
}
llply(metrics,get_diversity)
############################################
tcrs=all_TR %>% select(count,sample,aaSeqCDR3)
tcr_counts = aggregate(tcrs[,1],tcrs[,-1],sum)
colnames(tcr_counts)[3]="count"
tcr_counts=data.frame(tcr_counts[order(tcr_counts$sample),])

#get total counts (number of clones) per unique aaSeqCDR3 per sample
s_tcr=tcr_counts$sample
get_prob = function(sample){
  x=tcr_counts
  totalc=sum(x[x$sample==sample,]$count)
  return(totalc)
}
all_reads=unlist(llply(s_tcr,get_prob))

#get clonal fraction of those clones
tcr_counts$prob=(tcr_counts$count/all_reads)
tcr_counts$treads=all_reads

tcr_all_samples = tcr_counts %>% slice(rep(1:n(), tcr_counts$count))
count_50=filter(tcr_all_samples, treads >= 50)
#94 samples
count_100=filter(tcr_all_samples, treads >= 100)
#80 samples
count_150=filter(tcr_all_samples, treads >= 150)
#60 samples
counts=list("50"=count_50,"100"=count_100,"150"=count_150)

ID=unique(count_50$sample)
sample_cdr3 = function(ids){
  x=count_50
  y = sample_n(x[x$sample==ids,], size = 50, replace = FALSE)
  return(y)
}
samp50=as.data.table(ldply(llply(ID,sample_cdr3)))

x1=table(samp50$sample,samp50$aaSeqCDR3)
x2=table(samp100$sample,samp100$aaSeqCDR3)
x3=table(samp150$sample,samp150$aaSeqCDR3)

shannon_50=diversity(x1, index="shannon")
invsimpson_50=diversity(x1, index="invsimpson")
indices_50=data.frame(shannon_50,invsimpson_50)
indices_50$sample=rownames(indices_50)
indices_50=as.data.table(indices_50)

shannon_100=diversity(x2, index="shannon")
invsimpson_100=diversity(x2, index="invsimpson")
indices_100=data.frame(shannon_100,invsimpson_100)
indices_100$sample=rownames(indices_100)
indices_100=as.data.table(indices_100)


shannon_150=diversity(x3, index="shannon")
invsimpson_150=diversity(x3, index="invsimpson")
indices_150=data.frame(shannon_150,invsimpson_150)
indices_150$sample=rownames(indices_150)
indices_150=as.data.table(indices_150)

qc_all=as.data.frame(merge(qc,indices_50,by="sample",all=T))
qc_all=as.data.frame(merge(qc_all,indices_100,by="sample",all=T))
qc_all=as.data.frame(merge(qc_all,indices_150,by="sample",all=T))

metrics=c("shannon_50","invsimpson_50","shannon_100","invsimpson_100","shannon_150","invsimpson_150")
get_diversity = function(met,
                         samp){
  g1=ggboxplot(qc_all[qc_all$type=="FL",], x = "SITE_BIOPSY", y = met,
               title = " ", ylab = met, xlab = "Biopsy Site",
               color = "SITE_BIOPSY", palette = "jco")
  my_comparisons1 <- list( c("LN", "EN"))

  g2=ggboxplot(qc_all[qc_all$type %in% c("FL","DLBCL"),], x = "TYPE_BIOPSY", y = met,
               title = " ", ylab = met, xlab = "Biopsy Type",
               color = "TYPE_BIOPSY", palette = "jco", facet.by = "type")
  my_comparisons2 <- list( c("TISSUE", "CORE"))

  g3=ggboxplot(qc_all[qc_all$type=="FL" & qc_all$GRADE %in% c("1","2","3A"),], x = "GRADE", y = met,
               title = " ", ylab = met, xlab = "Sample Grade",
               color = "GRADE", palette = "jco")
  my_comparisons3 <- list( c("1", "2"),c("2","3A"),c("1","3A"))
  ######
    (g1 + stat_compare_means(comparisons = my_comparisons1) + grids(linetype = "dashed")) +
    (g2 + stat_compare_means(comparisons = my_comparisons2) + grids(linetype = "dashed")) +
    (g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed"))

  ggsave(filename = paste("/Users/sarahrussell/R/plots/",met, "diversity.pdf",sep="_"), width=10, height=5)
}
llply(metrics,get_diversity)

