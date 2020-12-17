#----------------------------------------------------------------------------------
#006_immunarch_diversity_indices.R
#----------------------------------------------------------------------------------

#Sarah Russell
#Date: Sept. 20, 2020
#This script takes raw MiXCR results files, loads into immunarch data structure
#and computes different diversity indices. See OneNote document in MiXCR Teams
#folder for explanation of different metrics.
#also added last part for getting vegan package diversity with added metadata

#ran locally, easier to download immunarch package locally (needs many dependencies)

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
library(vegan)
library(reshape2)
library(immunarch)

setwd("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Mixcr")

#----------------------------------------------------------------------------------
#DATA
#----------------------------------------------------------------------------------
#can load many seperate files if located in the same folder, will make a list of all datasets
immdata <- repLoad("/mixcr_data_immunarch/")
#passes through only coding clonotype, no filtering needed in this case

#filter for TRA/TRB V chains only, just looking at TCR set of data
tr_imm=immdata
x=''
for(i in 1:length(tr_imm$data)){
  tr_imm$data[[i]] <-tr_imm$data[[i]][grepl("TRA|TRB",tr_imm$data[[i]]$V.name),]
  x[i]=nrow(tr_imm$data[[i]]) > 0
  z=names(tr_imm$data[as.logical(x)])
  dat=tr_imm$data[as.logical(x)]
  met=filter(tr_imm$meta, Sample %in% z)
  repetoire=list(data=dat,meta=met)
}

#could also filter for IG chains only, just looking at BCR set of data
#br_imm=immdata
#x=''
#for(i in 1:length(br_imm$data)){
#  br_imm$data[[i]] <-br_imm$data[[i]][grepl("IG",br_imm$data[[i]]$V.name),]
#  x[i]=nrow(br_imm$data[[i]]) > 0
#  z=names(br_imm$data[as.logical(x)])
#  dat=br_imm$data[as.logical(x)]
#  met=(filter(bcr_metrics, sample %in% z)) %>% select(sample,stage,type,cluster,SITE_BIOPSY,TYPE_BIOPSY,GRADE) 
#  b_repetoire=list(data=dat,meta=met)
#}

#calculate diversity metrics
    div_chao <- as.data.frame(repDiversity(repetoire$data, "chao1"))
    div_chao$Sample=rownames(div_chao)
    div_chao=merge(div_chao,met,by="Sample")

    div_hill <- repDiversity(repetoire$data, "hill")

    div_div <- as.data.frame(repDiversity(repetoire$data, .col="aa", "div"))
    div_div=merge(div_div,met,by="Sample")

    div_gini.simp <- as.data.frame(repDiversity(repetoire$data, "gini.simp"))
    div_gini.simp=merge(div_gini.simp,met,by="Sample")

    div_gini <- as.data.frame(repDiversity(repetoire$data, .col="aa", "gini"))
    div_gini$Sample=rownames(div_gini)
    div_gini=merge(div_gini,met,by="Sample")

    div_inv.simp <- as.data.frame(repDiversity(repetoire$data, .col="aa", "inv.simp"))
    div_inv.simp=merge(div_inv.simp,met,by="Sample")

    div_d50 <- as.data.frame(repDiversity(repetoire$data, .col="aa", "d50"))
    div_d50$Sample=rownames(div_d50)
    div_d50=merge(div_d50,met,by="Sample")

all_div = list(div_chao=div_chao,div_div=div_div,div_gini.simp=div_gini.simp,div_gini=div_gini,
               div_inv.simp=div_inv.simp,div_d50=div_d50)

for(i in 1:length(all_div)){
  g1=ggboxplot(all_div[[i]], x = "TYPE", y = names(all_div[[i]][2]),
               title = "sample diversity estimation", ylab = names(all_div[i]), xlab = "Sample Type",
               color = "TYPE", palette = "jco")
  my_comparisons1 <- list( c("DLBCL", "FL"), c("FL", "RLN"), c("DLBCL", "RLN") )

  #######
  g2=ggboxplot(all_div[[i]][all_div[[i]]$TYPE=="FL",], x = "STAGE", y = names(all_div[[i]][2]),
               title = "sample diversity estimation", ylab = names(all_div[i]), xlab = "Sample Stage",
               color = "STAGE", palette = "jco")
  my_comparisons2 <- list( c("ADVANCED", "LIMITED"))

  ######
  g3=ggboxplot(all_div[[i]], x = "CLUSTER", y = names(all_div[[i]][2]),
               title = "sample diversity estimation", ylab = names(all_div[i]), xlab = "Sample Cluster",
               color = "CLUSTER", palette = "jco")
  my_comparisons3 <- list( c(1, 2))
  g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed")
  ######
  (g1 + stat_compare_means(comparisons = my_comparisons1) + grids(linetype = "dashed")) +
    (g2 + stat_compare_means(comparisons = my_comparisons2) + grids(linetype = "dashed")) +
    (g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed"))

  ggsave(filename = paste("/Users/sarahrussell/R/plots/",names(all_div[i]), "diversity.pdf",sep="_"), width=10, height=5)
}

p1 <- vis(div_hill, .by = c("TYPE"), .meta = repetoire$meta, .errorbars.off=T)
p2 <- vis(div_hill[div_hill$Sample %in% repetoire$meta[repetoire$meta$TYPE=="FL",]$Sample,], .by = c("STAGE"), .meta = repetoire$meta[repetoire$meta$TYPE=="FL",], .errorbars.off=T)
p3 <- vis(div_hill, .by = c("CLUSTER"), .meta = repetoire$meta, .errorbars.off=T)
p1 + p2 + p3
ggsave(filename = paste("/Users/sarahrussell/R/plots/", "hill_diversity.pdf",sep="_"), width=10, height=4)
################################################################################################################################################################################
#try sampling clones

#downsample would be similar to the sampling I had done with the vegan package
##basically expanding clone counts and choosing from a population of all clones, not
##all clonotypes - unclear if this is with or without replacement.. might sample with replacement

#sample chooses from clonotypes not clones, according to size if prob = T.

#should try downsampling, see what results look like, could also try sampling population similar to vegan package
#try n = 50,100,150

# Downsampling to 1000 clones (not clonotypes!)


x=''
y=''
z=''
for(i in 1:length(repetoire$data)){
  x[i]=sum(repetoire$data[[i]]$Clones) >= 50
  y[i]=sum(repetoire$data[[i]]$Clones) >= 100
  z[i]=sum(repetoire$data[[i]]$Clones) >= 150
}
pop=list(x,y,z)
get_samp_pops = function(pops){
  z=names(repetoire$data[as.logical(pops)])
  dat=repetoire$data[as.logical(pops)]
  met=filter(repetoire$meta, Sample %in% z)
  sample.pop=list(data=dat,meta=met)
  return(sample.pop)
}
sample_pops=llply(pop,get_samp_pops)

samp50 <- repSample(sample_pops[[1]]$data[c(1:length(sample_pops[[1]]$data))], .method = "downsample", .n=50)
#> length(sample_pops[[1]]$data) == length(samp50)  ##94
#TRUE
samp100 <- repSample(sample_pops[[2]]$data[c(1:length(sample_pops[[2]]$data))], .method = "downsample", .n=100)
#> length(sample_pops[[2]]$data) == length(samp100)   ##80
#[1] TRUE
samp150 <- repSample(sample_pops[[3]]$data[c(1:length(sample_pops[[3]]$data))], .method = "downsample", .n=150)
#> length(sample_pops[[3]]$data) == length(samp150)   ##60
#[1] TRUE
#is as expected based off previous sampling with vegan package


#calculate diversity metrics
get_sample_metrics = function(sampling,
                              metadata){

  div_chao <- as.data.frame(repDiversity(sampling, "chao1"))
  div_chao$Sample=rownames(div_chao)
  div_chao=merge(div_chao,metadata,by="Sample")

  #div_hill <- repDiversity(sampling, "hill")

  div_div <- as.data.frame(repDiversity(sampling, .col="aa", "div"))
  div_div=merge(div_div,metadata,by="Sample")

  div_gini.simp <- as.data.frame(repDiversity(sampling, "gini.simp"))
  div_gini.simp=merge(div_gini.simp,metadata,by="Sample")

  div_gini <- as.data.frame(repDiversity(sampling, .col="aa", "gini"))
  div_gini$Sample=rownames(div_gini)
  div_gini=merge(div_gini,metadata,by="Sample")

  div_inv.simp <- as.data.frame(repDiversity(sampling, .col="aa", "inv.simp"))
  div_inv.simp=merge(div_inv.simp,metadata,by="Sample")

  div_d50 <- as.data.frame(repDiversity(sampling, .col="aa", "d50"))
  div_d50$Sample=rownames(div_d50)
  div_d50=merge(div_d50,metadata,by="Sample")


  all_div = list(div_chao=div_chao,div_div=div_div,div_gini.simp=div_gini.simp,div_gini=div_gini,
                 div_inv.simp=div_inv.simp,div_d50=div_d50)
  return(all_div)
}
diversity.50 = get_sample_metrics(samp50,sample_pops[[1]]$meta)
diversity.100 = get_sample_metrics(samp100,sample_pops[[2]]$meta)
diversity.150 = get_sample_metrics(samp150,sample_pops[[3]]$meta)

#visualize sampling

get_sample_plots = function(samp_mets){

for(i in 1:length(samp_mets)){
  g1=ggboxplot(samp_mets[[i]], x = "TYPE", y = names(samp_mets[[i]][2]),
               title = "sample diversity estimation", ylab = names(samp_mets[i]), xlab = "Sample Type",
               color = "TYPE", palette = "jco")
  my_comparisons1 <- list( c("DLBCL", "FL"))

  #######
  g2=ggboxplot(samp_mets[[i]][samp_mets[[i]]$TYPE=="FL",], x = "STAGE", y = names(samp_mets[[i]][2]),
               title = "sample diversity estimation", ylab = names(samp_mets[i]), xlab = "Sample Stage",
               color = "STAGE", palette = "jco")
  my_comparisons2 <- list( c("ADVANCED", "LIMITED"))

  ######
  g3=ggboxplot(samp_mets[[i]], x = "CLUSTER", y = names(samp_mets[[i]][2]),
               title = "sample diversity estimation", ylab = names(samp_mets[i]), xlab = "Sample Cluster",
               color = "CLUSTER", palette = "jco")
  my_comparisons3 <- list( c(1, 2))
  ######
    (g1 + stat_compare_means(comparisons = my_comparisons1) + grids(linetype = "dashed")) +
    (g2 + stat_compare_means(comparisons = my_comparisons2) + grids(linetype = "dashed")) +
    (g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed"))

  ggsave(filename = paste("/Users/sarahrussell/R/plots/",names(samp_mets[i]), "diversity.pdf",sep="_"), width=10, height=5)
}
}
get_sample_plots(diversity.50)
get_sample_plots(diversity.100)
get_sample_plots(diversity.150)

get_hill_diversity = function(sampling,
                              metadata){
  div_hill <- repDiversity(sampling, "hill")
  p1 <- vis(div_hill, .by = c("TYPE"), .meta = repetoire$meta, .errorbars.off=T)
  p2 <- vis(div_hill[div_hill$Sample %in% repetoire$meta[repetoire$meta$TYPE=="FL",]$Sample,], .by = c("STAGE"), .meta = repetoire$meta[repetoire$meta$TYPE=="FL",], .errorbars.off=T)
  p3 <- vis(div_hill, .by = c("CLUSTER"), .meta = repetoire$meta, .errorbars.off=T)
  p1 + p2 + p3
  ggsave(filename = paste("/Users/sarahrussell/R/plots/", "hill_diversity.pdf",sep="_"), width=10, height=4)
}
get_hill_diversity(samp50,sample_pops[[1]]$meta)
get_hill_diversity(samp100,sample_pops[[2]]$meta)
get_hill_diversity(samp150,sample_pops[[3]]$meta)

################################################################################################################################################################################
#look at other metadata grouping


#qc=fread("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Mixcr/_2020-09-17_qc_annos.csv")

#format sample annotations
#sample_annos=fread("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/metadata/sample_annotations_rcd6Nov2019.csv")
#sample_annos=sample_annos %>% select(SAMPLE_ID,SITE_BIOPSY,TYPE_BIOPSY, LY_FL_ID)
#colnames(sample_annos)[1]="sample"
#a=filter(sample_annos,sample_annos$sample %in% qc$sample)
#b=filter(sample_annos, sample_annos$LY_FL_ID %in% qc$sample)
#colnames(b)[4]="sample"
#b[,1]=NULL
#a$LY_FL_ID=NULL
#c=rbind(a,b)

#clinical_annos=fread("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/metadata/clinical_data_rcd11Aug2020.csv")
#clinical_annos=clinical_annos %>% select(LY_FL_ID,GRADE)
#clinical_annos$sample=paste(clinical_annos$LY_FL_ID,"T1",sep="_")
#x=filter(clinical_annos,clinical_annos$sample %in% qc$sample)
#y=filter(clinical_annos, clinical_annos$LY_FL_ID %in% qc$sample)
#y[,3]=NULL
#colnames(y)[1]="sample"
#x[,1]=NULL
#z=rbind(x,y)

#qc_test=merge(qc,c,by="sample",all=T)
#qc_test=merge(qc_test,z,by="sample",all=T)
#qc=qc_test
#"SITE_BIOPSY"      "TYPE_BIOPSY"      "GRADE"


#immdata <- repLoad("/Users/sarahrussell/mixcr/")
#write.csv(qc_test,"_2020-09-18_qc_annos.csv")
qc=fread("/Users/sarahrussell/UHN/kridel-lab - FLOMICS/Analysis-Files/Mixcr/_2020-09-17_qc_annos.csv",header=T)
##################################

qc=qc %>% select(sample,SITE_BIOPSY,TYPE_BIOPSY,GRADE)
colnames(qc)[1]="Sample"
repetoire$meta=merge(repetoire$meta,qc,by="Sample")

#calculate diversity metrics
div_chao <- as.data.frame(repDiversity(repetoire$data, "chao1"))
div_chao$Sample=rownames(div_chao)
div_chao=merge(div_chao,repetoire$meta,by="Sample")

div_hill <- repDiversity(repetoire$data, "hill")

div_div <- as.data.frame(repDiversity(repetoire$data, .col="aa", "div"))
div_div=merge(div_div,repetoire$meta,by="Sample")

div_gini.simp <- as.data.frame(repDiversity(repetoire$data, "gini.simp"))
div_gini.simp=merge(div_gini.simp,repetoire$meta,by="Sample")

div_gini <- as.data.frame(repDiversity(repetoire$data, .col="aa", "gini"))
div_gini$Sample=rownames(div_gini)
div_gini=merge(div_gini,repetoire$meta,by="Sample")

div_inv.simp <- as.data.frame(repDiversity(repetoire$data, .col="aa", "inv.simp"))
div_inv.simp=merge(div_inv.simp,repetoire$meta,by="Sample")

div_d50 <- as.data.frame(repDiversity(repetoire$data, .col="aa", "d50"))
div_d50$Sample=rownames(div_d50)
div_d50=merge(div_d50,repetoire$meta,by="Sample")

all_div = list(div_chao=div_chao,div_div=div_div,div_gini.simp=div_gini.simp,div_gini=div_gini,
               div_inv.simp=div_inv.simp,div_d50=div_d50)

for(i in 1:length(all_div)){
  g1=ggboxplot(all_div[[i]][all_div[[i]]$TYPE=="FL",], x = "SITE_BIOPSY", y = names(all_div[[i]][2]),
               title = " ", ylab = names(all_div[i]), xlab = "Biopsy Site",
               color = "SITE_BIOPSY", palette = "jco")
  my_comparisons1 <- list( c("LN", "EN"))
  #######
  g2=ggboxplot(all_div[[i]][all_div[[i]]$TYPE %in% c("FL","DLBCL"),], x = "TYPE_BIOPSY", y = names(all_div[[i]][2]),
               title = " ", ylab = names(all_div[i]), xlab = "Biopsy Type",
               color = "TYPE_BIOPSY", palette = "jco", facet.by="TYPE")
  my_comparisons2 <- list( c("TISSUE", "CORE"))
  ######
  g3=ggboxplot(all_div[[i]][all_div[[i]]$TYPE=="FL" & all_div[[i]]$GRADE %in% c("1","2","3A"),], x = "GRADE", y = names(all_div[[i]][2]),
               title = " ", ylab = names(all_div[i]), xlab = "Sample Grade",
               color = "GRADE", palette = "jco")
  my_comparisons3 <- list( c("1", "2"),c("2","3A"),c("1","3A"))
  ######
  (g1 + stat_compare_means(comparisons = my_comparisons1) + grids(linetype = "dashed")) +
    (g2 + stat_compare_means(comparisons = my_comparisons2) + grids(linetype = "dashed")) +
    (g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed"))

  ggsave(filename = paste("/Users/sarahrussell/R/plots/",names(all_div[i]), "diversity.pdf",sep="_"), width=10, height=5)
}

#p1 <- vis(div_hill[div_hill$Sample %in% repetoire$meta[repetoire$meta$TYPE=="FL",]$Sample,], .by = c("SITE_BIOPSY"), .meta = repetoire$meta, .errorbars.off=T)
#p2 <- vis(div_hill, .by = c("TYPE_BIOPSY"), .meta = repetoire$meta[repetoire$meta$TYPE=="FL",], .errorbars.off=T)
#p3 <- vis(div_hill[div_hill$Sample %in% repetoire$meta[repetoire$meta$TYPE=="FL",]$Sample,], .by = c("GRADE"), .meta = repetoire$meta, .errorbars.off=T)
#p1 + p2 + p3
#ggsave(filename = paste("/Users/sarahrussell/R/plots/", "hill_diversity.pdf",sep="_"), width=10, height=4)
############################################
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
