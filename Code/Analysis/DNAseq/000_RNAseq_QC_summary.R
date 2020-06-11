#local exploration of mutations across clusters and stages
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors = F)
#avoid scientific notation
options(scipen=999)
#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr",
"mclust", "data.table", "plyr",
"ggrepel", "stringr", "maftools", "ggpubr", "readxl", "skimr",
 "edgeR", "annotables", "EnvStats")
library(gridExtra)
lapply(packages, require, character.only = TRUE)

#directory with FLOMICS related matrices
setwd("/Users/kisaev/github/FLOMICS/Data")

#----------------------------------------------------------------------
#data filtering
#----------------------------------------------------------------------

qc = as.data.table(read_excel("FL_TGL_STAR_logQC_2020-05-13_summary_KI_ClusterContamAdded.xlsx"))

#generate report summarizing key variables

vars = c("rrna_contam",
"Average.input.read.length",
"Average.mapped.length",
"Number.of.input.reads",
"X..of.reads.mapped.to.multiple.loci",
"Uniquely.mapped",
"Uniquely.mapped.reads..",
"Unmapped.reads", "X..of.reads.unmapped..too.short")

qc$STAGE[which(str_detect(qc$rna_seq_file_sample_ID, "RLN"))] = "RLN"
qc$STAGE[which(str_detect(qc$rna_seq_file_sample_ID, "DLC"))] = "DLBCL"
qc$CLUSTER_InfiniumClust = factor(qc$CLUSTER_InfiniumClust)

#for each variable generate summary stats
get_sum = function(var){
	print(var)
	z = which(colnames(qc) %in% c(var, "STAGE", "CLUSTER_InfiniumClust", "rna_seq_file_sample_ID"))
	dat = qc[,..z]
	colnames(dat)[which(colnames(dat)==var)] = "Variable"

	#histogram stage
	#hist1 = ggdensity(dat, x = "Variable", rug=TRUE, color="STAGE",
  #           add = "median", fill="STAGE", title=paste("Distribution of" , var))+
	#					 xlab(var)+theme_minimal()
	#print(hist1)


	#ggboxplot stage
	box1 = ggboxplot(dat, x = "STAGE", y="Variable", add ="jitter",
             fill="STAGE", title=paste("Distribution of" , var))+
						 ylab(var) + theme_minimal() + stat_n_text()+
						 stat_compare_means()
  print(box1)

	#histogram clusters
	#hist2 = ggdensity(dat, x = "Variable", rug=TRUE,color="CLUSTER_InfiniumClust",
  #           add = "median", fill="CLUSTER_InfiniumClust", title=paste("Distribution of" , var))+
	#					 xlab(var)+theme_minimal()
	#print(hist2)

	#ggboxplot clusters
	box2 = ggboxplot(dat, x = "CLUSTER_InfiniumClust", y="Variable", add="jitter",
             fill="CLUSTER_InfiniumClust", title=paste("Distribution of" , var))+
						 ylab(var) + theme_minimal() + stat_n_text()+
						 stat_compare_means()
  print(box2)

	#get summary stats
	summ = as.data.table(dat %>%
	dplyr::group_by(STAGE) %>%
	dplyr::summarise(
  	mean = mean(Variable, na.rm = TRUE),
		median = median(Variable, na.rm = TRUE),
  	sd = sd(Variable, na.rm = TRUE)))
	summ$variable = var
	p<-tableGrob(summ)
	print(grid.arrange(p))

}

pdf("FLOMICS_RNAseq_QC_STATS_by_STAGE_methyCLUSTER.pdf")
sapply(vars, get_sum)
dev.off()

#summarize how many patients lose from each filter
z1 = which(qc$rrna_contam > 40) # 3
z2 = which(qc$Uniquely.mapped < 10000000) # 22
z3 = which(qc$'Uniquely.mapped.reads..' < 60) # 60
z4 = which(qc$'X..of.reads.mapped.to.multiple.loci' > 20) # 15

total_lost = unique(c(z1, z2, z3, z4)) #63
#low unique mapped reads
low=qc[z3,]

pdf("FL_RNAseq_low_unique_mapped_reads_percentage.pdf")
hist(low$X..of.reads.unmapped..too.short)
dev.off()

qc <- as.data.table(filter(qc, rrna_contam < 40, Uniquely.mapped > 10000000))
z = which(qc$'X..of.reads.mapped.to.multiple.loci' < 20)
qc = qc[z,]
z = which(qc$'Uniquely.mapped.reads..' > 60)
qc = qc[z,] #69
print(dim(qc))
