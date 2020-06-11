#local exploration of mutations across clusters and stages
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
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

#date
date=Sys.Date()

#directory with FLOMICS related matrices
setwd("/Users/kisaev/github/FLOMICS/Data")

#----------------------------------------------------------------------
#data filtering
#----------------------------------------------------------------------

muts <- fread("2020-06-08_opossum_variant_FL_rna-seq_filtered.txt")
qc = as.data.table(read_excel("FL_TGL_STAR_logQC_2020-05-13_summary_KI_ClusterContamAdded.xlsx"))

length(unique(muts$SAMPLE_ID))
muts = as.data.table(filter(muts, Tumor_Sample_Barcode %in% qc$rna_seq_file_sample_ID))
length(unique(muts$SAMPLE_ID))

#remove bad quality rnaseq samples
#Uniquely mapped reads >= 70 [Keep]
#Uniquely.mapped.reads..
#reads mapped to multiple loci < 20 [Keep]
#X..of.reads.mapped.to.multiple.loci
#Uniquely mapped > 30000000 [Keep] here try 10000000
#rrna contamination > 40 [Remove]

muts <- as.data.table(filter(muts, rrna_contam < 40, Uniquely.mapped > 10000000))
z = which(muts$'X..of.reads.mapped.to.multiple.loci' < 20)
muts = muts[z,]
z = which(muts$'Uniquely.mapped.reads..' > 60)
muts = muts[z,]

#how many patients here?
length(unique(muts$SAMPLE_ID)) #68

#remove groups of mutations we are not currently interested in

muts <- as.data.table(filter(muts, !(Variant_Classification %in% c("UTR3 .",
"downstream .", "UTR5 .", "upstream .", "exonic unknown"))))

length(unique(muts$SAMPLE_ID)) #66 patient samples

#----------------------------------------------------------------------
#analysis of mutations
#----------------------------------------------------------------------

#[1] how many mutations in each gene
genes_sum <- as.data.table(table(muts$Hugo_Symbol))
genes_sum <- genes_sum[order(-N)]
genes_sum = as.data.table(filter(genes_sum, N > 5)) #genes wtih at least 5 mutations in them

#aside from TTN which is a super long gene
#SON has the most mutations, n=119 - artifact? real?

#[2] which genes are most often mutated in one of two clusters
clusters <- as.data.table(table(muts$Hugo_Symbol, muts$CLUSTER_InfiniumClust))
colnames(clusters) <- c("gene", "cluster", "num_muts")
clusters <- clusters[order(-num_muts)]

#[3] for each gene compare proportion of patients in each cluster that have mutation

genes = unique(genes_sum$V1)

muts$CLUSTER_InfiniumClust = factor(muts$CLUSTER_InfiniumClust, levels=c(1,2))
#get num unique patients in each cluster
clusts_tot = as.data.table(table(unique(muts[,c("SAMPLE_ID",
"CLUSTER_InfiniumClust")])$CLUSTER_InfiniumClust))
colnames(clusts_tot)[2] = "tot_pats"

get_comp = function(gene){
	print(gene)
	gene_dat = as.data.table(filter(muts, Hugo_Symbol == gene))
	#only count gene being mutated once per patient
	gene_dat = unique(gene_dat[,c("SAMPLE_ID", "STAGE", "CLUSTER_InfiniumClust")])

	gene_dat$STAGE[which(str_detect(gene_dat$SAMPLE_ID, "RLN"))] = "RLN"
	gene_dat$STAGE[which(str_detect(gene_dat$SAMPLE_ID, "DLC"))] = "DLBCL"
	gene_dat$CLUSTER_InfiniumClust = factor(gene_dat$CLUSTER_InfiniumClust)

	#plot summary
	plotsum = as.data.table(table(gene_dat$CLUSTER_InfiniumClust))
	g = ggbarplot(plotsum, x="V1", y="N", fill="grey", title=paste("mutations in", gene)) + xlab("cluster") +
	ylab("Num of patient with at least 1 mutation")
	print(g)

	#make matrix for fishers test
	plotsum = merge(plotsum, clusts_tot)
	mat = matrix(ncol=2, nrow=2)
	mat[1:2,1] = plotsum$N
	mat[1:2,2] = clusts_tot$tot_pats - plotsum$N
	p = fisher.test(mat)$p.value
	or = fisher.test(mat)$estimate
	row = c(gene, p, or)
	names(row) = c("gene", "fishers_p", "fishers_or")
	return(row)

}

pdf(paste(date, "Gene_based_rnaseq_variants_summary.pdf", sep="_"))
summ = as.data.table(ldply(llply(genes, get_comp)))
dev.off()

summ$fishers_p <- as.numeric(summ$fishers_p)
summ = summ[order(fishers_p)]
summ$fishers_fdr <- p.adjust(summ$fishers_p, method="fdr")
