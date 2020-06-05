#----------------------------------------------------------------------
#variants_004_summary_visualize.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
	"ggrepel", "stringr", "maftools", "VariantAnnotation", "ggpubr")
lapply(packages, require, character.only = TRUE)
library(Biobase)
# Convert the row names to entrez ids
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

date = Sys.Date()

#do this locally - transfer files from cluster
#/cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder/SNP_matrices_all_algortihms

setwd("/Users/kisaev/FLOMICS_Anjali/FLOMICS_Anajli/Data")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?


#----------------------------------------------------------------------
#data 
#----------------------------------------------------------------------

all_variants = readRDS(list.files(pattern="all_algos_merged")[length(list.files(pattern="all_algos_merged"))])

#gene annotations
genes = unique(fread("ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#[1.] remove likely SNPs based on rs id and population frequency when available 
all_variants$avsnp142 = as.character(all_variants$avsnp142)
z = which(str_detect(all_variants$avsnp142, "rs"))

#get list of all unique variants that are likely SNPs 
snps = all_variants[z,]
snps_sum = as.data.table(table(snps$avsnp142)) ; snps_sum = snps_sum[order(-N)]

#variants likely not to be SNPs
all_variants = all_variants[-z,]
algo_sum = as.data.table(table(all_variants$algo)) ; algo_sum = algo_sum[order(-N)]
all_variants$mut_id = paste(all_variants$Chromosome, all_variants$Start_Position, all_variants$End_Position, sep="_")
all_variants$combo = paste(all_variants$Chromosome, all_variants$Start_Position, all_variants$End_Position,  
	all_variants$patient, sep="_")
all_variants = unique(all_variants)

all_variants$combo_2 = paste(all_variants$Chromosome, all_variants$Start_Position, all_variants$End_Position,  
	all_variants$patient, all_variants$algo,  sep="_")
#if patitent has multiple mappings for same mutation keep only one mapping 
z = which(duplicated(all_variants$combo_2))
all_variants = all_variants[-z,]
all_variants$score=1

#should be in MAF format now 
summary_muts = dcast(all_variants, combo  ~ algo)
summary_muts[is.na(summary_muts)] <- 0
summary_muts$tot_algos = apply(summary_muts[,2:ncol(summary_muts)], 1, sum)

summary_muts = summary_muts[order(-tot_algos)]
write.table(summary_muts, file=paste("mutations_across_algorithms_marix", date, ".txt", sep="_"), quote=F, row.names=F, sep="\t")








