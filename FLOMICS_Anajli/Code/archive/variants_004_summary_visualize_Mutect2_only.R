#----------------------------------------------------------------------
#variants_004_summary_visualize_Mutect2_only.R
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
library(openxlsx)

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

all_variants = readRDS(list.files(pattern="all_mutect2_vcfs_merged")[length(list.files(pattern="all_mutect2_vcfs_merged"))])

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
#if patient has multiple mappings for same mutation keep only one mapping 
#z = which(duplicated(all_variants$combo_2))
#all_variants = all_variants[-z,]
#all_variants$score=1

#id mapping
part1 = read.xlsx("part1_mapping_bc_dna.xlsx")
part2 = read.xlsx("part2_mapping_bc_dna.xlsx")
all_parts= rbind(part1, part2)
colnames(all_parts)[8] = "patient"
all_parts = unique(all_parts[,c("patient", "External_ID")])
dim(all_variants)
all_variants = merge(all_variants, all_parts, by="patient")
all_variants = as.data.table(filter(all_variants))

summary_muts = summary_muts[order(-tot_algos)]
write.table(summary_muts, file=paste("mutations_across_algorithms_marix", date, ".txt", sep="_"), quote=F, row.names=F, sep="\t")

#summary of AF and depth 
all_variants$DP_other_numeric = as.numeric(all_variants$DP_other)
all_variants$AF_numeric = as.numeric(all_variants$AF)

#multiallelic clean up
all_variants$multiallelic = paste(all_variants$patient, all_variants$Gene.refGene, all_variants$GT, all_variants$AD, all_variants$AF, all_variants$DP_other, sep="_")
multis = unique(all_variants$multiallelic[which(is.na(all_variants$AF_numeric))])
all_variants$multiallelic_tag = ""

#clean up multtis 
clean_up = function(multi){
	dat = as.data.table(filter(all_variants, multiallelic == multi))
	dat = dat[nrow(dat),]
	dat$AF = as.numeric(unlist(strsplit(dat$AF, ","))[length(unlist(strsplit(dat$AF, ",")))])
	dat$multiallelic_tag = "yes"
	return(dat)
}
all_multis = as.data.table(ldply(llply(multis, clean_up)))
z = which(all_variants$multiallelic %in% all_multis$multiallelic)
all_variants = all_variants[-z,]
all_variants = rbind(all_variants, all_multis)

all_variants$score_new_new = 1

mutect2_matrix = dcast(all_variants, mut_id  ~ External_ID)
mutect2_matrix[is.na(mutect2_matrix)] <- 0

write.table(mutect2_matrix, file=paste("mutect2_matrix", date, ".txt", sep="_"), quote=F, row.names=F, sep="\t")
write.table(all_variants, file=paste("mutect2_full_data_w_AF_DP", date, ".txt", sep="_"), quote=F, row.names=F, sep="\t")

all_variants_filtered = as.data.table(filter(all_variants, AF > 0.05, DP_other_numeric > 20))
write.table(all_variants_filtered, file=paste("mutect2_full_data", date, ".maf", sep="_"), quote=F, row.names=F, sep="\t")

#read in maf object 
fl = read.maf(maf = paste("mutect2_full_data", date, ".maf", sep="_"))
pdf("maftools_summary_plots_SNVs_Indels_Mutect2_FL_min_AF_0.05.pdf", width=14)
plotmafSummary(maf = fl, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf = fl, top = 20)
laml.titv = titv(maf = fl, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)
lollipopPlot(maf = fl, gene = 'CREBBP', AACol = 'AAChange.refGene', showMutationRate = TRUE)
lollipopPlot(maf = fl, gene = 'EZH2', AACol = 'AAChange.refGene', showMutationRate = TRUE)
lollipopPlot(maf = fl, gene = 'KMT2D', AACol = 'AAChange.refGene', showMutationRate = TRUE)
dev.off()

all_variants_filtered = as.data.table(filter(all_variants, AF > 0.1, DP_other_numeric > 20))
write.table(all_variants_filtered, file=paste("mutect2_full_data", date, ".maf", sep="_"), quote=F, row.names=F, sep="\t")

#read in maf object 
fl = read.maf(maf = paste("mutect2_full_data", date, ".maf", sep="_"))
pdf("maftools_summary_plots_SNVs_Indels_Mutect2_FL_min_AF_0.1.pdf", width=14)
plotmafSummary(maf = fl, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf = fl, top = 20)
laml.titv = titv(maf = fl, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)
lollipopPlot(maf = fl, gene = 'CREBBP', AACol = 'AAChange.refGene', showMutationRate = TRUE)
lollipopPlot(maf = fl, gene = 'EZH2', AACol = 'AAChange.refGene', showMutationRate = TRUE)
lollipopPlot(maf = fl, gene = 'KMT2D', AACol = 'AAChange.refGene', showMutationRate = TRUE)
dev.off()



