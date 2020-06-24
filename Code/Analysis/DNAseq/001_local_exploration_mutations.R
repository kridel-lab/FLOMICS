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

muts <- fread("2020-06-23_opossum_variant_FL_rna-seq_filtered.txt")
qc = as.data.table(read_excel("FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.xlsx"))

muts = as.data.table(filter(muts, Tumor_Sample_Barcode %in% qc$rna_seq_file_sample_ID))
length(unique(muts$SAMPLE_ID)) #128

#filter out mutation types we dont want
muts <- as.data.table(filter(muts, !(Variant_Classification %in% c("UTR3 .",
"downstream .", "UTR5 .", "upstream .", "exonic unknown"))))

#-----------
#ANALYSIS
#-----------

#compare gene mutation burden between clusters

get_comp = function(gene, dat, clusts_tot){
	print(gene)
	gene_dat = as.data.table(filter(dat, Hugo_Symbol == gene))
	#only count gene being mutated once per patient
	gene_dat = unique(gene_dat[,c("SAMPLE_ID", "STAGE", "Cluster")])

	gene_dat$STAGE[which(str_detect(gene_dat$SAMPLE_ID, "RLN"))] = "RLN"
	gene_dat$STAGE[which(str_detect(gene_dat$SAMPLE_ID, "DLC"))] = "DLBCL"
	gene_dat$Cluster = factor(gene_dat$Cluster)

	#plot summary
	plotsum = as.data.table(table(gene_dat$Cluster))
	g = ggbarplot(plotsum, x="V1", y="N", fill="grey", title=paste("mutations in", gene)) + xlab("cluster") +
	ylab("# of patients with at least 1 mutation") 
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


#tiers
tiers=c("tier3", "tier2", "tier1")

#apply downstream code to each tier and then compare results
get_tier_mut_summary = function(tier){

  if(tier == "tier3"){
    dat=muts
    print(length(unique(dat$SAMPLE_ID)))}

  if(tier == "tier2"){
    dat=muts
    z1 = which(dat$rRNAcontam > 40)
    z2 = which(dat$Uniquely.mapped < 10000000)
    z3 = which(dat$'Uniquely.mapped.reads..' < 50)
    z4 = which(dat$'X..of.reads.mapped.to.multiple.loci' > 20)

    total_lost = unique(c(z1, z2, z3, z4))
    dat=muts[-total_lost]
    print(length(unique(dat$SAMPLE_ID)))}

  if(tier == "tier1"){
    dat=muts
    z1 = which(dat$rRNAcontam > 40)
    z2 = which(dat$Uniquely.mapped < 10000000)
    z3 = which(dat$'Uniquely.mapped.reads..' < 70)
    z4 = which(dat$'X..of.reads.mapped.to.multiple.loci' > 20)

    total_lost = unique(c(z1, z2, z3, z4))
    dat=muts[-total_lost]
    print(length(unique(dat$SAMPLE_ID)))}

    #----------------------------------------------------------------------
    #analysis of mutations
    #----------------------------------------------------------------------

    #[1] which genes are most often mutated in one of two clusters
    clusters <- as.data.table(table(dat$Hugo_Symbol, dat$Cluster))
    colnames(clusters) <- c("gene", "cluster", "num_muts")
    clusters <- clusters[order(-num_muts)]

    #[2] for each gene compare proportion of patients in each cluster that have mutation
    #[1] how many mutations in each gene
    genes_sum <- as.data.table(table(dat$Hugo_Symbol))
    genes_sum <- genes_sum[order(-N)]
    genes_sum = as.data.table(filter(genes_sum, N > 5)) #genes wtih at least 5 mutations in them

    genes = unique(genes_sum$V1)

    dat$Cluster = factor(dat$Cluster, levels=c(1,2))
    #get num unique patients in each cluster
    clusts_tot = as.data.table(table(unique(dat[,c("SAMPLE_ID", "STAGE",
    "Cluster")])$Cluster))
    colnames(clusts_tot)[2] = "tot_pats"

    pdf(paste(date, tier, "Gene_based_rnaseq_variants_summary.pdf", sep="_"))
    summ = as.data.table(ldply(llply(genes, get_comp, dat, clusts_tot)))
    dev.off()

    summ$fishers_p <- as.numeric(summ$fishers_p)
    summ = summ[order(fishers_p)]
    summ$fishers_fdr <- p.adjust(summ$fishers_p, method="fdr")
    summ$tier = tier
    return(summ)

    print("done tier analysis")

}

all_tiers = as.data.table(ldply(llply(tiers, get_tier_mut_summary)))
saveRDS(all_tiers, paste(date, "mutation_analysis_rnaseq_all_tiers.rds", sep="_"))
