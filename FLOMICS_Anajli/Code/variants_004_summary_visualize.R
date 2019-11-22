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

setwd("/Users/kisaev/FLOMICS_Anajli/Data")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?


#----------------------------------------------------------------------
#data 
#----------------------------------------------------------------------

variant_files = list.files(pattern=".rds")

#gene annotations
genes = unique(fread("ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#generate summary results 

summary_vars = function(variants){
  
  variants_callers=paste(unlist(strsplit(variants, "_"))[3:5], collapse=" ")
  
  var_dat = readRDS(variants)
  #remove population variants 
  z = which(str_detect(var_dat$avsnp142, "rs"))
  if(!(length(z)==0)){
    var_dat = var_dat[-z,]
  }

  #get gene names 
  #ensembl to gene ID
  
  #if multiple genes mapped to variant (usually if in between genes or upstream of genes)
  #keep id of first gene
  var_dat$Gene.ensGene = sapply(var_dat$Gene.ensGene, function(t){unlist(strsplit(t, '\\x', fixed=TRUE))[1]})
  var_dat = merge(var_dat, genes , by = "Gene.ensGene", allow.cartesian=TRUE)
  var_dat = as.data.table(filter(var_dat, !(is.na(hg19.ensemblToGeneName.value))))
  
  pdf(paste("/Users/kisaev/FLOMICS_Anajli/Plots/", variants_callers, "plots.pdf", sep="_"), width=15, height=10)
  
  #1. summary of # of muts per patient
  freq_pats = as.data.table(table(var_dat$patient))
  freq_pats = freq_pats[order(-N)]
  
  print(ggbarplot(freq_pats, x = "V1", y="N", fill="grey") +theme_bw() + 
    rotate_x_text(65) + ylab("Number of called mutations") + xlab("Patient") + ggtitle(paste(variants_callers, ",n=", dim(freq_pats)[1])))

  #2. summary genes  mutations 
  freq_genes = as.data.table(filter(as.data.table(table(var_dat$patient, var_dat$hg19.ensemblToGeneName.value)), N >=1))
  freq_genes = as.data.table(table(freq_genes$V2))
  freq_genes = freq_genes[order(-N)]
  
  print(ggbarplot(freq_genes, x = "V1", y="N", fill="grey") +theme_bw() + 
    rotate_x_text(65) + ylab("Number of patients with mutation") + xlab("Patient") + ggtitle(paste(variants_callers, ",n=", dim(freq_pats)[1])))
  
  #2. summary of actual mutations 
  freq_muts = as.data.table(filter(as.data.table(table(var_dat$patient, var_dat$id)), N >=1))
  freq_muts = as.data.table(table(freq_muts$V2))
  freq_muts = freq_muts[order(-N)]
  
  print(ggbarplot(freq_muts, x = "V1", y="N", fill="grey") +theme_bw() + 
    rotate_x_text(90) + ylab("Number of patients with mutation") + xlab("Patient") + ggtitle(paste(variants_callers, ",n=", dim(freq_pats)[1])))
 
  dev.off() 

}


#apply to all 
llply(variant_files, summary_vars)


















