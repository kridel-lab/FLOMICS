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

variant_files = list.files(pattern="_.rds")
z = which(str_detect(variant_files, "Anjali"))
variant_files=variant_files[-z]

#gene annotations
genes = unique(fread("ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#generate summary results 

summary_vars = function(variants){
  
  #extract name of tool that was used to call variants and whether they are INDELs or SNPs
  variants_callers=paste(unlist(strsplit(variants, "_"))[3:5], collapse=" ")
  print(variants_callers)
  mut_type = "SNV"
  
  indel1 = !(length(which(str_detect(variants, "indel")))==0)
  indel2 = !(length(which(str_detect(variants, "INDEL")))==0)
  if(indel1 | indel2){
    mut_type="INDEL"
  }
  
  #read in file 
  var_dat = readRDS(variants)
  
  #remove population variants via snpdb 
  var_dat$avsnp142 = as.character(var_dat$avsnp142)
  z = which(str_detect(var_dat$avsnp142, "rs"))
  if(!(length(z)==0)){
    var_dat = var_dat[-z,]
  }

  #remove potential population variants with high popluation allele frequencies 
  var_dat$AF_popmax = as.numeric(as.character(var_dat$AF_popmax))
  var_dat = as.data.table(filter(var_dat, ((is.na(AF_popmax) )| AF_popmax < 0.005)))
  var_dat$Gene.ensGene = as.character(var_dat$Gene.ensGene)
  #check if some mutations are potentially artifacts
  #summary of actual mutations 
  freq_muts = as.data.table(filter(as.data.table(table(var_dat$patient, var_dat$id, var_dat$Gene.ensGene)), N >=1))
  freq_muts$combo = paste(freq_muts$V2, freq_muts$V3)
  freq_muts = as.data.table(table(freq_muts$combo))
  freq_muts = freq_muts[order(-N)]
  freq_muts$tag[freq_muts$N > 0.05 * length(unique(var_dat$patient))] = "tagged"
  var_dat$combo = paste(var_dat$id, var_dat$Gene.ensGene)
  z = which(var_dat$combo %in% freq_muts$V1[freq_muts$tag=="tagged"])
  var_dat$tag = ""
  if(!(length(z)==0)){
    var_dat$tag[z] = "tagged"
  }
  
  #get gene names 
  #ensembl to gene ID
  
  #if multiple genes mapped to variant (usually if in between genes or upstream of genes, or if gene has multiple ENSG ids)
  #keep id of first gene
  var_dat$Gene.ensGene = sapply(var_dat$Gene.ensGene, function(t){unlist(strsplit(t, '\\x', fixed=TRUE))[1]})
  var_dat$Gene.ensGene = as.character(var_dat$Gene.ensGene)
  var_dat = merge(var_dat, genes , by = "Gene.ensGene", allow.cartesian=TRUE)
  var_dat = as.data.table(filter(var_dat, !(is.na(hg19.ensemblToGeneName.value))))
  
  pdf(paste("/Users/kisaev/FLOMICS_Anjali/FLOMICS_Anajli/Plots/", variants_callers, "plots.pdf", sep="_"), width=15, height=10)
  
  #1. summary of # of muts per patient
  freq_pats = as.data.table(table(var_dat$patient))
  freq_pats = freq_pats[order(-N)]
  
  print(ggbarplot(freq_pats, x = "V1", y="N", fill="grey") +theme_bw() + 
    rotate_x_text(65) + ylab("Number of called mutations") + xlab("Patient") + ggtitle(paste(variants_callers, ",n=", dim(freq_pats)[1], "at least 1 mutation")))

  #2. summary genes  mutations 
  freq_genes = as.data.table(filter(as.data.table(table(var_dat$patient, var_dat$hg19.ensemblToGeneName.value)), N >=1))
  freq_genes = as.data.table(table(freq_genes$V2))
  freq_genes = freq_genes[order(-N)]
  
  print(ggbarplot(freq_genes, x = "V1", y="N", fill="grey") +theme_bw() + 
    rotate_x_text(65) + ylab("Number of patients with mutation") + xlab("Patient") + ggtitle(paste(variants_callers, ",n=", dim(freq_pats)[1], "at least 1 mutation")))
  
  #3. summary of where in genome mutations fall 
  table(as.character(var_dat$Func.ensGene))
  
  
  #print(ggbarplot(freq_muts, x = "V1", y="N", fill="grey") +theme_bw() + 
  #  rotate_x_text(90) + ylab("Number of patients with mutation") + xlab("Patient") + ggtitle(paste(variants_callers, ",n=", dim(freq_pats)[1], "at least 1 mutation")))
 
  dev.off() 
  
  #finally save final RDS files with either INDELs or SNPs for Anjali to conduct further analysis..
  data_type = paste(unlist(strsplit(variants, "_"))[3:6], collapse="_")
  var_dat$data_type = data_type
  saveRDS(var_dat, file=paste(date, data_type, "final_version_for_Anjali.rds", sep="_"))
  var_dat$mut_type = mut_type
  return(var_dat)
}

#apply to all 
all_muts = as.data.table(ldply(llply(variant_files, summary_vars)))
all_muts$combo = paste(all_muts$id, all_muts$hg19.ensemblToGeneName.value)

#summary 
all_muts_sum = as.data.table(table(all_muts$combo, all_muts$data_type, all_muts$mut_type))
all_muts_sum = as.data.table(filter(all_muts_sum, N >0))
colnames(all_muts_sum) = c("mut_id", "mutation_caller", "type_mut", "num_patients")
















