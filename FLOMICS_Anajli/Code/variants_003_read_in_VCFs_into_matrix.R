#----------------------------------------------------------------------
#variants_003_read_in_VCFs_into_matrix.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
	"ggrepel", "stringr", "maftools", "VariantAnnotation")
lapply(packages, require, character.only = TRUE)
library(maftools)

date = Sys.Date()

print(date)
#args = commandArgs(trailingOnly = TRUE)
#index = args[1]
#print(index) #index is the name of the algorithm that was used, loFreq_VCFs for example

setwd("/cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?

#----------------------------------------------------------------------
#data 
#----------------------------------------------------------------------

vcfs = list.files(pattern="multianno.txt") #normalized annotated vcf files (n=131)

#gene annotations
genes = unique(fread("/cluster/projects/kridelgroup/paired_cns/ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. clean up individual paired vcf files

test=function(vcf){
  print(vcf)
  pat = unlist(strsplit(vcf, "\\."))[1]
  algo = paste(unlist(strsplit(vcf, "\\."))[c(3,4,5)], collapse="_")
  if(algo == "vcf"){
    algo = unlist(strsplit(vcf, "\\."))[4]
  }
  if(algo == "vcf_hg19_multianno_txt"){
    algo = unlist(strsplit(vcf, "\\."))[2]
  }
  if(algo == "txt"){
    algo = unlist(strsplit(vcf, "\\."))[2]
  }
  print(pat)
  print(algo)
}

check=vcfs[which(str_detect(vcfs, "B48364"))]
sapply(check, test)

#what is this? B48258.variant_calls.pass.vcf.snp_B48258.variant_calls.pass.vcf.snp.vcf.hg19_multianno.txt
#deleted it not sure why it came up like this 

clean_up_001 = function(vcf){

  print(vcf)
  pat = unlist(strsplit(vcf, "\\."))[1]
  algo = unlist(strsplit(vcf, "\\."))[5]
  print(pat)
  algo = paste(unlist(strsplit(vcf, "\\."))[c(3,4,5)], collapse="_")
  if(algo == "vcf"){
    algo = unlist(strsplit(vcf, "\\."))[4]
  }
  if(algo == "vcf_hg19_multianno_txt"){
    algo = unlist(strsplit(vcf, "\\."))[2]
  }
  if(algo == "txt"){
    algo = unlist(strsplit(vcf, "\\."))[2]
  }

  if(!(dim(fread(vcf))[1] == 0)){ #to prevent empty VCF files from throwing errors

  #transform into maf format for downstream analysis 
  var.annovar <- vcf
  var.annovar.maf <- annovarToMaf(annovar = var.annovar, refBuild = 'hg19', table = 'ensGene')

  var.annovar.maf$patient = pat
  var.annovar.maf$algo = algo
  var.annovar.maf = as.data.table(filter(var.annovar.maf, V39=="PASS"))

  #split V40 into DP and V41 GT, AD, AF, DP
  #var.annovar.maf = var.annovar.maf %>% separate(V40, c("DP"), sep=";") %>% 
  #  separate(V42, c("GT", "AD", "AF", "DP_other"), sep=":")

  #generating MAF files:
  #Mandatory fields: Hugo_Symbol, Chromosome, Start_Position, End_Position, 
  #Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, 
  #Variant_Type and Tumor_Sample_Barcode.

  cols_keep=c("patient", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
    "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "tx", "exon", "txChange", "aaChange", 
    "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene", "AF_popmax", 
    "hgnc_symbol", 
    "algo", "id", "seqnames", "start", "end", "REF", "ALT", "DP", "POP_AF", "Func.ensGene", 
    "Gene.ensGene", "GeneDetail.ensGene", "ExonicFunc.ensGene", "AAChange.ensGene", "cosmic68", "avsnp142") 
  z = which(colnames(var.annovar.maf) %in% cols_keep)
  var.annovar.maf = var.annovar.maf[,..z]
  print(head(var.annovar.maf))
  print("done")
  #retrun
  return(var.annovar.maf)
}
}#(summary(vcf_dat_info$id)[1] == 0)

all_vcfs_text = as.data.table(ldply(llply(vcfs, clean_up_001, .progress="text")))
print("done")
#print(index)
saveRDS(all_vcfs_text, file=paste("all_VCFs_all_algos_merged", date, ".rds", sep="_"))


