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

date = Sys.Date()

print(date)
args = commandArgs(trailingOnly = TRUE)
index = args[1]
print(index) #index is the name of the algorithm that was used, loFreq_VCFs for example

setwd(paste("/cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder/", index, sep=""))

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?

#----------------------------------------------------------------------
#data 
#----------------------------------------------------------------------

vcfs = list.files(pattern=".vcf") #normalized annotated vcf files (n=131)

#gene annotations
genes = unique(fread("/cluster/projects/kridelgroup/paired_cns/ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. clean up individual paired vcf files

clean_up_001 = function(vcf){

  pat = unlist(strsplit(vcf, "_"))[1]
  if(!(dim(fread(vcf))[1] == 3)){ #to prevent empty VCF files from throwing errors

  #read in VCF file 
  vcf_dat <- readVcf(vcf, "hg19")

  #extract genotype info into data table 
  vcf_dat_ranges = as.data.frame(rowRanges(vcf_dat))
  vcf_dat_ranges$id = rownames(vcf_dat_ranges)
  vcf_dat_info = info(vcf_dat)  
  vcf_dat_info$id = rownames(vcf_dat_info)
  vcf_dat = as.data.table(merge(vcf_dat_ranges, vcf_dat_info, by = "id"))

  #add patient ID 
  vcf_dat$patient = pat

  #retrun
  return(vcf_dat)
}
}

all_vcfs_text = as.data.table(ldply(llply(vcfs, clean_up_001, .progress="text")))
print("done")
print(index)
saveRDS(all_vcfs_text, file=paste("/cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder/SNP_matrices_all_algortihms/", date, index, ".rds", sep="_"))


