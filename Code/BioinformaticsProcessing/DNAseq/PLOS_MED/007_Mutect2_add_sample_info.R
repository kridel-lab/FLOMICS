#----------------------------------------------------------------------
#karin isaev
#post mutect2 and annovar soft filtering and matrix prep
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)
date=Sys.Date()

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)
library(readxl)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#mutect2 run on tumour only mode
#annotated by annovar

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#mutations
muts=fread(list.files(pattern="Mutect2_filtered_mutations_PLOS_MED")[length(list.files(pattern="Mutect2_filtered_mutations_PLOS_MED"))])

#samples
samp_info = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/sample-annot_PLOSMED.tsv")
colnames(samp_info)[3] = "sample_name"

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#remove multiallelic variants and keep (indels)
#muts = as.data.table(filter(muts, Tumor_Seq_Allele2 %in% c("A", "C", "G", "T"),
#Reference_Allele %in% c("A", "C", "G", "T")))

z = which(str_detect(muts$Tumor_Seq_Allele2, ","))
if(!(length(z)==0)){
muts = muts[-z,]}

muts$sample_mut = paste(muts$Chromosome, muts$Start_Position, muts$sample_name, sep="_")
all_pairs = unique(muts$sample_mut)

get_mut_info = function(pair){

	mut_dat = as.data.table(filter(muts, sample_mut == pair))
	pat_clean= unlist(strsplit(mut_dat$sample_name, "T1"))[1]

	#add sample information
	mut_dat$clean_sample = pat_clean

	return(mut_dat)
}


all_muts_checked = as.data.table(ldply(llply(all_pairs, get_mut_info, .progress="text")))

unique(all_muts_checked$sample_name[which(!(all_muts_checked$sample_name %in% samp_info$sample_name))])
#[1] "FL1176T2" "FL1177T2" "FL1185T2" "FL1216T2"

#merge with sample information
all_muts_checked = merge(all_muts_checked, samp_info, by="sample_name")
all_muts_checked$Tumor_Sample_Barcode = NULL

length(unique(all_muts_checked$patient_id))
#[1] 277

write.table(all_muts_checked, file=paste("/cluster/projects/kridelgroup/FLOMICS/DATA/",
date, "_Mutect2_filtered_mutations_PLOS_MED_mut_info.txt", sep=""), quote=F, row.names=F, sep=";")
