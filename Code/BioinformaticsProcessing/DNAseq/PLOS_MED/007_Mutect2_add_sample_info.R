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

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#remove multiallelic variants and keep (indels)
#muts = as.data.table(filter(muts, Tumor_Seq_Allele2 %in% c("A", "C", "G", "T"),
#Reference_Allele %in% c("A", "C", "G", "T")))

z = which(str_detect(muts$Tumor_Seq_Allele2, ","))
muts = muts[-z,]

muts$sample_mut = paste(muts$Chromosome, muts$Start_Position, muts$Tumor_Sample_Barcode, sep="_")
all_pairs = unique(muts$sample_mut)

get_mut_info = function(pair){

	mut_dat = as.data.table(filter(muts, sample_mut == pair))
	mut_dat$rnaseq_detect = ""
	mut_dat$old_muts = ""
	mut_dat$platypus = ""

	pat_clean= unlist(strsplit(mut_dat$Tumor_Sample_Barcode, "_T1"))[1]

	#[1] check if in old mutation data
	old = as.data.table(filter(old_muts, Chromosome == as.numeric(mut_dat$Chromosome),
		Start_Position == mut_dat$Start_Position ,
		External_ID %in% c(pat_clean, mut_dat$Tumor_Sample_Barcode)))
	if(!(dim(old)[1]) == 0){
			mut_dat$old_muts = "yes"
	}

	#[2] check if in rnaseq data
	rna = as.data.table(filter(rnaseq_muts, chr == as.numeric(mut_dat$Chromosome),
		Start_Position == mut_dat$Start_Position ,
		Tumor_Sample_Barcode %in% c(pat_clean, mut_dat$Tumor_Sample_Barcode)))
	if(!(dim(rna)[1]) == 0){
		mut_dat$rnaseq_detect = "yes"
	}

	#[3] check if in platypus data
	#plat = as.data.table(filter(platypus_muts, as.numeric(CHROM) == as.numeric(mut_dat$Chromosome),
	#	POS == mut_dat$Start_Position ,
	#	Indiv %in% c(pat_clean, mut_dat$Tumor_Sample_Barcode)))
	#if(!(dim(plat)[1]) == 0){
	#	mut_dat$platypus = "yes"
	#}

	return(mut_dat)
}


all_muts_checked = as.data.table(ldply(llply(all_pairs, get_mut_info, .progress="text")))
write.table(all_muts_checked, file=paste("/cluster/projects/kridelgroup/FLOMICS/DATA/",
date, "_Mutect2_filtered_mutations_FL_wRNASeq_mut_info.txt", sep=""), quote=F, row.names=F, sep=";")
