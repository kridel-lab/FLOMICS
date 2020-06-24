#----------------------------------------------------------------------
#processing_annovar_results.R
#karin isaev
#July 11th, 2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/annovar")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools", "readxl")
lapply(packages, require, character.only = TRUE)
library("readxl")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#mutect2 was run on paired mode compaing cns to diagnostic tumour
#now it's time to:
#summarize cns specific mutations
#but first should still filter out false positives (note, these are unfilitered variants)
#see how many appear in multiple comparisons (n=5 total)

#note these vcf files have been normalized and fed through annovar
#for annotations

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

date=Sys.Date()

paired = list.files(pattern="*final_vcf_file_filtered_maf_input.bed")
print(paired)

#gene annotations
genes = unique(fread("/cluster/home/kisaev/data/annotables_grch37.txt"))
genes = as.data.table(filter(genes, biotype == "protein_coding"))
genes = as.data.table(filter(genes, !(is.na(entrez))))

#bc mutation data
bc_dat = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/BC_mutation_data/BC_Cancer_capseq_data.csv")

#read in mutation files
read_f = function(f){
	ff = fread(f)
	return(ff)
}

all_dat = as.data.table(ldply(llply(paired, read_f)))
#write.table(all_dat, file="maftools_file_all_samples.txt", quote=F, row.names=F, sep="\t")

pats = as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.xlsx"))
colnames(pats)[3] = "Tumor_Sample_Barcode" #132 here

all_dat = merge(all_dat, pats, by = "Tumor_Sample_Barcode")

length(unique(all_dat$Tumor_Sample_Barcode)) #132

all_dat = as.data.table(filter(all_dat, Hugo_Symbol %in% genes$symbol, !(Variant_Classification %in% c("intronic .", "ncRNA_intronic .", "exonic synonymous_SNV",
"ncRNA_exonic .", "ncRNA_splicing ."))))

length(unique(all_dat$Tumor_Sample_Barcode)) #128

#add BC data - check if mutations present
bc_dat$pat_mut = paste(bc_dat$SAMPLE_ID, bc_dat$Hugo_Symbol, bc_dat$Chromosome, bc_dat$Start_Position, sep="_")
all_dat$chr = sapply(all_dat$Chromosome, function(x){unlist(strsplit(x, "chr"))[2]})
all_dat$pat_mut = paste(all_dat$SAMPLE_ID, all_dat$Hugo_Symbol, all_dat$chr, all_dat$Start_Position, sep="_")
z = which(bc_dat$pat_mut %in% all_dat$pat_mut)
bc_dat$mut_found_rnaseq = ""
bc_dat$mut_found_rnaseq[z] = "yes"

z = which(all_dat$pat_mut %in% bc_dat$pat_mut)
all_dat$mut_found_rnaseq = ""
all_dat$mut_found_rnaseq[z] = "yes"

write.table(all_dat, file=paste(date, "opossum_variant_FL_rna-seq_filtered.txt", sep="_"), sep="\t", quote=F, row.names=F)
