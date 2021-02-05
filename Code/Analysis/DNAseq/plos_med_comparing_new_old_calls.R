#check overlap of mutation calls for subset of samples from PLOS medicine

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
#setwd("FLOMICS") #From here /Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/
getwd()
#"/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS"

#----------------------------------------------------------------------
#data filtering
#----------------------------------------------------------------------

new_plos = fread("DNAseq/Mutation_calls_KI/PLOS_MED/2021-02-02_Mutect2_filtered_mutations_PLOS_MED_mut_info.txt")
colnames(new_plos)[which(colnames(new_plos) == "patient_id")] = "OTHER_ID"

old_plos = read.csv("DNAseq/BC_Cancer_capseq_data/BC_Cancer_capseq_data.csv", header = T) %>%
  mutate(Cohort = "PLOSMED") %>%
  select(SAMPLE_ID, Hugo_Symbol, Cohort, Chromosome, Start_Position, Reference_Allele,
  Variant_Allele, End_Position)

clin_data = fread("metadata/sample_annotations_rcd6Nov2019.csv")
clin_data = unique(clin_data[,c("SAMPLE_ID", "SAMPLE_ID_TGL", "RES_ID", "LY_FL_ID", "OTHER_ID", "INSTITUTION")])

#1. convert plos IDs to IDs we use mainly in our project
new_plos = merge(new_plos, clin_data, by = "OTHER_ID")
length(unique(new_plos$OTHER_ID))
#31

#2. keep only libraries from T1
new_plos = filter(new_plos, sample_type  == "T1")
length(unique(new_plos$OTHER_ID))
#31

#3. set up new data to look like old data

new_plos = new_plos %>%
  filter(avsnp142 == "." | (avsnp142 != "." & cosmic68 != ".")) %>%
  filter(Var_Freq > 0.1)

new_plos = new_plos %>%
  mutate(Cohort ="PLOSMED") %>%
    select(SAMPLE_ID, Hugo_Symbol, Cohort, Chromosome, Start_Position,
      Reference_Allele, Tumor_Seq_Allele2,
       End_Position)
colnames(new_plos)[7] = "Variant_Allele"

#4. check overlap

#for each patient check overlap
pats=unique(new_plos$SAMPLE_ID)

get_overlap = function(patient){
  muts_old = filter(old_plos, SAMPLE_ID == patient)
  muts_new = filter(new_plos, SAMPLE_ID == patient)
  old_muts = dim(muts_old)[1]
  new_muts = dim(muts_new)[1]
  muts_both = dim(merge(muts_old, muts_new, by=c("Hugo_Symbol", "Cohort", "Chromosome", "Start_Position")))[1]
  result = c(old_muts, new_muts, muts_both, patient)
  names(result) = c("num_old_muts", "num_new_muts", "num_both", "patient")
  return(result)
}

all_muts_summary = as.data.table(ldply(llply(pats, get_overlap)))
write.csv(all_muts_summary, file="DNAseq/Mutation_calls_KI/PLOS_MED/overlap_PLOS_medicine_calls_with_previous_calls.csv", quote=F, row.names=F)
write.csv(new_plos, file="DNAseq/Mutation_calls_KI/PLOS_MED/clean_up_PLOS_mutations_mutect2_Feb2021.csv", quote=F, row.names=F)
