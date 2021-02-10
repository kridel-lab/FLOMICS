#----------------------------------------------------------------------
#check overlap of mutation calls for subset of samples from PLOS medicine
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
length(unique(new_plos$OTHER_ID)) #31

#2. keep only libraries from T1
new_plos = filter(new_plos, sample_type  == "T1")
length(unique(new_plos$OTHER_ID)) #31
write.table(new_plos, file="DNAseq/Mutation_calls_KI/PLOS_MED/clean_up_PLOS_mutations_mutect2_unfiltered_Feb2021.txt", sep=";", quote=F, row.names=F)

#3. set up new data to look like old data

#rescue hotspot mutations that may have been missed
new_plos_rescue = new_plos
new_plos_rescue <- new_plos_rescue %>%
  mutate(Cohort ="PLOSMED") %>%
  filter(avsnp142 == "." | (avsnp142 != "." & cosmic68 != ".")) %>%
  filter(Var_Freq <= 0.1) %>%
  filter((Chromosome == "1" & Start_Position %in% c("150727482")) |                                                      # CTSS hotspot
         (Chromosome == "7" & Start_Position %in% c("148506437", "148506467", "148508727", "148508728")) |               # EZH2 hotspots
         (Chromosome == "8" & Start_Position %in% c("20074767", "20074768")) |                                           # ATP6V1B2 hotspots
         (Chromosome == "12" & Start_Position %in% c("57496660", "57496661", "57496662", "57498345", "57499079")))  %>%       # STAT6 hotspots
         select(SAMPLE_ID, Hugo_Symbol, Cohort, Chromosome, Start_Position,
           Reference_Allele, Tumor_Seq_Allele2,
            End_Position)
new_plos_rescue$type = "rescued_hotspots"
colnames(new_plos_rescue)[7] = "Variant_Allele"

#filter main dataset using same filters as applied to main FL dataset
new_plos = new_plos %>%
  filter(avsnp142 == "." | (avsnp142 != "." & cosmic68 != ".")) %>%
  filter(Var_Freq > 0.1)
new_plos = new_plos %>%
  mutate(Cohort ="PLOSMED") %>%
    select(SAMPLE_ID, Hugo_Symbol, Cohort, Chromosome, Start_Position,
      Reference_Allele, Tumor_Seq_Allele2,
       End_Position)
new_plos$type = "passed_filters"
colnames(new_plos)[7] = "Variant_Allele"

new_plos = rbind(new_plos, new_plos_rescue)
print(dim(new_plos)) #277   9

#keep only genes present in original plos call set
new_plos = filter(new_plos, Hugo_Symbol %in% old_plos$Hugo_Symbol)
print(dim(new_plos)) #266   9

#4. check overlap

#for each patient check overlap
pats=unique(new_plos$SAMPLE_ID)

get_overlap = function(patient){
  muts_old = filter(old_plos, SAMPLE_ID == patient)
  muts_new = filter(new_plos, SAMPLE_ID == patient)
  old_muts = dim(muts_old)[1]
  new_muts = dim(muts_new)[1]
  muts_both = dim(merge(muts_old, muts_new, by=c("Hugo_Symbol", "Cohort", "Chromosome", "Start_Position")))[1]
  perc_overlap_wnew_calls = muts_both/new_muts
  perc_overlap_wold_calls = muts_both/old_muts
  result = c(old_muts, new_muts, muts_both, patient, perc_overlap_wnew_calls, perc_overlap_wold_calls)
  names(result) = c("num_old_muts", "num_new_muts", "num_both", "patient", "perc_overlap_wnew_calls", "perc_overlap_wold_calls")
  return(result)
}

all_muts_summary = as.data.table(ldply(llply(pats, get_overlap)))
all_muts_summary$perc_overlap_wnew_calls = as.numeric(all_muts_summary$perc_overlap_wnew_calls)
all_muts_summary$perc_overlap_wold_calls = as.numeric(all_muts_summary$perc_overlap_wold_calls)

write.csv(all_muts_summary, file="DNAseq/Mutation_calls_KI/PLOS_MED/overlap_PLOS_medicine_calls_with_previous_calls.csv", quote=F, row.names=F)
write.csv(new_plos, file="DNAseq/Mutation_calls_KI/PLOS_MED/clean_up_PLOS_mutations_mutect2_Feb2021.csv", quote=F, row.names=F)
