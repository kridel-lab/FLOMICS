#local exploration of mutations across clusters and stages
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors = F)
#avoid scientific notation
options(scipen=999)
#load libraries
packages <- c("dplyr", "readr", "ggplot2", "tidyr",
"data.table", "plyr",
"ggrepel", "stringr", "maftools", "ggpubr", "readxl", "skimr",
 "edgeR", "annotables", "EnvStats", "maftools", "gridExtra")
lapply(packages, require, character.only = TRUE)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
library(NMF)
library(pheatmap)

#date
date=Sys.Date()

setwd("~/UHN/kridel-lab - Documents (1)/FLOMICS")

#directory with FLOMICS related matrices
#setwd("FLOMICS") should equal to this /Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/

#----------------------------------------------------------------------
#data filtering
#----------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#load data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. Gene panel
gene.panel <- read.csv("DNAseq/Final_target_files/coding_gene_panel_PLOSMED_FLOMICS.csv", header = TRUE)
genes.common <- gene.panel %>%
  filter(PLOS_MED_PANEL == "YES" & FLOMICS_PANEL == "YES") %>%
  .$Gene.Name %>% as.character() # n = 57

#2. All FLOMICS samples included
all.samples.DNAseq.FLOMICS <- read.csv("metadata/sample_annotations_rcd6Nov2019.csv", header = T) %>%
  filter(CAPSEQ_DATA == TRUE & CAPSEQ_DATA_YEAR == "2019" & CAPSEQ_INCLUDE == "YES") %>%
  .$SAMPLE_ID %>% as.character() # n=131 (n=123 T1, n=8 T2)

#3. All PLOSMED samples included
all.samples.DNAseq.PLOSMED <- read.csv("metadata/sample_annotations_rcd6Nov2019.csv", header = T) %>%
  filter(CAPSEQ_DATA == TRUE & CAPSEQ_DATA_YEAR == "2015" & CAPSEQ_INCLUDE == "YES") %>%
  .$SAMPLE_ID %>% as.character() # n=31 (only T1)

#4. Coverage FLOMICS DNAseq
coverage.samples <- read.csv("DNAseq/Targeted_Seq_Coverage/08-07-2020/2020-08-07_sample_based_coverage_summary.csv")

#5. Read in mutation calls for FLOMICS and filter
mut.FLOMICS <- fread("DNAseq/Mutation_calls_KI/08-13-2020/with_indels/2020-08-13_Mutect2_filtered_mutations_FL_wRNASeq_mut_info.txt") %>%
  mutate(Cohort = "FLOMICS") %>%
  filter(avsnp142 == "." | (avsnp142 != "." & cosmic68 != ".")) %>%
  filter(Var_Freq > 0.1)
dim(mut.FLOMICS)

#6. Hotspot mutations curated by RK
hotspots <- read.table("DNAseq/hotspot_mutations.txt", sep = "\t", header = T) %>%
  mutate(mutation_id = paste(Chromosome, Start_Position, sep = "_")) %>% .$mutation_id

#7. rescued FL mutations
mut.FLOMICS.rescued <- fread("DNAseq/Mutation_calls_KI/08-13-2020/with_indels/2020-08-13_Mutect2_filtered_mutations_FL_wRNASeq_mut_info.txt") %>%
  mutate(Cohort = "FLOMICS") %>%
  filter(Var_Freq <= 0.1 & Var_Freq > 0.05) %>%
  mutate(mutation_id = paste(Chromosome, Start_Position, sep = "_")) %>%
  filter(cosmic68 != "." | mutation_id %in% hotspots) %>%
  select(-mutation_id)
dim(mut.FLOMICS.rescued)

#8. combined FL mutations
mut.FLOMICS <- mut.FLOMICS %>% rbind(mut.FLOMICS.rescued)
dim(mut.FLOMICS)
length(unique(mut.FLOMICS$Tumor_Sample_Barcode)) # n = 818 rows & 131 unique patients

#9. read in mutation calls for PLOS MED and filter
mut.PLOSMED <- fread("DNAseq/Mutation_calls_KI/PLOS_MED/clean_up_PLOS_mutations_mutect2_unfiltered_Feb2021.txt") %>%
  mutate(Cohort = "PLOSMED") %>%
  filter(avsnp142 == "." | (avsnp142 != "." & cosmic68 != ".")) %>%
  filter(Var_Freq > 0.1)
dim(mut.PLOSMED)

#10. rescued PLOS mutations
mut.PLOSMED.rescued <- fread("DNAseq/Mutation_calls_KI/PLOS_MED/clean_up_PLOS_mutations_mutect2_unfiltered_Feb2021.txt") %>%
  mutate(Cohort = "PLOSMED") %>%
  filter(Var_Freq <= 0.1 & Var_Freq > 0.05) %>%
  mutate(mutation_id = paste(Chromosome, Start_Position, sep = "_")) %>%
  filter(cosmic68 != "." | mutation_id %in% hotspots) %>%
  select(-mutation_id)
dim(mut.PLOSMED.rescued)

#11. combined PLOS mutations
mut.PLOSMED <- mut.PLOSMED %>% rbind(mut.PLOSMED.rescued)
dim(mut.PLOSMED)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Analyze FLOMICS data in more detail
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Identify poor coverage samples and filter them out
poor.coverage.samples <- coverage.samples %>%
  filter(mean_cov < 50) %>% .$External_identifier %>% as.character() #18 samples (16=T1 and 2=T2)
length(poor.coverage.samples)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#combine mutations FLOMICs and PLOS MED
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Merge FLOMICS and PLOSMED
mut.FLOMICS = mut.FLOMICS %>% mutate(SAMPLE_ID = Tumor_Sample_Barcode)

#get common columns
cols = which(colnames(mut.FLOMICS) %in% colnames(mut.PLOSMED))
mut.FLOMICS = as.data.frame(mut.FLOMICS[,..cols])

cols = which(colnames(mut.PLOSMED) %in% colnames(mut.FLOMICS))
mut.PLOSMED = as.data.frame(mut.PLOSMED[,..cols])

mut.PLOSMED <- mut.PLOSMED[names(mut.FLOMICS)]

mut.merged <- as.data.table(mut.FLOMICS %>%
  rbind(mut.PLOSMED) %>%
  filter(Hugo_Symbol %in% genes.common)) # n = 954

# Percentage of mutations per gene in merged cohort (max 1 mutation per gene per sample counted)
mut.merged = mut.merged %>%
  mutate(TIME_POINT = substr(SAMPLE_ID, 11, 12)) %>%
  filter(TIME_POINT != "T2") %>%
  filter(!SAMPLE_ID %in% poor.coverage.samples) #816

#clinical data
clin = fread("metadata/clinical_data_rcd11Aug2020.csv")
colnames(clin)[1] = "Tumor_Sample_Barcode"
colnames(clin)[19] = "Overall_Survival_Status"
colnames(clin)[23] = "days_to_last_followup"

#SNF labels
labels = fread("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/Cluster Labels/InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv")
colnames(labels)[2] = "Tumor_Sample_Barcode"
labels$Tumor_Sample_Barcode = as.character(labels$Tumor_Sample_Barcode)
labels$status = sapply(labels$Tumor_Sample_Barcode, function(x){unlist(strsplit(x, "_"))[4]})
labels$Tumor_Sample_Barcode = sapply(labels$Tumor_Sample_Barcode, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
labels$SNFClust = labels$SNFClust10Feb2021
labels = as.data.table(filter(labels, status == "T1"))
labels = labels[,c("Tumor_Sample_Barcode", "SNFClust")]
labels$SNFClust[is.na(labels$SNFClust)] = "na"
print(table(labels$SNFClust))

#-----------
#ANALYSIS
#-----------

#maftools
#Mandatory fields: Hugo_Symbol, Chromosome, Start_Position, End_Position,
#Reference_Allele, Tumor_Seq_Allele2,
#Variant_Classification, Variant_Type and Tumor_Sample_Barcode.

#> unique(mut.merged$Variant_Classification)
# [1] "exonic frameshift_deletion"        "exonic frameshift_insertion"
# [3] "exonic nonsynonymous_SNV"          "exonic nonframeshift_deletion"
# [5] "exonic stopgain"                   "exonic nonframeshift_insertion"
# [7] "exonic nonframeshift_substitution" "exonic stoploss"
# [9] "Missense_Mutation"                 "Splice_Site"
#[11] "Frame_Shift_Ins"                   "Frame_Shift_Del"
#[13] "Nonsense_Mutation"                 "In_Frame_Del"

#change to match maftools requirements
mut.merged$Variant_Classification[mut.merged$Variant_Classification == "exonic frameshift_deletion"] = "Frame_Shift_Del"
mut.merged$Variant_Classification[mut.merged$Variant_Classification == "exonic frameshift_insertion"] = "Frame_Shift_Ins"
mut.merged$Variant_Classification[mut.merged$Variant_Classification == "exonic nonsynonymous_SNV"] = "Missense_Mutation"
mut.merged$Variant_Classification[mut.merged$Variant_Classification == "exonic nonframeshift_deletion"] = "In_Frame_Del"
mut.merged$Variant_Classification[mut.merged$Variant_Classification == "exonic stopgain"] = "Nonsense_Mutation"
mut.merged$Variant_Classification[mut.merged$Variant_Classification == "exonic nonframeshift_insertion"] = "In_Frame_Ins"
mut.merged$Variant_Classification[mut.merged$Variant_Classification == "exonic nonframeshift_substitution"] = "Inframe_INDEL"
mut.merged$Variant_Classification[mut.merged$Variant_Classification == "exonic stoploss"] = "Nonstop_Mutation"

mut.merged[,ref_alt_diff := nchar(Reference_Allele) - nchar(Tumor_Seq_Allele2)]
mut.merged[, Variant_Type := ifelse(ref_alt_diff == 0 , yes = "SNP", no = ifelse(ref_alt_diff < 0 , yes = "INS", no = "DEL"))]

colnames(mut.merged)[which(colnames(mut.merged) == "SAMPLE_ID")] = "Tumor_Sample_Barcode"

#seperate into T1 and T2 samples
mut.merged$status = sapply(mut.merged$Tumor_Sample_Barcode, function(x){unlist(strsplit(x, "_"))[4]})
mut.merged$Tumor_Sample_Barcode = sapply(mut.merged$Tumor_Sample_Barcode, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
mut.merged = as.data.table(filter(mut.merged, status=="T1"))
mut.merged = merge(mut.merged, labels, by="Tumor_Sample_Barcode")
mut.merged$Chromosome = paste("chr", mut.merged$Chromosome, sep="")

snf1 = filter(mut.merged, SNFClust=="1")
snf2 = filter(mut.merged, SNFClust=="2")

clin$PRIM_TX = NULL
clin=merge(clin, labels, by="Tumor_Sample_Barcode")
write.table(mut.merged, file="Analysis-Files/maftools/maftools_mutations_test.txt", quote=F, row.names=F, sep="\t")
write.table(snf1, file="Analysis-Files/maftools/maftools_mutations_snf1_only.txt", quote=F, row.names=F, sep="\t")
write.table(snf2, file="Analysis-Files/maftools/maftools_mutations_snf2_only.txt", quote=F, row.names=F, sep="\t")
write.table(clin, file="Analysis-Files/maftools/maftools_clinical_file.txt", quote=F, row.names=F, sep="\t")

#snf1 FL
snf1 = read.maf(maf = "Analysis-Files/maftools/maftools_mutations_snf1_only.txt", clinicalData="Analysis-Files/maftools/maftools_clinical_file.txt") #56 unique patients
#snf2 FL
snf2 = read.maf(maf = "Analysis-Files/maftools/maftools_mutations_snf2_only.txt", clinicalData="Analysis-Files/maftools/maftools_clinical_file.txt") #45 unqiue patients
#all samples
all = read.maf(maf = "Analysis-Files/maftools/maftools_mutations_test.txt", clinicalData="Analysis-Files/maftools/maftools_clinical_file.txt") #125 unique patients

get_maf_plots = function(maf_ob, type_analysis){

  pdf(paste("Analysis-Files/maftools/", type_analysis, "_", date, "_maftools_prelim_plots.pdf", sep=""))
  #global overview of mutations
  plotmafSummary(maf = maf_ob, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  #in depth summary of top 30 genes
  oncoplot(maf = maf_ob, fontSize=0.5, clinicalFeatures ='TYPE', sortByAnnotation = TRUE, top=30)
  #onco plot using SNF labels
  oncoplot(maf = maf_ob, fontSize=0.5, clinicalFeatures ='SNFClust', sortByAnnotation = TRUE, top=30)

  #somatic interactions
  ints_res = somaticInteractions(maf = maf_ob, top = 25, pvalue = c(0.05, 0.01), fontSize=0.5)
  #types of mutation in depth - snf1 only
  fl.titv = titv(maf = maf_ob, plot = FALSE, useSyn = TRUE)
  #plot titv summary
  plotTiTv(res = fl.titv)

  dgi = drugInteractions(maf = maf_ob, fontSize = 0.5)
  #OncogenicPathways(maf = maf_ob)

  #Using top 55 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
  prog_geneset = survGroup(maf = maf_ob, geneSetSize = 2, #to test all gene pairs
    time = "days_to_last_followup", Status = "Overall_Survival_Status", verbose = TRUE)
  print(prog_geneset)
  #multiple testing correction
  prog_geneset$fdr = p.adjust(prog_geneset$P_value, method="fdr")

  gene1=unlist(strsplit(as.character(prog_geneset[1,1]), "_"))[1]
  gene2=unlist(strsplit(as.character(prog_geneset[1,1]), "_"))[2]

  mafSurvGroup(maf = maf_ob, geneSet = c(gene1, gene2),
  time = "days_to_last_followup", Status = "Overall_Survival_Status")

  dev.off()

}

get_maf_plots(all, "all_patients")
get_maf_plots(snf1, "snf1")
get_maf_plots(snf2, "snf2")

pdf(paste("Analysis-Files/maftools/", date, "_SNF1_vs2_maftools_prelim_plots.pdf", sep=""))
#compare the two cohorts
pt.vs.rt <- mafCompare(m1 = snf1, m2 = snf2, m1Name = 'SNF1 FL', m2Name = 'SNF2 FL', minMut = 2)
print(pt.vs.rt)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
#coBarplot(m1 = snf1, m2 = snf2, m1Name = "SNF1 FL", m2Name = "SNF2 FL")
dev.off()

#mutation signature analysis

get_mut_signatures = function(maf_file, type_analysis, n_sigs){
  
  pdf(paste("Analysis-Files/maftools/", type_analysis, "_", date, "_maftools_mutation_signature_analysis.pdf", sep=""), width=10)
  
  snf=maf_file
  
  snf.tnm = trinucleotideMatrix(maf = snf, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
  
  snf.sig = extractSignatures(mat = snf.tnm, n = n_sigs, pConstant = 1e-9)
  
  #Compate against original 30 signatures 
  snf.og30.cosm = compareSignatures(nmfRes = snf.sig, sig_db = "legacy")
  #Compate against updated version3 60 signatures 
  snf.v3.cosm = compareSignatures(nmfRes = snf.sig, sig_db = "SBS")
  pheatmap::pheatmap(mat = snf.og30.cosm$cosine_similarities, 
                     cluster_rows = FALSE, main = "cosine similarity against validated signatures (old)")
  pheatmap::pheatmap(mat = snf.v3.cosm$cosine_similarities, 
                     cluster_rows = FALSE, main = "cosine similarity against validated signatures (new)")
  maftools::plotSignatures(nmfRes = snf.sig, title_size = 1.2, sig_db = "SBS")
  dev.off()
}

get_mut_signatures(all, "all_patients", 2)
get_mut_signatures(snf1, "snf1", 2)
get_mut_signatures(snf2, "snf2", 2)


