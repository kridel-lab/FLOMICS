#local exploration of mutations across clusters and stages
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
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
 "edgeR", "annotables", "EnvStats", "maftools", "gridExtra")
lapply(packages, require, character.only = TRUE)

#date
date=Sys.Date()

#directory with FLOMICS related matrices
#setwd("FLOMICS") should equal to this /Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/

#----------------------------------------------------------------------
#data filtering
#----------------------------------------------------------------------

# Gene panel
gene.panel <- read.csv("DNAseq/Final_target_files/coding_gene_panel_PLOSMED_FLOMICS.csv", header = TRUE)
genes.common <- gene.panel %>%
  filter(PLOS_MED_PANEL == "YES" & FLOMICS_PANEL == "YES") %>%
  .$Gene.Name %>% as.character() # n = 57

muts <- fread("RNAseq/variant_calls/2020-06-23_opossum_variant_FL_rna-seq_filtered.txt")
qc = fread("metadata/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")

muts = as.data.table(filter(muts, Tumor_Sample_Barcode %in% qc$rna_seq_file_sample_ID))
length(unique(muts$SAMPLE_ID)) #128

#filter out mutation types we dont want
muts <- as.data.table(filter(muts, !(Variant_Classification %in% c("UTR3 .",
"downstream .", "UTR5 .", "upstream .", "exonic unknown"))))

# Coverage FLOMICS DNAseq
coverage.samples <- read.csv("DNAseq/Targeted_Seq_Coverage/08-07-2020/2020-08-07_sample_based_coverage_summary.csv")

#mutations obtained via dna sequencing for 131 samples (some overlap)
#some patients have very low coverage
# Read in mutation calls for FLOMICS
mut.FLOMICS <- fread("DNAseq/Mutation_calls_KI/08-13-2020/with_indels/2020-08-13_Mutect2_filtered_mutations_FL_wRNASeq_mut_info.txt") %>%
  mutate(Cohort = "FLOMICS") %>%
  filter(avsnp142 == "." | (avsnp142 != "." & cosmic68 != ".")) %>%
  filter(Var_Freq > 0.1)
dim(mut.FLOMICS)
length(unique(mut.FLOMICS$Tumor_Sample_Barcode))

# Identify poor coverage samples and filter them out
poor.coverage.samples <- coverage.samples %>%
filter(mean_cov < 50) %>% .$External_identifier %>% as.character() #18 samples
length(poor.coverage.samples)

# Read in mutation calls for PLOSMED
mut.PLOSMED <- read.csv("DNAseq/BC_Cancer_capseq_data/BC_Cancer_capseq_data.csv", header = T) %>%
  mutate(Cohort = "PLOSMED") %>%
  select(SAMPLE_ID, Hugo_Symbol, Cohort, Chromosome, Start_Position, Reference_Allele,
  Tumor_Seq_Allele2=Variant_Allele, End_Position, Variant_Classification=Mutation_Type)

# Merge FLOMICS and PLOSMED
mut.merged <- mut.FLOMICS %>%
  select(SAMPLE_ID = Tumor_Sample_Barcode, Hugo_Symbol, Cohort, Chromosome,
  Start_Position, Reference_Allele, Tumor_Seq_Allele2,
  End_Position, Variant_Classification) %>%
  rbind(mut.PLOSMED) %>%
  mutate(mut_id = paste(.$Chromosome, .$Start_Position, sep="_")) %>%
  filter(Hugo_Symbol %in% genes.common) #55 unique genes

#clinical data
clin = fread("metadata/clinical_data_rcd11Aug2020.csv")
colnames(clin)[1] = "Tumor_Sample_Barcode"
colnames(clin)[19] = "Overall_Survival_Status"
colnames(clin)[23] = "days_to_last_followup"

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

colnames(mut.merged)[1] = "Tumor_Sample_Barcode"

#seperate into T1 and T2 samples
mut.merged$status = sapply(mut.merged$Tumor_Sample_Barcode, function(x){unlist(strsplit(x, "_"))[4]})
mut.merged$Tumor_Sample_Barcode = sapply(mut.merged$Tumor_Sample_Barcode, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
t1 = filter(mut.merged, status=="T1")
t2 = filter(mut.merged, status=="T2")
clin$PRIM_TX = NULL
write.table(mut.merged, file="Analysis-Files/maftools/maftools_mutations_test.txt", quote=F, row.names=F, sep="\t")
write.table(t1, file="Analysis-Files/maftools/maftools_mutations_T1_only.txt", quote=F, row.names=F, sep="\t")
write.table(t2, file="Analysis-Files/maftools/maftools_mutations_T2_only.txt", quote=F, row.names=F, sep="\t")
write.table(clin, file="Analysis-Files/maftools/maftools_clinical_file.txt", quote=F, row.names=F, sep="\t")

#T1 FL
t1 = read.maf(maf = "Analysis-Files/maftools/maftools_mutations_T1_only.txt", clinicalData="Analysis-Files/maftools/maftools_clinical_file.txt")
#T2 FL
t2 = read.maf(maf = "Analysis-Files/maftools/maftools_mutations_T2_only.txt", clinicalData="Analysis-Files/maftools/maftools_clinical_file.txt")

pdf("Analysis-Files/maftools/maftools_prelim_plots.pdf")

#global overview of mutations - t1 only
plotmafSummary(maf = t1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#in depth summary of top 30 genes -t1 only
oncoplot(maf = t1, fontSize=0.5, clinicalFeatures ='TYPE', sortByAnnotation = TRUE, top=30)
#somatic interactions - t1 only
ints_res = somaticInteractions(maf = t1, top = 55, pvalue = c(0.05, 0.01), fontSize=0.5)
ints_res$fdr = p.adjust(ints_res$pValue, method="fdr")
write.table(ints_res, file="Analysis-Files/maftools/somatic_interactions_results_T1_only.txt", sep=";", quote=F, row.names=F)

#types of mutation in depth - T1 only
fl.titv = titv(maf = t1, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = fl.titv)

#Using top 55 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
prog_geneset = survGroup(maf = t1, top = 55, geneSetSize = 2, #to test all gene pairs
  time = "days_to_last_followup", Status = "Overall_Survival_Status", verbose = FALSE)
print(prog_geneset)
#multiple testing correction
prog_geneset$fdr = p.adjust(prog_geneset$P_value, method="fdr")
write.csv(prog_geneset, file="Analysis-Files/maftools/prognostic_partners_results_T1_only.csv", quote=F, row.names=F)

mafSurvGroup(maf = t1, geneSet = c("CREBBP", "STAT6"),
time = "days_to_last_followup", Status = "Overall_Survival_Status")

mafSurvGroup(maf = t1, geneSet = c("FOXO1", "SOCS1"),
time = "days_to_last_followup", Status = "Overall_Survival_Status")

mafSurvGroup(maf = t1, geneSet = c("CREBBP", "FOXO1"),
time = "days_to_last_followup", Status = "Overall_Survival_Status")

mafSurvGroup(maf = t1, geneSet = c("TNFRSF14", "B2M"),
time = "days_to_last_followup", Status = "Overall_Survival_Status")

#compare the two cohorts
pt.vs.rt <- mafCompare(m1 = t1, m2 = t2, m1Name = 'T1 FL', m2Name = 'T2 FL', minMut = 2)
print(pt.vs.rt)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1, color = c('royalblue', 'maroon'), geneFontSize = 0.8)

#genes = c("CREBBP")
#coOncoplot(m1 = t1, m2 = t2,
#  m1Name = 'T1', m2Name = 'T2', genes=c("TP53", "EZH2"), removeNonMutated = TRUE)
coBarplot(m1 = t1, m2 = t2, m1Name = "T1 FL", m2Name = "T2 FL")

fab.ce = clinicalEnrichment(maf = t1, clinicalFeature = 'TYPE')
fab.ce$groupwise_comparision[p_value < 0.05]
plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05)
dgi = drugInteractions(maf = t1, fontSize = 0.5)

dnmt3a.dgi = drugInteractions(genes = "CREBBP", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

OncogenicPathways(maf = t1)

PlotOncogenicPathways(maf = t1, pathways = "NOTCH")

dev.off()
