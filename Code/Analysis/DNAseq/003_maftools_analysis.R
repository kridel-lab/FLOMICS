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
  Variant_Allele, End_Position)

# Merge FLOMICS and PLOSMED
mut.merged <- mut.FLOMICS %>%
  select(SAMPLE_ID = Tumor_Sample_Barcode, Hugo_Symbol, Cohort, Chromosome,
  Start_Position, Reference_Allele, Variant_Allele=Tumor_Seq_Allele2,
  End_Position) %>%
  rbind(mut.PLOSMED) %>%
  mutate(mut_id = paste(.$Chromosome, .$Start_Position, sep="_")) %>%
  filter(Hugo_Symbol %in% genes.common) #55 unique genes

#-----------
#ANALYSIS
#-----------

#maftools
#Mandatory fields: Hugo_Symbol, Chromosome, Start_Position, End_Position,
#Reference_Allele, Tumor_Seq_Allele2,
#Variant_Classification, Variant_Type and Tumor_Sample_Barcode.
mut.merged$Variant_Classification = "Missense_Mutation"
snp = which((mut.merged$Variant_Allele %in% c("A", "G", "C", "T")) &
(mut.merged$Reference_Allele %in% c("A", "T", "C", "G")))
mut.merged$Variant_Type = "SNP"
colnames(mut.merged)[1] = "Tumor_Sample_Barcode"
colnames(mut.merged)[7] = "Tumor_Seq_Allele2"
write.table(mut.merged, file="Analysis-Files/maftools/maftools_mutations_test.txt", quote=F, row.names=F, sep="\t")

maffile=read.maf("Analysis-Files/maftools/maftools_mutations_test.txt")
pdf("Analysis-Files/maftools/maftools_prelim_plots.pdf")
plotmafSummary(maf = maffile, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#oncoplot for top ten mutated genes.
oncoplot(maf = maffile, top = 10)
ints_res = somaticInteractions(maf = maffile, top = 55, pvalue = c(0.05, 0.01), fontSize=0.5)
ints_res$fdr = p.adjust(ints_res$pValue, method="fdr")
write.table(ints_res, file="Analysis-Files/maftools/somatic_interactions_results.txt", sep=";", quote=F, row.names=F)
dev.off()
