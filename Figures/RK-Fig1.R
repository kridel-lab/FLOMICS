#---
# This script reads in mutation (capseq) data and plots all panels from Figure 1
# Author: Robert Kridel
#---

packages <- c("dplyr", "ggplot2", "maftools", "tidyverse", "gridExtra", "Rediscover")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

setwd("~/github/FLOMICS/")

#--
# Read in data
#--

# Read in gene panels and define common genes
genes_PLOSMED <- read.csv("DNAseq/Final_target_files/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
    filter(PLOS_MED_PANEL == "YES") %>% .$Gene.Name # n=86
genes_UHN <- read.csv("DNAseq/Final_target_files/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
    filter(FLOMICS_PANEL == "YES") %>% .$Gene.Name # n=71
genes_common <- intersect(genes_PLOSMED, genes_UHN) # n=57

# Read in list of all mutations
capseq_T1 <- read.csv("DNAseq/CAPSEQ/Mutect2/all-cohorts/2022-05-24_AllCohorts_E2408-TO_filtered_MUTECT2_calls_BCL2_coding_and_noncoding.csv") %>%
  mutate(TIMEPOINT = str_sub(SAMPLE_ID, -2, -1)) %>%
  filter(Hugo_Symbol %in% genes_common) %>%
  filter(SAMPLE_ID != "LY_FL_571_T1") %>%
  filter(!SAMPLE_ID %in% c("LY_FL_1135_T1", "LY_FL_1156_T1")) %>%
  filter(TIMEPOINT == "T1")

# Read in mutation matrix
capseq_T1_mat <- read.table("DNAseq/CAPSEQ/SNV_Cluster_Analysis/2022-05-15_UHN_E4402_PLOSMED-T-O_E2408-T-O_codingBCL2only/2022-05-15_gene_vs_sample_SNV_matrix.txt",
                            sep = ";", header = TRUE)
capseq_T1_mat <- capseq_T1_mat[ , !(names(capseq_T1_mat) %in% c("LY_FL_1135_T1", "LY_FL_1156_T1"))]

# Read in samples with no mutations but with sufficient coverage
samples_nomuts_goodcoverage <- read.table("DNAseq/CAPSEQ/Mutect2/2022-05-25_samples_without_SNVs_with_high_coverage_in_gene_vs_sample_mat.txt", sep = "\t", header = TRUE) 
samples_nomuts <- names(capseq_T1_mat[,(colSums(capseq_T1_mat) == 0)])

# Read in sample annotations
sample_annot <- read.table("metadata/20221228_sample_annotations.txt", sep = "\t", header = TRUE) %>%
  filter(SAMPLE_ID %in% capseq_T1$SAMPLE_ID | SAMPLE_ID %in% samples_nomuts) %>%
  mutate(cohort = ifelse(INSTITUTION %in% c("AARHUS", "JGH", "KINGSTON", "OSLO", "PAH", "UHN"), "UHN",
                  ifelse(INSTITUTION == "E4402", "E4402",
                  ifelse(INSTITUTION == "E2408", "E2408",
                  ifelse(INSTITUTION == "BCCA", "PLOSMED", NA)))))

#--
# Create maf file
#--

capseq_T1_maf <- capseq_T1 %>%
  dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Protein_Change = AAChange.ensGene,
         Variant_Classification, Variant_Type, Tumor_Sample_Barcode = SAMPLE_ID, VAF = Var_Freq) %>%
  mutate(Protein_Change = gsub("(.*),.*", "\\1", Protein_Change)) %>%
  mutate(Protein_Change = sub(".*p.", "", Protein_Change)) %>%
  filter(!Variant_Classification %in% c("UTR3 .", "UTR5 .", "intronic .", "upstream .")) %>%
  mutate(Variant_Classification = ifelse(Variant_Classification == "exonic nonsynonymous_SNV", "Missense_Mutation",
                                  ifelse(Variant_Classification == "exonic stopgain", "Nonsense_Mutation",
                                  ifelse(Variant_Classification == "exonic frameshift_deletion", "Frame_Shift_Del",
                                  ifelse(Variant_Classification == "exonic frameshift_insertion", "Frame_Shift_Ins",
                                  ifelse(Variant_Classification == "exonic nonframeshift_deletion", "In_Frame_Del",
                                  ifelse(Variant_Classification == "exonic nonframeshift_insertion", "In_Frame_Ins",
                                  ifelse(Variant_Classification == "exonic nonframeshift_substitution", "Missense_Mutation",
                                  ifelse(Variant_Classification == "exonic stoploss", "Missense_Mutation",
                                  ifelse(Variant_Classification == "splicing .", "Splice_Site", NA)))))))))) %>%
  mutate(Variant_Type = ifelse(Variant_Classification == "Frame_Shift_Del", "DEL", Variant_Type)) %>%
  mutate(Variant_Type = ifelse(Variant_Classification == "Frame_Shift_Ins", "INS", Variant_Type))

#--
# Define cases in individual cohorts
#--

UHN.cases_T1 <- sample_annot %>% filter(cohort == "UHN") %>% .$SAMPLE_ID # n=208
length(UHN.cases_T1) #203
E4402.cases_T1 <- sample_annot %>% filter(cohort == "E4402") %>% .$SAMPLE_ID # n=185
length(E4402.cases_T1) #185
E2408.cases_T1 <- sample_annot %>% filter(cohort == "E2408") %>% .$SAMPLE_ID # n=85
length(E2408.cases_T1) #83
PLOSMED.cases_T1 <- sample_annot %>% filter(cohort == "PLOSMED") %>% .$SAMPLE_ID #n=386
length(PLOSMED.cases_T1) #242
length(c(UHN.cases_T1, E4402.cases_T1, E2408.cases_T1, PLOSMED.cases_T1)) #713
# total: n=713

#--
# Generate Fig 1 panel A
#--

# Plot mutation frequency in all cohorts
plot.mut.freq.all.cohorts <- capseq_T1_mat %>%
  mutate(Hugo_Symbol = row.names(.)) %>%
  pivot_longer(!Hugo_Symbol, names_to = "SAMPLE_ID") %>%
  filter(value == "1") %>%
  dplyr::count(Hugo_Symbol) %>%
  arrange(desc(n)) %>%
  mutate(proportion_mutated = n/ncol(capseq_T1_mat)) %>%
  filter(proportion_mutated > 0.02) %>%
  dplyr::select(Hugo_Symbol, n, proportion_mutated) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = proportion_mutated)) +
  geom_bar(stat = "identity") + 
  ylab("% of mutated samples") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic", size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_blank())

# Plot mutation frequency (all mutations) in all cohorts
order.genes <- capseq_T1_mat %>%
  mutate(Hugo_Symbol = row.names(.)) %>%
  pivot_longer(!Hugo_Symbol, names_to = "SAMPLE_ID") %>%
  filter(value == "1") %>%
  dplyr::count(Hugo_Symbol) %>%
  arrange(desc(n)) %>%
  mutate(proportion_mutated = n/ncol(capseq_T1_mat)) %>%
  filter(proportion_mutated > 0.02) %>%
  dplyr::select(Hugo_Symbol, n, proportion_mutated) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  .$Hugo_Symbol

plot.mut.freq.all.cohorts.all.muts <- capseq_T1_maf %>%
  group_by(Hugo_Symbol, Variant_Classification) %>%
  summarize(count = n()) %>%
  filter(Hugo_Symbol %in% order.genes) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = order.genes)) %>%
  mutate(Variant_Classification = ifelse(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins"), "Frame shift indel",
                                  ifelse(Variant_Classification %in% c("In_Frame_Del", "In_Frame_Ins"), "In frame indel", 
                                  ifelse(Variant_Classification == "Nonsense_Mutation", "Nonsense mutation",
                                  ifelse(Variant_Classification == "Missense_Mutation", "Missense mutation",
                                  ifelse(Variant_Classification == "Splice_Site", "Splice site", Variant_Classification)))))) %>%
  mutate(Variant_Classification = factor(Variant_Classification, 
                                         levels = c("Nonsense mutation", "Frame shift indel", "Splice site", "In frame indel", "Missense mutation"))) %>%
  ggplot(aes(x = Hugo_Symbol, y = count, fill = Variant_Classification)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Nonsense mutation" = "#b2182b", "Frame shift indel" = "#d6604d", "Splice site" = "#fddbc7",
                               "In frame indel" = "#bababa", "Missense mutation" = "#878787")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, face = "italic", size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position = c(0.85, 0.55),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2)) +
  ylab("Number of mutations")

g <- arrangeGrob(plot.mut.freq.all.cohorts, plot.mut.freq.all.cohorts.all.muts, nrow = 2, heights = c(0.8, 1))
ggsave(file = paste0("img/",  date, " Fig1A.pdf"), g, width = 16, height = 7, units = "cm")

#--
# Run Maftools
# (https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#1_Introduction)
#--

# Need to take into account those samples with no mutations
df_tmp <- cbind(Hugo_Symbol = rep(NA, length(samples_nomuts)),
                Chromosome = rep(NA, length(samples_nomuts)),
                Start_Position = rep(NA, length(samples_nomuts)),
                End_Position = rep(NA, length(samples_nomuts)),
                Reference_Allele = rep(NA, length(samples_nomuts)),
                Tumor_Seq_Allele2 = rep(NA, length(samples_nomuts)),
                Protein_Change = rep(NA, length(samples_nomuts)),
                Variant_Classification = rep(NA, length(samples_nomuts)),
                Variant_Type = rep(NA, length(samples_nomuts)),
                Tumor_Sample_Barcode = samples_nomuts,
                VAF = rep(NA, length(samples_nomuts)))

df_tmp <- data.frame(df_tmp)
capseq_T1_maf_ <- capseq_T1_maf %>% rbind(df_tmp) #4088 variants

#all cohorts
maf = read.maf(maf = capseq_T1_maf_)

#multicentre
uhn_capseq_T1_maf_ <- capseq_T1_maf_ %>% 
  dplyr::filter(Tumor_Sample_Barcode %in% UHN.cases_T1) #1045 variants
uhn_maf = read.maf(maf = uhn_capseq_T1_maf_)

#E4402
e4402_capseq_T1_maf_ <- capseq_T1_maf_ %>% 
  dplyr::filter(Tumor_Sample_Barcode %in% E4402.cases_T1) #1028 variants
e4402_maf = read.maf(maf = e4402_capseq_T1_maf_)

#E2408
e2408_capseq_T1_maf_ <- capseq_T1_maf_ %>% 
  dplyr::filter(Tumor_Sample_Barcode %in% E2408.cases_T1) #516 variants
e2408_maf = read.maf(maf = e2408_capseq_T1_maf_)

#PLOSMED
plosmed_capseq_T1_maf_ <- capseq_T1_maf_ %>% 
  dplyr::filter(Tumor_Sample_Barcode %in% PLOSMED.cases_T1) #1499 variants
plosmed_maf = read.maf(maf = plosmed_capseq_T1_maf_)

#--
# Generate Fig 1 panel B - somatic interactions
#--

#load in the modified funtion
source("~/Downloads/Re_ Fig 1/2023-03-11_modifed_discoversomaticInteractions_VS.R")

#please desginate which cohort we are plotting here,
 #and the directory for the output table
cohort <- "all-cohorts"
dir <- "~/github/FLOMICS/"

pdf(file = paste0("img/",  date, " Fig 1B.pdf"), width = 7, height = 7)
all_cohorts_somaticInteractions <- MODdiscoversomaticInteractions(maf = maf, top = 47,
                                                   pvalue = c(0.001, 0.001), # not show p value annotation in tiles
                                                   getMutexMixed = FALSE, fontSize = 0.4, 
                                                   returnAll = TRUE, showSum = FALSE, showCounts = FALSE)
all_cohorts_somaticInteractions
dev.off()

all_cohorts_somaticInteractions <- data.frame(all_cohorts_somaticInteractions)
all_cohorts_somaticInteractions$fdr = p.adjust(all_cohorts_somaticInteractions$pValue, method = "fdr")

#--
# Generate Fig 1 panel B - Odds ratios
#--

### collect odds ratios for the individual cohorts
cohort <- "Multicentre"
uhn_somaticInteractions <- MODdiscoversomaticInteractions(maf = uhn_maf, top = 47,
                                                                  pvalue = c(0.001, 0.001), # not show p value annotation in tiles
                                                                  getMutexMixed = FALSE, fontSize = 0.4, 
                                                                  returnAll = TRUE, showSum = FALSE, showCounts = FALSE)
uhn_somaticInteractions <- data.frame(uhn_somaticInteractions)
uhn_somaticInteractions$fdr = p.adjust(uhn_somaticInteractions$pValue, method = "fdr")

cohort <- "E4402"
e4402_somaticInteractions <- MODdiscoversomaticInteractions(maf = e4402_maf, top = 47,
                                                          pvalue = c(0.001, 0.001), # not show p value annotation in tiles
                                                          getMutexMixed = FALSE, fontSize = 0.4, 
                                                          returnAll = TRUE, showSum = FALSE, showCounts = FALSE)
e4402_somaticInteractions <- data.frame(e4402_somaticInteractions)
e4402_somaticInteractions$fdr = p.adjust(e4402_somaticInteractions$pValue, method = "fdr")

cohort <- "E2408"
e2408_somaticInteractions <- MODdiscoversomaticInteractions(maf = e2408_maf, top = 47,
                                                            pvalue = c(0.001, 0.001), # not show p value annotation in tiles
                                                            getMutexMixed = FALSE, fontSize = 0.4, 
                                                            returnAll = TRUE, showSum = FALSE, showCounts = FALSE)
e2408_somaticInteractions <- data.frame(e2408_somaticInteractions)
e2408_somaticInteractions$fdr = p.adjust(e2408_somaticInteractions$pValue, method = "fdr")

cohort <- "PLOSMED"
plosmed_somaticInteractions <- MODdiscoversomaticInteractions(maf = plosmed_maf, top = 47,
                                                            pvalue = c(0.001, 0.001), # not show p value annotation in tiles
                                                            getMutexMixed = FALSE, fontSize = 0.4, 
                                                            returnAll = TRUE, showSum = FALSE, showCounts = FALSE)
plosmed_somaticInteractions <- data.frame(plosmed_somaticInteractions)
plosmed_somaticInteractions$fdr = p.adjust(plosmed_somaticInteractions$pValue, method = "fdr")

### extract the odds ratios from each somatic interactions analysis
gna13_mef2b_dat <- rbind(
  (all_cohorts_somaticInteractions %>% 
  dplyr::filter(gene1 == 'GNA13', gene2 == 'MEF2B') %>%
    dplyr::mutate(cohort = "Combined", color = "combined", interaction = paste0(gene1, "-", gene2))),
  (uhn_somaticInteractions %>% 
     dplyr::filter(gene1 == 'GNA13', gene2 == 'MEF2B') %>%
    dplyr::mutate(cohort = "Multicentre",  color = "individual", interaction = paste0(gene1, "-", gene2))),
  (e4402_somaticInteractions %>% 
     dplyr::filter(gene1 == 'GNA13', gene2 == 'MEF2B') %>%
     dplyr::mutate(cohort = "E4402", color = "individual", interaction = paste0(gene1, "-", gene2))),
  (e2408_somaticInteractions %>% 
     dplyr::filter(gene1 == 'GNA13', gene2 == 'MEF2B') %>%
     dplyr::mutate(cohort = "E2408",  color = "individual", interaction = paste0(gene1, "-", gene2))),
  (plosmed_somaticInteractions %>% 
     dplyr::filter(gene1 == 'GNA13', gene2 == 'MEF2B') %>%
     dplyr::mutate(cohort = "PLOSMED",  color = "individual", interaction = paste0(gene1, "-", gene2)))
  )

crebbp_ep300_dat <- rbind(
  (all_cohorts_somaticInteractions %>% 
     dplyr::filter(gene1 == 'EP300', gene2 == 'CREBBP') %>%
     dplyr::mutate(cohort = "Combined", color = "combined", interaction = paste0(gene2, "-", gene1))),
  (uhn_somaticInteractions %>% 
     dplyr::filter(gene1 == 'EP300', gene2 == 'CREBBP') %>%
     dplyr::mutate(cohort = "Multicentre",  color = "individual", interaction = paste0(gene2, "-", gene1))),
  (e4402_somaticInteractions %>% 
     dplyr::filter(gene1 == 'EP300', gene2 == 'CREBBP') %>%
     dplyr::mutate(cohort = "E4402", color = "individual", interaction = paste0(gene2, "-", gene1))),
  (e2408_somaticInteractions %>% 
     dplyr::filter(gene1 == 'EP300', gene2 == 'CREBBP') %>%
     dplyr::mutate(cohort = "E2408",  color = "individual", interaction = paste0(gene2, "-", gene1))),
  (plosmed_somaticInteractions %>% 
     dplyr::filter(gene1 == 'EP300', gene2 == 'CREBBP') %>%
     dplyr::mutate(cohort = "PLOSMED",  color = "individual", interaction = paste0(gene2, "-", gene1)))
)

card11_ezh2_dat <- rbind(
  (all_cohorts_somaticInteractions %>% 
     dplyr::filter(gene1 == 'CARD11', gene2 == 'EZH2') %>%
     dplyr::mutate(cohort = "Combined", color = "combined", interaction = paste0(gene1, "-", gene2))),
  (uhn_somaticInteractions %>% 
     dplyr::filter(gene1 == 'CARD11', gene2 == 'EZH2') %>%
     dplyr::mutate(cohort = "Multicentre",  color = "individual", interaction = paste0(gene1, "-", gene2))),
  (e4402_somaticInteractions %>% 
     dplyr::filter(gene1 == 'CARD11', gene2 == 'EZH2') %>%
     dplyr::mutate(cohort = "E4402", color = "individual", interaction = paste0(gene1, "-", gene2))),
  (e2408_somaticInteractions %>% 
     dplyr::filter(gene1 == 'CARD11', gene2 == 'EZH2') %>%
     dplyr::mutate(cohort = "E2408",  color = "individual", interaction = paste0(gene1, "-", gene2))),
  (plosmed_somaticInteractions %>% 
     dplyr::filter(gene1 == 'CARD11', gene2 == 'EZH2') %>%
     dplyr::mutate(cohort = "PLOSMED",  color = "individual", interaction = paste0(gene1, "-", gene2)))
)

hvcn1_mef2b_dat <- rbind(
  (all_cohorts_somaticInteractions %>% 
     dplyr::filter(gene1 == 'HVCN1', gene2 == 'MEF2B') %>%
     dplyr::mutate(cohort = "Combined", color = "combined", interaction = paste0(gene1, "-", gene2))),
  (uhn_somaticInteractions %>% 
     dplyr::filter(gene1 == 'HVCN1', gene2 == 'MEF2B') %>%
     dplyr::mutate(cohort = "Multicentre",  color = "individual", interaction = paste0(gene1, "-", gene2))),
  (e4402_somaticInteractions %>% 
     dplyr::filter(gene1 == 'MEF2B', gene2 == 'HVCN1') %>%
     dplyr::mutate(cohort = "E4402", color = "individual", interaction = paste0(gene2, "-", gene1))),
  (e2408_somaticInteractions %>% 
     dplyr::filter(gene1 == 'MEF2B', gene2 == 'HVCN1') %>%
     dplyr::mutate(cohort = "E2408",  color = "individual", interaction = paste0(gene2, "-", gene1))),
  (plosmed_somaticInteractions %>% 
     dplyr::filter(gene1 == 'HVCN1', gene2 == 'MEF2B') %>%
     dplyr::mutate(cohort = "PLOSMED",  color = "individual", interaction = paste0(gene1, "-", gene2)))
)

data <- as.data.frame(rbind(gna13_mef2b_dat, crebbp_ep300_dat, card11_ezh2_dat, hvcn1_mef2b_dat))
data$OddsRatio <- as.numeric(data$OddsRatio)


data %>% mutate(cohort = factor(cohort, levels = c("Combined", "E4402", "E2408", "PLOSMED", "Multicentre"))) %>%
  mutate(interaction = factor(interaction, levels = c("CREBBP-EP300", "GNA13-MEF2B","CARD11-EZH2", "HVCN1-MEF2B"))) %>%
  ggplot(aes(x = cohort, y = OddsRatio, ymin = -8, ymax = 8, color = color)) +
  geom_point(aes(col = color)) +
  geom_hline(aes(yintercept = 1)) +
  coord_flip() + theme_bw() + theme(legend.position = "none") +
  scale_colour_manual(values = c("individual" = "#878787", "combined" = "#d6604d")) +
  #scale_y_continuous(trans = "log2", labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_discrete(position = "top") +
  facet_wrap(~interaction, strip.position = "left", nrow = 2) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 5),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 5)) +
  ylab("Odds Ratio")
ggsave(paste0("img/",  date, " Fig1B_OR.pdf"), width = 8, height = 6, units = "cm")

### ---

# 
# x1 <- fisher.test(table(capseq_T1_mat["GNA13",UHN.cases_T1], capseq_T1_mat["MEF2B",UHN.cases_T1], useNA = "ifany"))
# x2 <- fisher.test(table(capseq_T1_mat["GNA13",E4402.cases_T1], capseq_T1_mat["MEF2B",E4402.cases_T1], useNA = "ifany"))
# x3 <- fisher.test(table(capseq_T1_mat["GNA13",E2408.cases_T1], capseq_T1_mat["MEF2B",E2408.cases_T1], useNA = "ifany"))
# x4 <- fisher.test(table(capseq_T1_mat["GNA13",PLOSMED.cases_T1], capseq_T1_mat["MEF2B",PLOSMED.cases_T1], useNA = "ifany"))
# x5 <- fisher.test(table(capseq_T1_mat["GNA13",], capseq_T1_mat["MEF2B",], useNA = "ifany"))
# 
# y1 <- fisher.test(table(capseq_T1_mat["CREBBP",UHN.cases_T1], capseq_T1_mat["EP300",UHN.cases_T1], useNA = "ifany"))
# y2 <- fisher.test(table(capseq_T1_mat["CREBBP",E4402.cases_T1], capseq_T1_mat["EP300",E4402.cases_T1], useNA = "ifany"))
# y3 <- fisher.test(table(capseq_T1_mat["CREBBP",E2408.cases_T1], capseq_T1_mat["EP300",E2408.cases_T1], useNA = "ifany"))
# y4 <- fisher.test(table(capseq_T1_mat["CREBBP",PLOSMED.cases_T1], capseq_T1_mat["EP300",PLOSMED.cases_T1], useNA = "ifany"))
# y5 <- fisher.test(table(capseq_T1_mat["CREBBP",], capseq_T1_mat["EP300",], useNA = "ifany"))
# 
# z1 <- fisher.test(table(capseq_T1_mat["CARD11",UHN.cases_T1], capseq_T1_mat["EZH2",UHN.cases_T1], useNA = "ifany"))
# z2 <- fisher.test(table(capseq_T1_mat["CARD11",E4402.cases_T1], capseq_T1_mat["EZH2",E4402.cases_T1], useNA = "ifany"))
# z3 <- fisher.test(table(capseq_T1_mat["CARD11",E2408.cases_T1], capseq_T1_mat["EZH2",E2408.cases_T1], useNA = "ifany"))
# z4 <- fisher.test(table(capseq_T1_mat["CARD11",PLOSMED.cases_T1], capseq_T1_mat["EZH2",PLOSMED.cases_T1], useNA = "ifany"))
# z5 <- fisher.test(table(capseq_T1_mat["CARD11",], capseq_T1_mat["EZH2",], useNA = "ifany"))
# 
# u1 <- fisher.test(table(capseq_T1_mat["HVCN1",UHN.cases_T1], capseq_T1_mat["MEF2B",UHN.cases_T1], useNA = "ifany"))
# u2 <- fisher.test(table(capseq_T1_mat["HVCN1",E4402.cases_T1], capseq_T1_mat["MEF2B",E4402.cases_T1], useNA = "ifany"))
# u3 <- fisher.test(table(capseq_T1_mat["HVCN1",E2408.cases_T1], capseq_T1_mat["MEF2B",E2408.cases_T1], useNA = "ifany"))
# u4 <- fisher.test(table(capseq_T1_mat["HVCN1",PLOSMED.cases_T1], capseq_T1_mat["MEF2B",PLOSMED.cases_T1], useNA = "ifany"))
# u5 <- fisher.test(table(capseq_T1_mat["HVCN1",], capseq_T1_mat["MEF2B",], useNA = "ifany"))
# 
# dat <- rbind(c("Multicentre", x1$estimate, x1$conf.int, x1$p.value, "individual", "GNA13-MEF2B"),
#              c("E4402", x2$estimate, x2$conf.int, x2$p.value, "individual", "GNA13-MEF2B"),
#              c("E2408", x3$estimate, x3$conf.int, x3$p.value, "individual", "GNA13-MEF2B"),
#              c("PLOSMED", x4$estimate, x4$conf.int, x4$p.value, "individual", "GNA13-MEF2B"),
#              c("Combined", x5$estimate, x5$conf.int, x5$p.value, "combined", "GNA13-MEF2B"),
#              c("Multicentre", y1$estimate, y1$conf.int, y1$p.value, "individual", "CREBBP-EP300"),
#              c("E4402", y2$estimate, y2$conf.int, y2$p.value, "individual", "CREBBP-EP300"),
#              c("E2408", y3$estimate, y3$conf.int, y3$p.value, "individual", "CREBBP-EP300"),
#              c("PLOSMED", y4$estimate, y4$conf.int, y4$p.value, "individual", "CREBBP-EP300"),
#              c("Combined", y5$estimate, y5$conf.int, y5$p.value, "combined", "CREBBP-EP300"),
#              c("Multicentre", z1$estimate, z1$conf.int, z1$p.value, "individual", "CARD11-EZH2"),
#              c("E4402", z2$estimate, z2$conf.int, z2$p.value, "individual", "CARD11-EZH2"),
#              c("E2408", z3$estimate, z3$conf.int, z3$p.value, "individual", "CARD11-EZH2"),
#              c("PLOSMED", z4$estimate, z4$conf.int, z4$p.value, "individual", "CARD11-EZH2"),
#              c("Combined", z5$estimate, z5$conf.int, z5$p.value, "combined", "CARD11-EZH2"),
#              c("Multicentre", u1$estimate, u1$conf.int, u1$p.value, "individual", "HVCN1-MEF2B"),
#              c("E4402", u2$estimate, u2$conf.int, u2$p.value, "individual", "HVCN1-MEF2B"),
#              c("E2408", u3$estimate, u3$conf.int, u3$p.value, "individual", "HVCN1-MEF2B"),
#              c("PLOSMED", u4$estimate, u4$conf.int, u4$p.value, "individual", "HVCN1-MEF2B"),
#              c("Combined", u5$estimate, u5$conf.int, u5$p.value, "combined", "HVCN1-MEF2B"))
# colnames(dat) <- c("cohort", "OR", "lower", "higher", "P", "color", "interaction")
# dat <- data.frame(dat)
# dat$OR <- as.numeric(dat$OR)
# dat$lower <- as.numeric(dat$lower)
# dat$higher <- as.numeric(dat$higher)
# 
# dat %>% mutate(cohort = factor(cohort, levels = c("Combined", "E4402", "E2408", "PLOSMED", "Multicentre"))) %>%
#   mutate(interaction = factor(interaction, levels = c("CREBBP-EP300", "GNA13-MEF2B","CARD11-EZH2", "HVCN1-MEF2B"))) %>%
#   ggplot(aes(x = cohort, y = OR, ymin = lower, ymax = higher, color = color)) +
#   geom_pointrange(aes(col = color), lwd = 0.4) +
#   geom_hline(aes(fill = cohort), yintercept = 1, linetype = 2) +
#   coord_flip() + theme_bw() + theme(legend.position = "none") +
#   scale_colour_manual(values = c("individual" = "#878787", "combined" = "#d6604d")) +
#   scale_y_continuous(trans = "log2", labels = function(x) ifelse(x == 0, "0", x)) +
#   scale_x_discrete(position = "top") +
#   facet_wrap(~interaction, strip.position = "left", nrow = 2) +
#   theme(axis.title.y = element_blank(),
#         axis.title.x = element_text(size = 5),
#         axis.text.x = element_text(size = 5),
#         axis.text.y = element_text(size = 5),
#         strip.background = element_blank(),
#         strip.text = element_text(face = "italic", size = 5)) +
#   ylab("Odds Ratio (95% Confidence Interval)")
# ggsave(paste0("img/",  date, " Fig1B_OR.pdf"), width = 8, height = 6, units = "cm")
# 
# scale_fill_manual(values = c("Nonsense mutation" = "#b2182b", "Frame shift indel" = "#d6604d", "Splice site" = "#fddbc7",
#                              "In frame indel" = "#bababa", "Missense mutation" = "#878787")) +
