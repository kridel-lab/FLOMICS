#---
# This script reads in mutation (capseq) data and plots all panels from Figure 1
# Author: Robert Kridel
#---

packages <- c("dplyr", "ggplot2", "maftools", "tidyverse", "gridExtra",
 "Rediscover")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

setwd("~/your working directory/FLOMICS/")

#--
# Read in data
#--

# Read in gene panels and define common genes
genes_plosmed <- read.csv(
   "DNAseq/Final_target_files/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
    filter(PLOS_MED_PANEL == "YES") %>%
   .$Gene.Name
genes_uhn <- read.csv(
   "DNAseq/Final_target_files/coding_gene_panel_PLOSMED_FLOMICS.csv") %>%
    filter(FLOMICS_PANEL == "YES") %>%
     .$Gene.Name
genes_common <- intersect(genes_plosmed, genes_uhn)

# Read in list of all mutations
capseq_t1 <- read.csv("filtered_MUTECT2_calls.csv") %>%
  mutate(TIMEPOINT = str_sub(SAMPLE_ID, -2, -1)) %>%
  filter(Hugo_Symbol %in% genes_common) %>%
  filter(SAMPLE_ID != "LY_FL_571_T1") %>%
  filter(!SAMPLE_ID %in% c("LY_FL_1135_T1", "LY_FL_1156_T1")) %>%
  filter(TIMEPOINT == "T1")

# Read in mutation matrix
capseq_t1_mat <- read.table("gene_vs_sample_SNV_matrix.txt",
                            sep = ";", header = TRUE)
capseq_t1_mat <- capseq_t1_mat[,
 !(names(capseq_t1_mat) %in% c("LY_FL_1135_T1", "LY_FL_1156_T1"))]

# Read in samples with no mutations but with sufficient coverage
samples_nomuts_goodcoverage <- read.table(
   "samples_without_SNVs_with_high_coverage_in_gene_vs_sample_mat.txt",
    sep = "\t", header = TRUE)
samples_nomuts <- names(capseq_t1_mat[, (colSums(capseq_t1_mat) == 0)])

# Read in sample annotations
sample_annot <- read.table(
   "metadata/20221228_sample_annotations.txt", sep = "\t", header = TRUE) %>%
  filter(SAMPLE_ID %in% capseq_t1$SAMPLE_ID | SAMPLE_ID %in% samples_nomuts) %>%
  mutate(cohort = ifelse(INSTITUTION %in%
   c("AARHUS", "JGH", "KINGSTON", "OSLO", "PAH", "UHN"), "UHN",
                  ifelse(INSTITUTION == "E4402", "E4402",
                  ifelse(INSTITUTION == "E2408", "E2408",
                  ifelse(INSTITUTION == "BCCA", "PLOSMED", NA)))))

#--
# Create maf file
#--

capseq_t1_maf <- capseq_t1 %>%
  dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position,
   Reference_Allele, Tumor_Seq_Allele2, Protein_Change = AAChange.ensGene,
         Variant_Classification, Variant_Type,
          Tumor_Sample_Barcode = SAMPLE_ID, VAF = Var_Freq) %>%
  mutate(Protein_Change = gsub("(.*),.*", "\\1", Protein_Change)) %>%
  mutate(Protein_Change = sub(".*p.", "", Protein_Change)) %>%
  filter(!Variant_Classification %in%
   c("UTR3 .", "UTR5 .", "intronic .", "upstream .")) %>%
  mutate(Variant_Classification =
   ifelse(Variant_Classification == "exonic nonsynonymous_SNV",
    "Missense_Mutation",
   ifelse(Variant_Classification == "exonic stopgain", "Nonsense_Mutation",
   ifelse(Variant_Classification == "exonic frameshift_deletion",
    "Frame_Shift_Del",
   ifelse(Variant_Classification == "exonic frameshift_insertion",
    "Frame_Shift_Ins",
   ifelse(Variant_Classification == "exonic nonframeshift_deletion",
    "In_Frame_Del",
   ifelse(Variant_Classification == "exonic nonframeshift_insertion",
    "In_Frame_Ins",
   ifelse(Variant_Classification == "exonic nonframeshift_substitution",
   "Missense_Mutation",
   ifelse(Variant_Classification == "exonic stoploss", "Missense_Mutation",
   ifelse(Variant_Classification == "splicing .",
    "Splice_Site", NA)))))))))) %>%
  mutate(Variant_Type = ifelse(Variant_Classification == "Frame_Shift_Del",
   "DEL", Variant_Type)) %>%
  mutate(Variant_Type = ifelse(Variant_Classification == "Frame_Shift_Ins",
   "INS", Variant_Type))

#--
# Define cases in individual cohorts
#--

uhn_cases_t1 <- sample_annot %>% filter(cohort == "UHN") %>% .$SAMPLE_ID
length(uhn_cases_t1)
e4402_cases_t1 <- sample_annot %>% filter(cohort == "E4402") %>% .$SAMPLE_ID
length(e4402_cases_t1)
e2408_cases_t1 <- sample_annot %>% filter(cohort == "E2408") %>% .$SAMPLE_ID
length(e2408_cases_t1)
plosmed_cases_t1 <- sample_annot %>%
 filter(cohort == "PLOSMED") %>%
 .$SAMPLE_ID
length(plosmed_cases_t1) #242
length(c(uhn_cases_t1, e4402_cases_t1, e2408_cases_t1, plosmed_cases_t1))

#--
# Generate Fig 1 panel B
#--

# Plot mutation frequency in all cohorts
plot.mut.freq.all.cohorts <- capseq_t1_mat %>%
  mutate(Hugo_Symbol = row.names(.)) %>%
  pivot_longer(!Hugo_Symbol, names_to = "SAMPLE_ID") %>%
  filter(value == "1") %>%
  dplyr::count(Hugo_Symbol) %>%
  arrange(desc(n)) %>%
  mutate(proportion_mutated = n / ncol(capseq_t1_mat)) %>%
  filter(proportion_mutated > 0.02) %>%
  dplyr::select(Hugo_Symbol, n, proportion_mutated) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = proportion_mutated)) +
  geom_bar(stat = "identity") +
  ylab("% of mutated samples") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
   face = "italic", size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_blank())

# Plot mutation frequency (all mutations) in all cohorts
order_genes <- capseq_t1_mat %>%
  mutate(Hugo_Symbol = row.names(.)) %>%
  pivot_longer(!Hugo_Symbol, names_to = "SAMPLE_ID") %>%
  filter(value == "1") %>%
  dplyr::count(Hugo_Symbol) %>%
  arrange(desc(n)) %>%
  mutate(proportion_mutated = n / ncol(capseq_t1_mat)) %>%
  filter(proportion_mutated > 0.02) %>%
  dplyr::select(Hugo_Symbol, n, proportion_mutated) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  .$Hugo_Symbol

plot.mut.freq.all.cohorts.all.muts <- capseq_t1_maf %>%
  group_by(Hugo_Symbol, Variant_Classification) %>%
  summarize(count = n()) %>%
  filter(Hugo_Symbol %in% order_genes) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = order_genes)) %>%
  mutate(Variant_Classification =
   ifelse(Variant_Classification %in%
   c("Frame_Shift_Del", "Frame_Shift_Ins"), "Frame shift indel",
   ifelse(Variant_Classification %in%
    c("In_Frame_Del", "In_Frame_Ins"), "In frame indel",
   ifelse(Variant_Classification == "Nonsense_Mutation", "Nonsense mutation",
   ifelse(Variant_Classification == "Missense_Mutation", "Missense mutation",
   ifelse(Variant_Classification == "Splice_Site", "Splice site",
    Variant_Classification)))))) %>%
  mutate(Variant_Classification = factor(Variant_Classification,
  levels = c("Nonsense mutation", "Frame shift indel", "Splice site",
   "In frame indel", "Missense mutation"))) %>%
  ggplot(aes(x = Hugo_Symbol, y = count, fill = Variant_Classification)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Nonsense mutation" = "#b2182b",
   "Frame shift indel" = "#d6604d", "Splice site" = "#fddbc7",
   "In frame indel" = "#bababa", "Missense mutation" = "#878787")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, face = "italic",
   size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position = c(0.85, 0.55),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_rect(fill = "white", color = "black",
         size = 0.2)) +
  ylab("Number of mutations")

g <- arrangeGrob(plot.mut.freq.all.cohorts, plot.mut.freq.all.cohorts.all.muts,
 nrow = 2, heights = c(0.8, 1))
ggsave(file = paste0("img/",  date, " Fig1A.pdf"), g,
 width = 16, height = 7, units = "cm")

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
capseq_t1_maf_2 <- capseq_t1_maf %>% rbind(df_tmp)

#all cohorts
maf <- read.maf(maf = capseq_t1_maf_2)

#multicentre
uhn_capseq_t1_maf <- capseq_t1_maf_2 %>%
  dplyr::filter(Tumor_Sample_Barcode %in% uhn_cases_t1)
uhn_maf <- read.maf(maf = uhn_capseq_t1_maf)

#E4402
e4402_capseq_t1_maf <- capseq_t1_maf_2 %>%
  dplyr::filter(Tumor_Sample_Barcode %in% e4402_cases_t1)
e4402_maf <- read.maf(maf = e4402_capseq_t1_maf)

#E2408
e2408_capseq_t1_maf <- capseq_t1_maf_2 %>%
  dplyr::filter(Tumor_Sample_Barcode %in% e2408_cases_t1)
e2408_maf <- read.maf(maf = e2408_capseq_t1_maf)

#PLOSMED
plosmed_capseq_t1_maf <- capseq_t1_maf_2 %>%
  dplyr::filter(Tumor_Sample_Barcode %in% plosmed_cases_t1)
plosmed_maf <- read.maf(maf = plosmed_capseq_t1_maf)

#--
# Generate Fig 1 panel C - somatic interactions
#--

#load in the modified funtion
source("modifed_discoversomaticInteractions_VS.R")

#please desginate which cohort we are plotting here,
 #and the directory for the output table
cohort <- "all-cohorts"
dir <- "~/your working directory/FLOMICS/"

pdf(file = paste0("img/",  date, " Fig 1B.pdf"), width = 7, height = 7)
all_somatic_interactions <- MODdiscoversomaticInteractions(maf = maf, top = 47,
   pvalue = c(0.001, 0.001), # not show p value annotation in tiles
   getMutexMixed = FALSE, fontSize = 0.4,
   returnAll = TRUE, showSum = FALSE, showCounts = FALSE)
all_somatic_interactions
dev.off()

all_somatic_interactions <- data.frame(all_somatic_interactions)
all_somatic_interactions$fdr <- p.adjust(all_somatic_interactions$pValue,
 method = "fdr")

#--
# Generate Fig 1 panel C - Odds ratios
#--

### collect odds ratios for the individual cohorts
cohort <- "Multicentre"
uhn_somatic_interactions <- MODdiscoversomaticInteractions(maf = uhn_maf,
   top = 47,
   pvalue = c(0.001, 0.001), # not show p value annotation in tiles
   getMutexMixed = FALSE, fontSize = 0.4,
   returnAll = TRUE, showSum = FALSE, showCounts = FALSE)
uhn_somatic_interactions <- data.frame(uhn_somatic_interactions)
uhn_somatic_interactions$fdr <- p.adjust(uhn_somatic_interactions$pValue,
 method = "fdr")

cohort <- "E4402"
e4402_somatic_interactions <- MODdiscoversomaticInteractions(maf = e4402_maf,
   top = 47,
   pvalue = c(0.001, 0.001), # not show p value annotation in tiles
   getMutexMixed = FALSE, fontSize = 0.4,
   returnAll = TRUE, showSum = FALSE, showCounts = FALSE)
e4402_somatic_interactions <- data.frame(e4402_somatic_interactions)
e4402_somatic_interactions$fdr <- p.adjust(e4402_somatic_interactions$pValue,
 method = "fdr")

cohort <- "E2408"
e2408_somatic_interactions <- MODdiscoversomaticInteractions(maf = e2408_maf,
   top = 47,
   pvalue = c(0.001, 0.001), # not show p value annotation in tiles
   getMutexMixed = FALSE, fontSize = 0.4,
   returnAll = TRUE, showSum = FALSE, showCounts = FALSE)
e2408_somatic_interactions <- data.frame(e2408_somatic_interactions)
e2408_somatic_interactions$fdr <- p.adjust(e2408_somatic_interactions$pValue,
 method = "fdr")

cohort <- "PLOSMED"
plosmed_somatic_interactions <- MODdiscoversomaticInteractions(
   maf = plosmed_maf,
   top = 47,
   pvalue = c(0.001, 0.001), # not show p value annotation in tiles
   getMutexMixed = FALSE, fontSize = 0.4,
   returnAll = TRUE, showSum = FALSE, showCounts = FALSE)
plosmed_somatic_interactions <- data.frame(plosmed_somatic_interactions)
plosmed_somatic_interactions$fdr <- p.adjust(
   plosmed_somatic_interactions$pValue,
   method = "fdr")

### extract the odds ratios from each somatic interactions analysis
gna13_mef2b_dat <- rbind(
  (all_somatic_interactions %>%
  dplyr::filter(gene1 == "GNA13", gene2 == "MEF2B") %>%
      dplyr::mutate(cohort = "Combined", color = "combined",
         interaction = paste0(gene1, "-", gene2))),
         (uhn_somatic_interactions %>%
  dplyr::filter(gene1 == "GNA13", gene2 == "MEF2B") %>%
      dplyr::mutate(cohort = "Multicentre",  color = "individual",
         interaction = paste0(gene1, "-", gene2))),
         (e4402_somatic_interactions %>%
  dplyr::filter(gene1 == "GNA13", gene2 == "MEF2B") %>%
      dplyr::mutate(cohort = "E4402", color = "individual",
         interaction = paste0(gene1, "-", gene2))),
         (e2408_somatic_interactions %>%
  dplyr::filter(gene1 == "GNA13", gene2 == "MEF2B") %>%
      dplyr::mutate(cohort = "E2408",  color = "individual",
         interaction = paste0(gene1, "-", gene2))),
         (plosmed_somatic_interactions %>%
  dplyr::filter(gene1 == "GNA13", gene2 == "MEF2B") %>%
      dplyr::mutate(cohort = "PLOSMED",  color = "individual",
         interaction = paste0(gene1, "-", gene2)))
   )

crebbp_ep300_dat <- rbind(
  (all_somatic_interactions %>%
     dplyr::filter(gene1 == "EP300", gene2 == "CREBBP") %>%
     dplyr::mutate(cohort = "Combined", color = "combined",
         interaction = paste0(gene2, "-", gene1))),
         (uhn_somatic_interactions %>%
     dplyr::filter(gene1 == "EP300", gene2 == "CREBBP") %>%
     dplyr::mutate(cohort = "Multicentre",  color = "individual",
         interaction = paste0(gene2, "-", gene1))),
         (e4402_somatic_interactions %>%
     dplyr::filter(gene1 == "EP300", gene2 == "CREBBP") %>%
     dplyr::mutate(cohort = "E4402", color = "individual",
         interaction = paste0(gene2, "-", gene1))),
         (e2408_somatic_interactions %>%
     dplyr::filter(gene1 == "EP300", gene2 == "CREBBP") %>%
     dplyr::mutate(cohort = "E2408",  color = "individual",
         interaction = paste0(gene2, "-", gene1))),
         (plosmed_somatic_interactions %>%
     dplyr::filter(gene1 == "EP300", gene2 == "CREBBP") %>%
     dplyr::mutate(cohort = "PLOSMED",  color = "individual",
         interaction = paste0(gene2, "-", gene1)))
   )

card11_ezh2_dat <- rbind(
  (all_somatic_interactions %>%
     dplyr::filter(gene1 == "CARD11", gene2 == "EZH2") %>%
     dplyr::mutate(cohort = "Combined", color = "combined",
         interaction = paste0(gene1, "-", gene2))),
         (uhn_somatic_interactions %>%
     dplyr::filter(gene1 == "CARD11", gene2 == "EZH2") %>%
     dplyr::mutate(cohort = "Multicentre",  color = "individual",
         interaction = paste0(gene1, "-", gene2))),
         (e4402_somatic_interactions %>%
     dplyr::filter(gene1 == "CARD11", gene2 == "EZH2") %>%
     dplyr::mutate(cohort = "E4402", color = "individual",
         interaction = paste0(gene1, "-", gene2))),
         (e2408_somatic_interactions %>%
     dplyr::filter(gene1 == "CARD11", gene2 == "EZH2") %>%
     dplyr::mutate(cohort = "E2408",  color = "individual",
         interaction = paste0(gene1, "-", gene2))),
         (plosmed_somatic_interactions %>%
     dplyr::filter(gene1 == "CARD11", gene2 == "EZH2") %>%
     dplyr::mutate(cohort = "PLOSMED",  color = "individual",
         interaction = paste0(gene1, "-", gene2)))
   )

hvcn1_mef2b_dat <- rbind(
  (all_somatic_interactions %>%
     dplyr::filter(gene1 == "HVCN1", gene2 == "MEF2B") %>%
     dplyr::mutate(cohort = "Combined", color = "combined",
         interaction = paste0(gene1, "-", gene2))),
         (uhn_somatic_interactions %>%
     dplyr::filter(gene1 == "HVCN1", gene2 == "MEF2B") %>%
     dplyr::mutate(cohort = "Multicentre",  color = "individual",
         interaction = paste0(gene1, "-", gene2))),
         (e4402_somatic_interactions %>%
     dplyr::filter(gene1 == "MEF2B", gene2 == "HVCN1") %>%
     dplyr::mutate(cohort = "E4402", color = "individual",
         interaction = paste0(gene2, "-", gene1))),
         (e2408_somatic_interactions %>%
     dplyr::filter(gene1 == "MEF2B", gene2 == "HVCN1") %>%
     dplyr::mutate(cohort = "E2408",  color = "individual",
         interaction = paste0(gene2, "-", gene1))),
         (plosmed_somatic_interactions %>%
     dplyr::filter(gene1 == "HVCN1", gene2 == "MEF2B") %>%
     dplyr::mutate(cohort = "PLOSMED",  color = "individual",
         interaction = paste0(gene1, "-", gene2)))
   )

data <- as.data.frame(rbind(gna13_mef2b_dat, crebbp_ep300_dat, card11_ezh2_dat,
 hvcn1_mef2b_dat))
data$OddsRatio <- as.numeric(data$OddsRatio)


data %>%
 mutate(cohort = factor(cohort, levels = c("Combined", "E4402", "E2408",
  "PLOSMED", "Multicentre"))) %>%
 mutate(interaction = factor(interaction, levels = c("CREBBP-EP300",
  "GNA13-MEF2B", "CARD11-EZH2", "HVCN1-MEF2B"))) %>%
 ggplot(aes(x = cohort, y = OddsRatio, ymin = -8, ymax = 8, color = color)) +
  geom_point(aes(col = color)) +
  geom_hline(aes(yintercept = 1)) +
  coord_flip() + theme_bw() + theme(legend.position = "none") +
  scale_colour_manual(values = c(
   "individual" = "#878787", "combined" = "#d6604d")) +
  scale_x_discrete(position = "top") +
  facet_wrap(~interaction, strip.position = "left", nrow = 2) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 5),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 5)) +
  ylab("Odds Ratio")
ggsave(paste0("img/",  date, " Fig1B_OR.pdf"),
 width = 8, height = 6, units = "cm")
