
packages <- c("dplyr", "ggplot2", "data.table", "EnvStats", "ggpubr")
lapply(packages, require, character.only = TRUE)

# setwd("~/FLOMICS/") <- FLOMICS teams folder

# This script explores mutation data

# FLOMICS refers to current, UHN-led project.
# PLOSMED refers to samples sequenced as part of Kridel and Chan et al, Plos Med, 2016

# Gene panel
gene.panel <- read.csv("DNAseq/Final_target_files/coding_gene_panel_PLOSMED_FLOMICS.csv", header = TRUE)
genes.common <- gene.panel %>%
  filter(PLOS_MED_PANEL == "YES" & FLOMICS_PANEL == "YES") %>%
  .$Gene.Name %>% as.character() # n = 57

# All FLOMICS samples included
all.samples.DNAseq.FLOMICS <- read.csv("metadata/sample_annotations_rcd6Nov2019.csv", header = T) %>%
  filter(CAPSEQ_DATA == TRUE & CAPSEQ_DATA_YEAR == "2019" & CAPSEQ_INCLUDE == "YES") %>%
  .$SAMPLE_ID %>% as.character() # n=131 (n=123 T1, n=8 T2)

# All PLOSMED samples included
all.samples.DNAseq.PLOSMED <- read.csv("metadata/sample_annotations_rcd6Nov2019.csv", header = T) %>%
  filter(CAPSEQ_DATA == TRUE & CAPSEQ_DATA_YEAR == "2015" & CAPSEQ_INCLUDE == "YES") %>%
  .$SAMPLE_ID %>% as.character() # n=31 (only T1)

# Coverage FLOMICS DNAseq
coverage.samples <- read.csv("DNAseq/Targeted_Seq_Coverage/08-07-2020/2020-08-07_sample_based_coverage_summary.csv")

# Read in mutation calls for FLOMICS
mut.FLOMICS <- fread("DNAseq/Mutation_calls_KI/08-13-2020/with_indels/2020-08-13_Mutect2_filtered_mutations_FL_wRNASeq_mut_info.txt") %>%
  mutate(Cohort = "FLOMICS") %>%
  filter(avsnp142 == "." | (avsnp142 != "." & cosmic68 != ".")) %>%
  filter(Var_Freq > 0.1)
dim(mut.FLOMICS)
length(unique(mut.FLOMICS$Tumor_Sample_Barcode))
# n = 679 rows & 129 unique patients

# Plot nb of mutations per sample in FLOMICS
mut.FLOMICS %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  full_join(coverage.samples, by = c("Tumor_Sample_Barcode" = "External_identifier")) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  ggplot(aes(x = count)) +
  geom_histogram(bins = 30)

# Plot nb of mutations vs. coverage in FLOMICS
mut.FLOMICS %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  full_join(coverage.samples, by = c("Tumor_Sample_Barcode" = "External_identifier")) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  ggplot(aes(x = mean_cov, y = count)) + geom_point() + geom_smooth(method = lm)
# Conclusion: no correlation between coverage and nb of mutations
# some of the very low coverage cases have high nb of mutations, however
# this means we should probably filter such cases out.

# Identify poor coverage samples and filter them out
poor.coverage.samples <- coverage.samples %>%
  filter(mean_cov < 50) %>% .$External_identifier %>% as.character() #18 samples
length(poor.coverage.samples)

# Nb of mutations per gene in FLOMICS (max 1 mutation per gene per sample counted)
mut.FLOMICS %>%
  group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  group_by(Hugo_Symbol) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = count)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Read in mutation calls for PLOSMED
mut.PLOSMED <- read.csv("DNAseq/BC_Cancer_capseq_data/BC_Cancer_capseq_data.csv", header = T) %>%
  mutate(Cohort = "PLOSMED") %>%
  select(SAMPLE_ID, Hugo_Symbol, Cohort)
# n = 380 rows

# Merge FLOMICS and PLOSMED
mut.merged <- mut.FLOMICS %>%
  select(SAMPLE_ID = Tumor_Sample_Barcode, Hugo_Symbol, Cohort) %>%
  rbind(mut.PLOSMED) %>%
  filter(Hugo_Symbol %in% genes.common) # n = 932

# Compare average nb of mutations by cohort (max 1 mutation per gene per sample counted)
mut.merged %>%
  group_by(SAMPLE_ID, Hugo_Symbol, Cohort) %>%
  summarize(count = n()) %>%
  group_by(Cohort) %>%
  summarize(count = n()) %>%
  mutate(mean.nb.mut.by.sample = ifelse(Cohort == "FLOMICS", count/131, count/length(all.samples.DNAseq.PLOSMED)))
# on average, PLOSMED cases have 7.32 mutations per sample, whereas FLOMICS have 4.21 mutations per sample
# possible reasons:
# - no min VAF filter in PLOSMED
# - adv st cases have likely more mutations on average as tumour content is higher

# Compare mutations by stage
# Need sample annotation
sample.annotation <- read.csv("metadata/sample_annotations_rcd6Nov2019.csv", header = T) %>%
  select(SAMPLE_ID, TIME_POINT, STAGE) %>%
  filter(SAMPLE_ID %in% c(all.samples.DNAseq.FLOMICS, all.samples.DNAseq.PLOSMED))

mut.merged %>%
  distinct(SAMPLE_ID, Hugo_Symbol, Cohort) %>%
  group_by(SAMPLE_ID, Cohort) %>%
  summarize(count = n()) %>%
  right_join(sample.annotation) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  filter(TIME_POINT == "T1") %>%
  ggpubr::ggboxplot(x = "STAGE", y = "count", color = "STAGE", add = "jitter")+
  stat_compare_means() + stat_n_text()
# higher nb of mutations in advanced vs limited

# subset on only FLOMICS cases
mut.merged %>%
  distinct(SAMPLE_ID, Hugo_Symbol, Cohort) %>%
  group_by(SAMPLE_ID, Cohort) %>%
  summarize(count = n()) %>%
  right_join(sample.annotation) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  filter(TIME_POINT == "T1") %>%
  filter(Cohort == "FLOMICS") %>%
  ggpubr::ggboxplot(x = "STAGE", y = "count", color = "STAGE", add = "jitter")+
  stat_compare_means() + stat_n_text()
# still slight increase of nb of mutations in advanced vs limited, although no statistically significant (wilcoxon p-value = 0.056)
# i.e. difference previously seen is in part artifact of differences between FLOMICS and PLOSMED
# but also smaller sample size now in advanced

# Explore differences between cohorts, gene by gene, for advanced stage cases
mut.merged %>%
  distinct(SAMPLE_ID, Hugo_Symbol, Cohort) %>%
  left_join(sample.annotation) %>%
  filter(TIME_POINT == "T1" & STAGE == "ADVANCED") %>%
  group_by(Hugo_Symbol, Cohort) %>%
  summarize(count = n()) %>%
  reshape2::dcast(Hugo_Symbol ~ Cohort) %>%
  mutate(percentage.FLOMICS = (FLOMICS/45)*100, # n=44 if filtered based on poor coverage
         percentage.PLOSMED = (PLOSMED/31)*100)

# Explore differences between stages, gene by gene
n.lim <- sample.annotation %>%
  filter(SAMPLE_ID %in% all.samples.DNAseq.FLOMICS & TIME_POINT == "T1" & STAGE == "LIMITED") %>%
  nrow()

n.adv <- sample.annotation %>%
  filter(SAMPLE_ID %in% c(all.samples.DNAseq.FLOMICS, all.samples.DNAseq.PLOSMED) & TIME_POINT == "T1" & STAGE == "ADVANCED") %>%
  nrow()

fisher <- function(a, b, c, d) {
  data <- matrix(c(a,b,c,d), ncol = 2)
  c(P = fisher.test(data)$p.value,
    OR = fisher.test(data)$estimate,
    CI = fisher.test(data)$conf.int)
}

mut.lim.vs.adv <- mut.merged %>%
  distinct(SAMPLE_ID, Hugo_Symbol, Cohort) %>%
  left_join(sample.annotation) %>%
  filter(TIME_POINT == "T1") %>%
  group_by(Hugo_Symbol, STAGE) %>%
  summarize(count = n()) %>%
  reshape2::dcast(Hugo_Symbol ~ STAGE) %>%
  mutate(ADVANCED = ifelse(is.na(ADVANCED), 0, ADVANCED)) %>%
  mutate(LIMITED = ifelse(is.na(LIMITED), 0, LIMITED)) %>%
  mutate(not.ADVANCED = n.adv - ADVANCED) %>%
  mutate(not.LIMITED = n.lim - LIMITED) %>%
  mutate(percentage.ADVANCED = (ADVANCED/n.adv)*100) %>%
  mutate(percentage.LIMITED = (LIMITED/n.lim)*100) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(p = fisher(ADVANCED, not.ADVANCED, LIMITED, not.LIMITED)[[1]],
         OR = fisher(ADVANCED, not.ADVANCED, LIMITED, not.LIMITED)[[2]],
         CI1 = fisher(ADVANCED, not.ADVANCED, LIMITED, not.LIMITED)[[3]],
         CI2 = fisher(ADVANCED, not.ADVANCED, LIMITED, not.LIMITED)[[4]])

mut.lim.vs.adv$fdr <- p.adjust(mut.lim.vs.adv$p, method = "fdr")

# Mutation matrix
no.mut.cases <- setdiff(sample.annotation$SAMPLE_ID, mut.merged$SAMPLE_ID)
# n = 3 without mutation calls "LY_FL_179_T1" "LY_FL_181_T1" "LY_FL_399_T2"

mut.merged.df.incomplete <- mut.merged %>%
  distinct(SAMPLE_ID, Hugo_Symbol) %>%
  left_join(sample.annotation) %>%
  select(SAMPLE_ID, Hugo_Symbol) %>%
  mutate(var = 1) %>%
  reshape2::dcast(Hugo_Symbol ~ SAMPLE_ID) %>%
  replace(is.na(.), 0)

no.mut.cases.df <-  data.frame(matrix(vector(), nrow(mut.merged.matrix.incomplete), 3,
                               dimnames = list(c(), no.mut.cases)),
                               stringsAsFactors = F)

mut.merged.df.T1.T2 <- cbind(mut.merged.df.incomplete, no.mut.cases.df)
mut.merged.df.T1.T2 <- mut.merged.df.T1.T2[,order(colnames(mut.merged.df.T1.T2))]

mut.merged.df.T1.T2.poor.cov.excl <- mut.merged.df.T1.T2[, c("Hugo_Symbol", setdiff(sample.annotation$SAMPLE_ID, poor.coverage.samples))]

T1.samples <- sample.annotation %>% filter(TIME_POINT == "T1") %>% .$SAMPLE_ID %>% as.character() # n = 154
T1.samples.poor.cov.excl <- setdiff(T1.samples, poor.coverage.samples) # n = 138

mut.merged.df.T1 <- mut.merged.df.T1.T2[,c("Hugo_Symbol", T1.samples)]
mut.merged.df.T1.poor.cov.excl <- mut.merged.df.T1.T2[,c("Hugo_Symbol", T1.samples.poor.cov.excl)]

write.csv(mut.merged.df.T1.T2, file = "mut.merged.df.T1.T2.csv", row.names = FALSE)
write.csv(mut.merged.df.T1.T2.poor.cov.excl, file = "mut.merged.df.T1.T2.poor.cov.excl.csv", row.names = FALSE)
write.csv(mut.merged.df.T1, file = "mut.merged.df.T1.csv", row.names = FALSE)
write.csv(mut.merged.df.T1.poor.cov.excl, file = "mut.merged.df.T1.poor.cov.excl.csv", row.names = FALSE)
