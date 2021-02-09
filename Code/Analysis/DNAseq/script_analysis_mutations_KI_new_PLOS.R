packages <- c("dplyr", "ggplot2", "data.table", "EnvStats", "ggpubr")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

# setwd("~/github/FLOMICS/") <- FLOMICS teams folder

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

mut.FLOMICS.rescued <- fread("DNAseq/Mutation_calls_KI/08-13-2020/with_indels/2020-08-13_Mutect2_filtered_mutations_FL_wRNASeq_mut_info.txt") %>%
  mutate(Cohort = "FLOMICS") %>%
  filter(avsnp142 == "." | (avsnp142 != "." & cosmic68 != ".")) %>%
  filter(Var_Freq <= 0.1) %>%
  filter((Chromosome == "1" & Start_Position %in% c("150727482")) |                                                      # CTSS hotspot
         (Chromosome == "7" & Start_Position %in% c("148506437", "148506467", "148508727", "148508728")) |               # EZH2 hotspots
         (Chromosome == "8" & Start_Position %in% c("20074767", "20074768")) |                                           # ATP6V1B2 hotspots
         (Chromosome == "12" & Start_Position %in% c("57496660", "57496661", "57496662", "57498345", "57499079")))       # STAT6 hotspots

dim(mut.FLOMICS.rescued)

mut.FLOMICS <- mut.FLOMICS %>% rbind(mut.FLOMICS.rescued)

dim(mut.FLOMICS)
length(unique(mut.FLOMICS$Tumor_Sample_Barcode))
# n = 692 rows & 130 unique patients

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
  filter(mean_cov < 50) %>% .$External_identifier %>% as.character() #18 samples (16=T1 and 2=T2)
length(poor.coverage.samples)

# Mean coverage of remaining samples
coverage.samples %>%
  filter(!External_identifier %in% poor.coverage.samples) %>%
  summarize(mean = mean(mean_cov))

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
mut.PLOSMED <- read.csv("DNAseq/Mutation_calls_KI/PLOS_MED/clean_up_PLOS_mutations_mutect2_Feb2021.csv")
# n = 266 rows
mut.PLOSMED = mut.PLOSMED[,c("SAMPLE_ID", "Hugo_Symbol", "Cohort")]

# Merge FLOMICS and PLOSMED
mut.merged <- mut.FLOMICS %>%
  select(SAMPLE_ID = Tumor_Sample_Barcode, Hugo_Symbol, Cohort) %>%
  rbind(mut.PLOSMED) %>%
  filter(Hugo_Symbol %in% genes.common) # n = 867

# Percentage of mutations per gene in merged cohort (max 1 mutation per gene per sample counted)
mut.merged %>%
  mutate(TIME_POINT = substr(SAMPLE_ID, 11, 12)) %>%
  filter(TIME_POINT != "T2") %>%
  group_by(Hugo_Symbol, SAMPLE_ID) %>%
  summarize(count = n()) %>%
  group_by(Hugo_Symbol) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(percentage = count/136*100) %>%
  filter(percentage > 5) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = percentage)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic", angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#ggsave(paste0("img/",  date, " Mutation_percentage_FLOMICS_PLOSMED_merged.pdf"), width = 14, height = 8, units = "cm")
ggsave(paste0("Analysis-Files/",  date, " Mutation_percentage_FLOMICS_PLOSMED_merged.pdf"), width = 14, height = 8, units = "cm")

# Compare average nb of mutations by cohort (max 1 mutation per gene per sample counted)
mut.merged %>%
  group_by(SAMPLE_ID, Hugo_Symbol, Cohort) %>%
  summarize(count = n()) %>%
  group_by(Cohort) %>%
  summarize(count = n()) %>%
  mutate(mean.nb.mut.by.sample = ifelse(Cohort == "FLOMICS", count/131, count/length(all.samples.DNAseq.PLOSMED)))
# on average, PLOSMED cases have 6.06 mutations per sample, whereas FLOMICS have 4.31 mutations per sample
# possible reasons:
# - no min VAF filter in PLOSMED
# - adv st cases have likely more mutations on average as tumour content is higher

# Compare mutations by stage
# Need sample annotation
sample.annotation <- read.csv("metadata/sample_annotations_rcd6Nov2019.csv", header = T) %>%
  select(SAMPLE_ID, TIME_POINT, STAGE) %>%
  filter(SAMPLE_ID %in% c(all.samples.DNAseq.FLOMICS, all.samples.DNAseq.PLOSMED))

sample.no.mut <- data.frame(SAMPLE_ID  = "LY_FL_179_T1", Cohort = "FLOMICS", count = 0, TIME_POINT = "T1", STAGE = "LIMITED")

mut.merged %>%
  distinct(SAMPLE_ID, Hugo_Symbol, Cohort) %>%
  group_by(SAMPLE_ID, Cohort) %>%
  summarize(count = n()) %>%
  right_join(sample.annotation) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  filter(TIME_POINT == "T1") %>%
  data.frame() %>%
  rbind(sample.no.mut) %>%
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
  data.frame() %>%
  rbind(sample.no.mut) %>%
  ggpubr::ggboxplot(x = "STAGE", y = "count", color = "STAGE", add = "jitter")+
  stat_compare_means() + stat_n_text()
# still slight increase of nb of mutations in advanced vs limited, statistically significant (wilcoxon p-value = 0.038)

# Explore differences between cohorts, gene by gene, for advanced stage cases
tmp <- mut.merged %>%
  distinct(SAMPLE_ID, Hugo_Symbol, Cohort) %>%
  left_join(sample.annotation) %>%
  filter(TIME_POINT == "T1" & STAGE == "ADVANCED")

n.adv.FLOMICS <- tmp %>%
  filter(Cohort == "FLOMICS") %>%
  .$SAMPLE_ID %>% unique() %>% length()

n.adv.PLOSMED <- tmp %>%
  filter(Cohort == "PLOSMED") %>%
  .$SAMPLE_ID %>% unique() %>% length()

tmp %>%
  group_by(Hugo_Symbol, Cohort) %>%
  summarize(count = n()) %>%
  reshape2::dcast(Hugo_Symbol ~ Cohort) %>%
  mutate(percentage.FLOMICS = (FLOMICS/n.adv.FLOMICS)*100,
         percentage.PLOSMED = (PLOSMED/n.adv.PLOSMED)*100)

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

# Nb mutations by SNF cluster
# Compare mutations by stage
# Need sample annotation
cluster.annotation <- read.csv("Cluster Labels/InfiniumClust_SNF_tSeq_Labels_18Nov2020.csv", header = T) %>%
  select(SAMPLE_ID = ID, SNFClust) %>%
  filter(SAMPLE_ID %in% c(all.samples.DNAseq.FLOMICS, all.samples.DNAseq.PLOSMED))

sample.no.mut <- data.frame(SAMPLE_ID  = "LY_FL_179_T1", Cohort = "FLOMICS", count = 0, SNFClust = "1", TIME_POINT = "T1")

mut.merged %>%
  distinct(SAMPLE_ID, Hugo_Symbol, Cohort) %>%
  group_by(SAMPLE_ID, Cohort) %>%
  summarize(count = n()) %>%
  right_join(cluster.annotation) %>%
  filter(!is.na(SNFClust)) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  mutate(TIME_POINT = substr(SAMPLE_ID, 11, 12)) %>%
  filter(TIME_POINT == "T1") %>%
  data.frame() %>%
  rbind(sample.no.mut) %>%
  ggpubr::ggboxplot(x = "SNFClust", y = "count", color = "SNFClust", add = "jitter")+
  stat_compare_means() + stat_n_text()
#ggsave(paste0("img/",  date, " Mutation_count_by_SNFcluster.pdf"), width = 10, height = 8, units = "cm")
ggsave(paste0("Analysis-Files/",  date, " Mutation_percentage_FLOMICS_PLOSMED_merged.pdf"), width = 14, height = 8, units = "cm")

# higher nb of mutations in C2 vs C1

# Mutation matrix
no.mut.cases <- setdiff(sample.annotation$SAMPLE_ID, mut.merged$SAMPLE_ID)
# n = 2 without mutation calls "LY_FL_179_T1" "LY_FL_399_T2"

mut.merged.df.incomplete <- mut.merged %>%
  distinct(SAMPLE_ID, Hugo_Symbol) %>%
  left_join(sample.annotation) %>%
  select(SAMPLE_ID, Hugo_Symbol) %>%
  mutate(var = 1) %>%
  reshape2::dcast(Hugo_Symbol ~ SAMPLE_ID) %>%
  replace(is.na(.), 0)

no.mut.cases.df <-  data.frame(matrix(vector(), nrow(mut.merged.df.incomplete), 2,
                               dimnames = list(c(), no.mut.cases)),
                               stringsAsFactors = F)
no.mut.cases.df[is.na(no.mut.cases.df)] <- 0

mut.merged.df.T1.T2 <- cbind(mut.merged.df.incomplete, no.mut.cases.df)
mut.merged.df.T1.T2 <- mut.merged.df.T1.T2[,order(colnames(mut.merged.df.T1.T2))]

mut.merged.df.T1.T2.poor.cov.excl <- mut.merged.df.T1.T2[, c("Hugo_Symbol", setdiff(sample.annotation$SAMPLE_ID, poor.coverage.samples))]

T1.samples <- sample.annotation %>% filter(TIME_POINT == "T1") %>% .$SAMPLE_ID %>% as.character() # n = 154
T1.samples.poor.cov.excl <- setdiff(T1.samples, poor.coverage.samples) # n = 138

mut.merged.df.T1 <- mut.merged.df.T1.T2[,c("Hugo_Symbol", T1.samples)]
mut.merged.df.T1.poor.cov.excl <- mut.merged.df.T1.T2[,c("Hugo_Symbol", T1.samples.poor.cov.excl)]

write.csv(mut.merged.df.T1.T2, file = "DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.T2.csv", row.names = FALSE)
write.csv(mut.merged.df.T1.T2.poor.cov.excl, file = "DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.T2.poor.cov.excl.csv", row.names = FALSE)
write.csv(mut.merged.df.T1, file = "DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.csv", row.names = FALSE)
write.csv(mut.merged.df.T1.poor.cov.excl, file = "DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl.csv", row.names = FALSE)

# Added by Anjali - generate counts by stage
# Nb of mutations per gene in FLOMICS (max 1 mutation per gene per sample counted)
length(unique(mut.FLOMICS$sample_mut)) # 128

mut.FLOMICS$sample_mut <- gsub(".*LY", "LY", mut.FLOMICS$sample_mut)
matchingEntries <- match(mut.FLOMICS$sample_mut, ClinicalFile_T1$SAMPLE_ID)
ClinicalFile_T1$SAMPLE_ID[matchingEntries[! is.na(matchingEntries)]]
length(ClinicalFile_T1$SAMPLE_ID[matchingEntries[! is.na(matchingEntries)]]) # 440 entries
length(unique(ClinicalFile_T1$SAMPLE_ID[matchingEntries[! is.na(matchingEntries)]])) # 105

mut.FLOMICS.Stage <- data.frame(mut.FLOMICS[match(ClinicalFile_T1$SAMPLE_ID[matchingEntries[! is.na(matchingEntries)]], mut.FLOMICS$sample_mut), ],
                                STAGE = ClinicalFile_T1$STAGE[matchingEntries[! is.na(matchingEntries)]],
                                TYPE = ClinicalFile_T1$TYPE[matchingEntries[! is.na(matchingEntries)]],
                                TRANSLOCATION1418 = ClinicalFile_T1$TRANSLOC_14_18[matchingEntries[! is.na(matchingEntries)]])
mut.FLOMICS.Stage
dim(mut.FLOMICS.Stage) # 440  19
length(unique(mut.FLOMICS.Stage$sample_mut)) # 105


# All cases
# Of 547 total cases, 440 entries when matched with 170 cases from methylation
# Of these only 105 out of 440 are unique patient IDs;
# This agrees with FLOMICs only cohort of Robert.
mut.FLOMICS.Stage %>%
  group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  group_by(Hugo_Symbol) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Advanced stage patients
mut.FLOMICS.Stage %>%
  filter(STAGE == "ADVANCED") %>%
  group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  group_by(Hugo_Symbol) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Limited stage patients
mut.FLOMICS.Stage %>%
  filter(STAGE == "LIMITED") %>%
  group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  group_by(Hugo_Symbol) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# TYPE
mut.FLOMICS.Stage %>%
  filter(TYPE == "FL") %>%
  group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  group_by(Hugo_Symbol) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


mut.FLOMICS.Stage %>%
  filter(TRANSLOCATION1418 == "1") %>%
  group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  group_by(Hugo_Symbol) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mut.FLOMICS.Stage %>%
  filter(TRANSLOCATION1418 == "0") %>%
  group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  summarize(count = n()) %>%
  group_by(Hugo_Symbol) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
