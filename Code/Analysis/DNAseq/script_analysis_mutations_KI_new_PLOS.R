#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#script_analysis_mutations_KI_new_PLOS.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#this script takes in mutations from our FL cohort as well as mutation
#calls from BC PLOS Med paper n=31 patient samples run through Mutect2 by KI
#mutations from both sources are further filtered and combined

# FLOMICS refers to current, UHN-led project.
# PLOSMED refers to samples sequenced as part of Kridel and Chan et al, Plos Med, 2016
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#load packages
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

packages <- c("dplyr", "ggplot2", "data.table", "EnvStats", "ggpubr")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()
# setwd("~/github/FLOMICS/") <- FLOMICS teams folder

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

mut.PLOSMED = mut.PLOSMED[,c("SAMPLE_ID", "Hugo_Symbol", "Cohort")]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Analyze FLOMICS data in more detail
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
# Conclusion: no positive correlation between coverage and nb of mutations
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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#combine mutations FLOMICs and PLOS MED
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Merge FLOMICS and PLOSMED
mut.merged <- mut.FLOMICS %>%
  select(SAMPLE_ID = Tumor_Sample_Barcode, Hugo_Symbol, Cohort) %>%
  rbind(mut.PLOSMED) %>%
  filter(Hugo_Symbol %in% genes.common) # n = 949

# Percentage of mutations per gene in merged cohort (max 1 mutation per gene per sample counted)
mut.merged %>%
  mutate(TIME_POINT = substr(SAMPLE_ID, 11, 12)) %>%
  filter(TIME_POINT != "T2") %>%
  filter(!SAMPLE_ID %in% poor.coverage.samples) %>%
  group_by(Hugo_Symbol, SAMPLE_ID) %>%
  summarize(count = n()) %>%
  group_by(Hugo_Symbol) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(percentage = count/138*100) %>%
  filter(percentage > 5) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = percentage)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic", angle = 90, vjust = 0.5, hjust = 1, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black"),
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
# on average, PLOSMED cases have 6.19 mutations per sample, whereas FLOMICS have 4.80 mutations per sample
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
# still slight increase of nb of mutations in advanced vs limited, statistically significant (wilcoxon p-value = 0.045)

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
cluster.annotation <- read.csv("Cluster Labels/InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv", header = T) %>%
  select(SAMPLE_ID = ID, SNFClust = SNFClust10Feb2021) %>%
  filter(SAMPLE_ID %in% c(all.samples.DNAseq.FLOMICS, all.samples.DNAseq.PLOSMED))

# sample.no.mut <- data.frame(SAMPLE_ID  = "LY_FL_179_T1", Cohort = "FLOMICS", count = 0, SNFClust = "1", TIME_POINT = "T1")

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
  # rbind(sample.no.mut) %>%
  ggpubr::ggboxplot(x = "SNFClust", y = "count", color = "SNFClust", add = "jitter")+
  stat_compare_means() + stat_n_text()
#ggsave(paste0("img/",  date, " Mutation_count_by_SNFcluster.pdf"), width = 10, height = 8, units = "cm")
ggsave(paste0("Analysis-Files/",  date, " Mutation_percentage_FLOMICS_PLOSMED_merged.pdf"), width = 14, height = 8, units = "cm")

# higher nb of mutations in C1 vs C2

# Calculate mean nb of genes mutated in C1 and C2
n.C1 <- cluster.annotation %>% filter(SNFClust == "1") %>% nrow()
n.C2 <- cluster.annotation %>% filter(SNFClust == "2") %>% nrow()

mut.merged %>%
  distinct(SAMPLE_ID, Hugo_Symbol, Cohort) %>%
  right_join(cluster.annotation) %>%
  filter(!is.na(SNFClust)) %>%
  group_by(SNFClust) %>%
  summarize(count = n()) %>%
  mutate(mean.nb.genes.mut = ifelse(SNFClust == "1", count/n.C1, count/n.C2))

# Mutation matrix
no.mut.cases <- setdiff(sample.annotation$SAMPLE_ID, unique(mut.merged$SAMPLE_ID))
# n = 0

mut.merged.df.T1.T2 <- mut.merged %>%
# mut.merged.df.incomplete <- mut.merged %>%
  distinct(SAMPLE_ID, Hugo_Symbol) %>%
  left_join(sample.annotation) %>%
  select(SAMPLE_ID, Hugo_Symbol) %>%
  mutate(var = 1) %>%
  reshape2::dcast(Hugo_Symbol ~ SAMPLE_ID) %>%
  replace(is.na(.), 0)

# no.mut.cases.df <-  data.frame(matrix(vector(), nrow(mut.merged.df.incomplete), 2,
#                                dimnames = list(c(), no.mut.cases)),
#                                stringsAsFactors = F)
# no.mut.cases.df[is.na(no.mut.cases.df)] <- 0

# mut.merged.df.T1.T2 <- cbind(mut.merged.df.incomplete, no.mut.cases.df)
# mut.merged.df.T1.T2 <- mut.merged.df.T1.T2[,order(colnames(mut.merged.df.T1.T2))]
# mut.merged.df.T1.T2 <- mut.merged.df.T1.T2[,order(colnames(mut.merged.df.T1.T2))]

#add 4 genes that were removed from matrix because no mutations
#in them that passed filters
genes_missing = genes.common[which(!(genes.common %in% mut.merged.df.T1.T2$Hugo_Symbol))]
genes_missing_df = mut.merged.df.T1.T2[1:4,]
genes_missing_df$Hugo_Symbol = genes_missing
genes_missing_df[,2:ncol(genes_missing_df)] = 0

mut.merged.df.T1.T2 = rbind(mut.merged.df.T1.T2, genes_missing_df)

mut.merged.df.T1.T2.poor.cov.excl <- mut.merged.df.T1.T2[, c("Hugo_Symbol", setdiff(sample.annotation$SAMPLE_ID, poor.coverage.samples))]

T1.samples <- sample.annotation %>% filter(TIME_POINT == "T1") %>% .$SAMPLE_ID %>% as.character() # n = 154
T1.samples.poor.cov.excl <- setdiff(T1.samples, poor.coverage.samples) # n = 138

mut.merged.df.T1 <- mut.merged.df.T1.T2[,c("Hugo_Symbol", T1.samples)]
mut.merged.df.T1.poor.cov.excl <- mut.merged.df.T1.T2[,c("Hugo_Symbol", T1.samples.poor.cov.excl)]

write.csv(mut.merged.df.T1.T2, file = "DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.T2.csv", row.names = FALSE)
write.csv(mut.merged.df.T1.T2.poor.cov.excl, file = "DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.T2.poor.cov.excl.csv", row.names = FALSE)
write.csv(mut.merged.df.T1, file = "DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.csv", row.names = FALSE)
write.csv(mut.merged.df.T1.poor.cov.excl, file = "DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl.csv", row.names = FALSE)

# Create file for mutation mapper with FLOMICS calls (can be uploaded here: https://www.cbioportal.org/mutation_mapper):
mut.FLOMICS.mutation.mapper <- mut.FLOMICS %>%
  select(Sample_ID = Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, Reference_Allele = Reference_Allele, Variant_Allele = Tumor_Seq_Allele2)
write.table(mut.FLOMICS.mutation.mapper, file = "DNAseq/Mutation_and_BA_matrices/mut.FLOMICS.mutation.mapper.txt", sep = "\t", row.names = FALSE, quote = FALSE)
