
library(dplyr)
library(data.table)
library(ggpubr)
library(ggmosaic)

date <- Sys.Date()

setwd("~/github/FLOMICS/")

#--
# some RNAseq QC data
#--

tiers <- read.csv("RNAseq/qc/Tiers1to3RNAseqFLOMICS_15Jan2021ASilva.csv")
tier1 <- tiers$T1[!is.na(tiers$T1)] # n = 81
tier2 <- tiers$T2[!is.na(tiers$T2)] # n = 104
tier3 <- tiers$T3[!is.na(tiers$T3)] # n = 132

# Tier 1: 81 samples
# RNAseqSampleCufoffUniqMapReadCount = 10000000,
# RNAseqSampleCufoffUniqMapReadPercentage = 70,
# RNAseqSampleCufoffReadsMappedMultipleLoci = 20,
# RNAseqSampleCutoffRRNAcontam = 40
# 
# Tier 2: 104 samples
# RNAseqSampleCufoffUniqMapReadCount = 10000000,
# RNAseqSampleCufoffUniqMapReadPercentage = 50,
# RNAseqSampleCufoffReadsMappedMultipleLoci = 20,
# RNAseqSampleCutoffRRNAcontam = 40
# 
# Tier 3: 132 samples
# RNAseqSampleCufoffUniqMapReadCount = 10000000,
# RNAseqSampleCufoffUniqMapReadPercentage = 0,
# RNAseqSampleCufoffReadsMappedMultipleLoci = 100,
# RNAseqSampleCutoffRRNAcontam = 100

qc <- read.csv("RNAseq/qc/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")
hist(qc$Number.of.input.reads)
hist(qc$Uniquely.mapped)
hist(qc$Uniquely.mapped.reads..)
hist(qc$rRNAcontam)
plot(qc$Number.of.input.reads, qc$Uniquely.mapped)
plot(qc$Uniquely.mapped, qc$rRNAcontam)

p <- qc %>%
  filter(TYPE == "FL") %>%
  ggboxplot(x = "STAGE", y = "Uniquely.mapped", color = "STAGE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- qc %>%
  filter(STAGE %in% c("ADVANCED", "LIMITED")) %>%
  filter(SAMPLE_ID %in% tier3) %>%
  ggboxplot(x = "STAGE", y = "rRNAcontam", color = "STAGE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

#--
# Read in clonotypes
#--

clonotypes <- fread("Analysis-Files/Mixcr/_2020-07-30_MIXCR_OUTPUT_WITH_SAMPLE_ANNOTATION.csv", sep = "\t") %>%
  mutate(x = nchar(sample)) %>%
  mutate(sample = ifelse(x == 9, paste0(sample, "_T1"), sample)) %>%
  select(-x, -CLUSTER) 

clonotypes$LY_FL_ID <- stringr::str_replace(clonotypes$sample, "_T1", "")

clin <- read.csv("metadata/clinical_data_rcd11Aug2020.csv")

clonotypes <- clonotypes %>%
  left_join(clin[,c("LY_FL_ID", "ANN_ARBOR_STAGE")])

SNF.clust <- read.csv("Cluster Labels/InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv") %>%
  select(ID, SNFClust = SNFClust10Feb2021)

clonotypes <- clonotypes %>%
  left_join(SNF.clust[,c("ID", "SNFClust")], by = c("sample" = "ID"))

BCL2_rearrangements <- read.csv("DNAseq/Mutation_and_BA_matrices/BA.results.T1.csv")

clonotypes <- clonotypes %>%
  left_join(BCL2_rearrangements[,c("SAMPLE_ID", "BCL2_BA_consensus", "BCL6_BA_consensus")], by = c("sample" = "SAMPLE_ID"))

length(unique(clonotypes$sample))
# n = 130 samples

#--
# Plot relative abundance of dominant Ig heavy chain clonotypes
#--

clonotypes.rel.abundance <- list()
for (i in unique(clonotypes$sample)) {
  tmp <- clonotypes %>% filter(sample == i) %>% filter(grepl("IGH", allCHitsWithScore))
  sums <- sum(tmp$cloneCount)
  dominant <- tmp %>%  top_n(n = 1, wt = cloneCount) %>% .$cloneCount
  dominant <- head(dominant, n = 1)
  clonotypes.rel.abundance[[i]] <- c(sample = i, dominant.count = dominant, total.count = sums)
}

clonotypes.rel.abundance <- do.call(rbind, clonotypes.rel.abundance)
clonotypes.rel.abundance <- data.frame(clonotypes.rel.abundance)
clonotypes.rel.abundance <- clonotypes.rel.abundance %>%
  mutate(total.count = ifelse(dominant.count == 0, 0, total.count)) %>%
  mutate(total.count = as.numeric(total.count)) %>%
  mutate(dominant.count = as.numeric(dominant.count))
hist(clonotypes.rel.abundance$total.count, breaks = 100)
hist(clonotypes.rel.abundance$total.count, breaks = 1000, xlim = c(0,1000))
mean(clonotypes.rel.abundance$total.count)

clonotypes.rel.abundance <- clonotypes.rel.abundance %>%
  mutate(dominant.fraction = dominant.count/total.count)

plot(clonotypes.rel.abundance$dominant.fraction, clonotypes.rel.abundance$dominant.count)
plot(clonotypes.rel.abundance$dominant.fraction, clonotypes.rel.abundance$dominant.count, ylim = c(0,1000))

clonotypes.IGH.confident <- clonotypes.rel.abundance %>%
  # filter(dominant.fraction > 0.2) %>%
  filter(dominant.count > 10) %>%
  .$sample

#--
# Filter
#--

clonotypes <- clonotypes %>%
 filter(sample %in% clonotypes.IGH.confident) %>%
 filter(sample %in% tier3)

length(unique(clonotypes$sample)) # n = 110 samples

#--
# Pull out IGH sequences
#--

clonotypes.IGH <- clonotypes %>% filter(grepl("IGH", allCHitsWithScore)) # 5660
clonotypes.IGKL <- clonotypes %>% filter(grepl("IGK|IGL", allCHitsWithScore)) # 16115

#--
# Select dominant IGH clonotype for each sample
#--

clonotypes.IGH <- clonotypes.IGH %>% group_by(sample) %>%
  top_n(n = 1, wt = cloneCount)

clonotypes.IGKL <- clonotypes.IGKL %>% group_by(sample) %>%
  top_n(n = 1, wt = cloneCount)

#--
# Merge dominant IGH and dominant 
#--

clonotypes <- clonotypes.IGH %>% rbind(clonotypes.IGKL)

#--
# Annotate aaSeqCDR3 with glycosylation site info (Asn-X-Ser/Thr)
#--

clonotypes <- clonotypes %>% 
  mutate(tmp1 = ifelse(stringr::str_detect(aaSeqCDR3, "N.{1}[ST]"), "YES", "NO")) %>%
  mutate(tmp2 = ifelse(stringr::str_detect(aaSeqCDR3, "N[PDE][ST]"), "YES", "NO")) %>%
  mutate(aaSeqCDR3_glycosyl = ifelse(tmp1 == "YES" & tmp2 == "NO", "YES", "NO")) %>%
  select(-tmp1, -tmp2)

table(clonotypes$aaSeqCDR3_glycosyl, clonotypes$TYPE)
fisher.test(clonotypes$aaSeqCDR3_glycosyl, clonotypes$TYPE)

table(clonotypes$aaSeqCDR3_glycosyl, clonotypes$STAGE)
fisher.test(clonotypes$aaSeqCDR3_glycosyl, clonotypes$STAGE)

table(clonotypes$aaSeqCDR3_glycosyl, clonotypes$SNFClust)
fisher.test(clonotypes$aaSeqCDR3_glycosyl, clonotypes$SNFClust)

table(clonotypes$aaSeqCDR3_glycosyl, clonotypes$BCL2_BA_consensus)
fisher.test(clonotypes$aaSeqCDR3_glycosyl, clonotypes$BCL2_BA_consensus)

#--
# IgM vs others, by stage
#--

df <- clonotypes.IGH %>%
  mutate(isotype = substring(allCHitsWithScore, 1, 4)) %>%
  filter(!is.na(STAGE)) %>%
  group_by(STAGE, isotype) %>%
  dplyr::summarize(n = n()) %>%
  mutate(isotype = factor(isotype, levels = c("IGHM", "IGHG", "IGHD", "IGHA", "IGHE")))

df.wide <- df %>%
  tidyr::spread(STAGE, n) %>%
  replace(is.na(.), 0)

chisq.test(df.wide[,2:3])

g <- ggbarplot(df, x = "isotype", y = "n", label = df$n , lab.size = 2,
               fill = "isotype", facet.by = c("STAGE"), palette = c("jco")) +
  theme_bw() +
  theme(text = element_text(size = 10)) + ylab("Number of samples with given isotype") + xlab("Stage")
ggpar(g, legend = "bottom", legend.title = "isotype") + geom_hline(yintercept = 0)

##

df <- clonotypes.IGH %>%
  mutate(isotype = substring(allCHitsWithScore, 1, 4)) %>%
  mutate(isotype = ifelse(isotype == "IGHM", "IGHM", "Not-IGHM")) %>%
  filter(!is.na(STAGE)) %>%
  group_by(STAGE, isotype) %>%
  dplyr::summarize(n = n())

df.wide <- df %>%
  tidyr::spread(STAGE, n) %>%
  replace(is.na(.), 0)

fisher.test(df.wide[,2:3])

g <- ggbarplot(df, x = "isotype", y = "n", label = df$n , lab.size = 2,
               fill = "isotype", facet.by = c("STAGE"), palette = c("jco")) +
  theme_bw()+
  theme(text = element_text(size=10)) + ylab("Number of samples with given isotype") + xlab("Stage")
ggpar(g, legend = "bottom", legend.title = "isotype") + geom_hline(yintercept = 0)

#--
# IgM vs others, by Ann Arbor stage
#--

df <- clonotypes.IGH %>%
  group_by(sample) %>%
  top_n(n = 1, wt = cloneCount) %>%
  mutate(isotype = substring(allCHitsWithScore, 1, 4)) %>%
  mutate(isotype = ifelse(isotype == "IGHM", "IGHM", "Not-IGHM")) %>%
  filter(!is.na(ANN_ARBOR_STAGE)) %>%
  group_by(ANN_ARBOR_STAGE, isotype) %>%
  dplyr::summarize(n = n())

df.wide <- df %>%
  tidyr::spread(ANN_ARBOR_STAGE, n) %>%
  replace(is.na(.), 0)

chisq.test(df.wide[,2:5])

g <- ggbarplot(df, x = "isotype", y = "n", label = df$n , lab.size = 2,
               fill = "isotype", facet.by = c("ANN_ARBOR_STAGE"), palette = c("jco")) +
  theme_bw()+
  theme(text = element_text(size=10)) + ylab("Number of samples with given isotype") + xlab("ANN_ARBOR_STAGE")
ggpar(g, legend = "bottom", legend.title = "isotype") + geom_hline(yintercept = 0)
  
#--
# IgM vs others, by SNF cluster
#--
  
df <- clonotypes.IGH %>%
  group_by(sample) %>%
  top_n(n = 1, wt = cloneCount) %>%
  mutate(isotype = substring(allCHitsWithScore, 1, 4)) %>%
  mutate(isotype = ifelse(isotype == "IGHM", "IGHM", "Not-IGHM")) %>%
  filter(!is.na(SNFClust)) %>%
  group_by(SNFClust, isotype) %>%
  dplyr::summarize(n = n())

df.wide <- df %>%
  tidyr::spread(SNFClust, n) %>%
  replace(is.na(.), 0)

fisher.test(df.wide[,2:3])

g <- ggbarplot(df, x = "isotype", y = "n", label = df$n , lab.size = 2,
               fill = "isotype", facet.by = c("SNFClust"), palette = c("jco")) +
  theme_bw()+
  theme(text = element_text(size=10)) + ylab("Number of samples with given isotype") + xlab("SNFClust")
ggpar(g, legend = "bottom", legend.title = "isotype") + geom_hline(yintercept = 0)

#--
# IgM vs others, by BCL2 translocation
#--

df <- clonotypes.IGH %>%
  group_by(sample) %>%
  top_n(n = 1, wt = cloneCount) %>%
  mutate(isotype = substring(allCHitsWithScore, 1, 4)) %>%
  mutate(isotype = ifelse(isotype == "IGHM", "IGHM", "Not-IGHM")) %>%
  filter(!is.na(BCL2_BA_consensus)) %>%
  group_by(BCL2_BA_consensus, isotype) %>%
  dplyr::summarize(n = n())

df.wide <- df %>%
  tidyr::spread(BCL2_BA_consensus, n) %>%
  replace(is.na(.), 0)

fisher.test(df.wide[,2:3])

g <- ggbarplot(df, x = "isotype", y = "n", label = df$n , lab.size = 2,
               fill = "isotype", facet.by = c("BCL2_BA_consensus"), palette = c("jco")) +
  theme_bw()+
  theme(text = element_text(size=10)) + ylab("Number of samples with given isotype") + xlab("BCL2_BA_consensus")
ggpar(g, legend = "bottom", legend.title = "isotype") + geom_hline(yintercept = 0)

#--
# IgM vs others, by BCL6 translocation
#--

df <- clonotypes.IGH %>%
  group_by(sample) %>%
  top_n(n = 1, wt = cloneCount) %>%
  mutate(isotype = substring(allCHitsWithScore, 1, 4)) %>%
  mutate(isotype = ifelse(isotype == "IGHM", "IGHM", "Not-IGHM")) %>%
  filter(!is.na(BCL6_BA_consensus)) %>%
  group_by(BCL6_BA_consensus, isotype) %>%
  dplyr::summarize(n = n())

df.wide <- df %>%
  tidyr::spread(BCL6_BA_consensus, n) %>%
  replace(is.na(.), 0)

fisher.test(df.wide[,2:3])

g <- ggbarplot(df, x = "isotype", y = "n", label = df$n , lab.size = 2,
               fill = "isotype", facet.by = c("BCL6_BA_consensus"), palette = c("jco")) +
  theme_bw()+
  theme(text = element_text(size=10)) + ylab("Number of samples with given isotype") + xlab("BCL6_BA_consensus")
ggpar(g, legend = "bottom", legend.title = "isotype") + geom_hline(yintercept = 0)

#--
# Mosaic plot
#--

df <- clonotypes.IGH %>%
  group_by(sample) %>%
  top_n(n = 1, wt = cloneCount) %>%
  mutate(isotype = substring(allCHitsWithScore, 1, 4)) %>%
  mutate(isotype = ifelse(isotype == "IGHM", "IGHM", "Not-IGHM")) %>%
  filter(TYPE == "FL") %>%
  select(BCL2_BA_consensus, isotype, STAGE) %>%
  mutate(BCL2_BA_consensus = as.factor(BCL2_BA_consensus))

ggplot(data = df) +
  geom_mosaic(aes(x = product(BCL2_BA_consensus, isotype),
                  fill = BCL2_BA_consensus), na.rm = TRUE) + 
  facet_grid(STAGE~.) +
  theme_bw()

ggplot(data = df) +
  geom_mosaic(aes(x = product(isotype, BCL2_BA_consensus),
                  fill = isotype), na.rm = TRUE) + 
  facet_grid(STAGE~.) +
  theme_bw()





df <- clonotypes.IGH %>%
  group_by(sample) %>%
  top_n(n = 1, wt = cloneCount) %>%
  mutate(isotype = substring(allCHitsWithScore, 1, 4)) %>%
  mutate(isotype = ifelse(isotype == "IGHM", "IGHM", "Not-IGHM")) %>%
  filter(TYPE == "FL") %>%
  select(BCL2_BA_consensus, isotype, SNFClust, STAGE,TYPE, ANN_ARBOR_STAGE, BCL6_BA_consensus) %>%
  mutate(BCL2_BA_consensus = as.factor(BCL2_BA_consensus)) %>%
  mutate(BCL6_BA_consensus = as.factor(BCL6_BA_consensus))

ggplot(data = df) +
  geom_mosaic(aes(x = product(BCL2_BA_consensus, BCL6_BA_consensus),
                  fill = BCL2_BA_consensus), na.rm = TRUE) + 
  facet_grid(ANN_ARBOR_STAGE~.) +
  theme_bw()

