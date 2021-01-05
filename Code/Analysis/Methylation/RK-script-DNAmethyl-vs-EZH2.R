
#---
# This script analyzes impact of EZH2 mutations on DNA methylation
# This script needs to be run from /FLOMICS/ folder
#---
  
packages <- c("dplyr", "ggplot2", "limma", "reshape2", "ggpubr", "DMRcate", "GenomicRanges")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

# setwd("~/github/FLOMICS/")

# Read in AnnotationFile
AnnotationFile <- read.csv(file = "methylation/Ann_IlluminaHumanMethylationEPICanno.ilm10b2.hg19.csv")
AnnotationFile$X <- NULL

# Read in methylation data
M.values <- data.table::fread("methylation/2_MvalueMatrix_updSamples_Ordered_T1_FilteredProbes.csv")
M.values <- data.frame(M.values)
row.names(M.values) <- M.values$V1
M.values$V1 <- NULL

beta.values <- data.table::fread("methylation/2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.csv")
beta.values <- data.frame(beta.values)
row.names(beta.values) <- beta.values$V1
beta.values$V1 <- NULL

# Read in mutation data
mutations <- read.csv("DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.csv") %>%
  filter(Hugo_Symbol == "EZH2")
row.names(mutations) <- mutations$Hugo_Symbol
mutations <- mutations[,2:ncol(mutations)]
mutations <- t(mutations)
mutations <- data.frame(mutations)
mutations$SAMPLE_ID <- row.names(mutations)
mutations <- mutations[!is.na(mutations$EZH2),]
samples.mut.and.methylation <- intersect(row.names(mutations), colnames(M.values))
mutations <- mutations %>% filter(SAMPLE_ID %in% samples.mut.and.methylation) %>%
  mutate(EZH2 = ifelse(EZH2 == "0", "wt",
                       ifelse(EZH2 == "1", "mut", NA)))

EZH2.wt.cases <- mutations %>% data.frame() %>% filter(EZH2 == "wt") %>% .$SAMPLE_ID
EZH2.mut.cases <- mutations %>% data.frame() %>% filter(EZH2 == "mut") %>% .$SAMPLE_ID

# Subset M and beta values data frame based on common samples
M.values <- M.values[,samples.mut.and.methylation]
beta.values <- beta.values[,samples.mut.and.methylation]

# EZH2 mutation vs. purity
read.csv("methylation/Purity_281probes_10Jan2020.csv") %>%
  mutate(EZH2 = ifelse(X %in% EZH2.mut.cases, "MUT",
                       ifelse(X %in% EZH2.wt.cases, "WT", NA))) %>%
  filter(!is.na(EZH2)) %>%
  ggboxplot(x = "EZH2", y = "purity", width = 0.8) +
  stat_compare_means()
# no difference in purity

# Differential methylation by EZH2 mutation status
design <- model.matrix(~0 + EZH2, data = mutations[,c("SAMPLE_ID", "EZH2")])
colnames(design) <- c("mut", "wt")
fit <- lmFit(M.values, design)
contMatrix <- makeContrasts(mut-wt, levels = design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2, adjust.method = "BH", p.value = 0.05, lfc = 2)) # DM CpGs at FDR < 0.05

diff_meth_mut_wt <- topTable(fit2, coef = "mut - wt", adjust.method = "BH", number = Inf)
diff_meth_mut_wt$Probe <- row.names(diff_meth_mut_wt)
# minimal if no CpG methylation changes in EZH2 mut vs wt cases

# Differentially methylated regions using DMRcate
myAnnotation <- cpg.annotate(object = as.matrix(M.values), datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "mut - wt", arraytype = "EPIC")
DMRs <- dmrcate(myAnnotation, lambda = 1000, C = 2)
DMRs_ranges <- extractRanges(DMRs, genome = "hg19")
DMRs_results <- data.frame(DMRs_ranges)
# similar findings, only 1 hit (DNAJA4) with weak evidence

# Determine mean methylation in EZH2 mut vs wt
melted.beta.values <- beta.values %>%
  mutate(probe = row.names(.)) %>%
  melt()

melted.beta.values %>%
  mutate(EZH2 = ifelse(variable %in% EZH2.mut.cases, "MUT",
                ifelse(variable %in% EZH2.wt.cases, "WT", NA))) %>%
  group_by(EZH2) %>%
  summarize(mean = mean(value))

melted.beta.values %>%
  left_join(AnnotationFile[,c("V1", "Regulatory_Feature_Group")], by = c("probe" = "V1")) %>%
  mutate(EZH2 = ifelse(variable %in% EZH2.mut.cases, "MUT",
                       ifelse(variable %in% EZH2.wt.cases, "WT", NA))) %>%
  group_by(probe, EZH2, Regulatory_Feature_Group) %>%
  summarize(mean = mean(value)) %>%
  ggboxplot(x = "EZH2", y = "mean", width = 0.8) +
  facet_grid(. ~ Regulatory_Feature_Group) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means()

melted.beta.values %>%
  left_join(AnnotationFile[,c("V1", "Relation_to_Island")], by = c("probe" = "V1")) %>%
  mutate(EZH2 = ifelse(variable %in% EZH2.mut.cases, "MUT",
                       ifelse(variable %in% EZH2.wt.cases, "WT", NA))) %>%
  filter(!is.na(EZH2)) %>%
  group_by(EZH2, Relation_to_Island) %>%
  summarize(mean = mean(value))

melted.beta.values %>%
  left_join(AnnotationFile[,c("V1", "Relation_to_Island")], by = c("probe" = "V1")) %>%
  mutate(EZH2 = ifelse(variable %in% EZH2.mut.cases, "MUT",
                       ifelse(variable %in% EZH2.wt.cases, "WT", NA))) %>%
  filter(!is.na(EZH2)) %>%
  group_by(probe, EZH2, Relation_to_Island) %>%
  summarize(mean = mean(value)) %>%
  ggviolin(x = "EZH2", y = "mean", width = 0.8, add = "boxplot", color = "EZH2") +
  facet_grid(. ~ Relation_to_Island) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means()

# TAD analysis
TADs <- read.table("methylation/41588_2018_338_MOESM4_ESM_tab1.txt", sep = "\t", header = T)
table(TADs$Class)

TADs %>% # hg19?
  group_by(Class) %>%
  summarize(mean = mean(Mean.H3K27me3.logFC))

TADs_GRanges <- TADs %>%
  select(chr = Chr, start = Start, end = End, Class)
TADs_GRanges <- makeGRangesFromDataFrame(TADs_GRanges, keep.extra.columns = TRUE)

AnnotationFile_GRanges <- AnnotationFile %>%
  select(chr = chr, start = pos, end = pos, strand = strand, V1)
AnnotationFile_GRanges <- makeGRangesFromDataFrame(AnnotationFile_GRanges, keep.extra.columns = TRUE)

m <- findOverlaps(AnnotationFile_GRanges, TADs_GRanges)
AnnotationFile_GRanges.matched <- AnnotationFile_GRanges[queryHits(m)]
mcols(AnnotationFile_GRanges.matched) <- cbind.data.frame(
  mcols(AnnotationFile_GRanges.matched),
  mcols(TADs_GRanges[subjectHits(m)]))
AnnotationFile_GRanges.matched <- data.frame(AnnotationFile_GRanges.matched)

melted.beta.values.TAD <- melted.beta.values %>%
  mutate(EZH2 = ifelse(variable %in% EZH2.mut.cases, "MUT",
                       ifelse(variable %in% EZH2.wt.cases, "WT", NA))) %>%
  left_join(AnnotationFile_GRanges.matched[,c("V1", "Class")], by = c("probe" = "V1")) %>%
  mutate(Class = as.character(Class)) %>%
  mutate(Class = ifelse(is.na(Class), "not.annotated", Class))

melted.beta.values.TAD %>%
  group_by(Class, EZH2) %>%
  summarize(mean = mean(value))

# Overall conclusion: despite that diff methylation in FL is characterized by hypermethylation, specifically in areas targeted by EZH2,
# and that EZH2 has been shown in other contexts to affect DNA methylation,
# we do not see meaningful association between EZH2 mutations and DNA methylation.
# In keeping with DNA methylation resulting from EZH2 not being direct effect of catalytic site,
# but from recruitment of DNA methyltransferases.
# Alternatively, EZH2 may not be implicated in DNA methylation in FL context.

# library(regioneR)

suffix <- seq(1:22)
chr.to.keep <- paste0("chr", suffix)

H3K27me3 <- read.table("~/FLOMICS/Git/Methylation/Pipeline/merged_peaks_h3k27me3_karpas_wsudlcl2_ocily.bed", sep = "\t", header = FALSE) %>%
  select(chr = V1, start = V2, end = V3) #%>%
  #filter(chr %in% chr.to.keep)
H3K27me3_GRanges <- makeGRangesFromDataFrame(H3K27me3)

EZH2.peaks <- read.table("~/EZH2_peaks.bed", sep = "\t", header = FALSE) %>%
  select(chr = V1, start = V2, end = V3) #%>%
  #filter(chr %in% chr.to.keep)
EZH2.peaks_GRanges <- makeGRangesFromDataFrame(EZH2.peaks)

DNMT1 <- read.table("~/DNMT1_peaks.bed", sep = "\t", header = FALSE) %>%
  select(chr = V1, start = V2, end = V3) %>%
  filter(chr %in% chr.to.keep)
DNMT1_GRanges <- makeGRangesFromDataFrame(DNMT1)

DNMT3B <- read.table("~/DNMT3A_peaks.bed", sep = "\t", header = FALSE) %>%
  select(chr = V1, start = V2, end = V3) %>%
  filter(chr %in% chr.to.keep)
DNMT3B_GRanges <- makeGRangesFromDataFrame(DNMT3B)

peaks_45870 <- read.table("~/45870_peaks.bed", sep = "\t", header = FALSE) %>%
  select(chr = V1, start = V2, end = V3) %>%
  filter(chr %in% chr.to.keep)
peaks_45870_GRanges <- makeGRangesFromDataFrame(peaks_45870)

library(regioneR)

A.input <- H3K27me3_GRanges 
B.input <- EZH2.peaks_GRanges
# hg 19???

#pt <-
overlapPermTest(A = A.input, B = B.input, ntimes = 100)
#pt <-
permTest(A = A.input, B = B.input, ntimes = 100,
               randomize.function =  randomizeRegions,
               non.overlapping = FALSE,  per.chromosome = TRUE,
               evaluate.function = numOverlaps, count.once = TRUE,
               genome = "hg19", mc.set.seed = FALSE, verbose = TRUE)

permTest(A=HepG2_Rad21_5K, B=HepG2_Ctcf, ntimes=1000,
         randomize.function=circularRandomizeRegions,
         evaluate.function=numOverlaps, count.once=TRUE,
         genome="hg19", mc.set.seed=FALSE, mc.cores=4)

pt
plot(pt)
overlapGraphicalSummary(A = A.input, B = B.input)
lz <- localZScore(pt = pt, A = A.input, B = B.input)
plot(lz)



pt <- overlapPermTest(A = H3K27me3_GRanges, B = EZH2.peaks_GRanges, ntimes=100)
pt <- overlapPermTest(A = H3K27me3_GRanges, B = DNMT1_GRanges, ntimes=100)
pt <- overlapPermTest(A = H3K27me3_GRanges, B = DNMT3B_GRanges, ntimes=100)
pt <- overlapPermTest(A = EZH2.peaks_GRanges, B = DNMT1_GRanges, ntimes=100)


pt <- overlapPermTest(A = H3K27me3_GRanges, B = DNMT1_GRanges, ntimes=100)
overlapGraphicalSummary(A = H3K27me3_GRanges, B = DNMT1_GRanges)
pt

pt <- overlapPermTest(A = EZH2.peaks_GRanges, B = H3K27me3_GRanges, ntimes=100)
pt

pt <- overlapPermTest(A = EZH2.peaks_GRanges, B = peaks_45870_GRanges, ntimes=100)
pt

genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
A <- createRandomRegions(nregions=20, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE) 
B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))
pt <- overlapPermTest(A=A, B=B, ntimes=100, genome=genome, non.overlapping=FALSE)
plot(pt)


pt
plot(pt)
lz <- localZScore(pt = pt, A = H3K27me3_GRanges, B = DNMT3B_GRanges)
lz <- localZScore(pt = pt, A = EZH2.peaks_GRanges, B = DNMT1_GRanges)
pt <- overlapPermTest(A = EZH2.peaks_GRanges, B = peaks_45870_GRanges, ntimes=100)
plot(lz)


x <- read.table("~/41590_2018_181_MOESM3_ESM.txt", sep = "\t", header = TRUE)

x %>%
  filter(QCfilters == "pass") %>%
  group_by(SimpleSortPheno) %>%
  summarize(count = n())
  
x %>%
  filter(QCfilters == "pass") %>%  
  # filter(MKI67 < 5) %>%
  ggscatter(x = "EZH2", y = "DNMT1",
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE) + # Add confidence interval
  facet_grid(. ~ SimpleSortPheno) +
  ylim(0,30) +
  stat_cor(method = "pearson", label.x = 3, label.y = 30)

ggsave("tmp.pdf", width = 40, height = 80, units = "cm")

# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 3, label.y = 30)
#> `geom_smooth()` using formula 'y ~ x'
  
  
  ggplot(aes(x = EZH2, y = DNMT1)) +
  facet_grid(. ~ SimpleSortPheno) +
  geom_point() +
  geom_smooth(method = lm, se = TRUE) +
  theme_bw()

ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point()


plot(x$EZH2, x$DNMT1)







#--
# Analysis by EZH2 mutation status
#---

mutations <- read.csv("DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.csv") %>%
  filter(Hugo_Symbol == "EZH2")
row.names(mutations) <- mutations$Hugo_Symbol
mutations <- mutations[,2:ncol(mutations)]
mutations <- t(mutations)
mutations <- data.frame(mutations)
mutations$SAMPLE_ID <- row.names(mutations)
mutations <- mutations[!is.na(mutations$EZH2),]
samples.mut.and.methylation <- intersect(row.names(mutations), colnames(Combined_M_valueMatrix))
mutations <- mutations %>% filter(SAMPLE_ID %in% samples.mut.and.methylation) %>%
  mutate(EZH2 = ifelse(EZH2 == "0", "wt",
                       ifelse(EZH2 == "1", "mut", NA)))
Combined_M_valueMatrix_EZH2 <- Combined_M_valueMatrix[,mutations$SAMPLE_ID]

design <- model.matrix(~0 + EZH2, data = mutations[,c("SAMPLE_ID", "EZH2")])
colnames(design) <- c("mut", "wt")
fit <- lmFit(Combined_M_valueMatrix_EZH2, design)
contMatrix <- makeContrasts(mut-wt, levels = design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2, adjust.method = "BH", p.value = 0.05, lfc = 2)) # DM CpGs at FDR < 0.05

diff_meth_mut_wt <- topTable(fit2, coef = "mut - wt", adjust.method = "BH", number = Inf)
diff_meth_mut_wt$Probe <- row.names(diff_meth_mut_wt)

ggplot(data = diff_meth_mut_wt, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 2) +
  theme(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  scale_color_manual(values = c("up" = "#d6604d", "down" = "#4393c3", not.sig = "grey")) +
  geom_vline(xintercept = c(-2,2), linetype = "dashed") +
  geom_hline(yintercept = 1.30103, linetype = "dashed") +
  guides(color = FALSE) +
  theme_bw()
# minimal if no CpG methylation changes in EZH2 mut vs wt cases

EZH2.wt.cases <- mutations %>%
  data.frame() %>%
  filter(EZH2 == "wt") %>%
  .$SAMPLE_ID

EZH2.mut.cases <- mutations %>%
  data.frame() %>%
  filter(EZH2 == "mut") %>%
  .$SAMPLE_ID

# EZH2 mutation vs. purity
data.table::fread("2020-02-16 InfininiumPurify.txt") %>%
  mutate(EZH2 = ifelse(V1 %in% EZH2.mut.cases, "MUT",
                       ifelse(V1 %in% EZH2.wt.cases, "WT", NA))) %>%
  filter(!is.na(EZH2)) %>%
  ggboxplot(x = "EZH2", y = "x", width = 0.8) +
  stat_compare_means()
# no difference in purity
ggsave(paste0("img/",  date, " RK-purity-by-EZH2.pdf"), width = 7, height = 7, units = "cm")

# Mean methylation by EZH2
Combined_Beta_valueMatrix <- getBeta(Combined_mSet_Sw)

df1 <- Combined_Beta_valueMatrix %>%
  data.frame() %>%
  mutate(probe = row.names(.)) %>%
  melt() %>%
  mutate(variable = gsub("\\.", "-", variable)) %>%
  left_join(pData, by = c("variable" = "SAMPLE_ID"))

df1 %>%
  filter(variable %in% c(EZH2.mut.cases, EZH2.wt.cases)) %>%
  mutate(EZH2 = ifelse(variable %in% EZH2.mut.cases, "MUT",
                       ifelse(variable %in% EZH2.wt.cases, "WT", NA))) %>%
  group_by(EZH2) %>%
  summarize(mean = mean(value))

df1 %>%
  filter(variable %in% c(EZH2.mut.cases, EZH2.wt.cases)) %>%
  mutate(EZH2 = ifelse(variable %in% EZH2.mut.cases, "MUT",
                       ifelse(variable %in% EZH2.wt.cases, "WT", NA))) %>%
  ggplot(aes(x = value, color = EZH2)) +
  geom_density()

df1 %>%
  filter(variable %in% c(EZH2.mut.cases, EZH2.wt.cases)) %>%
  left_join(AnnotationFile[,c("V1", "Relation_to_Island")], by = c("probe" = "V1")) %>%
  mutate(EZH2 = ifelse(variable %in% EZH2.mut.cases, "MUT",
                       ifelse(variable %in% EZH2.wt.cases, "WT", NA))) %>%
  group_by(probe, EZH2, Relation_to_Island) %>%
  summarize(mean = mean(value)) %>%
  ggviolin(x = "EZH2", y = "mean", width = 0.8, add = "boxplot", color = "EZH2") +
  facet_grid(. ~ Relation_to_Island) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means()

df1 %>%
  filter(variable %in% c(EZH2.mut.cases, EZH2.wt.cases)) %>%
  left_join(AnnotationFile[,c("V1", "Regulatory_Feature_Group")], by = c("probe" = "V1")) %>%
  mutate(EZH2 = ifelse(variable %in% EZH2.mut.cases, "MUT",
                       ifelse(variable %in% EZH2.wt.cases, "WT", NA))) %>%
  group_by(probe, EZH2, Regulatory_Feature_Group) %>%
  summarize(mean = mean(value)) %>%
  ggboxplot(x = "EZH2", y = "mean", width = 0.8) +
  facet_grid(. ~ Regulatory_Feature_Group) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means()

#---
# TAD analysis
#---

TADs <- read.table("methylation/41588_2018_338_MOESM4_ESM_tab1.txt", sep = "\t", header = T)
table(TADs$Class)

TADs %>% # hg19?
  group_by(Class) %>%
  summarize(mean = mean(Mean.H3K27me3.logFC))

TADs_GRanges <- TADs %>%
  select(chr = Chr, start = Start, end = End, Class)
TADs_GRanges <- makeGRangesFromDataFrame(TADs_GRanges, keep.extra.columns = TRUE)

AnnotationFile_GRanges <- AnnotationFile %>%
  select(chr = chr, start = pos, end = pos, strand = strand, V1) %>%
  filter(V1 %in% row.names(Combined_M_valueMatrix))
AnnotationFile_GRanges <- makeGRangesFromDataFrame(AnnotationFile_GRanges, keep.extra.columns = TRUE)

m <- findOverlaps(AnnotationFile_GRanges, TADs_GRanges)
AnnotationFile_GRanges.matched <- AnnotationFile_GRanges[queryHits(m)]
mcols(AnnotationFile_GRanges.matched) <- cbind.data.frame(
  mcols(AnnotationFile_GRanges.matched),
  mcols(TADs_GRanges[subjectHits(m)]))
AnnotationFile_GRanges.matched <- data.frame(AnnotationFile_GRanges.matched)



df2 <- df1 %>%
  group_by(probe, TYPE) %>%
  summarize(mean = mean(value)) %>%
  left_join(AnnotationFile_GRanges.matched[,c("V1", "Class")], by = c("probe" = "V1")) %>%
  mutate(Class = as.character(Class)) %>%
  mutate(Class = ifelse(is.na(Class), "not.annotated", Class))

# df2 %>%
#   filter(probe %in% hypermeth_FL_intersect) %>%
#   ggboxplot(x = "Class", y = "mean", width = 0.8) +
#   facet_grid(. ~ TYPE) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   stat_compare_means()
# 
# ggsave(paste0("img/",  date, " RK-methylation-TADs.pdf"), width = 40, height = 10, units = "cm")

# need to compare methylation within TADs between EZH2 mut and wt cases

mutations <- read.csv("~/FLOMICS/DNAseq/Mutation_and_BA_results/mut.merged.df.T1.csv") %>%
  filter(Hugo_Symbol == "EZH2")
row.names(mutations) <- mutations$Hugo_Symbol
mutations <- mutations[,2:ncol(mutations)]
mutations <- t(mutations)





uu <- AnnotationFile %>%
  # filter(Regulatory_Feature_Group == "Promoter_Associated") %>%
  # filter(Relation_to_Island == "Island") %>%
  filter(Relation_to_Island == "OpenSea") %>%
  # filter(UCSC_RefGene_Group == "Island") %>%
  .$V1




df1 %>%
  filter(variable %in% c(EZH2.mut.cases, EZH2.wt.cases)) %>%
  mutate(EZH2 = ifelse(variable %in% EZH2.mut.cases, "MUT",
                       ifelse(variable %in% EZH2.wt.cases, "WT", NA))) %>%
  left_join(AnnotationFile_GRanges.matched[,c("V1", "Class")], by = c("probe" = "V1")) %>%
  filter(probe %in% hypermeth_FL_intersect) %>%
  mutate(Class = as.character(Class)) %>%
  mutate(Class = ifelse(is.na(Class), "not.annotated", Class)) %>%
  group_by(Class, EZH2) %>%
  summarize(mean = mean(value))

df3 <- df1 %>%
  filter(variable %in% c(EZH2.mut.cases, EZH2.wt.cases)) %>%
  mutate(EZH2 = ifelse(variable %in% EZH2.mut.cases, "MUT",
                       ifelse(variable %in% EZH2.wt.cases, "WT", NA))) %>%
  group_by(probe, EZH2) %>%
  summarize(mean = mean(value)) %>%
  left_join(AnnotationFile_GRanges.matched[,c("V1", "Class")], by = c("probe" = "V1"))



uu <- AnnotationFile %>%
  # filter(Regulatory_Feature_Group == "Promoter_Associated") %>%
  # filter(Relation_to_Island == "Island") %>%
  filter(Relation_to_Island == "OpenSea") %>%
  # filter(UCSC_RefGene_Group == "Island") %>%
  .$V1

df4 <- df3 %>%
  filter(probe %in% hypermeth_FL_intersect) %>%
  mutate(Class = as.character(Class)) %>%
  mutate(Class = ifelse(is.na(Class), "not.annotated", Class)) %>%
  filter(!Class == "not.annotated")
dim(df4)

df5 <- df4 %>%
  filter(probe %in% uu)

dim(df5)
table(df5$Class)

df5 %>% 
  ggboxplot(x = "EZH2", y = "mean", width = 0.8) +
  # facet_grid(. ~ EZH2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means()

ggsave(paste0("img/",  date, " RK-mean-meth-byEZH2-TAD.pdf"), width = 14, height = 14, units = "cm")



  
