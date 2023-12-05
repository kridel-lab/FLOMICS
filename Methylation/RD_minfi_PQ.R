# R Packages
# ==========
library(minfi)
library(minfiData)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(limma)
library(RColorBrewer)
library(maxprobes)
library(sva)
library(lumi)
library(dplyr)
library(stringr)

# Read IDAT Files
# ===============
# set the working directory
baseDir <- "/cluster/home/t110989uhn/kridelgroup/rajesh/01_DNA_Methylome_Analysis/00_minfi/02_FL_Samples/"
list.files(baseDir)

# Read a sample sheet
targets <- read.metharray.sheet(baseDir)

# Remove 2 machine control samples
targets <- targets %>% filter(Sample_Name != 'HCT116_DKO_methylated')

# Remove 2 samples for clinical reason
targets <- targets %>% filter(Sample_Name != 'LY_FL_311_T1' & Sample_Name != 'LY_FL_159_T1')

# Rename replicate of sample "LY_FL_159_T1"
targets$Sample_Name <- str_replace(targets$Sample_Name, "LY_FL_159_T1_rep", "LY_FL_159_T1")

# Remove 1 TFL sample
targets <- targets %>% filter(Sample_Name != 'LY_FL_127_T2')

# Remove 2 suspicious repeated samples
targets <- targets %>% filter(Sample_Name != 'LY_FL_1156_T1' & Sample_Name != 'LY_FL_571_T1')

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,539 features
RGset <- read.metharray.exp(targets = targets, force=TRUE)
dim(RGset)

# Give the samples descriptive names
targets$ID <- gsub(" ", "", paste(targets$Sample_Name,".",targets$Array))
sampleNames(RGset) <- targets$ID
RGset
dim(RGset)

# See which annotation package being used by
annotation(RGset)

# Read in the sample annotation
sample_ann <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/01_DNA_Methylome_Analysis/00_minfi/db/20221228_sample_annotations.txt",
                         sep = "\t",
                         header = TRUE,
                         na.strings=c("","NA"))

# Left outer join sample sheet and sample annotation
targets <- merge(x = targets,
              y = sample_ann,
              by.x = "Sample_Name",
              by.y = "SAMPLE_ID",
              all.x = TRUE)


## Set Working Current Directory
## =============================
setwd("/cluster/home/t110989uhn/kridelgroup/rajesh/01_DNA_Methylome_Analysis/00_minfi/02_FL_Samples/")

# Filter Samples By p-value
# =========================
# Calculate the detection p-values
detP <- detectionP(RGset)
head(detP)

# Examine mean detection p-values across all samples to identify any failed samples
pval <- (colMeans(detP)) %>% as.data.frame.list()
pval <- as.data.frame(t(pval))
pval <- cbind(rownames(pval), data.frame(pval, row.names=NULL))
colnames(pval) <- c("Sample_Name", "mean_pval")
write.table(pval, file = "pval.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# Examine mean detection p-values across all samples to identify any failed samples
pdf("pval.pdf", height = 40)
par(mar=c(6,6,2,2)) # Increase margin size
barplot(pval$mean_pval,
        names=pval$Sample_Name,
        col='#454545',
        xlab='Mean detection p-values',
        horiz = TRUE,
        las=1,
        cex.names=0.4,
        xlim = c(0, 0.6),
        border = NA
)
abline(v = 0.01, col = "red")
dev.off()

# remove poor quality samples
keep <- colMeans(detP) < 0.01
RGset <- RGset[,keep]
dim(RGset)

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

# Re-Examine mean detection p-values across all samples to identify any failed samples
pval <- (colMeans(detP)) %>% as.data.frame.list()
pval <- as.data.frame(t(pval))
pval <- cbind(rownames(pval), data.frame(pval, row.names=NULL))
colnames(pval) <- c("Sample_Name", "mean_pval")
write.table(pval, file = "pval_updated.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# Re-Examine mean detection p-values across all samples to identify any failed samples
pdf("pval_updated.pdf", height = 40)
par(mar=c(6,6,2,2)) # Increase margin size
barplot(pval$mean_pval,
        names=pval$Sample_Name,
        col='#454545',
        xlab='Mean detection p-values',
        horiz = TRUE,
        las=1,
        cex.names=0.4,
        xlim = c(0, 0.6),
        border = NA
)
abline(v = 0.01, col = "red")
dev.off()

# Filter known repeated samples by a # of significant probes
# ==========================================================
# Examine significant number of probes in each samples
sig_probes <- (colSums(detP <= 0.01)) %>% as.data.frame.list()
sig_probes <- as.data.frame(t(sig_probes))
sig_probes <- cbind(rownames(sig_probes), data.frame(sig_probes, row.names=NULL))
colnames(sig_probes) <- c("Sample_Name", "significant_probes")
write.table(sig_probes, file = "significant_probes.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# Examine significant number of probes in each samples
pdf("significant_probes.pdf", height = 40)
par(mar=c(6,6,2,2)) # Increase margin size
barplot(sig_probes$significant_probes,
        names=sig_probes$Sample_Name,
        col='#454545',
        xlab = "Significant number of probes",
        horiz = TRUE,
        las=1,
        cex.names=0.4,
        xlim = c(0,865859),
        border = NA
)
abline(v = 780000, col = "red")
dev.off()

# Identify known repeated sample with minimum # of significant probes
sig_probes_t <- data.frame(sig_probes[,-1], row.names=sig_probes[,1])
colnames(sig_probes_t) <- c("significant_probes")

# Compare and find known repeated sample names with minimum # of significant probes
sample_1 <- c("LY_FL_529_T1.R01C01",
              "LY_FL_527_T1.R01C01",
              "LY_FL_525_T1.R03C01",
              "LY_FL_524_T1.R02C01",
              "LY_FL_523_T1.R01C01",
              "LY_FL_488_T1.R06C01",
              "LY_FL_479_T1.R02C01",
              "LY_FL_158_T1.R01C01",
              "LY_FL_498_T1.R06C01")

sample_2 <- c("LY_FL_529_T1.R05C01",
              "LY_FL_527_T1.R04C01",
              "LY_FL_525_T1.R02C01",
              "LY_FL_524_T1.R03C01",
              "LY_FL_523_T1.R04C01",
              "LY_FL_488_T1.R04C01",
              "LY_FL_479_T1.R07C01",
              "LY_FL_158_T1.R08C01",
              "LY_FL_498_T1.R05C01")

rep_samples <- data.frame(sample_1, sample_2)

sample_min <- c()

for (i in 1:9) {
  df <- subset(sig_probes_t, rownames(sig_probes_t) %in% c(rep_samples$sample_1[i], rep_samples$sample_2[i]))
  sample_min <- append(rownames(df)[which(df == min(df), arr.ind = TRUE)[ , 1]], sample_min)
}

sample_min
# [1] "LY_FL_498_T1.R06C01" "LY_FL_158_T1.R01C01" "LY_FL_479_T1.R02C01" "LY_FL_488_T1.R06C01" "LY_FL_523_T1.R01C01"
# [6] "LY_FL_524_T1.R02C01" "LY_FL_525_T1.R03C01" "LY_FL_527_T1.R04C01" "LY_FL_529_T1.R01C01"

# remove repeated 9 samples from mean detection p-values
pval <- pval[!(pval$Sample_Name %in% sample_min),]

# filter sample sheet with final remaining samples
targets <- targets[(targets$ID %in% pval$Sample_Name),]

# Re-generate Final RGset, detP, pval, sig_probes with final unique samples
# =========================================================================
# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,539 features
RGset <- read.metharray.exp(targets = targets, force=TRUE)
RGset

# Give the samples original and unique descriptive names
sampleNames(RGset) <- targets$Sample_Name
RGset

# Calculate the detection p-values
detP <- detectionP(RGset)
head(detP)

# Examine mean detection p-values across all samples to identify any failed samples
pval <- (colMeans(detP)) %>% as.data.frame.list()
pval <- as.data.frame(t(pval))
pval <- cbind(rownames(pval), data.frame(pval, row.names=NULL))
colnames(pval) <- c("Sample_Name", "mean_pval")
write.table(pval, file = "pval_final.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# Examine mean detection p-values across all samples to identify any failed samples
pdf("pval_final.pdf", height = 40)
par(mar=c(6,6,2,2)) # Increase margin size
# pal <- brewer.pal(8,"Dark2")
pal <- c("#fc1417", "#6dcc04", "#8b25fa")
barplot(pval$mean_pval,
        names=pval$Sample_Name,
        col=pal[factor(targets$TYPE)],
        xlab='Mean detection p-values',
        horiz = TRUE,
        las=1,
        cex.names=0.4,
        xlim = c(0, 0.6),
        border = NA
)
abline(v = 0.01, col = "red")
dev.off()

# Examine significant number of probes in each samples
sig_probes <- (colSums(detP <= 0.01)) %>% as.data.frame.list()
sig_probes <- as.data.frame(t(sig_probes))
sig_probes <- cbind(rownames(sig_probes), data.frame(sig_probes, row.names=NULL))
colnames(sig_probes) <- c("Sample_Name", "significant_probes")
write.table(sig_probes, file = "significant_probes_final.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# Examine significant number of probes in each samples
pdf("significant_probes_final.pdf", height = 40)
par(mar=c(6,6,2,2)) # Increase margin size
# pal <- brewer.pal(8,"Dark2")
pal <- c("#fc1417", "#6dcc04", "#8b25fa")
barplot(sig_probes$significant_probes,
        names=sig_probes$Sample_Name,
        col=pal[factor(targets$TYPE)],
        xlab = "Significant number of probes",
        horiz = TRUE,
        las=1,
        cex.names=0.4,
        xlim = c(0,865859),
        border = NA
)
abline(v = 780000, col = "red")
dev.off()

# Also take a look, If any samples got significant probes <= 780000
filter(sig_probes, significant_probes <= 780000)

# Pre-Processing And Filtering
# ============================
# create a MethylSet object from the raw data for plotting
# Mset: MethylSet
Mset <- preprocessQuantile(RGset)
Mset
dim(Mset)

# ensure probes are in the same order in the Mset and detP objects
detP <- detP[match(featureNames(Mset),rownames(detP)),]

# keep only probes that have passed in at least 90% of the samples
per_samples <- round(((length(colnames(Mset)) * 90) / 100))
keep <- rowSums(detP < 0.01) >= per_samples
table(keep)

MsetFlt <- Mset[keep,]
message ("After keeping only probes that have passed in at least 90% of the samples")
dim(MsetFlt)

# get the 850k annotation data
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(MsetFlt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
table(keep)
MsetFlt <- MsetFlt[keep,]
message ("After removing probes on the sex chromosomes")
dim(MsetFlt)

# remove probes with SNPs at CpG site
MsetFlt <- dropLociWithSnps(MsetFlt)
message ("After removing probes with SNPs at CpG site")
dim(MsetFlt)

# exclude cross reactive probes
MsetFlt <- maxprobes::dropXreactiveLoci(MsetFlt)
message ("After excluding cross-reactive probes")
dim(MsetFlt)

# Extracting Beta and M-values
beta <- getBeta(MsetFlt)
mval <- getM(MsetFlt)

# Batch effects correction with SVA
batch <- targets$Batch
pheno <- targets %>% select(Sample_Name, SEX, TYPE, STAGE) # No need to include batch here
modcombat = model.matrix(~1, data=pheno)
combat_mval <- sva::ComBat(dat = mval, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
combat_beta <- lumi::m2beta(combat_mval)

# Writing out M and beta matrices
saveRDS(mval, "mval.rds")
saveRDS(beta, "beta.rds")

saveRDS(combat_mval, "mval_combat.rds")
saveRDS(combat_beta, "beta_combat.rds")

# Data exploration 1 - Before and After normalization
# ===================================================
# factor:Type, Batch - visualize what the data looks like before and after normalization
pdf("Density-Raw-Beta-RGset-Type.pdf")
pal <- c("#fc1417", "#6dcc04", "#8b25fa")
densityPlot(RGset,sampGroups=targets$TYPE, main="Beta-Value - Raw Data - Type", legend=FALSE, pal=pal)
legend("top", legend = levels(factor(targets$TYPE)),
       text.col=pal)
dev.off()

pdf("Density-Raw-Beta-RGset-Batch.pdf")
pal <- c("#FF7F00", "#6dcc04", "#fc51a8")
densityPlot(RGset, sampGroups=targets$Batch, main="Beta-Value - Raw Data - Batch", legend=FALSE, pal=pal)
legend("top", legend = levels(factor(targets$Batch)),
       text.col=pal)
dev.off()

pdf("Density-Norm-Beta-MsetFlt-Type.pdf")
pal <- c("#fc1417", "#6dcc04", "#8b25fa")
densityPlot(beta, sampGroups=targets$TYPE, main="Beta-Value - Normalized & Filtered Data - Type", legend = FALSE, pal=pal)
legend("top", legend = levels(factor(targets$TYPE)), 
       text.col=pal)
dev.off()

pdf("Density-Norm-Beta-MsetFlt-Batch.pdf")
pal <- c("#FF7F00", "#6dcc04", "#fc51a8")
densityPlot(beta, sampGroups=targets$Batch, main="Beta-Value - Normalized and Filtered - Batch", legend = FALSE, pal=pal)
legend("top", legend = levels(factor(targets$Batch)),
       text.col=pal)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Beta-MsetFlt-Type.pdf")
pal <- c("#fc1417", "#6dcc04", "#8b25fa")
plotMDS(beta,
        # top=1000, gene.selection="common",
        col=pal[factor(targets$TYPE)],
        pch=19,
        main="Beta-Value - Normalized & Filtered Data - Type"
)
legend("topleft",
       legend=levels(factor(targets$TYPE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Beta-MsetFlt-Batch.pdf")
pal <- c("#FF7F00", "#6dcc04", "#fc51a8")
plotMDS(beta,
        # top=1000, gene.selection="common",
        col=pal[factor(targets$Batch)],
        pch=19,
        main="Beta-Value - Normalized & Filtered Data - Batch"
)
legend("topleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-M-MsetFlt-Type.pdf")
pal <- c("#fc1417", "#6dcc04", "#8b25fa")
plotMDS(mval,
        # top=1000, gene.selection="common",
        col=pal[factor(targets$TYPE)],
        pch=19,
        main="M-Value - Normalized & Filtered Data - Type"
)
legend("topleft",
       legend=levels(factor(targets$TYPE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-M-MsetFlt-Batch.pdf")
pal <- c("#FF7F00", "#6dcc04", "#fc51a8")
plotMDS(mval,
        # top=1000, gene.selection="common",
        col=pal[factor(targets$Batch)],
        pch=19,
        main="M-Value - Normalized & Filtered Data - Batch"
)
legend("topleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# Data exploration 2 - after batch effects correction 
# ===================================================
# factor: Type, Batch - Visualize what the data looks like after batch effects correction
pdf("Density-Norm-Combat-Beta-Type.pdf")
pal <- c("#fc1417", "#6dcc04", "#8b25fa")
densityPlot(combat_beta, sampGroups=targets$TYPE, main="Beta-Value - Normalized, Filtered, Batch Corrected Data - Type", legend = FALSE, pal=pal)
legend("top", legend = levels(factor(targets$TYPE)), 
       text.col=pal)
dev.off()

pdf("Density-Norm-Combat-Beta-Batch.pdf")
pal <- c("#FF7F00", "#6dcc04", "#fc51a8")
densityPlot(combat_beta, sampGroups=targets$Batch, main="Beta-Value - Normalized, Filtered, Batch Corrected Data - Batch", legend = FALSE, pal=pal)
legend("top", legend = levels(factor(targets$Batch)), 
       text.col=pal)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-Beta-Type.pdf")
pal <- c("#fc1417", "#6dcc04", "#8b25fa")
plotMDS(combat_beta,
        # top=1000, gene.selection="common",
        col=pal[factor(targets$TYPE)],
        pch=19,
        main="Beta-Value - Normalized, Filtered, Batch Corrected Data - Type"
)
legend("topleft",
       legend=levels(factor(targets$TYPE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-Beta-Batch.pdf")
pal <- c("#FF7F00", "#6dcc04", "#fc51a8")
plotMDS(combat_beta,
        # top=1000, gene.selection="common",
        col=pal[factor(targets$Batch)],
        pch=19,
        main="Beta-Value - Normalized, Filtered, Batch Corrected Data - Batch"
)
legend("topleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-Beta-Stage.pdf")
pal <- c("#FF7F00", "#6dcc04")
plotMDS(combat_beta,
        # top=1000, gene.selection="common",
        col=pal[factor(targets$STAGE)],
        pch=19,
        main="Beta-Value - Normalized, Filtered, Batch Corrected Data - Stage"
)
legend("topleft",
       legend=levels(factor(targets$STAGE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-M-Type.pdf")
pal <- c("#fc1417", "#6dcc04", "#8b25fa")
plotMDS(combat_mval,
        # top=1000, gene.selection="common",
        col=pal[factor(targets$TYPE)],
        pch=19,
        main="M-Value - Normalized, Filtered, Batch Corrected Data - Type"
)
legend("topleft",
       legend=levels(factor(targets$TYPE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-M-Batch.pdf")
pal <- c("#FF7F00", "#6dcc04", "#fc51a8")
plotMDS(combat_mval,
        # top=1000, gene.selection="common",
        col=pal[factor(targets$Batch)],
        pch=19,
        main="M-Value - Normalized, Filtered, Batch Corrected Data - Batch"
)
legend("topleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-M-Stage.pdf")
pal <- c("#FF7F00", "#6dcc04")
plotMDS(combat_mval,
        # top=1000, gene.selection="common",
        col=pal[factor(targets$STAGE)],
        pch=19,
        main="M-Value - Normalized, Filtered, Batch Corrected Data - Stage"
)
legend("topleft",
       legend=levels(factor(targets$STAGE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()
