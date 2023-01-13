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
baseDir <- setwd("/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/01_DNA_Methylome_Analysis/02_minfi/02_FL_Samples/")
list.files(baseDir)

# Read a sample sheet
targets <- read.metharray.sheet(baseDir)

# Remove 2 machine control samples
targets <- targets %>% filter(Sample_Name != 'HCT116_DKO_methylated')

# Remove 2 samples for clinical reason
targets <- targets %>% filter(Sample_Name != 'LY_FL_311_T1' & Sample_Name != 'LY_FL_159_T1')

# Rename replicate of sample "LY_FL_159_T1"
targets$Sample_Name <- str_replace(targets$Sample_Name, "LY_FL_159_T1_rep", "LY_FL_159_T1")

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
sample_ann <- read.table(file = "20221228_sample_annotations.txt",
                         sep = "\t",
                         header = TRUE,
                         na.strings=c("","NA"))

# Left outer join sample sheet and sample annotation
targets <- merge(x = targets,
              y = sample_ann,
              by.x = "Sample_Name",
              by.y = "SAMPLE_ID",
              all.x = TRUE)

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
# [1] "LY_FL_158_T1.R01C01" "LY_FL_479_T1.R02C01" "LY_FL_488_T1.R06C01"
# [4] "LY_FL_498_T1.R06C01" "LY_FL_523_T1.R01C01" "LY_FL_524_T1.R02C01"
# [7] "LY_FL_525_T1.R03C01" "LY_FL_527_T1.R04C01" "LY_FL_529_T1.R01C01"

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
pal <- brewer.pal(8,"Dark2")
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
pal <- brewer.pal(8,"Dark2")
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
# per_samples <- round(((length(colnames(Mset)) * 90) / 100))
# keep <- rowSums(detP < 0.01) >= per_samples
# table(keep)

# MsetFlt <- Mset[keep,]
# message ("After keeping only probes that have passed in at least 90% of the samples")
# dim(MsetFlt)

# get the 850k annotation data
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)

# if your data includes males and females, remove probes on the sex chromosomes
# keep <- !(featureNames(MsetFlt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
# table(keep)
# MsetFlt <- MsetFlt[keep,]
# message ("After removing probes on the sex chromosomes")
# dim(MsetFlt)

# remove probes with SNPs at CpG site
# MsetFlt <- dropLociWithSnps(MsetFlt)
# message ("After removing probes with SNPs at CpG site")
# dim(MsetFlt)

# exclude cross reactive probes
# MsetFlt <- maxprobes::dropXreactiveLoci(MsetFlt)
# message ("After excluding cross-reactive probes")
# dim(MsetFlt)

MsetFlt <- Mset
dim(MsetFlt)

# Data exploration 1 - Before and After normalization
# ===================================================
# factor:Type, visualize what the data looks like before and after normalization
pdf("betaval-RGset-Type.pdf")
densityPlot(RGset,sampGroups=targets$TYPE, main="Beta-Value Raw Data", legend=FALSE)
legend("top", legend = levels(factor(targets$TYPE)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

pdf("betaval-MsetFlt-Type.pdf")
densityPlot(getBeta(MsetFlt),sampGroups=targets$TYPE, main="Beta-Value Normalized (preprocessQuantile) & Filtered Data", legend=FALSE)
legend("top", legend = levels(factor(targets$TYPE)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

# factor:STAGE, visualize what the data looks like before and after normalization
pdf("betaval-RGset-STAGE.pdf")
densityPlot(RGset,sampGroups=targets$STAGE, main="Beta-Value Raw Data", legend=FALSE)
legend("top", legend = levels(factor(targets$STAGE)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

pdf("betaval-MsetFlt-STAGE.pdf")
densityPlot(getBeta(MsetFlt),sampGroups=targets$STAGE, main="Beta-Value Normalized (preprocessQuantile) & Filtered Data", legend=FALSE)
legend("top", legend = levels(factor(targets$STAGE)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

# MDS plots to look at largest sources of variation
pdf("mds-plot-MsetFlt-type.pdf")
plotMDS(getM(MsetFlt),
        col=pal[factor(targets$TYPE)],
        pch=19,
        main="Sample Type Variation Normalized & Filtered Data"
)
legend("topleft",
       legend=levels(factor(targets$TYPE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()


pdf("mds-plot-MsetFlt-stage.pdf")
plotMDS(getM(MsetFlt),
        col=pal[factor(targets$STAGE)],
        pch=19,
        main="Sample Stage Variation Normalized & Filtered Data"
)
legend("topleft",
       legend=levels(factor(targets$STAGE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# Extracting Beta and M-values
beta <- minfi::getBeta(MsetFlt)
mval <- minfi::getM(MsetFlt)

# Batch effects correction with SVA
batch <- targets$Batch
pheno <- targets %>% select(Sample_Name, SEX, TYPE, STAGE) # No need to include batch here
modcombat = model.matrix(~1, data=pheno)
combat_mval <- sva::ComBat(dat = mval, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
combat_beta <- lumi::m2beta(combat_mval)

# Remove suspicious repeated samples
# ==================================
beta <- subset(beta, select = -c(LY_FL_1156_T1,
                                 LY_FL_571_T1
))

mval <- subset(mval, select = -c(LY_FL_1156_T1,
                                 LY_FL_571_T1
))

combat_mval <- subset(combat_mval, select = -c(LY_FL_1156_T1,
                                               LY_FL_571_T1
))

combat_beta <- subset(combat_beta, select = -c(LY_FL_1156_T1,
                                               LY_FL_571_T1
))

# Writing out M and beta matrices
saveRDS(mval, "mval.rds")
saveRDS(beta, "beta.rds")

saveRDS(combat_mval, "mval_combat.rds")
saveRDS(combat_beta, "beta_combat.rds")

# Data exploration 2 - Before and after batch effects correction 
# ==============================================================

# factor: sample plate, visualize what the data looks like before and after batch effects correction 
pdf("betaval-before-batch.pdf")
densityPlot(beta,sampGroups=targets$Batch, main="Beta-Value Before Batch Effects Correction", legend=FALSE)
legend("top", legend = levels(factor(targets$Batch)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

pdf("betaval-after-batch.pdf")
densityPlot(combat_beta,sampGroups=targets$Batch, main="Beta-Value After Batch Effects Correction", legend=FALSE)
legend("top", legend = levels(factor(targets$Batch)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

pdf("mval-before-batch.pdf")
densityPlot(mval,sampGroups=targets$Batch, main="M-Value Before Batch Effects Correction", legend=FALSE)
legend("topleft", legend = levels(factor(targets$Batch)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

pdf("mval-after-batch.pdf")
densityPlot(combat_mval,sampGroups=targets$Batch, main="M-Value After Batch Effects Correction", legend=FALSE)
legend("topleft", legend = levels(factor(targets$Batch)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

# MDS plots to look at largest sources of variation

# Beta-Value before and after Batch Effects Correction
pdf("betavalue-mds-plot-before-batch.pdf")
plotMDS(beta,
        col=pal[factor(targets$Batch)],
        pch=19,
        main="Beta Value Variation Before Batch Effects Correction"
)
legend("bottomleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

pdf("betavalue-mds-plot-after-batch.pdf")
plotMDS(combat_beta,
        col=pal[factor(targets$Batch)],
        pch=19,
        main="Beta Value Variation After Batch Effects Correction"
)
legend("bottomleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# M-Value before and after Batch Effects Correction
pdf("mval-mds-plot-before-batch.pdf")
plotMDS(mval,
        col=pal[factor(targets$Batch)],
        pch=19,
        main="M Value Variation Before Batch Effects Correction"
)
legend("bottomleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()


pdf("mval-mds-plot-after-batch.pdf")
plotMDS(combat_mval,
        col=pal[factor(targets$Batch)],
        pch=19,
        main="M Value Variation After Batch Effects Correction"
)
legend("bottomleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation - after batch affect correction
pdf("mds-plot-MsetFlt-type-after-batch.pdf")
plotMDS(combat_mval,
        col=pal[factor(targets$TYPE)],
        pch=19,
        main="Sample Type Variation After Batch Affect Correction"
)
legend("topleft",
       legend=levels(factor(targets$TYPE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()


pdf("mds-plot-MsetFlt-stage-after-batch.pdf")
plotMDS(combat_mval,
        col=pal[factor(targets$STAGE)],
        pch=19,
        main="Sample Stage Variation After Batch Affect Correction"
)
legend("topleft",
       legend=levels(factor(targets$STAGE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()
