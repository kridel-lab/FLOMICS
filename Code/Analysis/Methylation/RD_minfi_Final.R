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

# Read IDAT Files
# ===============
# set the working directory
baseDir <- setwd("/cluster/home/t110989uhn/kridelgroup/rajesh/01_DNA_Methylome_Analysis/02_minfi/02_FL_Samples")
list.files(baseDir)

# Read a sample sheet
targets <- read.metharray.sheet(baseDir)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet
RGset <- read.metharray.exp(targets = targets, force=TRUE)
RGset

# Give the samples descriptive names
targets$ID <- paste(targets$Sample_Name)
sampleNames(RGset) <- targets$ID
RGset

# See which annotation package being used by
annotation(RGset)

# get the 850k EPIC annotation data
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


# Quality control
# ===============
# Calculate the detection p-values
detP <- detectionP(RGset)
head(detP)

# Examine mean detection p-values across all samples to identify any failed samples
pval <- (colMeans(detP))
write.table(pval, file = "pval.tsv", sep="\t", quote = FALSE)

# Examine mean detection p-values across all samples to identify any failed samples
pdf("pval.pdf")
pal <- brewer.pal(8,"Dark2")
barplot(colMeans(detP),
        col=pal[factor(targets$Type)],
        las=1,
        cex.names=0.1,
        xlab = "Mean detection p-values",
        xlim = c(0,0.6),
        horiz = TRUE
)
abline(v = 0.01, col = "red")
dev.off()

# remove poor quality samples
keep <- colMeans(detP) < 0.01
RGset <- RGset[,keep]
RGset

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

# Re-Examine mean detection p-values across all samples to identify any failed samples
pval <- (colMeans(detP))
write.table(pval, file = "pval_updated.tsv", sep="\t", quote = FALSE)

# Re-Examine mean detection p-values across all samples to identify any failed samples
pdf("pval_updated.pdf")
pal <- brewer.pal(8,"Dark2")
barplot(colMeans(detP),
        col=pal[factor(targets$Type)],
        las=1,
        cex.names=0.1,
        xlab = "Mean detection p-values",
        xlim = c(0,0.6),
        horiz = TRUE
)
abline(v = 0.01, col = "red")
dev.off()

# qcReport
# qcReport(RGset, sampNames=targets$Unique_ID, sampGroups=targets$Type, pdf="qcReport.pdf")
# Error: subscript contains out-of-bounds indices


#   Pre-Processing And Filtering
#   ============================
# create a MethylSet object from the raw data for plotting
# Mset: MethylSet, 865,859 features, 381 arrays, 381 samples
Mset <- preprocessQuantile(RGset)
Mset

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(Mset),rownames(detP)),]

# Examine significant number of probes in each samples
Sample_Probes <- (colSums(detP <= 0.01))
write.table(Sample_Probes, file = "Sample_Probes.tsv", sep="\t", quote = FALSE)

pdf("Sample_Probes.pdf")
pal <- brewer.pal(8,"Dark2")
barplot(Sample_Probes,
        col=pal[factor(targets$Type)],
        las=1,
        cex.names=0.1,
        xlab = "Significant number of probes in each samples",
        xlim = c(0,865859),
        horiz = TRUE
)
abline(v = 780000, col = "red")
dev.off()

# Examine significant number of probes detected across all samples
All_Sample_Probes <- (rowSums(detP <= 0.01))
write.table(All_Sample_Probes, file = "All_Sample_Probes.tsv", sep="\t", quote = FALSE)

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(Mset)
table(keep)
write.table(keep, file = "Passed_Failed_Probes.tsv", sep="\t", quote = FALSE)

MsetFlt <- Mset[keep,]
message ("After removing probes that have failed in one or more samples")
MsetFlt

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(MsetFlt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
table(keep)
MsetFlt <- MsetFlt[keep,]
message ("After removing probes on the sex chromosomes")
MsetFlt

# remove probes with SNPs at CpG site
MsetFlt <- dropLociWithSnps(MsetFlt)
message ("After removing probes with SNPs at CpG site")
MsetFlt

# exclude cross reactive probes
MsetFlt <- maxprobes::dropXreactiveLoci(MsetFlt)
message ("After excluding cross-reactive probes")
MsetFlt


# Data exploration 1
# ===================

# factor:Type, visualize what the data looks like before and after normalization
pdf("betaval-RGset-Type.pdf")
densityPlot(RGset,sampGroups=targets$Type, main="Beta-Value Raw Data", legend=FALSE)
legend("top", legend = levels(factor(targets$Type)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

pdf("betaval-MsetFlt-Type.pdf")
densityPlot(getBeta(MsetFlt),sampGroups=targets$Type, main="Beta-Value Normalized (preprocessQuantile) & Filtered Data", legend=FALSE)
legend("top", legend = levels(factor(targets$Type)), 
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
        col=pal[factor(targets$Type)],
        pch=19,
        main="Sample Type Variation Normalized & Filtered Data"
)
legend("topleft",
       legend=levels(factor(targets$Type)),
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
mval <- minfi::getM(MsetFlt)
beta <- minfi::getBeta(MsetFlt)

# Batch effects correction with SVA
batch <- targets$Batch
pheno <- targets %>% select(Unique_ID, Sex, Type, STAGE) # No need to include batch here
modcombat = model.matrix(~1, data=pheno)
combat_mval <- sva::ComBat(dat = mval, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
combat_beta <- lumi::m2beta(combat_mval)

# Writing out M and beta matrices
saveRDS(mval, "mval.rds")
saveRDS(beta, "beta.rds")

saveRDS(combat_mval, "mval_combat.rds")
saveRDS(combat_beta, "beta_combat.rds")


# Data exploration 2
# ===================

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
