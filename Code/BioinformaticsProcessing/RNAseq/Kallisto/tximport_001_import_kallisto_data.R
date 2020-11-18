#----------------------------------------------------------------------------------
#tximport_001_import_kallisto_data.R
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
#load packages
#----------------------------------------------------------------------------------

date = Sys.Date()

library(data.table)
library(dplyr)
library(plyr)
library(reshape2)
library(edgeR)
library(tidyverse)
library(readxl)
library(tximportData)
library(tximport)

#prepare gene annotation
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(readr)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
dir <- system.file("extdata", package = "tximportData")
library(readr)

library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file="/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v27lift37.annotation.gtf")
saveDb(x=txdb, file = "gencode.v27.annotation.TxDb")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)
dim(tx2gene)

#tx2gene = read.csv("/cluster/projects/kridelgroup/FLOMICS/genome_files/tx2gene.gencode.v27.lift37.csv")
#tx2gene$X = NULL
tx2gene$TXNAME = sapply(tx2gene$TXNAME, function(x){unlist(strsplit(x, "_"))[1]})
tx2gene$GENEID = sapply(tx2gene$GENEID, function(x){unlist(strsplit(x, "_"))[1]})
head(tx2gene)

#----------------------------------------------------------------------------------
#data
#----------------------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/KALLISTO")

#----------------------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------------------

#import kallisto results for each sample
samples = as.data.frame(list.files())
samples$type = "kallisto"
colnames(samples)[1] = "run"

#define all abundance h5 files for each sample
files <- file.path(getwd(), samples$run, "abundance.h5")
names(files) <- samples$run

#test
files = files[1:2]

#read in all files into matrix
#txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
head(txi.kallisto$counts)

txi.kallisto_t <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto_t$counts)

tpm = txi.kallisto$counts
#remove the "." in the gene id
rownames(tpm) = sapply(rownames(tpm), function(x){unlist(strsplit(x, "\\."))[1]})

#save results
write.table(tpm, file=paste("/cluster/projects/kridelgroup/FLOMICS/DATA/", date,
"_kallisto_gene_based_counts.txt", sep=""), quote=F,
row.names=T, sep=";")
