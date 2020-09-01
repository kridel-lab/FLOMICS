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
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
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
#files = files[1:5]

#read in all files into matrix
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
head(txi.kallisto$counts)

tpm = txi.kallisto$counts
#remove the "." in the gene id
rownames(tpm) = sapply(rownames(tpm), function(x){unlist(strsplit(x, "\\."))[1]})

#save results
write.table(tpm, file=paste("/cluster/projects/kridelgroup/FLOMICS/DATA/", date,
"_kallisto_gene_based_counts.txt", sep=""), quote=F,
row.names=F, sep=";")
