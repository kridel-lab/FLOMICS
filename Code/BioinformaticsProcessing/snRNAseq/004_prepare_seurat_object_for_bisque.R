#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

options(stringsAsFactors=F)

#make sure loading R/3.6.1 before running this script

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA")

#load functions to analyze seurat
source("/cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/doSeuratProc.R")

#output directory
setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/")

#load libraries
library(ellipsis)
library(tidyverse)
library(splitstackshape)
packages <- c("readr", "data.table", "plyr",
	"stringr",
  "Seurat",
  "cowplot",
	"patchwork", "Biobase", "BisqueRNA")

lapply(packages, require, character.only = TRUE)

#-------------------------------------------------------------------------------
#purpose
#-------------------------------------------------------------------------------

#use seurat object to create single cell input data for bisque

#-------------------------------------------------------------------------------
#data
#-------------------------------------------------------------------------------

#1. read in processed seurat object
#combined = readRDS("seurat_integrated_dim_10_2000_samples_clusters.rds")
output="/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/Feb2020/"
combined=readRDS(paste(output, "pc_genes_only_no_seurat_integrated_dim_20_2000_2021-02-23_samples_clusters.rds", sep=""))

#-------------------------------------------------------------------------------
#analysis
#-------------------------------------------------------------------------------

#prepare for bisque
seurat.object = combined
delimiter= "_"
position = 2

get.cell.names <- function(obj) base::colnames(obj)
get.ident <- function(obj) Seurat::Idents(object=obj)
get.raw.data <- function(obj) Seurat::GetAssayData(object = obj[["RNA"]], slot = "counts")

individual.ids <- base::sapply(base::strsplit(get.cell.names(seurat.object),
                                                delimiter), `[[`, position)

base::names(individual.ids) <- get.cell.names(seurat.object)
individual.ids <- base::factor(individual.ids)
n.individuals <- base::length(base::levels(individual.ids))
base::message(base::sprintf("Split sample names by \"%s\"", delimiter),
                base::sprintf(" and checked position %i.", position),
                base::sprintf(" Found %i individuals.", n.individuals))
base::message(base::sprintf("Example: \"%s\" corresponds to individual \"%s\".",
                              get.cell.names(seurat.object)[1], individual.ids[1]))
sample.ids <- base::names(get.ident(seurat.object))
sc.pheno <- base::data.frame(check.names=F, check.rows=F,
                               stringsAsFactors=F,
                               row.names=sample.ids,
                               SubjectName=individual.ids,
                               cellType=get.ident(seurat.object))
sc.meta <- base::data.frame(labelDescription=base::c("SubjectName",
                                                       "cellType"),
                              row.names=base::c("SubjectName",
                                                "cellType"))
sc.pdata <- methods::new("AnnotatedDataFrame",
                           data=sc.pheno,
                           varMetadata=sc.meta)

sc.data <- base::as.matrix(get.raw.data(seurat.object)[,sample.ids,drop=F])
sc.eset <- Biobase::ExpressionSet(assayData=sc.data,
                                    phenoData=sc.pdata)

#save final object ready for BISQUE, right now the sample names are 1,2,3,4
saveRDS(sc.eset, "sc_eset_for_bisque.rds")
