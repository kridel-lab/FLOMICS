#----------------------------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
  "ggrepel", "stringr", "maftools", "VariantAnnotation")
lapply(packages, require, character.only = TRUE)
library(GenomicRanges)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
library(biomaRt)
library(readxl)

date = Sys.Date()

print(date)

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MANTA/PROCESSED/")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#gene annotations
genes = unique(fread("/cluster/projects/kridelgroup/paired_cns/ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#ensembl gene annotations with coordaintes
all_genes = readRDS("/cluster/home/kisaev/data/ensembl_biomaRt_coordinate_data.rds")
all_genes$chromosome_name = paste("chr", all_genes$chromosome_name, sep="")
all_genes$geneid = paste(all_genes$chromosome_name, all_genes$start_position, all_genes$end_position, sep="_")
all_genes = unique(all_genes[,c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "geneid")])

genes = unique(fread("/cluster/home/kisaev/data/annotables_grch37.txt"))
genes = as.data.table(filter(genes, biotype == "protein_coding"))
genes = as.data.table(filter(genes, !(is.na(entrez))))

#add extra window to start and end of gene to capture gene even if not directly overlapping
all_genes$start_position = all_genes$start_position - 500
all_genes$end_position = all_genes$end_position + 500
all_genes = as.data.table(filter(all_genes, chromosome_name %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
  "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
  "chr20", "chr21", "chr22", "chrX", "chrY", "chrMT")))
all_genes = as.data.table(filter(all_genes, !(hgnc_symbol == "")))

#Sample info (originally provided by BC team in first data upload)
samp_info = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/library_mapping_BC.csv")

#detailed sample info (provided by Anjali and Robert)
more_samp_info = as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/sample_annotations_rcd6Nov2019.xlsx"))

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

all_svs = list.files(pattern="_2020-08-13_")

get_res = function(file){
  f = readRDS(file)
  if(!(dim(f)[1]==0)){
  pat = f$pat[1]
  sample = unlist(strsplit(pat, "_"))[1]
  sample_check = unlist(strsplit(pat, "-"))[1]
  if(!(sample_check == pat)){
      sample = sample_check
    }
  print(sample)
  f$Library = sample
  return(f)
}}

all_res = as.data.table(ldply(llply(all_svs, get_res)))
all_res = merge(all_res, samp_info, by="Library")
#merge with detaild sample information
colnames(more_samp_info)[1] = "External_identifier"
all_res = merge(all_res, more_samp_info, by = "External_identifier")
all_res$RNAseq_COMMENT=NULL

write.table(all_res, file=paste("/cluster/projects/kridelgroup/FLOMICS/DATA/",
date, "_Manta_SVs_pass.txt", sep=""), quote=F, row.names=F, sep=";")
