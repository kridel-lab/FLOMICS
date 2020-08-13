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

date = Sys.Date()

print(date)
args = commandArgs(trailingOnly = TRUE) #patient ID
index = args[1]
print(index)

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA")
setwd(paste("MANTA_WORKDIR_", index, "/", "results/variants", sep=""))

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

vcf = "tumorSV.vcf.gz" #normalized annotated vcf files (n=131)

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

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. clean up individual paired vcf files

clean_up_001 = function(vcf){

  pat = index
  if(!(dim(fread(vcf))[1] == 3)){ #to prevent empty VCF files from throwing errors

  #read in VCF file
  vcf_dat <- readVcf(vcf, "hg19")

  #extract genotype info into data table
  vcf_dat_ranges = as.data.frame(rowRanges(vcf_dat))
  vcf_dat_ranges$id = rownames(vcf_dat_ranges)
  vcf_dat_info = info(vcf_dat)
  vcf_dat_info$id = rownames(vcf_dat_info)
  vcf_dat = as.data.table(merge(vcf_dat_ranges, vcf_dat_info, by = "id"))

  #add patient ID
  vcf_dat$patient = pat

  #filter to include SVs that passed
  vcf_dat = as.data.table(filter(vcf_dat, IMPRECISE == "FALSE"))

  #turn SV coordinates into Granges object and intersect with genes
  vcf_dat$MATEID = as.character(vcf_dat$MATEID)
  vcf_dat$ALT = as.character(vcf_dat$ALT)
  vcf_dat$SVLEN = as.character(vcf_dat$SVLEN)
  vcf_dat$HOMLEN = as.character(vcf_dat$HOMLEN)
  vcf_dat$HOMSEQ = as.character(vcf_dat$HOMSEQ)
  vcf_dat$SVINSLEN = as.character(vcf_dat$SVINSLEN)
  vcf_dat$SVINSSEQ = as.character(vcf_dat$SVINSSEQ)

  vcf_dat = unique(vcf_dat[,c("seqnames", "start", "end", "MATEID", "REF", "ALT", "id", "SVTYPE", "SVLEN", "IMPRECISE", "BND_DEPTH", "MATE_BND_DEPTH", "HOMLEN", "HOMSEQ", "SVINSLEN",
    "SVINSSEQ")])
  vcf_dat$seqnames = paste("chr", vcf_dat$seqnames, sep="")
  vcf_dat_coords = vcf_dat
  vcf_dat_coords = makeGRangesFromDataFrame(vcf_dat_coords)

  all_genes_coords = unique(all_genes[,c("chromosome_name", "start_position", "end_position", "hgnc_symbol")])
  colnames(all_genes_coords) = c("seqnames", "start", "end", "gene")
  all_genes_coords = makeGRangesFromDataFrame(all_genes_coords, keep.extra.columns =TRUE)

  hits <- findOverlaps(vcf_dat_coords, all_genes_coords, ignore.strand=TRUE)
  hits_overlap = cbind(as.data.table(vcf_dat[queryHits(hits),]), as.data.table(all_genes_coords)[subjectHits(hits),])
  print(head(hits_overlap))
  colnames(hits_overlap)[1:3] = c("SV_CHR", "SV_start", "SV_end")

  #if translocation (BND) - need to make sure we preserved the mate
  #means region didnt map to a gene
  #bring it back to dataset though otherwise will miss it

  z = which(hits_overlap$SVTYPE == "BND")
  BND = hits_overlap[z,]
  nonBND = hits_overlap[-z,]
  bnds_ids = unique(BND$id)

  get_BND_mate = function(bnd){
    print(bnd)
    get_mate = filter(BND, id == bnd)$MATEID
    #check if in dataset
    check = (dim(filter(BND, id == get_mate))[1] >=1)
    print(check)
    if(!(check)){
      dat = as.data.table(filter(vcf_dat, id == get_mate))
      dat_old = dat
      #add colnames
      cols_add = colnames(BND)[which(!(colnames(BND) %in% colnames(dat)))]
      new_dat = as.data.table(matrix(ncol=length(cols_add), nrow=1))
      colnames(new_dat) = cols_add
      dat = cbind(dat, new_dat)
      dat$SV_CHR = dat$seqnames
      dat$SV_start = dat$start
      dat$SV_end = dat$end

      return(dat)
    }
    }

  missing_BNDs = as.data.table(ldply(llply(bnds_ids, get_BND_mate)))

  #add back all pieces
  vcf_dat = rbind(BND, missing_BNDs, nonBND)
  vcf_dat$pat = pat

  #retrun
  return(vcf_dat)
}
}

processed_vcf = clean_up_001(vcf)
print("done")
print(index)
saveRDS(processed_vcf, file=paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA/PROCESSED/", date, index, ".rds", sep="_"))
