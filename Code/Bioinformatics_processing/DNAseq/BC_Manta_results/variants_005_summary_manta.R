#----------------------------------------------------------------------
#variants_003_read_in_VCFs_into_matrix.R
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
setwd("/cluster/projects/kridelgroup/FLOMICS")

print(date)
args = commandArgs(trailingOnly = TRUE) #sample ID 
index = args[1] 
print(index) 

setwd(paste("TargetedDNAseq/DATA-157/structural_variants/", index, "/", "manta-1.4.0/bwa-mem-0.7.6a-sb/125nt/hg19a/", sep=""))
setwd(list.files())

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?

#----------------------------------------------------------------------
#data 
#----------------------------------------------------------------------

vcf = "tumorSV.vcf" 

#gene annotations
genes = unique(fread("/cluster/projects/kridelgroup/paired_cns/ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#ensembl gene annotations with coordaintes 
all_genes = readRDS("/cluster/home/kisaev/data/ensembl_biomaRt_coordinate_data.rds")
all_genes$chromosome_name = paste("chr", all_genes$chromosome_name, sep="")
all_genes$geneid = paste(all_genes$chromosome_name, all_genes$start_position, all_genes$end_position, sep="_")
all_genes = unique(all_genes[,c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "geneid")])

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
  vcf_dat = as.data.table(filter(vcf_dat, FILTER == "PASS", IMPRECISE == "FALSE"))

  #turn SV coordinates into Granges object and intersect with genes 
  vcf_dat$MATEID = as.character(vcf_dat$MATEID)
  vcf_dat$ALT = as.character(vcf_dat$ALT)
  vcf_dat = unique(vcf_dat[,c("seqnames", "start", "end", "MATEID", "REF", "ALT", "id", "SVTYPE", "IMPRECISE", "BND_DEPTH", "MATE_BND_DEPTH")])
  vcf_dat$seqnames = paste("chr", vcf_dat$seqnames, sep="")
  vcf_dat_coords = vcf_dat
  vcf_dat_coords = makeGRangesFromDataFrame(vcf_dat_coords)
   
  all_genes_coords = unique(all_genes[,c("chromosome_name", "start_position", "end_position")])
  colnames(all_genes_coords) = c("seqnames", "start", "end")
  all_genes_coords = makeGRangesFromDataFrame(all_genes_coords)

  hits <- findOverlaps(vcf_dat_coords, all_genes_coords, ignore.strand=TRUE)
  hits_overlap = cbind(as.data.table(vcf_dat[queryHits(hits),]), as.data.table(all_genes_coords)[subjectHits(hits),])
  print(head(hits_overlap))
  colnames(hits_overlap)[1:3] = c("SV_CHR", "SV_start", "SV_end")
  hits_overlap$geneid = paste(hits_overlap$seqnames, hits_overlap$start, hits_overlap$end, sep="_")

  hits_overlap = merge(hits_overlap, all_genes, by = c("geneid"))

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
      #add colnames
      cols_add = colnames(BND)[which(!(colnames(BND) %in% colnames(dat)))]
      new_dat = as.data.table(matrix(ncol=length(cols_add), nrow=1))
      colnames(new_dat) = cols_add
      dat = cbind(dat, new_dat)
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
saveRDS(processed_vcf, file=paste("/cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder/MANTA_SVs/", date, index, ".rds", sep="_"))
