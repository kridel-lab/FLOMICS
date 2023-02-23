
#--
# This script executes first pass filtering upon variant calls
# post Mutect2 and Annovar analysis and prepares a matrix of
# variants calls for further filtering downstream.
# Author: Victoria Shelton
# Last modified: February 23rd, 2023
# Tested on R/4.1.0
#--

#--
# Load in packages
#--

packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust",
 "data.table", "plyr", "ggrepel", "stringr", "maftools", "magrittr")
lapply(packages, require, character.only = TRUE)

#--
# Setup working space
#--

options(stringsAsFactors = FALSE)
date <- Sys.Date()

# set directory to where the annotated vcf files of the capseq pipeline
## were deposited
setwd("/your working directory/")
dir <- "/your working directory/"

#--
# Read in data
#--

#gene annotations
genes <- unique(fread("ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] <- "Gene.ensGene"

files <- list.files(pattern = "_annovar.vcf.gz.hg19_multianno.vcf")

#--
# Filter calls and prepare matrix
#--

###clean up individual paired vcf files
clean_up_001 <- function(paired_vcf) {

  mutations_t1 <- read.vcfR(paired_vcf)
  mutations_t1 <- vcfR2tidy(mutations_t1)
  meta <- as.data.table(mutations_t1$meta)

  vcf <- as.data.table(mutations_t1$fix)
  chr_conv <- unique(vcf[, c("ChromKey", "CHROM")])

  gt <- as.data.table(mutations_t1$gt)
  gt <- merge(gt, chr_conv, by = "ChromKey", allow.cartesian = TRUE)

  gt$mut_id <- paste(gt$CHROM, gt$POS, sep = "_")
  vcf$mut_id <- paste(vcf$CHROM, vcf$POS, sep = "_")

  #1. keep only variants that have filter annotated as 'pass'
  vcf <- as.data.table(filter(vcf, FILTER == "PASS"))
  print(paste("number of variants that passed filtering=", dim(vcf)[1]))

  #2. combine vcf and gt info
  cols <- colnames(gt)[which(colnames(gt) %in% colnames(vcf))]
  gt <- merge(gt, vcf, by = cols, allow.cartesian = TRUE)

  #3. get hugo gene names
  gt <- merge(gt, genes, by = "Gene.ensGene", allow.cartesian = TRUE)
  print(paste("number of variants that passed gene merge", dim(gt)[1]))

  #4. filter by coverage (depth)
  vcf <- as.data.table(filter(vcf, DP >= 10))
  print(paste("number of variants that passed coverage=", dim(vcf)[1]))

  #5. keep only those with population allele frequency < 0.001
  gt$controls_AF_popmax = as.numeric(gt$controls_AF_popmax)
  gt <- as.data.table(filter(gt,
   (controls_AF_popmax < 0.001 | is.na(controls_AF_popmax))))
  print(paste("number of variants that passed controls_AF_popmax=", dim(gt)[1]))

  #6. keep only mutations where VAF > 0.1
  gt$gt_AF <- as.numeric(gt$gt_AF)
  gt <- as.data.table(filter(gt, gt_AF >= 0.1))
  print(paste("number of variants that passed vaf >= 0.1=", dim(gt)[1]))

  #8. generate bed file -
  # summary of mutation and coordinates
  pat <- unlist(strsplit(paired_vcf, "_annovar"))[1]
  gt$sample <- pat

  colnames(gt)[which(colnames(gt) == "Gene.ensGene")] <- "hg19.ensGene.name2"
  gt$Hugo_Symbol <- 	gt$hg19.ensemblToGeneName.value
  gt$End_Position <- gt$POS
  gt$Start_Position <- gt$POS
  gt$Chromosome <- gt$CHROM
  gt$Reference_Allele <- gt$REF
  gt$Tumor_Seq_Allele2 <- gt$ALT
  gt$Variant_Classification <- paste(gt$Func.ensGene, gt$ExonicFunc.ensGene)
  gt$Variant_Type <- "SNP"
  gt$Tumor_Sample_Barcode <- gt$sample
  gt$Var_Freq <- gt$gt_AF
  gt$Allelic_Depth <- gt$gt_AD
  gt$Read_Depth <- gt$gt_DP

  gt <- gt[, c("hg19.ensGene.name2", "Hugo_Symbol", "Chromosome",
   "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
    "avsnp142", "cosmic68", "AAChange.ensGene", "Variant_Classification",
    "Variant_Type", "Tumor_Sample_Barcode", "Var_Freq", "Allelic_Depth",
     "Read_Depth", "FILTER")]

  return(gt)

  print("done")
}

all_muts <- as.data.table(ldply(llply(files, clean_up_001, .progress = "text")))
print(all_muts$Allelic_Depth)

all_depth <- str_split_fixed(all_muts$Allelic_Depth, ",", 2)
all_muts <- cbind(all_muts, all_depth)
all_muts <- rename(all_muts, c("V1" = "Ref_counts", "V2" = "alt_counts"))
all_muts <- select(all_muts, -15)
print("all_muts made")

###write out matrix
table_path = paste0(dir, date, "variant_calls.txt")
write.table(all_muts, file = table_path, quote = F, row.names = F, sep = ";")
