#local exploration of mutations across clusters and stages
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors = F)
#avoid scientific notation
options(scipen=999)
#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr",
"mclust", "data.table", "plyr",
"ggrepel", "stringr", "maftools", "ggpubr", "readxl", "skimr",
 "edgeR", "annotables", "EnvStats", "MutationalPatterns", "BSgenome")
library(gridExtra)
lapply(packages, require, character.only = TRUE)

#date
date=Sys.Date()

#directory with FLOMICS related matrices
setwd("/Users/kisaev/github/FLOMICS/Data")

#----------------------------------------------------------------------
#data filtering
#----------------------------------------------------------------------

muts <- fread("2020-06-22_opossum_variant_FL_rna-seq_filtered.txt")
qc = as.data.table(read_excel("FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.xlsx"))

muts = as.data.table(filter(muts, Tumor_Sample_Barcode %in% qc$rna_seq_file_sample_ID))
length(unique(muts$SAMPLE_ID)) #128

#filter out mutation types we dont want
muts <- as.data.table(filter(muts, !(Variant_Classification %in% c("UTR3 .",
"downstream .", "UTR5 .", "upstream .", "exonic unknown"))))

muts$id = paste(muts$chr, muts$Start_Position, sep=":")
muts$bases = paste(muts$Reference_Allele, muts$Tumor_Seq_Allele2, sep="/")
muts$mut_id = paste(muts$id, muts$bases, sep="_")
muts$mut_gene = paste(muts$Hugo_Symbol, muts$mut_id)

#potential artificats
arts = as.data.table(table(muts$mut_gene))
arts = arts[order(-N)]
arts = as.data.table(filter(arts, N > (0.1*132)))
arts = as.data.table(filter(muts, mut_gene %in% arts$V1))

#-----------
#ANALYSIS
#-----------

#tiers
tiers=c("tier3", "tier2", "tier1")

#apply downstream code to each tier and then compare results
get_tier_input_dat = function(tier){

  #set up ref genome
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
  library(ref_genome, character.only = TRUE)

  if(tier == "tier3"){
    dat=muts
    print(length(unique(dat$SAMPLE_ID)))}

  if(tier == "tier2"){
    dat=muts
    z1 = which(dat$rRNAcontam > 40)
    z2 = which(dat$Uniquely.mapped < 10000000)
    z3 = which(dat$'Uniquely.mapped.reads..' < 50)
    z4 = which(dat$'X..of.reads.mapped.to.multiple.loci' > 20)

    total_lost = unique(c(z1, z2, z3, z4))
    dat=muts[-total_lost]
    print(length(unique(dat$SAMPLE_ID)))}

  if(tier == "tier1"){
    dat=muts
    z1 = which(dat$rRNAcontam > 40)
    z2 = which(dat$Uniquely.mapped < 10000000)
    z3 = which(dat$'Uniquely.mapped.reads..' < 70)
    z4 = which(dat$'X..of.reads.mapped.to.multiple.loci' > 20)

    total_lost = unique(c(z1, z2, z3, z4))
    dat=muts[-total_lost]
    print(length(unique(dat$SAMPLE_ID)))}

    #----------------------------------------------------------------------
    #reorganize into GRangesList
    #----------------------------------------------------------------------

    #granges object list
    #one data table for each patient
    #names of list are patient ID
    #seqnames, ranges, strand, paramRangeID, REF, ALT, QUAL, FILTER
    dat$strand = "*"
    dat$paramRangeID = NA
    dat$FILTER = "PASS"
    dat$id = paste(dat$chr, dat$Start_Position, sep=":")
    dat$bases = paste(dat$Reference_Allele, dat$Tumor_Seq_Allele2, sep="/")
    dat$mut_id = paste(dat$id, dat$bases, sep="_")
    dat = unique(dat[,c("Chromosome", "Start_Position", "End_Position",
    "strand", "paramRangeID",
    "Reference_Allele", "Tumor_Seq_Allele2", "Var_Freq", "FILTER", "mut_id",
    "SAMPLE_ID")])
    all_sampls=split(dat, by="SAMPLE_ID")

    get_df_granges = function(datt){
      datt = as.data.frame(datt)
      rownames(datt) = datt$mut_id
      datt$SAMPLE_ID = NULL
      datt$mut_id = NULL
      colnames(datt) = c("seqnames", "start", "end" ,"strand", "paramRangeID",
      "REF", "ALT", "QUAL", "FILTER")
      gr = makeGRangesFromDataFrame(datt, keep.extra.columns=T)
      return(gr)
    }

    #successful gr object lets see if it works as input for mutationmapper
    all_sampls_gr = llply(all_sampls, get_df_granges)
    type_occurrences <- mut_type_occurrences(all_sampls_gr, ref_genome)
    mut_mat <- mut_matrix(vcf_list = all_sampls_gr, ref_genome = ref_genome)
    #Mutation spectrum
    p1 <- plot_spectrum(type_occurrences)
    p2 <- plot_spectrum(type_occurrences, CT = TRUE)
    library("gridExtra")
    grid.arrange(p1, p2, ncol=2)

}
