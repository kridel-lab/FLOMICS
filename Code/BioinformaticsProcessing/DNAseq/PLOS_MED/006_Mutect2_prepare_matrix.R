#----------------------------------------------------------------------
#karin isaev
#post mutect2 and annovar soft filtering and matrix prep
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)
date=Sys.Date()

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/EGA_calls/MUTECT2/annovar")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#mutect2 run on tumour only mode
#annotated by annovar

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#gene annotations
genes = unique(fread("/cluster/projects/kridelgroup/paired_cns/ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#genes = fread("/cluster/home/kisaev/data/annotables_grch37.txt")
#genes = as.data.table(filter(genes, biotype == "protein_coding"))
#genes = as.data.table(filter(genes, !(is.na(entrez))))

files = list.files(pattern=".vcf.gz_annovar.vcf.gz.hg19_multianno.vcf")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. clean up individual paired vcf files

clean_up_001 = function(paired_vcf){

  mutations_T1 =read.vcfR(paired_vcf)
  mutations_T1 = vcfR2tidy(mutations_T1)
  meta = as.data.table(mutations_T1$meta)

  vcf = as.data.table(mutations_T1$fix)
  chr_conv = unique(vcf[,c("ChromKey", "CHROM")])

  gt = as.data.table(mutations_T1$gt)
  gt = merge(gt, chr_conv, by="ChromKey")

  gt$mut_id = paste(gt$CHROM, gt$POS, sep="_")
  vcf$mut_id = paste(vcf$CHROM, vcf$POS, sep="_")

  #1. keep only the ones that passed default mutect2 filter
  vcf = as.data.table(filter(vcf, FILTER=="PASS"))
  print(paste("number of variants that passed filtering=", dim(vcf)[1]))

	#2. combine vcf and gt info
  cols = colnames(gt)[which(colnames(gt) %in% colnames(vcf))]
  gt = merge(gt, vcf, by= cols)

	#3. get hugo gene names
  gt = merge(gt, genes, by= "Gene.ensGene")
  print(paste("number of variants that passed gene merge", dim(gt)[1]))

  #4. filter by coverage
  #try 30X
  vcf = as.data.table(filter(vcf, DP >=10))
  print(paste("number of variants that passed coverage=", dim(vcf)[1]))

  #5. keep only those with population allele frequency < 0.001
  gt$controls_AF_popmax = as.numeric(gt$controls_AF_popmax)
  gt = as.data.table(filter(gt, (controls_AF_popmax < 0.001 | is.na(controls_AF_popmax))))
  print(paste("number of variants that passed controls_AF_popmax=", dim(gt)[1]))

  #6. keep only mutations where t2 VAF > 0.1
	#gt$gt_AF = as.numeric(gt$gt_AF)
	#gt = as.data.table(filter(gt, gt_AF >= 0.1))
	print(paste("number of variants that passed vaf >= 0.1=", dim(gt)[1]))

	#7. keep only muts affecting pcg regions
	gt = as.data.table(filter(gt, !(ExonicFunc.ensGene == "."),
	!(ExonicFunc.ensGene == "synonymous_SNV")))
	print(paste("number of variants that passed pcg muts only", dim(gt)[1]))

  #8. generate bed file - summary of mutation and coordinates to intersect with cnvkit output
  pat = unlist(strsplit(paired_vcf, ".filter."))[1]
	gt$sample=pat

	colnames(gt)[which(colnames(gt)=="hg19.ensemblToGeneName.value")] = "Hugo_Symbol"
	gt$End_Position = gt$POS
	gt$Start_Position = gt$POS
	gt$Chromosome = gt$CHROM
	gt$Reference_Allele = gt$REF
	gt$Tumor_Seq_Allele2 = gt$ALT
	gt$Variant_Classification = paste(gt$Func.ensGene, gt$ExonicFunc.ensGene)
	gt$Variant_Type = "SNP"
	gt$Tumor_Sample_Barcode = gt$sample
	gt$Var_Freq = gt$gt_AF

	gt = gt[,c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
	"Tumor_Seq_Allele2", "avsnp142", "cosmic68", "AAChange.ensGene",
	"Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "Var_Freq")]

	return(gt)

  print("done")
}

all_muts = as.data.table(ldply(llply(files, clean_up_001, .progress="text")))
write.table(all_muts, file=paste("/cluster/projects/kridelgroup/FLOMICS/DATA/",
date, "_Mutect2_filtered_mutations_PLOS_MED.txt", sep=""), quote=F, row.names=F, sep=";")
