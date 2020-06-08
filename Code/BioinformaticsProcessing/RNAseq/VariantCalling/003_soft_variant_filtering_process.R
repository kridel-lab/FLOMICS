#----------------------------------------------------------------------
#processing_annovar_results.R
#karin isaev
#July 11th, 2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/annovar")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#mutect2 was run on paired mode compaing cns to diagnostic tumour
#now it's time to:
#summarize cns specific mutations
#but first should still filter out false positives (note, these are unfilitered variants)
#see how many appear in multiple comparisons (n=5 total)

#note these vcf files have been normalized and fed through annovar
#for annotations

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE) #patient ID
index = args[1]
print(index)

#gene annotations
genes = unique(fread("/cluster/projects/kridelgroup/paired_cns/ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#genes = fread("/cluster/home/kisaev/data/annotables_grch37.txt")
#genes = as.data.table(filter(genes, biotype == "protein_coding"))
#genes = as.data.table(filter(genes, !(is.na(entrez))))

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

  #2. filter by coverage
  #try 60X
  vcf = as.data.table(filter(vcf, TC >=10))
  print(paste("number of variants that passed coverage=", dim(vcf)[1]))

  #3. combine vcf and gt info
  cols = colnames(gt)[which(colnames(gt) %in% colnames(vcf))]
  gt = merge(gt, vcf, by= cols)

  #4. remove potential snps annotated by gnomad
  #z = which(str_detect(gt$avsnp142, "rs"))
  #if(!(length(z)==0)){
  #gt = gt[-z,]}
  #print(paste("number of variants that passed avsnp142=", dim(gt)[1]))

  #6. keep only those with population allele frequency < 0.001
  gt$controls_AF_popmax = as.numeric(gt$controls_AF_popmax)
  gt = as.data.table(filter(gt, (controls_AF_popmax < 0.001 | is.na(controls_AF_popmax))))
  print(paste("number of variants that passed controls_AF_popmax=", dim(gt)[1]))

  #7. keep only mutations where t2 VAF > 0.1
#  gt$gt_AF = as.numeric(gt$gt_AF)
#  gt = as.data.table(filter(gt, gt_AF >= 0.1))
#  print(paste("number of variants that passed vaf >= 0.1=", dim(gt)[1]))

  #8. remove variants from chromosome X and Y
  gt = as.data.table(filter(gt, !(CHROM %in% c("chrX", "chrY", "chrM", "chrUn_gl000220"))))
	z = which(str_detect(gt$CHROM, "chrUn"))
	if(!(length(z) == 0)){
		gt=gt[-z,]
	}
	z = which(str_detect(gt$CHROM, "gl"))
	if(!(length(z) == 0)){
		gt=gt[-z,]
	}
	z = which(str_detect(gt$CHROM, "mcf"))
	if(!(length(z) == 0)){
		gt=gt[-z,]
	}

  print(paste("number of variants that passed X Y=", dim(gt)[1]))

  #9. get hugo gene names
  gt = merge(gt, genes, by= "Gene.ensGene")
  print(paste("number of variants that passed gene merge", dim(gt)[1]))

  #10. generate bed file - summary of mutation and coordinates to intersect with cnvkit output
  pat = unlist(strsplit(paired_vcf, ".clean."))[1]
	gt$sample=pat
  file = paste(pat, "_final_vcf_file_filtered_maf_input.bed", collapse="_", sep="")

	#11. finally edit file so it is maf compatabile
	#Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele,
	#Tumor_Seq_Allele2, Variant_Classification, Variant_Type and Tumor_Sample_Barcode.

	#also keep only SNPs
	#gt = as.data.table(filter(gt, REF %in% c("A", "C", "G", "T"), ALT %in% c("A", "C", "G", "T")))
	#print(paste("number of variants that passed SNPs only", dim(gt)[1]))

	colnames(gt)[which(colnames(gt)=="hg19.ensemblToGeneName.value")] = "Hugo_Symbol"
	gt$End_Position = gt$POS
	gt$Start_Position = gt$POS
	gt$Chromosome = gt$CHROM
	gt$Reference_Allele = gt$REF
	gt$Tumor_Seq_Allele2 = gt$ALT
	gt$Variant_Classification = paste(gt$Func.ensGene, gt$ExonicFunc.ensGene)
	gt$Variant_Type = "SNP"
	gt$Tumor_Sample_Barcode = gt$sample
	gt$Var_Freq = gt$TC

	gt = gt[,c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
	"Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "Var_Freq")]
	write.table(gt, file, quote=F, row.names=F, sep="\t", col.names=T)

  print("done")
}

clean_up_001(index)
