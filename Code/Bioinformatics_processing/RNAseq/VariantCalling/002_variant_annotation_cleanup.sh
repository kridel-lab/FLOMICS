#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=20000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J annovar_clean
#SBATCH -o %j.out			# File to which standard out will be written
#SBATCH -e %j.err
#SBATCH --array=0-135 # job array index - number of jobs = numb of unique samples with top up runs

module load annovar
module load bam-readcount
module load gcc
module load STAR
module load rsem
module load perl
module load python/2.7
module load gatk
module load tabix
module load vcftools

cd /cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/opossum_platypus

#ls *.vcf > all_vcf_files

samples=all_vcf_files
names=($(cat $samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
sample=${MYVAR%%.vcf*}
echo $sample

gtf_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v19.annotation.gtf
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/ucsc.hg19.fasta
out_folder=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/annovar

#first keep only variants that passed baseline encoded variants already done
awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' ${sample}.vcf > ${sample}.clean.vcf

#remove all the weird chromosome names - already done
grep -Ev '^(chr17_ctg5_hap1|chrUn_gl000211|chrUn_gl000220|chrUn_gl000224|chr1_gl000192_random|chr6_apd_hap1|chr6_cox_hap2|chr6_dbb_hap3|chr6_mann_hap4|chr6_qbl_hap6|chr6_ssto_hap7|chrY|chrX|chrM)' ${sample}.clean.vcf > ${sample}.clean.chrs.vcf

#run annovar already done
table_annovar.pl --buildver hg19 ${sample}.clean.chrs.vcf /cluster/tools/software/annovar/humandb --protocol ensGene,gnomad211_genome,cosmic68,avsnp142 --operation g,f,f,f --vcfinput

#gzip, index and sort VCF file
bgzip ${sample}.clean.chrs.vcf.hg19_multianno.vcf #- already done on all 136 files

#remove two lines below:
#bcftools index ${sample}.clean.chrs.vcf.hg19_multianno.vcf.gz
#bcftools sort -Oz ${sample}.clean.chrs.vcf.hg19_multianno.vcf.gz -o ${sample}.clean.chrs.vcf.hg19_multianno.sorted.vcf.gz
vcf-sort ${sample}.clean.chrs.vcf.hg19_multianno.vcf.gz > ${sample}.clean.chrs.vcf.hg19_multianno.sorted.vcf.gz

#add cotig lengths to header and index
bcftools reheader ${sample}.clean.chrs.vcf.hg19_multianno.sorted.vcf.gz --fai /cluster/projects/kridelgroup/FLOMICS/genome_files/ucsc.hg19.fasta.fai > ${sample}.clean.chrs.vcf.hg19_multianno.sorted.reheader.vcf.gz

#Unzipping the VCFs and zipping them again with bgzip and later indexing with IndexFeatureFile resolved the issue.
bgzip ${sample}.clean.chrs.vcf.hg19_multianno.sorted.reheader.vcf.gz
bcftools index ${sample}.clean.chrs.vcf.hg19_multianno.sorted.reheader.vcf.gz.gz

#get proper GATK VCF index file
gatk IndexFeatureFile -F ${sample}.clean.chrs.vcf.hg19_multianno.sorted.reheader.vcf.gz.gz

#use GATK to filter VCF based on filters recommended for RNA-seq variant calling
#gatk VariantFiltration \
#   -R $fasta_file \
#   -V ${sample}.clean.chrs.vcf.hg19_multianno.sorted.reheader.vcf.gz.gz  \
#   -O ${sample}.clean.chrs.vcf.hg19_multianno.sorted.reheader.filtered.vcf.gz \
#   -window 35 -cluster 3 --filter-name "FS" --filter-expression "FS > 30.0" --filter-name "QD" --filter-expression "QD < 2.0"

#move VCF file to annovar results folder
mv ${sample}.clean.chrs.vcf.hg19_multianno.sorted.reheader.vcf.gz.gz /cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/annovar
