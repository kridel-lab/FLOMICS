#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=60000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J opossum_job

module load java/8  #8
module load samtools
module load gatk
module load annovar

cd /cluster/projects/kridelgroup/FLOMICS/DATA/TGL_BAM_RNASEQ_sorted_FASTQ

#submit patient specific job
#here the input to the job is the patient name...
index=$1
echo $index

MYVAR=$index
sample=${MYVAR%%Ali*}
echo $sample

out=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/opossum_platypus

snakemake -s /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/RNAseq/VariantCalling/001_call_variants_snakemake.snakefile $out/${sample}.vcf --rerun-incomplete --nolock

echo $index    #print folder or patient name...
