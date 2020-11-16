#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=60440M
#SBATCH -t 2-00:00 # Runtime in D-HH:MM

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar
module load bam-readcount
module load kallisto

cd /cluster/projects/burst2/TGL_BAM_RNASEQ_sorted_FASTQ
index_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/transcripts.idx
output=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/KALLISTO
gtf=/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v27lift37.annotation.gtf
chrs=/cluster/projects/kridelgroup/FLOMICS/genome_files/hg19.chrlist

index=$1 #this refers to the sample ID
echo $index

mkdir ${output}/${index}

#need to use sorted FASTQ files
kallisto quant -i $index_file -o ${output}/${index} \
-b 100 ${index}_R1_sorted.fastq.gz ${index}_R2_sorted.fastq.gz --genomebam --gtf $gtf --chromosomes $chrs
