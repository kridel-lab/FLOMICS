#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=40440M
#SBATCH -t 2-00:00 # Runtime in D-HH:MM

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar
module load bam-readcount
module load kallisto

cd /cluster/projects/kridelgroup/FLOMICS/DATA/TGL_FASTQ_RNASEQ
index_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/transcripts.idx
output=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/KALLISTO
gtf=/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v19.annotation.gtf
chrs=/cluster/projects/kridelgroup/FLOMICS/genome_files/hg19.chrlist

index=$1 #this refers to the sample ID
echo $index

mkdir ${output}/${index}

kallisto quant -i $index_file -o ${output}/${index} \
-b 100 ${index}_R1.fastq ${index}_R2.fastq --genomebam --gtf $gtf --chromosomes $chrs
