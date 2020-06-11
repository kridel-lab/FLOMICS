#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=60440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J index

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar
module load bam-readcount
module load gcc
module load STAR
module load rsem
module load perl

#pwd
cd /cluster/projects/kridelgroup/FLOMICS/DATA/TGL_FASTQ_RNASEQ

index=$1 #this refers to the patient sample
echo $index

path1=${index}_R1.fastq
path2=${index}_R2.fastq
genome=/cluster/projects/kridelgroup/FLOMICS/genome_files
out=/cluster/projects/kridelgroup/FLOMICS/DATA/TGL_BAM_RNASEQ_sorted_FASTQ

#sort FASTQ files by their sequence identifier
zcat $path1 \
| paste - - - - \
| sort -k1,1 -S 5G \
| tr '\t' '\n' \
| gzip > ${index}_R1_sorted.fastq.gz

zcat $path2 \
| paste - - - - \
| sort -k1,1 -S 5G \
| tr '\t' '\n' \
| gzip > ${index}_R2_sorted.fastq.gz

#zcat $path1 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${index}_R1_sorted.fastq
#zcat $path2 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${index}_R2_sorted.fastq

#two pass mapping + get gene counts
STAR --genomeDir /cluster/projects/kridelgroup/FLOMICS/genome_files/ \
--readFilesIn ${index}_R1_sorted.fastq.gz ${index}_R2_sorted.fastq.gz --outSAMtype BAM SortedByCoordinate \
--outSAMunmapped None --outFileNamePrefix $out/${index} --quantMode GeneCounts \
--twopassMode Basic --readFilesCommand zcat
