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
out=/cluster/projects/kridelgroup/FLOMICS/DATA/TGL_BAM_RNASEQ

#two pass mapping + get gene counts
STAR --genomeDir /cluster/projects/kridelgroup/FLOMICS/genome_files/ \
--readFilesIn $path1 $path2 --outSAMtype BAM SortedByCoordinate \
--outSAMunmapped None --outFileNamePrefix $out/${index} --quantMode GeneCounts \
--twopassMode Basic --readFilesCommand zcat