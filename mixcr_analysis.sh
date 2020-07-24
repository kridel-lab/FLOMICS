#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=60440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J mixcr

module load java/8  #8
module load python
module load gcc
module load perl
module load mixcr

#pwd
cd /cluster/projects/kridelgroup/FLOMICS/DATA/TGL_FASTQ_RNASEQ
sorted=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS
out=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS

export TMPDIR=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS
export TEMP=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS
export TMP=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS

index=$1 #this refers to the patient sample
echo $index

path1=${index}_R1.fastq
path2=${index}_R2.fastq

#sort files
zcat $path1 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" | gzip > $sorted/${index}_R1_sorted.fastq.gz
zcat $path2 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" | gzip > $sorted/${index}_R2_sorted.fastq.gz

#use MiXCR ANALYZE SHOTGUN command for entire pipeline
mixcr analyze shotgun \
  -s hsa \
  --starting-material rna \
  $sorted/${index}_R1_sorted.fastq.gz $sorted/${index}_R2_sorted.fastq.gz $out/${index}
