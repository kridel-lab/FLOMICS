#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=60440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J mixcr
#SBATCH --array=0-135 # job array index

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

names=($(cat all_fastq_files_cleaned.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

path1=${names[${SLURM_ARRAY_TASK_ID}]}_R1.fastq
path2=${names[${SLURM_ARRAY_TASK_ID}]}_R2.fastq

#sort files
zcat $path1 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" | gzip > $sorted/${names[${SLURM_ARRAY_TASK_ID}]}_R1_sorted.fastq.gz
zcat $path2 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" | gzip > $sorted/${names[${SLURM_ARRAY_TASK_ID}]}_R2_sorted.fastq.gz

#use MiXCR ANALYZE SHOTGUN command for entire pipeline
mixcr analyze shotgun \
  -s hsa \
  --starting-material rna \
  $sorted/${names[${SLURM_ARRAY_TASK_ID}]}_R1_sorted.fastq.gz $sorted/${names[${SLURM_ARRAY_TASK_ID}]}_R2_sorted.fastq.gz $out/${names[${SLURM_ARRAY_TASK_ID}]}
