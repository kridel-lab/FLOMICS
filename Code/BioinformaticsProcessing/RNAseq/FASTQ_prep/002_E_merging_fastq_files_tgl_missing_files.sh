#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=10440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J merging_fastq
#SBATCH --array=0-2 # job array index - number of jobs = numb of unique samples with top up runs

#need to run Mutect2 on all bam files from tumour samples
#author: Karin Isaev
#date started: June 25, 2019

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar
module load bam-readcount

#pwd
cd /cluster/projects/kridelgroup/TGL_transfers/TGL13_missing_data/fastq

#---------------------------------------------------------------------

#pwd
samples=missing_fastq.txt
names=($(cat $samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

find ./ -name "*${names[${SLURM_ARRAY_TASK_ID}]}*" -and -name "*R1*" -exec cat {} + > /cluster/projects/kridelgroup/TGL_transfers/TGL13_transfer/data/fastq/topups_new_fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R1.fastq
find ./ -name "*${names[${SLURM_ARRAY_TASK_ID}]}*" -and -name "*R2*" -exec cat {} + > /cluster/projects/kridelgroup/TGL_transfers/TGL13_transfer/data/fastq/topups_new_fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R2.fastq
