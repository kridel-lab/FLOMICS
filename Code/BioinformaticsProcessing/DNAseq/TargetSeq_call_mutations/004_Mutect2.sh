#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J picard
#SBATCH --array=0-130 # job array index - number of jobs = numb of unique samples with top up runs

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

cd /cluster/projects/kridelgroup/GSC-1741

samples=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/target_seq_samples_bam_locations.txt
names=($(cat $samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]}
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

name=${sample#*hg19a/}
echo $name

gtf_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v19.annotation.gtf
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/ucsc.hg19.fasta
out_folder=/cluster/projects/kridelgroup/FLOMICS/DATA/BC_TargetSeq_Calls
