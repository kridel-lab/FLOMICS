#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=21440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J Platypus
#SBATCH --array=0-130 # job array index

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
module load vt

cd /cluster/projects/kridelgroup/FLOMICS/ANALYSIS/PLATYPUS

#ls *.vcf > all_vcf_files

samples=all_vcf_files
names=($(cat $samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
sample=${MYVAR%%.platypus*}
echo $sample

out_folder=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/PLATYPUS/annovar
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/human_g1k_v37.decompressed.fasta

#normalize file
input_vcf=${names[${SLURM_ARRAY_TASK_ID}]}

#normalize variants, send to standard out and remove duplicates.
vt normalize $input_vcf -r $fasta_file -o ${sample}.normalized.vcf
