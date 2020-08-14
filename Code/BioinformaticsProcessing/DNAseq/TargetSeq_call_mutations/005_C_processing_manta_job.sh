#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J manta
#SBATCH --array=0-130 # job array index - number of jobs = numb of unique samples with top up runs

module load strelka/2.9.10
module load python
module load manta/1.6.0
module load tabix

cd /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/BC_TargetSeq_Aug2020
#ls *.bam > all_bam_files_FL

samples=all_bam_files_FL
names=($(cat $samples))
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

script=/cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/DNAseq/TargetSeq_call_mutations/005_B_processing_manta.R

Rscript $script $sample
