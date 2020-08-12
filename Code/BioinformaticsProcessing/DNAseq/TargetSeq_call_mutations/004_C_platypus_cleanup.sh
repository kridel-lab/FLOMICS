#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=31440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J Platypus

module load R/3.6.1
Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/DNAseq/TargetSeq_call_mutations/004_B_platypus_cleanup.R
