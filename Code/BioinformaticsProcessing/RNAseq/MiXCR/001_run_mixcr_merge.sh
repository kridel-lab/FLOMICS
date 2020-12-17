#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=60440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J mixcr_matrix

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load bedtools
module load R/3.6.1

Rscript /cluster/home/srussell/github/FLOMICS/Code/BioinformaticsProcessing/RNAseq/MiXCR/002_mixcr_concatenate_all_clonotypes.R
