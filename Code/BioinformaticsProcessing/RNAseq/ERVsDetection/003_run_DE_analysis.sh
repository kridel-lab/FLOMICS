#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=60440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J telescope_counts_matrix

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load bedtools
module load R/3.6.1

Rscript FLOMICS/Code/BioinformaticsProcessing/RNAseq/ERVsDetection/003_B_Clusters_differential_expression_analysis.R
Rscript FLOMICS/Code/BioinformaticsProcessing/RNAseq/ERVsDetection/003_B_Stages_differential_expression_analysis.R
