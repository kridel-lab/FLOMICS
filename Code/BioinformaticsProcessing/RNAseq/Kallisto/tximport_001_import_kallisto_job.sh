#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=60440M
#SBATCH -t 2-00:00 # Runtime in D-HH:MM

cd /cluster/projects/kridelgroup/FLOMICS

Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/RNAseq/Kallisto/tximport_001_import_kallisto_data.R
