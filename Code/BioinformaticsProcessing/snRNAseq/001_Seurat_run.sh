#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J seurat_run

module load R/4.0.0

#pwd
cd /cluster/projects/kridelgroup/FLOMICS

Rscript /cluster/home/kisaev/
