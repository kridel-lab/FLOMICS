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

#run once without mitochondrial genes and protein coding genes only
nc_genes_rm=yes
Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/001_Seurat_sarah_final_code_get_clusters.R $nc_genes_rm

#run once without mitochondrial genes and all genes in matrix
nc_genes_rm=no
Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/001_Seurat_sarah_final_code_get_clusters.R $nc_genes_rm
