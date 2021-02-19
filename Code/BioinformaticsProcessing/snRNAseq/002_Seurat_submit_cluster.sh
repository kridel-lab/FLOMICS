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

#SC norm PC only
#file_name=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/pc_genes_only_yes_seurat_integrated_SCnorm_dim_20_2000_2021-02-17_samples_clusters.rds
#Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/002_Seurat_sarah_final_code_visualize_genes_across_clusters.R $file_name

#SC norm all genes
#file_name=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/pc_genes_only_no_seurat_integrated_SCnorm_dim_20_2000_2021-02-17_samples_clusters.rds
#Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/002_Seurat_sarah_final_code_visualize_genes_across_clusters.R $file_name

#standard workflow PC only
file_name=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/Feb2020/pc_genes_only_yes_seurat_integrated_dim_20_2000_2021-02-19_samples_clusters.rds
Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/002_Seurat_sarah_final_code_visualize_genes_across_clusters.R $file_name

#standard workflow all genes
#file_name=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/snRNAseq/seurat/pc_genes_only_no_seurat_integrated_dim_20_2000_2021-02-17_samples_clusters.rds
#Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/002_Seurat_sarah_final_code_visualize_genes_across_clusters.R $file_name
