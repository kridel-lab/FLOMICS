#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Seurat/Bisque workflow
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cd
module load R/4.0.0

#1. run seurat on gene count data processed through cellranger (3 patient samples)
sbatch /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/001_Seurat_submit_cluster.sh

#2. get gene marker plots across seurat clusters (final seurat object with 20 dimensions)
sbatch /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/002_Seurat_submit_cluster.sh

#3. rename clusters after reviewing gene marker plots (with manual review by Robert)
Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/003_Seurat_rename_clusters_to_cell_types.R

#4. prepare seurat object for bisque analysis
Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/004_prepare_seurat_object_for_bisque.R

#5. perform DE analysis between clusters of interest using seurat object
Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/005_Seurat_differential_expression_analysis_clusters.R

#6. cell cycle regression from seurat object
Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/snRNAseq/006_Seurat_sarah_final_code_cell_cycle_analysis.R

#7. represent boxplots of bisque cell fractions
  #FLOMICS/Code/Analysis/RNAseq/RNAseq-immune-deconvolution-bisque.R

#8. generate descriptive table for bisque boxplots
  #FLOMICS/Code/Analysis/RNAseq/RNAseq-immune-deconvolution-bisque-table.R
