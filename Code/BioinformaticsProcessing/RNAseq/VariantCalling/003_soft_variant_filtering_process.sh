#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=5000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J annovar_sum

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar

cd /cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/annovar

index=$1 #this refers to the patient sample
echo $index

#submit patient specific job
#here the input to the job is the patient name...

Rscript /cluster/home/kisaev/bioinformatics_misc/RNASEQ_TGL_ALIGNMENT_VARIANT_CALLING/opossum_platypus_snakemake/variant_annotation_process.R  $index
