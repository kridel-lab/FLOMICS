#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=5000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J opossum_main

#python 2.7 needs to be automatically loaded

module load java/8  #8
module load samtools
module load gatk
module load annovar

cd /cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STAR/
#ls  *Aligned.sortedByCoord.out.bam  > all_star_files

for sample in $(cat all_star_files)    #for each patient sample...
do
    #submit patient specific job
    #here the input to the job is the patient name...
    sbatch /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/RNAseq/VariantCalling/001_call_variants_opossum_start_snake.sh  $sample
    echo $sample    #print folder or patient name...

done
