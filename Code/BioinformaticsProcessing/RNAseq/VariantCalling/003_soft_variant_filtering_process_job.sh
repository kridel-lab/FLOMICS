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
#ls *.clean.chrs.vcf.hg19_multianno.sorted.reheader.vcf.gz > all_vcf_files

for sample in $(cat all_vcf_files)    #for each patient sample...
do
    #submit patient specific job
    #here the input to the job is the patient name...
    sbatch /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/RNAseq/VariantCalling/003_soft_variant_filtering_process.sh  $sample
    echo $sample    #print folder or patient name...

done
