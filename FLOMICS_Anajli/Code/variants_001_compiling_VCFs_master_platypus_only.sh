#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=25440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J capture_merge

#load all potentially required tools 
module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load bedtools
module load tabix

wd=/cluster/projects/kridelgroup/FLOMICS/TargetedDNAseq/DATA-157/variants
cd $wd 

#-------
#ls */*/*.vcf* > /cluster/projects/kridelgroup/FLOMICS/all_VCFs_samples.txt 
#ls */*/*variant_calls.pass.vcf* > /cluster/projects/kridelgroup/FLOMICS/all_VCFs_platypus_samples.txt 
#-------

#generate text file with names of all samples 
#-------
#ls > /cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder/samples_names.txt
#-------

#cd /cluster/projects/kridelgroup/FLOMICS
#sed '/md5/d' ./all_VCFs_samples.txt > nomd5_all_VCFs_samples.txt
#sed '/md5/d' ./all_VCFs_platypus_samples.txt > nomd5_all_VCFs_platypus_samples.txt

#sed '/Recal/d' ./nomd5_all_VCFs_samples.txt > no_recal.txt
#sed '/Scylla/d' ./no_recal.txt > no_scylla.txt
#sed '/Pisces/d' ./no_scylla.txt > no_pisces_nomd5_all_VCFs_samples.txt

vcfs=/cluster/projects/kridelgroup/FLOMICS/nomd5_all_VCFs_platypus_samples.txt

#go into analysis folder where we have writing permissions
cd /cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder

#In the working directory listed above, there are 131 folders (one for each sample)
#This script will submit a job on each sample
#This script will obtain the VCF file 
#This VCF will then be normalized (vt tool) and annotated (annovar tool)
#Finally this VCF will then be saved in a new directory: /cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder

for f in $(cat $vcfs)    #for each VCF file (for each sample)... 
do     
    #submit patient specific job 
    #here the input to the job is the VCF name 
    sbatch /cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder/scripts/variants_002_helper_script_master_platypus_only.sh $f 
    echo $f    #print folder or patient name... 

done




