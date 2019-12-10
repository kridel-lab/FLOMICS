#!/bin/bash
#

#load all potentially required tools 
module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load bedtools
module load tabix

index=$1 #this refers to the algorithms used to call variants 
echo $index

wd=/cluster/projects/kridelgroup/FLOMICS/TargetedDNAseq/DATA-157/variants
cd $wd 

#generate text file with names of all samples 
#ls > /cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder/samples_names.txt

#go into analysis folder where we have writing permissions
cd /cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder

#In the working directory listed above, there are 131 folders (one for each sample)
#This script will submit a job on each sample
#This script will obtain the VCF file with SNPs obtained from merging Playpus and LoFreq 
#This VCF will then be normalized (vt tool) and annotated (annovar tool)
#Finally this VCF will then be saved in a new directory: /cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder

for dir in $(cat samples_names.txt)    #for each patient folder... 
do     
    #submit patient specific job 
    #here the input to the job is the patient name... 
    sbatch /cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder/scripts/variants_002_helper_script_lofreq_Indels_only.sh $dir $index
    echo $dir    #print folder or patient name... 

done




