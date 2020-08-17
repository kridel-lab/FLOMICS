#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=20440M
#SBATCH -t 2-00:00 # Runtime in D-HH:MM

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar
module load bam-readcount
module load kallisto

#pwd
cd /cluster/projects/kridelgroup/FLOMICS

#ls *_R1.fastq > all_fastq_files.txt
#remove _R1.fastq from all lines
#sed -e 's!_R1.fastq!!' all_fastq_files.txt > all_fastq_files_cleaned.txt

for sample in $(cat  /cluster/projects/kridelgroup/FLOMICS/DATA/TGL_FASTQ_RNASEQ/all_fastq_files_cleaned.txt)    #for each patient folder...
do
    #submit patient specific job
    #here the input to the job is the patient name...
    sbatch /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/RNAseq/Kallisto/kallisto_002_quantification_script.sh $sample
    echo $sample    #print folder or patient name...

done
