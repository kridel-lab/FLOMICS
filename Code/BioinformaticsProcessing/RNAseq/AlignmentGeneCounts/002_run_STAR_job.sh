#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=10440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J STAR_index

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar
module load bam-readcount
module load gcc
module load STAR

#pwd
cd /cluster/projects/kridelgroup/FLOMICS/DATA/TGL_FASTQ_RNASEQ

#ls *_R1.fastq > all_fastq_files.txt
#remove _R1.fastq from all lines
#sed -e 's!_R1.fastq!!' all_fastq_files.txt > all_fastq_files_cleaned.txt

for sample in $(cat all_fastq_files_cleaned.txt)    #for each patient sample...
do
    #submit patient specific job
    #here the input to the job is the patient name...
    sbatch /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/RNAseq/AlignmentGeneCounts/002_run_STAR_script.sh $sample
    echo $sample    #print folder or patient name...

done
