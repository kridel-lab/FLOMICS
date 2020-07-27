#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MIXCRJ

module load java/8  #8
module load python
module load gcc
module load perl
module load mixcr

cd /cluster/projects/kridelgroup/FLOMICS/DATA/TGL_FASTQ_RNASEQ

export TMPDIR=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS
export TEMP=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS
export TMP=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MIXCR_ANALYSIS

for sample in $(cat all_fastq_files_cleaned.txt)    #for each patient sample
do
    #submit patient specific job
    #input to the job is sample ID
    sbatch /cluster/home/srussell/mixcr_analysis.sh $sample
    echo $sample    #print folder or saple ID

done
