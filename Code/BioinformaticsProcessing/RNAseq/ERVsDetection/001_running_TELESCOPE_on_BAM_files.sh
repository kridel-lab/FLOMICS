#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J telescope
#SBATCH --array=0-135 # job array index

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load bedtools

#the two commands below need to be run in terminal before
#submitting as job
conda init bash
conda activate telescope_env

#Retrieve and print stats in the index file.
#The output is TAB-delimited with each line consisting of reference sequence name,
#sequence length, # mapped reads and # unmapped reads.

cd /cluster/projects/kridelgroup/FLOMICS/DATA/TGL_BAM_RNASEQ_sorted_FASTQ
#ls *Aligned.sortedByCoord.out.bam  > summarize_bams_files

output=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS

names=($(cat summarize_bams_files))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

annotation=/cluster/projects/kridelgroup/FLOMICS/DATA/GRCh37_GENCODE_rmsk_TE.gtf #downloaded from ucsc table browser

mkdir $output/${names[${SLURM_ARRAY_TASK_ID}]}_telescope_results #make directory to store report for each sample, later will concattenate them all somehow

#originally was getting an error because apparently bam files had some mispaired reads
#to overcome this, run this command
#samtools view -f 2 -o ${names[${SLURM_ARRAY_TASK_ID}]}_fixed_reads.bam ${names[${SLURM_ARRAY_TASK_ID}]} #already done DO NOT RUN

#make temp directory for collated files
export TMPDIR=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS
export TEMP=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS
export TMP=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS

#collate BAM files and read pairs will be sequential
#keeps redirecting to /tmp - do not run!
#samtools collate ${names[${SLURM_ARRAY_TASK_ID}]} -o ${names[${SLURM_ARRAY_TASK_ID}]}.collate.bam


#***fragment count estimates are available in the report***
telescope assign ${names[${SLURM_ARRAY_TASK_ID}]} $annotation --outdir $output/${names[${SLURM_ARRAY_TASK_ID}]}_telescope_results --tempdir $TMPDIR --attribute transcript_id

echo done
