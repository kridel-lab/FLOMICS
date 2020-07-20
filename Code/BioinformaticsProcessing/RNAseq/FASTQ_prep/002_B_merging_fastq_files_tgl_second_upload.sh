#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=10440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J merging_fastq
#SBATCH --array=0-5 # job array index - number of jobs = numb of unique samples with top up runs

#need to run Mutect2 on all bam files from tumour samples
#author: Karin Isaev
#date started: June 25, 2019

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar
module load bam-readcount

#pwd
cd /cluster/projects/kridelgroup/TGL_transfers/TGL13_transfer_part2/fastq

#---------------------------------------------------------------------

#DO NOT RUN (R) - generate file with samples that need lane merging
#module load R/3.5.0

#setwd("/cluster/projects/kridelgroup/TGL13_transfer_part2")
#library(stringr)
#library(data.table)

#qc = fread("QC.txt")
#qc$samples_edited = paste(qc$samples, "topup", sep="")
#z1 = which(str_detect(qc$samples_edited, "-1topup"))
#z2 = which(str_detect(qc$samples_edited, "-2topup"))
#topups = qc[c(z1,z2),]
#nontopup = qc[-c(z1,z2),]

#topups$new_samples = ""
#topups$new_samples = topups$samples
#topups$clean_id = sapply(topups$new_samples, function(x){unlist(strsplit(x, "-"))[1]})

#topups = (as.data.frame(unique(topups$clean_id)))
#nontopup = (as.data.frame(unique(nontopup$samples)))

#write.table(topups, file="/cluster/projects/kridelgroup/FLOMICS/TOPUP_ANALYSIS_RNASEQ/topups_unique_sample_IDs.txt", quote=F, row.names=F, col.names=F)
#write.table(nontopup, file="/cluster/projects/kridelgroup/FLOMICS/TOPUP_ANALYSIS_RNASEQ/topups_unique_sample_IDs_notopups.txt", quote=F, row.names=F, col.names=F)

#---------------------------------------------------------------------

#pwd
samples=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TOPUP_ANALYSIS_RNASEQ/topups_unique_sample_IDs.txt #manually added patient LY_DLC_001 because they do have top-up
names=($(cat $samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

find ./ -name "*${names[${SLURM_ARRAY_TASK_ID}]}*" -and -name "*R1*" -exec cat {} + > /cluster/projects/kridelgroup/TGL_transfers/TGL13_transfer/data/fastq/topups_new_fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R1.fastq
find ./ -name "*${names[${SLURM_ARRAY_TASK_ID}]}*" -and -name "*R2*" -exec cat {} + > /cluster/projects/kridelgroup/TGL_transfers/TGL13_transfer/data/fastq/topups_new_fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R2.fastq

#pwd
#samples=/cluster/projects/kridelgroup/FLOMICS/TOPUP_ANALYSIS_RNASEQ/topups_unique_sample_IDs_notopups.txt
#names=($(cat $samples))
#echo ${names[${SLURM_ARRAY_TASK_ID}]}

#find ./ -name "*${names[${SLURM_ARRAY_TASK_ID}]}*" -and -name "*R1*" -exec cat {} + > /cluster/projects/kridelgroup/TGL13_transfer/data/fastq/topups_new_fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R1.fastq
#find ./ -name "*${names[${SLURM_ARRAY_TASK_ID}]}*" -and -name "*R2*" -exec cat {} + > /cluster/projects/kridelgroup/TGL13_transfer/data/fastq/topups_new_fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R2.fastq
