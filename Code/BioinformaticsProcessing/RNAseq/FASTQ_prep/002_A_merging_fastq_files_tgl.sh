#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=10440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J merging_fastq
#SBATCH --array=0-54 # job array index - number of jobs = numb of unique samples with top up runs

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
cd /cluster/projects/kridelgroup/TGL_transfers/TGL13_transfer/data/fastq

#---------------------------------------------------------------------

#DO NOT RUN (R) - generate file with samples that need lane merging
#module load R/3.5.0

#setwd("/cluster/projects/kridelgroup/TGL13_transfer/QC")
#library(stringr)
#library(data.table)

#qc = fread("TGL13_run_lane_QC.txt")
#qc$samples_edited = paste(qc$samples, "topup", sep="")
#z1 = which(str_detect(qc$samples_edited, "_1topup"))
#z2 = which(str_detect(qc$samples_edited, "_2topup"))
#topups = qc[c(z1,z2),]

#nontopup = qc[-c(z1,z2),]

#topups$new_samples = ""
#z = which(str_detect(topups$samples, "T1"))
#topups$new_samples[z] = sapply(topups$samples[z], function(x){paste(unlist(strsplit(x, "_"))[1:4], collapse="_")})
#z = which(str_detect(topups$samples, "T2"))
#topups$new_samples[z] = sapply(topups$samples[z], function(x){paste(unlist(strsplit(x, "_"))[1:4], collapse="_")})
#z = which(!(str_detect(topups$samples, "T")))
#topups$new_samples[z] = sapply(topups$samples[z], function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})

#topups = (as.data.frame(unique(topups$new_samples)))
#nontopup = (as.data.frame(unique(nontopup$samples)))

#write.table(topups, file="/cluster/projects/kridelgroup/TGL13_transfer/data/fastq/topups_unique_sample_IDs.txt", quote=F, row.names=F, col.names=F)
#write.table(nontopup, file="/cluster/projects/kridelgroup/TGL13_transfer/data/fastq/nontopups_unique_sample_IDs.txt", quote=F, row.names=F, col.names=F)

#---------------------------------------------------------------------

#pwd
names=($(cat topups_unique_sample_IDs.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

find ./ -name "*${names[${SLURM_ARRAY_TASK_ID}]}*" -and -name "*R1*" -exec cat {} + > topups_new_fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R1.fastq
find ./ -name "*${names[${SLURM_ARRAY_TASK_ID}]}*" -and -name "*R2*" -exec cat {} + > topups_new_fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R2.fastq

#pwd
#names=($(cat nontopups_unique_sample_IDs.txt))
#echo ${names[${SLURM_ARRAY_TASK_ID}]}

#find ./ -name "*${names[${SLURM_ARRAY_TASK_ID}]}*" -and -name "*R1*" -exec cat {} + > topups_new_fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R1.fastq
#find ./ -name "*${names[${SLURM_ARRAY_TASK_ID}]}*" -and -name "*R2*" -exec cat {} + > topups_new_fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R2.fastq
