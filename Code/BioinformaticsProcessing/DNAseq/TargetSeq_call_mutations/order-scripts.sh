#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J run_all_scripts

#----navigate to directory with BAM files---------------------------------------
cd /cluster/projects/kridelgroup/GSC-1741

#FASTA file downloaded from 1000 genomes

#save all BAM files in text file that will be used as array input for indices
find . -type f -name '*.bam' > /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/target_seq_samples_bam_locations.txt

#----run PICARD collect metrics-------------------------------------------------

cd /cluster/projects/kridelgroup/FLOMICS

#1. prepare interval and amplicon files
Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/DNAseq/TargetSeq_call_mutations/001_Prepare_Amplicons_Targets.R

#2. run collect targeted pcr metrics
sbatch /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/DNAseq/TargetSeq_call_mutations/002_CollectTargetedPcrMetrics.sh

#3. process coverage results/merge with probe info which genes cover and convert
#sample id
Rscript /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/DNAseq/TargetSeq_call_mutations/003_summarize_probe_coverage.R

#----run PLATYPUS---------------------------------------------------------------

#mutation calling

#----run MUTECT2----------------------------------------------------------------

#----run strelka and manta----------------------------------------------------------------

#make file gzi for fasta
sbatch /cluster/home/kisaev/FLOMICS/Code/BioinformaticsProcessing/DNAseq/TargetSeq_call_mutations/005_strelka.sh
