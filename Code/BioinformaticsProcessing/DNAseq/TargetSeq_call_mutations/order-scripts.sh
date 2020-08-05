#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J run_all_scripts

#----navigate to directory with BAM files---------------------------------------
cd /cluster/projects/kridelgroup/GSC-1741

#save all BAM files in text file that will be used as array input for indices
find . -type f -name '*.bam' > /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/target_seq_samples_bam_locations.txt

#----run PICARD collect metrics-------------------------------------------------

#1. prepare interval and amplicon files


#----run PLATYPUS---------------------------------------------------------------

#mutation calling

#----run MUTECT2----------------------------------------------------------------
