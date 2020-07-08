#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=21440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J telescope_summary
#SBATCH --array=0-135 # job array index

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load bedtools

#conda activate telescope_env

#Retrieve and print stats in the index file.
#The output is TAB-delimited with each line consisting of reference sequence name,
#sequence length, # mapped reads and # unmapped reads.

output=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS

cd $output

new_output=/cluster/projects/kridelgroup/FLOMICS/TELESCOPE_ANALYSIS/concatenated_results

ls -d  *sortedByCoord.out.bam_telescope_results > indiv_telescope_results

names=($(cat indiv_telescope_results))
echo ${names[${SLURM_ARRAY_TASK_ID}]}
cd ${names[${SLURM_ARRAY_TASK_ID}]}
scp *.tsv ${names[${SLURM_ARRAY_TASK_ID}]}.tsv
mv ${names[${SLURM_ARRAY_TASK_ID}]}.tsv $new_output

#cd ..
#rm -r ${names[${SLURM_ARRAY_TASK_ID}]}
