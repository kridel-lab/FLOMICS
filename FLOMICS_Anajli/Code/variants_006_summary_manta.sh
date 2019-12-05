#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J process_manta
#SBATCH -c 8
	
module load R/3.5.0
cd /cluster/projects/kridelgroup/FLOMICS/

#define list of samples

for sample in $(cat samples_names.txt) 
do
	echo $sample
	Rscript /cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder/scripts/variants_005_summary_manta.R $sample 

done