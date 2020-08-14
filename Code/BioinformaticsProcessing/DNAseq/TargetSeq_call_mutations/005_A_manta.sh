#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J manta
#SBATCH --array=0-130 # job array index - number of jobs = numb of unique samples with top up runs

module load strelka/2.9.10
module load python
module load manta/1.6.0
module load tabix

cd /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/BC_TargetSeq_Aug2020
#ls *.bam > all_bam_files_FL

samples=all_bam_files_FL
names=($(cat $samples))
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

gtf_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v19.annotation.gtf
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/human_g1k_v37.decompressed.fasta
out_folder=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA
#targets_interval_list=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/picard_tools_targets_input.bed
targets_interval_list=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/picard_tools_amps_input.bed

#MANTA
MANTA_INSTALL_PATH=/cluster/tools/software/centos7/manta/1.6.0

echo ${names[${SLURM_ARRAY_TASK_ID}]}

mkdir ${out_folder}/MANTA_WORKDIR_nointervals_${names[${SLURM_ARRAY_TASK_ID}]}
MANTA_ANALYSIS_PATH=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA/MANTA_WORKDIR_nointervals_${names[${SLURM_ARRAY_TASK_ID}]}

#index bed file
sort -k 1,1 -k 2,2n -k 3,3n $targets_interval_list | bgzip -c > ${targets_interval_list}.gz
tabix -pbed $targets_interval_list.gz

${MANTA_INSTALL_PATH}/bin/configManta.py \
--tumorBam ${names[${SLURM_ARRAY_TASK_ID}]} \
--referenceFasta $fasta_file --exome \
--runDir ${MANTA_ANALYSIS_PATH} #\
#--callRegions ${targets_interval_list}.gz

#After succesfful configuration run the following:
/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA/MANTA_WORKDIR_nointervals_${names[${SLURM_ARRAY_TASK_ID}]}/runWorkflow.py -j 20
