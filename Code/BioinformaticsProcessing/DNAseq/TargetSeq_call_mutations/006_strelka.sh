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
#index bed file
sort -k 1,1 -k 2,2n -k 3,3n $targets_interval_list | bgzip -c > ${targets_interval_list}.gz
tabix -pbed $targets_interval_list.gz

mkdir ${out_folder}/STRELKA_WORKDIR_${names[${SLURM_ARRAY_TASK_ID}]}
MANTA_ANALYSIS_PATH=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA/MANTA_WORKDIR_${names[${SLURM_ARRAY_TASK_ID}]}

#STRELKA
STRELKA_INSTALL_PATH=/cluster/tools/software/centos7/strelka/2.9.10

STRELKA_ANALYSIS_PATH=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA/STRELKA_WORKDIR_${names[${SLURM_ARRAY_TASK_ID}]}
MANTA_ANALYSIS_PATH=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA/MANTA_WORKDIR_${names[${SLURM_ARRAY_TASK_ID}]}

${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
--tumorBam ${names[${SLURM_ARRAY_TASK_ID}]} \
--referenceFasta $fasta_file \
--indelCandidates ${MANTA_ANALYSIS_PATH}/results/variants/candidateSmallIndels.vcf.gz \
--runDir ${STRELKA_ANALYSIS_PATH}

#After succesfful configuration run the following:
/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA/STRELKA_WORKDIR_${sample}/runWorkflow.py -j 8 -m local
