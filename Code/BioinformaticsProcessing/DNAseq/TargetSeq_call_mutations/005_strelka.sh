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

cd /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/BC_TargetSeq_Aug2020
ls *.bam > all_bam_files_FL

samples=all_bam_files_FL
names=($(cat $samples))
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

gtf_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v19.annotation.gtf
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/human_g1k_v37.fasta
out_folder=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA

#MANTA
MANTA_INSTALL_PATH=/cluster/tools/software/centos7/manta/1.6.0

echo ${names[${SLURM_ARRAY_TASK_ID}]}

mkdir ${out_folder}/MANTA_WORKDIR_${names[${SLURM_ARRAY_TASK_ID}]}
MANTA_ANALYSIS_PATH=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA/MANTA_WORKDIR_${names[${SLURM_ARRAY_TASK_ID}]}

${MANTA_INSTALL_PATH}/bin/configManta.py \
--tumorBam ${names[${SLURM_ARRAY_TASK_ID}]} \
--referenceFasta $fasta_file \
--runDir ${MANTA_ANALYSIS_PATH}

#After succesfful configuration run the following:
/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/STRELKA_MANTA/MANTA_WORKDIR_${names[${SLURM_ARRAY_TASK_ID}]}/runWorkflow.py -j 20
