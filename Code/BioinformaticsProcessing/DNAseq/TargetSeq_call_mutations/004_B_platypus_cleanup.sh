#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J Platypus
#SBATCH --array=0-130 # job array index

module load annovar
module load bam-readcount
module load gcc
module load STAR
module load rsem
module load perl
module load python/2.7
module load gatk
module load tabix
module load vcftools

cd /cluster/projects/kridelgroup/FLOMICS/ANALYSIS/PLATYPUS

ls *.vcf > all_vcf_files

samples=all_vcf_files
names=($(cat $samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
sample=${MYVAR%%.vcf*}
echo $sample



samples=all_bam_files_FL
names=($(cat $samples))
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

gtf_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v19.annotation.gtf
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/human_g1k_v37.decompressed.fasta
out_folder=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/PLATYPUS
#targets_interval_list=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/picard_tools_targets_input.bed
targets_interval_list=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/picard_tools_amps_input.bed
ints=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/${sample}_targets.interval_list

tum=($(samtools view -H ${names[${SLURM_ARRAY_TASK_ID}]} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
echo "${tum}"
export tum

python /cluster/tools/software/platypus/0.8.1/Platypus.py callVariants \
--bamFiles ${names[${SLURM_ARRAY_TASK_ID}]} --refFile $fasta_file \
-o ${out_folder}/${tum}.platypus.vcf --filterDuplicates=0
