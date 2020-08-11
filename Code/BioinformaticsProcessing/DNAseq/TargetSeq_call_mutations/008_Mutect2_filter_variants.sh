#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2
#SBATCH --array=0-130 # job array index

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar

cd /cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MUTECT2
ls *.vcf.gz > all_vcf_files_FL

samples=all_vcf_files_FL
names=($(cat $samples))
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

gtf_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v19.annotation.gtf
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/human_g1k_v37.decompressed.fasta
out_folder=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MUTECT2
#targets_interval_list=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/picard_tools_targets_input.bed
targets_interval_list=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/picard_tools_amps_input.bed
ints=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/${sample}_targets.interval_list

gatk Mutect2 \
-R $fasta_file \
-I ${names[${SLURM_ARRAY_TASK_ID}]} \
-tumor ${tum} \
-L $ints \
-O $out_folder/${tum}.vcf.gz \
--germline-resource /cluster/projects/kridelgroup/RAP_ANALYSIS/af-only-gnomad.raw.sites.b37.vcf.gz

gatk FilterMutectCalls -R $fasta_file \
-V $sample -O $sample_filtered.vcf.gz
