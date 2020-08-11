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

cd /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/BC_TargetSeq_Aug2020
#ls *.bam > all_bam_files_FL

samples=all_bam_files_FL
names=($(cat $samples))
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

gtf_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v19.annotation.gtf
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/human_g1k_v37.decompressed.fasta
var_folder=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MUTECT2
out_folder=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MUTECT2/pileup_summaries
#targets_interval_list=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/picard_tools_targets_input.bed
targets_interval_list=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/picard_tools_amps_input.bed
ints=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/${sample}_targets.interval_list

tum=($(samtools view -H ${names[${SLURM_ARRAY_TASK_ID}]} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
echo "${tum}"
export tum

gatk GetPileupSummaries \
  -I ${sample} \
  -V /cluster/projects/kridelgroup/RAP_ANALYSIS/af-only-gnomad.raw.sites.b37.vcf.gz \
  -L $ints \
  -O $out_folder/${tum}_getpileupsummaries.table

gatk CalculateContamination \
  -I $out_folder/${tum}_getpileupsummaries.table \
  -O $out_folder/${tum}_calculatecontamination.table

gatk FilterMutectCalls -R $fasta_file \
  -V $var_folder/${tum}.vcf.gz -O $var_folder/${tum}_filtered.vcf.gz
