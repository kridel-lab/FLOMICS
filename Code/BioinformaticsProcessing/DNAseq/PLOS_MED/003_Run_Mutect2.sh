#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2
#SBATCH --array=0-483 # job array index

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar

cd /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/EGA_data_processing
#ls *.bam > all_bam_files_PLOS

samples=all_bam_files_PLOS
names=($(cat $samples))
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

gtf_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v19.annotation.gtf
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/human_g1k_v37.decompressed.fasta
out_folder=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/EGA_calls

ints=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/EGA_data_processing/interval_lists/${sample}_amplicon.interval_list

tum=($(samtools view -H ${names[${SLURM_ARRAY_TASK_ID}]} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
echo "${tum}"
export tum

gatk Mutect2 \
-R $fasta_file \
-I ${names[${SLURM_ARRAY_TASK_ID}]} \
-tumor ${tum} \
-L $ints \
-O $out_folder/${tum}.vcf.gz \
--germline-resource /cluster/projects/kridelgroup/RAP_ANALYSIS/af-only-gnomad.raw.sites.b37.vcf.gz