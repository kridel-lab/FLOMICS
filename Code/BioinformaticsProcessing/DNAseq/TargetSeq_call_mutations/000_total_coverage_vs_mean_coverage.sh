#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -J picard
#SBATCH --array=0-130 # job array index - number of jobs = numb of unique samples with top up runs

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
module load samtools

cd /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/RAW_DATA_UPLOADS/GSC-1741

#ls */*/*/*/*.bam > /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/FL_Aug2020_library_list.txt

samples=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/FL_Aug2020_library_list.txt
names=($(cat $samples))

#get sample file
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

#extract just the bam file name which includes sample name
name=${sample#*hg19a/}
echo $name

#prepare files and folders required for analysis and storage of downstream files
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/human_g1k_v37.fasta
#where to store output (summary of coverage and pcr metrics from each bam file)
out_folder=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/SYMLINKS_to_BAMs/GSC-1741-BCAug2020

#first create symbolic link for bam file so that can actually modify it and index it
cd /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/RAW_DATA_UPLOADS/
ln -s /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/RAW_DATA_UPLOADS/GSC-1741/${sample} ${out_folder}/${name}

cd ${out_folder}
sample=$name #linked bam file that we can use

#index BAM file
samtools index $sample

#get tot number of reads
samtools view -c $sample > /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/PROCESSING/library_sizes/${name}_num_reads_in_bam.txt
