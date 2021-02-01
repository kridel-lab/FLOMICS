#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J picard
#SBATCH --array=0-483 # job array index - number of jobs = numb of unique samples with top up runs

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

cd /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/EGAD00001002898

samples=/cluster/projects/kridelgroup/FLOMICS/DATA/EGA_samples_list.txt
names=($(cat $samples))

processing_folder=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/EGA_data_processing

#get sample file
sample=${names[${SLURM_ARRAY_TASK_ID}]}
cd $sample

#extract just the bam file name which includes sample name
name="$(ls  *bam)"
echo $name

#prepare files and folders required for analysis and storage of downstream files
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/human_g1k_v37.fasta

#where to store output (summary of coverage and pcr metrics from each bam file)
out_folder=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/EGA_calls

probe_cords=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/PLOS_MED_probes_input_picard.bed

#first create symbolic link for bam file so that can actually modify it and index it
ln -s /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/EGAD00001002898/${sample}/${name} ${processing_folder}/${name}

cd ${processing_folder}
sample=$name #linked bam file that we can use

#index BAM file
samtools index $sample

#make sample specific interval list files required for picard as it doesn't work
#with bed files
#SD is dictionary which in this case is just the BAM file

gatk BedToIntervalList \
      -I $probe_cords \
      -O /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/EGA_data_processing/interval_lists/${sample}_amplicon.interval_list \
      -SD $sample

#save new interval lists as variables which will be used as input for final picard function
ints=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/EGA_data_processing/interval_lists/${sample}_amplicon.interval_list
