#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
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

cd /cluster/projects/kridelgroup/GSC-1741

samples=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/target_seq_samples_bam_locations.txt
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
out_folder=/cluster/projects/kridelgroup/FLOMICS/DATA/BC_TargetSeq_Calls
#coordinates of probes (coding and non-coding) used as provided by IDT and Robert
amplicon_interval_list=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/picard_tools_amps_input.bed
#coordinates of targets (coding and non-coding) used as provided by IDT and Robert
targets_interval_list=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/picard_tools_targets_input.bed

#first create symbolic link for bam file so that can actually modify it and index it
cd /cluster/projects/kridelgroup
ln -s /cluster/projects/kridelgroup/GSC-1741/$sample /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/BC_TargetSeq_Aug2020/$name
cd /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/BC_TargetSeq_Aug2020
sample=$name #linked bam file that we can use

#index BAM file
samtools index $sample

#make sample specific interval list files required for picard as it doesn't work
#with bed files
#SD is dictionary which in this case is just the BAM file

gatk BedToIntervalList \
      -I $amplicon_interval_list \
      -O /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/${sample}_amplicon.interval_list \
      -SD $sample

gatk BedToIntervalList \
      -I $targets_interval_list \
      -O /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/${sample}_targets.interval_list \
      -SD $sample

#save new interval lists as variables which will be used as input for final picard function
amps=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/${sample}_amplicon.interval_list
ints=/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/${sample}_targets.interval_list

gatk CollectTargetedPcrMetrics \
       -I $sample \
       -O ${out_folder}/${name}.output_pcr_metrics.txt \
       -R $fasta_file \
       --PER_TARGET_COVERAGE ${out_folder}/${name}.per_target_coverage.txt \
       --AMPLICON_INTERVALS $amps \
       --TARGET_INTERVALS $ints
