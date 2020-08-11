#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=31440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J ANNOVAR
#SBATCH --array=0-130 # job array index

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load tabix
module load vt

cd /cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MUTECT2
ls *filtered.vcf.gz > all_vcf_files_FL

samples=all_vcf_files_FL
names=($(cat $samples))
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/human_g1k_v37.decompressed.fasta

#normalize file
input_vcf=$sample

#normalize variants, send to standard out and remove duplicates.
vt normalize $sample -r $fasta_file -o ${sample}.normalized.vcf.gz

#RUN ANNOVAR
anno_input=${sample}.normalized.vcf.gz
out_folder=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/MUTECT2

table_annovar.pl --buildver hg19 ${anno_input} /cluster/tools/software/annovar/humandb \
--protocol ensGene,gnomad211_genome,cosmic68,avsnp142 --operation g,f,f,f \
--outfile ${out_folder}/${sample}_annovar.vcf.gz --vcfinput

mv ${sample}_annovar.vcf.gz.hg19_multianno.vcf annovar
rm *${sample}_annovar*
