#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=31440M
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
module load vt

cd /cluster/projects/kridelgroup/FLOMICS/ANALYSIS/PLATYPUS

#ls *.vcf > all_vcf_files

samples=all_vcf_files
names=($(cat $samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
sample=${MYVAR%%.platypus*}
echo $sample

out_folder=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/PLATYPUS/annovar
fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/human_g1k_v37.decompressed.fasta

#normalize file
input_vcf=${names[${SLURM_ARRAY_TASK_ID}]}

#keep only passed variants
#first keep only variants that passed baseline encoded variants already done
awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $input_vcf > ${sample}.clean.vcf

#normalize variants, send to standard out and remove duplicates.
vt normalize ${sample}.clean.vcf -r $fasta_file -o ${sample}.normalized.vcf

#RUN ANNOVAR
anno_input=${sample}.normalized.vcf

table_annovar.pl --buildver hg19 $anno_input /cluster/tools/software/annovar/humandb \
--protocol ensGene,gnomad211_genome,cosmic68,avsnp142 --operation g,f,f,f \
--outfile ${out_folder}/${sample}_annovar.vcf --vcfinput

mv ${sample}_annovar.vcf.gz.hg19_multianno.vcf annovar
rm *${sample}_annovar*
