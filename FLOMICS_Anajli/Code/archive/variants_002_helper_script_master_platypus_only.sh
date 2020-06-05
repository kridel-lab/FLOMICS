#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=25440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J capture_merge
#date stamp = November 21, 2019 Karin Isaev 

#load all potentially required tools 
module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load bedtools
module load tabix
module load vt

index=$1 #this refers to the VCF file
echo $index

#[1]

#first we will enter the folder with name of patient 
cd /cluster/projects/kridelgroup/FLOMICS/TargetedDNAseq/DATA-157/variants/
pwd

#[2]

#Establish outdirectory for where the annotated VCFs will be saved to 
out_dir=/cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder

#[3]

#now normalize VCF coordinates 
input_vcf=$index
MYVAR=$index
tum_loc=${MYVAR%%/*}

#test "${INDEL#*$tum_loc}" != "$INDEL" && echo "$tum_loc found in $INDEL"
tum_name=${MYVAR##*/}
patient=${MYVAR%/*/*}
algo_type="$(cut -d'/' -f2 <<<$MYVAR)"
echo $patient
echo $algo_type

#[4]

#normalize variants, send to standard out and remove duplicates.
ref=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #fasta file for hg19 from previous work 
vt normalize $input_vcf -r $ref -o $out_dir/normalized_${tum_name}_${algo_type}.vcf

cd $out_dir
 
#[5]

#set up annovar input file whihc is the output file from the normalization command above 
anno_input=normalized_${tum_name}_${algo_type}.vcf

#[6]

#re-name AF column in VCF file so it's not confused later with the annovar AF tag 
#bgzip -c $anno_input > $anno_input.gz
#tabix -p vcf $anno_input.gz
#bcftools annotate -a $anno_input.gz -c INFO/AF_SEQ:=INFO/AF $anno_input.gz > retagged_$anno_input
#bgzip -c retagged_$anno_input.gz > retagged_$anno_input.gz

#[7]

#Run annovar 
table_annovar.pl --buildver hg19 $anno_input /cluster/tools/software/annovar/humandb --protocol ensGene,gnomad211_genome,cosmic68,avsnp142 --operation g,f,f,f --outfile ${tum_name}_${algo_type}.vcf --vcfinput

#[7] 

#save final VCF file into clean new directory and delete all other files produced since we don't need them 

#remove original AF column
vt rminfo ${tum_name}_${algo_type}.vcf.hg19_multianno.vcf -t AF -o ${tum_name}_${algo_type}.vcf.hg19_multianno.vcf
#mkdir final_VCFs
mv ${tum_name}_${algo_type}.vcf.hg19_multianno.vcf final_VCFs
#rm *${patient}*

#DONE 
