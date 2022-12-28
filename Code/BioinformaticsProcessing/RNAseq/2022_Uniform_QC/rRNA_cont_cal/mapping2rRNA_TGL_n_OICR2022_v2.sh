#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 2
#SBATCH --mem=61440M
#SBATCH -t 00-6:00 # Runtime in D-HH:MM
#SBATCH -J map_rrna_TGL
#SBATCH --array=0-154 

work_path=$SLURM_SUBMIT_DIR
cd ${work_path}


module load bwa/0.7.15
module load samtools/1.9
files_location=${work_path}/clean_TGL_n_OICR2022
out_path=${work_path}/TGL_n_OICR2022_rrna_mapped_v2
ref_bundle=${work_path}/rrna_ref

samples=${work_path}/fq1_files_list_TGL.txt
echo $sample
prefixs=($(cat $samples))
echo $prefixs
input_fq1=${prefixs[${SLURM_ARRAY_TASK_ID}]}

if [[ "$input_fq1" =~ ${files_location}/(.*).r1_paired.fq.gz ]]
 then
        echo "$input_fq1"
        prefix="${BASH_REMATCH[1]}"
        echo "$prefix"
cat "${files_location}/${prefix}.r1_paired.fq.gz" "${files_location}/${prefix}.r2_paired.fq.gz" >"${out_path}/${prefix}_concat.fq.gz"  && echo "** $prefix concat done **"
bwa mem -M -p -t 4 -R "@RG\tID:$prefix\tPL:ILLUMINA\tSM:$prefix" "${ref_bundle}/human_all_rRNA_mod.fasta" "${out_path}/${prefix}_concat.fq.gz" | samtools view -bS - > "${out_path}/${prefix}.bam" && echo "** $prefix bwa mem done **"  ## -p here is for smart pairing, get away with the error message "[mem_sam_pe] paired reads have different names:" for some samples
samtools flagstat "${out_path}/${prefix}.bam" >${out_path}/${prefix}_flagstat.tsv && echo "** $prefix flagstat done **"

fi
