#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 2
#SBATCH --mem=61440M
#SBATCH -t 00-6:00 # Runtime in D-HH:MM
#SBATCH -J fastp
#SBATCH --array=0-209 # job array index - number of jobs = numb of unique samples in file_path folder; e.g. 6samples then --array=0-5

module load fastp/0.23.1
module load fastqc/0.11.5

work_path=$SLURM_SUBMIT_DIR
cd ${work_path}

if [ ! -d "${work_path}/clean_fq" ]; then
        mkdir -p "${work_path}/clean_fq"
fi

files_location=${work_path}/rawdata_symlinks
out_path=${work_path}/clean_fq

#ls ${file_path}/*_R1.fastq.gz  > ${work_path}/fq1_files_list.txt
#file=${files_location}/*.r1_paired.fq.gz

samples=${work_path}/fq1_files_list.txt
echo $sample
prefixs=($(cat $samples))
echo $prefixs
input_fq1=${prefixs[${SLURM_ARRAY_TASK_ID}]}

if [[ "$input_fq1" =~ ${files_location}/(.*)_R1.fastq.gz ]]
 then
	echo "$input_fq1"
	prefix="${BASH_REMATCH[1]}"
	echo "$prefix"
fastp --in1 $input_fq1 --in2 ${files_location}/${prefix}_R2.fastq.gz --out1 ${out_path}/${prefix}.r1_paired.fq.gz --out2 ${out_path}/${prefix}.r2_paired.fq.gz --unpaired1 ${out_path}/${prefix}.r1_unpaired.fq.gz --unpaired2 ${out_path}/${prefix}.r2_unpaired.fq.gz -h "${out_path}/${prefix}.html" && echo "$num fastp reads qc done"
fastqc ${out_path}/${prefix}.r1_paired.fq.gz -o ${out_path} && echo "${prefix} fastqc done"
fastqc ${out_path}/${prefix}.r2_paired.fq.gz -o ${out_path} && echo "${prefix} fastqc done"
fi
