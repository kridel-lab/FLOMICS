#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 2
#SBATCH --mem=61440M
#SBATCH -t 00-6:00 # Runtime in D-HH:MM
#SBATCH -J trimmomatic
#SBATCH --array=0-n # job array index - number of jobs = numb of unique samples in file_path folder; e.g. 6samples then --array=0-5


##activate the conda env before submitting the jobs
#conda activate reads_QC_fastp_trimmomatic

#creaate the fq1 file
#ls $PWD/rawfq/*_R1.fastq.gz>fq1_files_list.txt
work_path=$SLURM_SUBMIT_DIR
cd ${work_path}


if [ ! -d "${work_path}/cleanfq" ]; then
        mkdir -p "${work_path}/cleanfq"
fi

files_location=${work_path}/rawfq
out_path=${work_path}/cleanfq

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

${miniconda}/envs/reads_QC_fastp_trimmomatic/bin/trimmomatic PE -threads 3 -phred33  $input_fq1 ${files_location}/${prefix}_R2.fastq.gz  ${out_path}/${prefix}.r1_paired.fq.gz ${out_path}/${prefix}.r1_unpaired.fq.gz \
${out_path}/${prefix}.r2_paired.fq.gz ${out_path}/${prefix}.r2_unpaired.fq.gz ILLUMINACLIP:${miniconda}/envs/reads_QC_fastp_trimmomatic/share/trimmomatic-0.39-2/adapters/all_adapter.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:5:15 MINLEN:50 && echo "** ${prefix} trimmomatic done **"
fastqc ${out_path}/${prefix}.r1_paired.fq.gz -o ${out_path} && echo "${prefix} fastqc done"
fastqc ${out_path}/${prefix}.r2_paired.fq.gz -o ${out_path} && echo "${prefix} fastqc done"
fi
