#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 4
#SBATCH --mem=61440M
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J star
#SBATCH --array=0-n # job array index - number of jobs = numb of unique samples in file_path folder; e.g. 6samples then --array=0-5

module load STAR/2.7.9a

work_path=$SLURM_SUBMIT_DIR
cd ${work_path}

if [ ! -d  $work_path/STAR_out ]
then mkdir -p $work_path/STAR_out
fi

files_location=${work_path}/cleanfq
out_path=${work_path}/STAR_out

#ls ${file_path}/*.r1_paired.fq.gz  > ${work_path}/fq1_files_list.txt

samples=${work_path}/fq1_files_clean_list.txt
echo $sample
prefixs=($(cat $samples))
echo $prefixs
input_fq1=${prefixs[${SLURM_ARRAY_TASK_ID}]}

if [[ "$input_fq1" =~ ${files_location}/(.*).r1_paired.fq.gz ]]
 then
	echo "$input_fq1"
	prefix="${BASH_REMATCH[1]}"
	echo "$prefix"
  STAR --genomeDir ${resources}/Gencode_human/hg37/STAR_genome_Files/ \
  --readFilesIn ${input_fq1} ${files_location}/${prefix}.r2_paired.fq.gz --outSAMtype BAM SortedByCoordinate \
  --readFilesCommand zcat \
  --outSAMunmapped Within --outFileNamePrefix ${out_path}/${prefix} --quantMode GeneCounts \
  --twopassMode Basic
fi
