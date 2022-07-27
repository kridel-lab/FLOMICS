#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 4
#SBATCH --mem=61440M
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J htseq
#SBATCH --array=0-154 # job array index - number of jobs = numb of unique samples in file_path folder; e.g. 6samples then --array=0-5

module load HTSeq/0.11.0


work_path=$SLURM_SUBMIT_DIR
cd ${work_path}

if [ ! -d  $work_path/count_out ]
then mkdir -p $work_path/count_out
fi

files_location=${work_path}/input_bams
out_path=${work_path}/count_out
cd $out_path

genomic_ref="/cluster/projects/kridelgroup/resources/Gencode_human/hg37/STAR_genome_Files"
#ls ${file_path}/*.r1_paired.fq.gz  > ${work_path}/fq1_files_list.txt
#file=${files_location}/*.r1_paired.fq.gz

samples=${work_path}/bam_files_list.txt
echo $sample
prefixs=($(cat $samples))
echo $prefixs
input_bam=${prefixs[${SLURM_ARRAY_TASK_ID}]}

if [[ "$input_bam" =~ ${files_location}/(.*)Aligned.sortedByCoord.out.bam ]]
 then
	echo "$input_bam"
	prefix="${BASH_REMATCH[1]}"
	echo "$prefix"
htseq-count -f bam -r pos -m union -s reverse ${input_bam} ${genomic_ref}/gencode.v37lift37.annotation.gtf >$prefix  && echo "$prefix htSeq done"

#the "-r pos" is necessary here if your bam/sam file is not sorted by name the read mates in two consecutive lines, otherwise you will have 2x read count per gene
#another option is sorting the bam/sam files, then you can use "-r name' option.

fi
