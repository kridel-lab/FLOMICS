#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p veryhimem
#SBATCH -c 2
#SBATCH --mem=991440M
#SBATCH -t 03-10:00 # Runtime in D-HH:MM
#SBATCH -J htseq

module load HTSeq/0.11.0
module load samtools/1.10

work_path=$SLURM_SUBMIT_DIR
cd ${work_path}

if [ ! -d  $work_path/count_out_pending ]
then mkdir -p $work_path/count_out_pending
fi

files_location=${work_path}/input_bams
out_path=${work_path}/count_out_pending
cd $out_path

genomic_ref="/cluster/projects/kridelgroup/resources/Gencode_human/hg37/STAR_genome_Files"
#ls ${file_path}/*.r1_paired.fq.gz  > ${work_path}/fq1_files_list.txt
#file=${files_location}/*.r1_paired.fq.gz

samples="${work_path}/unfinished_bam_list.txt"

for i in `cat ${samples}`
do
if [[ "$i" =~ ${files_location}/(.*)Aligned.sortedByCoord.out.bam ]]
 then
 prefix="${BASH_REMATCH[1]}"
htseq-count -f bam -r pos -m union -s reverse ${i} ${genomic_ref}/gencode.v37lift37.annotation.gtf >$prefix  && echo "$prefix htSeq done"
#the "-r pos" is necessary here if your bam/sam file is not sorted by name the read mates in two consecutive lines, otherwise you will have 2x read count per gene
#another option is sorting the bam/sam files, then you can use "-r name' option.
fi
done
