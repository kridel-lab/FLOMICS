#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 4
#SBATCH --mem=61440M
#SBATCH -t 01-10:00 # Runtime in D-HH:MM
#SBATCH -J RSeQC
#SBATCH --array=0-364 # job array index - number of jobs = numb of unique samples in file_path folder; e.g. 6samples then --array=0-5

work_path=$SLURM_SUBMIT_DIR
cd ${work_path}

module load samtools/1.10
module load RSeQC/3.0.1

if [ ! -d  $work_path/RSeQC ]
then mkdir -p $work_path/RSeQC
fi

files_location=${work_path}/input_bams
out_path=${work_path}/RSeQC
#cd $out_path

genomic_ref="/cluster/projects/kridelgroup/resources/Gencode_human/hg37/STAR_genome_Files"
#ls ${files_location}/*.bam > ${work_path}/bam_list.txt
samples="${work_path}/bam_list.txt"
prefixs=($(cat $samples))
#echo $prefixs
input_bam=${prefixs[${SLURM_ARRAY_TASK_ID}]}

if [[ "$input_bam" =~ ${files_location}/(.*)Aligned.sortedByCoord.out.bam ]]
 then
 echo "$input_bam"
 prefix="${BASH_REMATCH[1]}"
 echo "$prefix"
#qualimap rnaseq -bam ${input_bam} -gtf ${genomic_ref}/gencode.v37lift37.annotation.gtf -outdir $work_path/out_path/${prefix} --java-mem-size=41440M && echo "$prefix rnaseq done"
samtools index ${input_bam}
/cluster/tools/software/centos7/RSeQC/3.0.1/bin/split_bam.py -i ${input_bam} -r /cluster/projects/kridelgroup/resources/downloads/hg19_rRNA.bed -o ${out_path}/${prefix} >${prefix}.rrna.txt && echo "$prefix RSeQC done"

fi


