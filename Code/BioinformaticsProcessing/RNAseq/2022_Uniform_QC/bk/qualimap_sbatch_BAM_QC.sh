#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 4
#SBATCH --mem=61440M
#SBATCH -t 01-10:00 # Runtime in D-HH:MM
#SBATCH -J RNA-seq_QC
#SBATCH --array=0-364 # job array index - number of jobs = numb of unique samples in file_path folder; e.g. 6samples then --array=0-5

work_path=$SLURM_SUBMIT_DIR
cd ${work_path}

module load qualimap/2.2

if [ ! -d  $work_path/bamQC_out_path ]
then mkdir -p $work_path/bamQC_out_path
fi

files_location=${work_path}/input_bams
out_path=${work_path}/bamQC_out_path
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
 if [ ! -d  $work_path/bamQC_out_path/${prefix} ]
  then mkdir -p $work_path/bamQC_out_path/${prefix}
 fi
 qualimap bamqc -bam ${input_bam} -gff ${genomic_ref}/gencode.v37lift37.annotation.gtf -outdir $out_path/${prefix} --java-mem-size=41440M && echo "$prefix bamseq done"
fi


