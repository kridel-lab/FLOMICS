#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 2
#SBATCH --mem=31440M
#SBATCH -t 01-10:00 # Runtime in D-HH:MM
#SBATCH -J Picard_RNAQC
#SBATCH --array=0-364 # job array index - number of jobs = numb of unique samples in file_path folder; e.g. 6samples then --array=0-5

work_path=$SLURM_SUBMIT_DIR
cd ${work_path}

module load picard

if [ ! -d  $work_path/picard_RnaMetrics ]
then mkdir -p $work_path/picard_RnaMetrics
fi

files_location=${work_path}/input_bams
out_path=${work_path}/picard_RnaMetrics
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
java -jar $picard_dir/picard.jar CollectRnaSeqMetrics I=${input_bam} O=$out_path/${prefix} STRAND_SPECIFICITY=NONE REFERENCE_SEQUENCE=/cluster/projects/kridelgroup/resources/Gencode_human/hg37/STAR_genome_Files/GRCh37.primary_assembly.genome.fa REF_FLAT=/cluster/projects/kridelgroup/resources/Gencode_human/hg37/supp/refFlat.txt && echo "$prefix picard rnaseq done"

fi


