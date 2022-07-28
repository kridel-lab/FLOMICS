#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 4
#SBATCH --mem=61440M
#SBATCH -t 03-10:00 # Runtime in D-HH:MM
#SBATCH -J RNA-seq_QC

work_path=$SLURM_SUBMIT_DIR
cd ${work_path}

module load qualimap/2.2

if [ ! -d  $work_path/out_path ]
then mkdir -p $work_path/out_path
fi

files_location=${work_path}/input_bams
out_path=${work_path}/out_path
#cd $out_path

genomic_ref="/cluster/projects/kridelgroup/resources/Gencode_human/hg37/STAR_genome_Files"
ls ${files_location}/*.bam > ${work_path}/bam_list.txt
samples="${work_path}/bam_list.txt"

for i in `cat ${samples}`
do
if [[ "$i" =~ ${files_location}/(.*)Aligned.sortedByCoord.out.bam ]]
 then
 prefix="${BASH_REMATCH[1]}"

if [ ! -d  $work_path/out_path/${prefix} ]
then mkdir -p $work_path/out_path/${prefix}
fi
#htseq-count -f bam -r pos -m union -s reverse ${i} ${genomic_ref}/gencode.v37lift37.annotation.gtf >$prefix  && echo "$prefix htSeq done"
qualimap rnaseq -bam ${i} -gtf ${genomic_ref}/gencode.v37lift37.annotation.gtf -outdir $work_path/out_path/${prefix} --java-mem-size=41440M
qualimap bamqc -bam ${i} -gtf ${genomic_ref}/gencode.v37lift37.annotation.gtf -c -outdir $work_path/out_path/${prefix} --java-mem-size=41440M

fi
done


