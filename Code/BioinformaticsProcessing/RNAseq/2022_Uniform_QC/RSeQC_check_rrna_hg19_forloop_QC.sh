#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 4
#SBATCH --mem=61440M
#SBATCH -t 03-10:00 # Runtime in D-HH:MM
#SBATCH -J RSeQC_check_rrna

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

for i in `cat ${samples}`
do
if [[ "$i" =~ ${files_location}/(.*)Aligned.sortedByCoord.out.bam ]]
 then
 prefix="${BASH_REMATCH[1]}"

if [ ! -d  $work_path/out_path/${prefix} ]
then mkdir -p $work_path/out_path/${prefix}
fi
#htseq-count -f bam -r pos -m union -s reverse ${i} ${genomic_ref}/gencode.v37lift37.annotation.gtf >$prefix  && echo "$prefix htSeq done"
#qualimap rnaseq -bam ${i} -gtf ${genomic_ref}/gencode.v37lift37.annotation.gtf -outdir $work_path/out_path/${prefix} --java-mem-size=41440M
##BAM files should be sorted and indexed
#samtools sort -n -O bam -o ${out_path}/$prefix.sorted.bam ${i}
samtools index ${i}
/cluster/tools/software/centos7/RSeQC/3.0.1/bin/split_bam.py -i ${i} -r /cluster/projects/kridelgroup/resources/downloads/hg19_rRNA.bed -o ${out_path}/${prefix}  && echo "$prefix RSeQC done"

fi
done


