#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=60440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J index

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar
module load bam-readcount
module load gcc
module load STAR
module load rsem
module load perl


STAR --runMode genomeGenerate --genomeDir /cluster/projects/kridelgroup/FLOMICS/genome_files/ \
--genomeFastaFiles /cluster/projects/kridelgroup/FLOMICS/genome_files/ucsc.hg19.fasta --sjdbGTFfile /cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v19.annotation.gtf --runThreadN 10
