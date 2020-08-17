#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=30440M
#SBATCH -t 2-00:00 # Runtime in D-HH:MM

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar
module load bam-readcount
module load kallisto

#pwd
cd /cluster/projects/kridelgroup/genome_files

#building an index
fasta_file=Homo_sapiens.GRCh37.cdna.all.fa.gz
kallisto index -i transcripts.idx $fasta_file
