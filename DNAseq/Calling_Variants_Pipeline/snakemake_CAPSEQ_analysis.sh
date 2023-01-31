#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=20000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J ppE2408-CAPSEQ
#SBATCH -o /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_OUTPUT/E2408-Tumor-Only/run001/slurm-outs/%x-%j.out #redirect job output (both stdout and stderr) to a file called “<job name>-<job id>.out”

#load required modules
module load picard/2.10.9
module load bwa/0.7.15
module load python3/3.7.2
module load samtools
module load snakemake/5.20.1

cd /cluster/home/vshelton/CapSeq_pipeline/TargetedDNAseq_pipeline

snakemake -s /cluster/home/vshelton/CapSeq_pipeline/TargetedDNAseq_pipeline/CAPSEQ_Pipeline_BC.smk \
    -j157 \
    --latency-wait 432000 \
    --cluster-config /cluster/home/vshelton/CapSeq_pipeline/TargetedDNAseq_pipeline/cluster.yaml \
    --cluster-sync "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -o {cluster.out}" \
    --rerun-incomplete \
    --nolock
