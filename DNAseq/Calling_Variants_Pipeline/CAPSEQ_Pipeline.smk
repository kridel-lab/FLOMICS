#Victoria Shelton  - Last update: August 2021
#This CAPSEQ pipeline is intended to run upon the BCAugust2020, BCJuly2021 and E2408 Tumor BAM Files.
#The Analysis performed: coverage summary, Manta SV calling, CNVkit CNV calling,
    #Mutect2 SNV calling & Annovar SNV annotation

import pandas as pd
from pathlib import Path
import subprocess
from os.path import join

configfile: "/cluster/home/vshelton/CapSeq_pipeline/TargetedDNAseq_pipeline/cluster.yaml"
configfile: "/cluster/home/vshelton/CapSeq_pipeline/TargetedDNAseq_pipeline/config_2.yaml"

### Globals ---------------------------------------------------------------------
# A Snakemake regular expression matching BAM_symlinks files.
SAMPLES, = glob_wildcards(join(config["BAMsymlinks"], "{sample}.bam"))
print(SAMPLES)

wildcard_constraints:
    sample = "\w+"

### Rules -----------------------------------------------------------------------
# Pipeline output files.
rule all:
    input:
        expand(join(config["outputDIR"], "MANTA/MANTA_WORKDIR_nointervals_{sample}/runWorkflow.py"), sample=SAMPLES),
        expand(join(config["outputDIR"], "CNVkit/CNVkit_WORKDIR_{sample}/CNVkit_amps_input_{sample}.bed"), sample=SAMPLES),
        expand(join(config["outputDIR"], "MUTECT2/{sample}.vcf.gz"), sample=SAMPLES),
        expand(join(config["outputDIR"], "MUTECT2/pileup_summaries/{sample}_getpileupsummaries.table"), sample=SAMPLES),
        expand(join(config["outputDIR"], "MUTECT2/pileup_summaries/{sample}_calculatecontamination.table"), sample=SAMPLES),
        expand(join(config["outputDIR"], "MUTECT2/{sample}_normalized.vcf.gz"), sample=SAMPLES),
        expand(join(config["outputDIR"], "MUTECT2/{sample}_decomposed.vcf"), sample=SAMPLES),
        expand(join(config["outputDIR"], "ANNOVAR/{sample}_annovar.vcf.gz.hg19_multianno.vcf"), sample=SAMPLES),

#Collet target PCR metrics from Picard
rule coverage_summary:
    input:
        sample_bam = join(config["BAMsymlinks"], "{sample}.bam"),
        fasta_file = config["hg19Fasta"],
        amps = config["amplicon_bed"],
        targets = config["target_bed"],
    params:
        sample_ID = "{sample}.bam",
    output:
        amplicon_interval_list = join(config["outputDIR"], "amplicon_interval_list/{sample}_amplicon.interval_list"),
        targets_interval_list = join(config["outputDIR"], "target_interval_list/{sample}_targets.interval_list"),
        pcr_metrics = join(config["outputDIR"], "METRICS/{sample}.output_pcr_metrics.txt"),
        per_target_coverage = join(config["outputDIR"], "COVERAGE/{sample}.per_target_coverage.txt"),
    wildcard_constraints:
        sample = "\w+"
    message:
        "Indexing BAM files & collecting the probe target PCR metrics (coverage) from Targeted DNA-seq BAM files"
    shell:
            """
                #load modules
                module load annovar/20180416
                module load bam-readcount/0.7.4
                module load gcc/6.2.0
                module load STAR/2.7.3a
                module load rsem/1.3.0
                module load perl/5.30.0
                module load python/2.7
                module load gatk/4.0.5.1
                module load tabix/0.2.6
                module load vcftools/0.1.15
                module load samtools/1.10

                #index BAM file
                samtools index {input.sample_bam}

                gatk BedToIntervalList \
                        -I {input.amps} \
                        -O {output.amplicon_interval_list} \
                        -SD {input.sample_bam}

                gatk BedToIntervalList \
                        -I {input.targets} \
                        -O {output.targets_interval_list} \
                        -SD {input.sample_bam}

                gatk CollectTargetedPcrMetrics \
                        -I {input.sample_bam} \
                        -O {output.pcr_metrics} \
                        -R {input.fasta_file} \
                        --PER_TARGET_COVERAGE {output.per_target_coverage} \
                        --AMPLICON_INTERVALS {output.amplicon_interval_list} \
                        --TARGET_INTERVALS {output.targets_interval_list}
            """

#Run Manta to obtain SVs and indels
rule run_Manta:
    input:
        sample_bam = rules.coverage_summary.input.sample_bam,
        gtf_file = config["gencodehg19GTF"],
        decom_fasta_file = config["Decomhg19Fasta"],
        amps = config["amplicon_bed"],
        MantaPath = config["MantaPath"],
        MantaDir = config["MantaDir"],
        rule1 = rules.coverage_summary.output.targets_interval_list,
    output:
        manta_outpath = directory(join(config["outputDIR"], "MANTA/MANTA_WORKDIR_nointervals_{sample}/")),
        bed_index = join(config["outputDIR"], "MANTA/AMP_BED_ZIP/amp_bed_{sample}.gz"),
        manta_out1 = join(config["outputDIR"], "MANTA/MANTA_WORKDIR_nointervals_{sample}/runWorkflow.py"),
    wildcard_constraints:
        sample = "\w+"
    message:
        "Configuring and Running Manta."
    shell:
            """
            #load modules
            module load strelka/2.9.10
            module load python/3.4.3
            module load manta/1.6.0
            module load tabix/0.2.6

            #index bed file
            sort -k 1,1 -k 2,2n -k 3,3n {input.amps} | bgzip -c > {output.bed_index}
            tabix -pbed -f {output.bed_index}

            #configuration
            {input.MantaPath} \
            --tumorBam {input.sample_bam} \
            --referenceFasta {input.decom_fasta_file} --exome \
            --runDir {output.manta_outpath}

            #After succesfful configuration run the following:
            {output.manta_out1} -j 20

            """

#Run CNVkit to obtain copy number vraiants
rule run_CNVkit:
    input:
        sample_bam = rules.coverage_summary.input.sample_bam,
        gtf_file = config["gencodehg19GTF"],
        decom_fasta_file = config["Decomhg19Fasta"],
        amps = config["amplicon_bed"],
        CNVkitPath = config["CNVkitPath"],
        CNVkitDir = config["CNVkitDir"],
        gene_anno = config["Gene_Anno"],
        rule1 = rules.coverage_summary.output.targets_interval_list,
    output:
        CNVkit_outpath = directory(join(config["outputDIR"], "CNVkit/CNVkit_WORKDIR_{sample}/")),
        bed_index = join(config["outputDIR"], "CNVkit/CNVkit_WORKDIR_{sample}/CNVkit_amps_input_{sample}.bed"),
        CNVkit_flatref = join(config["outputDIR"], "CNVkit/CNVkit_WORKDIR_{sample}/flat_reference.cnn"),
    wildcard_constraints:
        sample = "\w+"
    message:
        "Running CNVkit."
    shell:
            """
            #load modules
            module load strelka/2.9.10
            module load python/3.4.3
            module load R/4.0.0
            module load CNVkit/0.9.3

            #preparing BED file of baited regions for use with CNVkit
            {input.CNVkitPath} \
            target  {input.amps} \
            --annotate {input.gene_anno} \
            -o {output.bed_index}

            #running CNVkit
            {input.CNVkitPath} \
            batch {input.sample_bam} \
            --normal \
            -m amplicon \
            --targets {output.bed_index} \
            --fasta {input.decom_fasta_file} \
            --output-reference {output.CNVkit_flatref} \
            --output-dir {output.CNVkit_outpath}

            """


#Run MuTect2 to obtain SNV calls
rule SNVs_mutect2_tumour_only:
    input:
        sample_bam = rules.coverage_summary.input.sample_bam,
        gtf_file = config["gencodehg19GTF"],
        decom_fasta_file = config["Decomhg19Fasta"],
        amps = config["amplicon_bed"],
        targets_interval_list = rules.coverage_summary.output.targets_interval_list,
        raw_sites = config["gnomadRawSites"],
    output:
        vcf = join(config["outputDIR"], "MUTECT2/{sample}.vcf.gz"),
    wildcard_constraints:
        sample = "\w+"
    message:
        "Calling SNVs with Mutect2 in tumour only mode."
    shell:
            """
                #load modules
                module load java/8
                module load samtools/1.10
                module load python/3.4.3
                module load gatk/4.1.8.1
                module load annovar/20180416

                gatk Mutect2 \
                        -R {input.decom_fasta_file} \
                        -I {input.sample_bam} \
                        -L {input.targets_interval_list} \
                        -O {output.vcf} \
                        --germline-resource {input.raw_sites}


                        """

#Retrieve contamination and filter SNVs called by MuTect2
rule filter_mutect2_SNVs:
    input:
        sample_bam = rules.coverage_summary.input.sample_bam,
        gtf = config["gencodehg19GTF"],
        decom_fasta_file = config["Decomhg19Fasta"],
        amps = config["amplicon_bed"],
        targets_interval_list = rules.coverage_summary.output.targets_interval_list,
        raw_sites = config["gnomadRawSites"],
        vcf = rules.SNVs_mutect2_tumour_only.output.vcf,
    output:
        pile_sum = join(config["outputDIR"], "MUTECT2/pileup_summaries/{sample}_getpileupsummaries.table"),
        cal_con = join(config["outputDIR"], "MUTECT2/pileup_summaries/{sample}_calculatecontamination.table"),
        vcffil = join(config["outputDIR"], "MUTECT2/{sample}_filtered.vcf.gz"),
    wildcard_constraints:
        sample = "\w+"
    message:
        "Retrieving contamination and filtering SNVs called by Mutect2."
    shell:
            """
                #load modules
                module load java/8
                module load samtools/1.10
                module load python/3.4.3
                module load gatk/4.1.8.1
                module load annovar/20180416

                gatk GetPileupSummaries \
                        -I {input.sample_bam} \
                        -V {input.raw_sites} \
                        -L {input.targets_interval_list} \
                        -O {output.pile_sum}

                gatk CalculateContamination \
                        -I {output.pile_sum} \
                        -O {output.cal_con}

                gatk FilterMutectCalls \
                        -R {input.decom_fasta_file} \
                        -V {input.vcf} \
                        -O {output.vcffil}

            """

#normalize VCF files
rule normalize_mutect2_SNVs:
    input:
        vcffil = rules.filter_mutect2_SNVs.output.vcffil,
        decom_fasta_file = config["Decomhg19Fasta"],
    output:
        vcf_norm = join(config["outputDIR"], "MUTECT2/{sample}_normalized.vcf.gz"),
        vcf_decom = join(config["outputDIR"], "MUTECT2/{sample}_decomposed.vcf"),
    wildcard_constraints:
        sample = "\w+"
    message:
        "Normalizing and decomposing filtered VCF files."
    shell:
            """
                #load modules
                module load java/8
                module load samtools/1.10
                module load python/3.4.3
                module load gatk/4.1.8.1
                module load tabix/0.2.6
                module load vt/0.577

                #normalize variants, send to standard out and remove duplicates.
                vt normalize {input.vcffil} -r {input.decom_fasta_file} -o {output.vcf_norm}

                #----------------------
                #split multiallelic variants to biallelic
                vt decompose -s {output.vcf_norm} -o {output.vcf_decom}

            """

#run Annovar to annotate SNVs
rule VCF_annovar_annotation:
    input:
        vcf_decom = rules.normalize_mutect2_SNVs.output.vcf_decom,
        decom_fasta_file = config["Decomhg19Fasta"],
        humandb =config["Annohumandb"],
    params:
        outfile = join(config["outputDIR"], "ANNOVAR/{sample}_annovar.vcf.gz"),
    output:
        anno = join(config["outputDIR"], "ANNOVAR/{sample}_annovar.vcf.gz.hg19_multianno.vcf"),
    wildcard_constraints:
        sample = "\w+",
    message:
        "Annotating decomposed VCF files with Annovar."
    shell:
            """
                #load modules
                module load java/8
                module load samtools/1.10
                module load python/3.4.3
                module load gatk/4.1.8.1
                module load annovar/20180416
                module load tabix/0.2.6
                module load vt/0.577

                #RUN ANNOVAR
                table_annovar.pl --buildver hg19 {input.vcf_decom} {input.humandb} \
                --protocol ensGene,gnomad211_genome,cosmic68,avsnp142 --operation g,f,f,f \
                --outfile {params.outfile} --vcfinput

            """
