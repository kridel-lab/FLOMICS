rule all:
    input:
        "/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/opossum_platypus/{sample}.clean.vcf"

rule samtools_md:
    input:
        "{sample}Aligned.sortedByCoord.out.bam"
    output:
        "{sample}.calmd.bam"
    shell:
        "samtools calmd {input} /cluster/projects/kridelgroup/FLOMICS/genome_files/ucsc.hg19.fasta > {output}"

rule opossum_run:
    input:
        "{sample}.calmd.bam"
    output:
        "{sample}.opossum.bam"
    shell:
        """
        module load python/3.4.3
        python /cluster/home/kisaev/Opossum/Opossum.py --BamFile={input} --OutFile={output} --SoftClipsExist True
        """

rule opossum_index:
    input:
        "{sample}.opossum.bam"
    shell:
        "samtools index {input}"

rule check_bam_index:
    input:
        "{sample}.opossum.bam.bai"
    shell:
        "ls -lt {input}"

rule platypus:
    input:
        "{sample}.opossum.bam"
    output:
        "/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/opossum_platypus/{sample}.vcf"
    shell:
        """
        module load platypus/0.8.1
        python /cluster/tools/software/platypus/0.8.1/Platypus.py callVariants --bamFiles {input} --refFile /cluster/projects/kridelgroup/FLOMICS/genome_files/ucsc.hg19.fasta  --filterDuplicates 0 --minMapQual 0 --minFlank 0 --maxReadLength 500 --minGoodQualBases 10 --minBaseQual 20 -o {output}
        """

#rule clean_vcf:
#    input:
#        "/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/opossum_platypus/{sample}.vcf"
#    output:
#        "/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/opossum_platypus/{sample}.clean.vcf"
#    shell:
#        """
#        awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' {input} > {output}
#        """

#rule annovar:#
#    input:#
#        "/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/opossum_platypus/{sample}.clean.vcf"
#    shell:
#        "table_annovar.pl --buildver hg19 {input} /cluster/tools/software/annovar/humandb --protocol ensGene,gnomad211_genome,cosmic68,avsnp142 --operation g,f,f,f --vcfinput"

#rule filter_annovar:
    #input:
#        "{sample}.clean.vcf.hg19_multianno.vcf"
    #output:
#        "{sample}.clean.vcf.hg19_multianno.filtered.vcf"

#add these rules below
#table_annovar.pl --buildver hg19 LY_FL_013_T1.clean.vcf /cluster/tools/software/annovar/humandb --protocol ensGene,gnomad211_genome,cosmic68,avsnp142 --operation g,f,f,f --vcfinput

#gtf_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/gencode.v19.annotation.gtf
#fasta_file=/cluster/projects/kridelgroup/FLOMICS/genome_files/ucsc.hg19.fasta
#out_folder=/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/RNAseq_variants/annovar

#awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' LY_FL_470_T2.clean.vcf.hg19_multianno.vcf > LY_FL_470_T2.clean.vcf.hg19_multianno.clean.vcf
#grep -Ev '^(GL|HAP)' LY_FL_470_T2.clean.vcf.hg19_multianno.vcf > LY_FL_470_T2.clean.vcf.hg19_multianno.clean.vcf

#grep -Ev '^(chr17_ctg5_hap1|chrUn_gl000211|chrUn_gl000220|chrUn_gl000224|chr1_gl000192_random|chr6_apd_hap1|chr6_cox_hap2|chr6_dbb_hap3|chr6_mann_hap4|chr6_qbl_hap6|chr6_ssto_hap7)' LY_FL_013_T1.clean.vcf.hg19_multianno.vcf > LY_FL_013_T1.clean.vcf.hg19_multianno.clean.chrs.vcf

#module load tabix
#bgzip LY_FL_013_T1.clean.vcf.hg19_multianno.clean.chrs.vcf
#bcftools index LY_FL_013_T1.clean.vcf.hg19_multianno.clean.chrs.vcf.gz
#bcftools sort -Oz LY_FL_013_T1.clean.vcf.hg19_multianno.clean.chrs.vcf.gz -o LY_FL_013_T1.clean.vcf.hg19_multianno.clean.chrs.sorted.vcf.gz

#bcftools reheader LY_FL_013_T1.clean.vcf.hg19_multianno.clean.chrs.sorted.vcf.gz --fai /cluster/projects/kridelgroup/FLOMICS/genome_files/ucsc.hg19.fasta.fai > LY_FL_013_T1.clean.vcf.hg19_multianno.clean.chrs.sorted.reheader.vcf.gz
#bcftools index LY_FL_013_T1.clean.vcf.hg19_multianno.clean.chrs.sorted.reheader.vcf.gz

#grep -v "^#" LY_FL_470_T2.clean.vcf.hg19_multianno.clean.sorted.vcf.gz | cut -f 1 | sort | uniq -c
#bcftools view ${sample}.pass.vcf.gz --regions chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22

#gatk IndexFeatureFile -F LY_FL_013_T1.clean.vcf.hg19_multianno.clean.chrs.sorted.reheader.vcf.gz

#gatk VariantFiltration \
#   -R $fasta_file \
#   -V LY_FL_013_T1.clean.vcf.hg19_multianno.clean.chrs.sorted.reheader.vcf.gz  \
#   -O LY_FL_013_T1.clean.vcf.hg19_multianno.clean.chrs.sorted.filtered.vcf.gz \
#   -window 35 -cluster 3 --filter-name "FS" --filter-expression "FS > 30.0" --filter-name "QD" --filter-expression "QD < 2.0"
