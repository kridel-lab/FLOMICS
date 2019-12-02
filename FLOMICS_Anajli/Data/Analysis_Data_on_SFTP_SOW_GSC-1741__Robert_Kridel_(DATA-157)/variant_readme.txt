Nov. 1 2019,

Some notes on running Platypus and LoFrec on Kridel special capture libraries

PLATYPUS:
This version was used : Platypus_0.8.1/Platypus.py

platypus allows for "regions of interest" to be provided and in this case that would be the specialCapture regions provided.  The appropriate ref. still also needs to be provided.  Otherwise only required parameters were provided.  Here is an example of a platypus command used for these runs:

 /gsc/software/linux-x86_64-centos7/python-2.7.11/bin/python2  Platypus.py callVariants --bamFiles=/projects/analysis/analysis33/B48323/merge47568_bwa-mem-0.7.6a-sb/125bp/hg19a/B48323_1_lane_dupsFlagged.bam --refFile=hg19a.fa --output=platypus/B48323.variant_calls.vcf --regions=IDT-specialCapture_sansCHR.bed --verbosity=2 --logFileName platypus_Log/B48323.txt --nCPU 2

Platypus automatically filters the variants into either a "PASS" or into about a half dozen fail modes (or combinations of ) from the inputs provided.  Eg.:

Q20
Q20;alleleBias
Q20;badReads
Q20;badReads;QD
Q20;badReads;QD;alleleBias
Q20;QD;alleleBias
SC;badReads
SC;badReads;alleleBias
SC;badReads;QD

The badReads filter is of particlar concern here as about 92% of potential variants were eliminated with this failure:

ls platypus/B48*.variant_calls.vcf|xargs grep -v ^#|wc -l
57073  (TOTAL VARIANTS ACROSS ALL 131 SAMPLES)

ls platypus/B48*.variant_calls.vcf|xargs grep badReads|wc -l
53000  (THOSE WHICH WERE badReads FAILS)

This filter triggers when across reads supporting a variant, the median of the minimum base quality close to the focal site (default 7 bp either side, configurable using ‐‐badReadsWindow) is too low (default 15 or less, configurable using ‐‐badReadsThreshold); this identifies systematic local sequencing issues causing an excess of read errors, which are not always accurately reflected in the read quality scores.

After collecting only PASSes,  the SNPs were seperated from INDELs as platypus outputs these as a single file.  These are provided as two files in the dissemination 

LOFREQ:
the version used:  /LoFreq/lofreq_star-2.1.3.1/bin/lofreq

Because we are intersted in INDEL variants as well a two step process was needed.
We added indel qualities to each BAM with a command like this for each sample:

/LoFreq/lofreq_star-2.1.3.1/bin/lofreq  indelqual --dindel --ref hg19a.fa --out DIndel/B.bam  /projects/analysis/analysis33/B48323/merge47568_bwa-mem-0.7.6a-sb/125bp/hg19a/B48323_1_lane_dupsFlagged.bam

then after indexing those the main loFreq executable was run to generate  the variant callset:

/LoFreq/lofreq_star-2.1.3.1/bin/lofreq call-parallel --pp-threads 3 -f hg19a.fa --call-indels  --bed IDT-specialCapture_sansCHR.bed -o loFreq/B48322.vars.vcf DIndel/B48322.dindel.bam

Lofreq onlys outputs PASSed variants.  The INDELs are annotated as such but outputted along with SNPs into a single VCF.  We have separated the two and provided them in separate files.

The intersection of the two sets of SNPs has been generated with bedtools as well as the two sets of INDELs. A typical command would look like:
 
/bedtools-2.27.1/bin/bedtools intersect -a loFreq/B48323.vars.vcf.pass.indel -b platypus/B48323.variant_calls.pass.vcf.indel -u -header


