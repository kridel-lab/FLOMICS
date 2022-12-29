# Follicular Lymphoma Multiomics Project (FLOMICS) - integrating data from targeted DNA sequencing, methylome profiling, bulk RNAseq and single nuclei RNAseq to discover new subtypes

Molecular classification or subtyping of cancers is becoming an essential step toward clinical practice of precision medicine. In follicular lymphoma (FL), molecular classification or subtyping has emerged as a major unmet need. This project focuses on the analysis of multi-omics data to subclassify FL. The project aims to study gene mutations, the transcriptome and the methylome of FL, to rigorously define and validate molecular subtypes by unravelling inter- and intra-patient heterogeneity. Our hypothesis is that FL is not just one disease, but that it can be classified into biologically distinct subgroups.

## Patient cohorts

We adopted an inclusive approach, aiming to enrol primary FL patient samples from both population-based cohorts and from clinical trial series, and from patients with diverse clinical presentations. Overall, we included 4 cohorts:
1. A retrospective, multicentre cohort of patients presenting with limited-stage disease and treated with radiation or with advanced-stage disease, and treated with immunochemotherapy, accrued from Toronto, Montreal, Kingston, Aarhus, Oslo and Brisbane (referred to as “retrospective” cohort);
2. A retrospective cohort of existing data from a previously published study in which treatments were heterogeneous (referred to as “PLOSMED” cohort);
3. The E2408 trial, a randomized phase II trial in which patients were treated with bendamustine and rituximab, with or without bortezomib, followed by maintenance with rituximab, with or without lenalidomide;
4. The E4402 trial, a randomized phase III trial in which patients with low-tumor burden FL were treated with single-agent rituximab, followed by rituximab maintenance or re-treatment, as needed.

## Sample overview

All samples included, post QC, are shown below, as of 15 Dec 2022:
(image to be updated)
<img src="AllSamples_22July2020.png" alt="AllSamples_22July2020" width="400"/>

### Methylome: DNA methylation EPIC Microarray Data Analysis

In February 2019, DNA methylation data for another 147 patients using Illumina MethylationEpic BeadChip were received. The breakdown of these samples were as follows:
- 10 diffuse large B-cell lymphoma (DLBCL) = “aggressive” control
- 5 reactive lymph node (RLN) = “normal” control
- 132 FL = cases
With the 30 FL samples from pilot dataset, a total of 177 samples and 162 FL cases. Of these, only 170 were used for downstream analysis based on QC. Copy number data was derived for 165 of these samples (170 - 5 RLN “normal” controls).

#### Main Scripts

- Preprocessing/QC [[Code]](https://github.com/kridel-lab/FLOMICS/blob/master/Code/Analysis/Methylation/1_QCRemoveSamples.R)
- Differential methylation via probes [[Code]](https://github.com/kridel-lab/FLOMICS/blob/master/Code/Analysis/Methylation/7_DifferentialMethylation.R)
- Differential methyaltion via regions [[Code]](https://github.com/kridel-lab/FLOMICS/blob/master/Code/Analysis/Methylation/9_DifferentiallyMethylatedRegions.R)
- Differential variability via probes [[Code]](https://github.com/kridel-lab/FLOMICS/blob/master/Code/Analysis/Methylation/8_DifferentialVariability.R)
- Tumor purity [[Code]](https://github.com/kridel-lab/FLOMICS/blob/master/Code/Analysis/Methylation/24_Tumor_purity_check.R)
- Clustering via Kmeans, Medoids, Hierarchical, RPMM and InfiniumClust [[Code]](https://github.com/kridel-lab/FLOMICS/blob/master/Code/Analysis/Methylation/14_Clustering_Kmeans_Medoids_Hierarchical_RPMM_InfiniumClust.R)
- tSNE [[Code]](https://github.com/kridel-lab/FLOMICS/blob/master/Code/Analysis/Methylation/13_tSNEPlot.R)
- Estimate cell counts [[Code]](https://github.com/kridel-lab/FLOMICS/blob/master/Code/Analysis/Methylation/12_EstimateCellCountsMethylation.R)
- Other scripts [[Code]](https://github.com/kridel-lab/FLOMICS/tree/master/Code/Analysis/Methylation)

### Transcriptome: RNA-seq Data Analysis

In February 2020, RNAseq data for 136 samples were obtained. After removing T2 samples (3) and an unmatching sample with methyaltion data (1), 132 samples were available.  The breakdown of these samples were as follows:
- 10 diffuse large B-cell lymphoma (DLBCL) = “aggressive” control
- 1 reactive lymph node (RLN) = “normal” control
- 121 FL = cases

OICR returned 19 samples that passed their QC criteria in June 2022
- 19 FL samples returned (rawdata:/cluster/projects/kridelgroup/FLOMICS/2022_combined_RNA-Seq_analysis/rawdata/2022_OICR/)

E4402 samples that were sequenced at BC Cancer in 2017 were also included in the 2022 uniform QC analyses
- 210 FL samples (rawdata:/cluster/projects/kridelgroup/FLOMICS/DATA/E4402/RNAseq/GSC-1464_fastq/)

#### RNA-seq workflow

<img src="RNAseq_2022_QC.png" alt="RNAseq_2022_QC" width="600"/>

#### Main QC steps

- Screen for rRNA contamination [[Code]](Code/BioinformaticsProcessing/RNAseq/2022_Uniform_QC/rRNA_cont_cal)
- Calculate the percentage of aligned coding bases [[Code]](Code/BioinformaticsProcessing/RNAseq/2022_Uniform_QC/rRNA_cont_cal)
- Running qualimap bamqc to calculate the insert size [[Code]](Code/BioinformaticsProcessing/RNAseq/2022_Uniform_QC/coding_bases_collectRnaSeqMetrics)
- Collect the STAR log files [[Code]](Code/BioinformaticsProcessing/RNAseq/2022_Uniform_QC/STAR_QC)

Tier2: 290 sample passed
rrna_contam_perct<=35 &&
picard_RnaMetrics_perct>=5 (PF_BASES/PF_ BASES)


#### RNAseq data processing

- Pre-processing: merging or renaming the samples (TGL 136, OICR 19, E4402 210); remove adapters and low-quality bases (trimmomatic-0.39) [[Code]](https://github.com/kridel-lab/FLOMICS/Code/BioinformaticsProcessing/RNAseq/E4402/trimmomatic-0.39-2_conda_QC_parallel.sh)
- mapping: mapping against refence genome – STAR/2.7.9a (Spliced Transcripts Alignment to a Reference), which a splice-aware alignment tool with two-step process: [[Code]](https://github.com/kridel-lab/FLOMICS/Code/BioinformaticsProcessing/RNAseq/E4402/STAR_parallel_sbatch_v37.sh)

  - create a genome index (consistent with the software version)
human genome build- “GRCh37.primary_assembly.genome.fa”
annotation file  - “gencode.v37lift37.annotation.gtf”

  - map reads to the genome
[[Code]](Code/BioinformaticsProcessing/RNAseq/AlignmentGeneCounts/)
STAR_log files per sample were collected as well to evaluate the mapping quality
- counting:use the resulting BAM files as input to count tools htseq-count /0.11.0 to obtain the raw counts per gene per sample, then merge into the final expression matrix [[Code]](https://github.com/kridel-lab/FLOMICS/Code/BioinformaticsProcessing/RNAseq/E4402/htseq_parallel_sbatch_v3_grch37.sh)


#### investigate and adjust the Batch-effect:
- BactchQC was used to investigate the batch effect:
Running BatchQC, you will need two files:
  - A gene by sample matrix with gene IDs in the first column and sample IDs as column headers. The cells contain quantile normalized expression values.
  - A metadata file with sample IDs in the first column and information about the samples in the remainder It should include the suspected batch variables, such as Sequencing Platform, Data, Biopsy Site, etc., as well as your classifier (e.g. tumor type).
- ComBat-seq was used to adjust the batch effect:
ComBat-seq takes untransformed, raw count matrix as input, and it requires a known batch variable.
- filter out the low-exp genes (optional): filterByExpr function from edgeR automatically filter low exps genes

#### Mutation profiling

- Variant calling from RNA-seq aligned bam files [[Code]](Code/BioinformaticsProcessing/RNAseq/VariantCalling/
)
- Mutation association [[Code]](https://github.com/kridel-lab/FLOMICS/blob/master/Code/Analysis/RNAseq/38_RNAseqToMutationCalls01.R)
- Visualization of mutations across clusters and stages [Code]
- Extact BCL2 and BCL6 translocation info from Manta predictions and merge with previous data from BC [[Code]](Code/Analysis/DNAseq/xxx_script_extract_BCL2_BCL6_translocations_from_Manta.R)


### Genome: Targeted DNAseq Data Analysis

Intial data for 133 samples submitted for hybridization-based capture sequencing were received after 22 September 2019 from Centre for Lymphoid Cancer at British Columbia Cancer Agency (BCCA). Final raw data were obtained in July 2020. Of these, x were used for downstream analysis based on QC.
- diffuse large B-cell lymphoma (DLBCL) = “aggressive” control
- reactive lymph node (RLN) = “normal” control
- FL = cases
With the 30 patients from pilot dataset, a total of x samples.

#### Scripts

#### Mutation profiling

- Variant calling and annotation from RNA-seq alignments [[Code]](Code/BioinformaticsProcessing/RNAseq/VariantCalling/
)
- Mutation calling and annotation from targeted DNA sequencing [Code]
- Visualization of mutations across clusters and stages [Code]
- Extact BCL2 and BCL6 translocation info from Manta predictions and merge with previous data from BC [[Code]](Code/Analysis/DNAseq/xxx_script_extract_BCL2_BCL6_translocations_from_Manta.R)

### Immune deconvolution
- Seurat analysis using snRNAseq data [[Code]](Code/BioinformaticsProcessing/snRNAseq/)
- Bisque analysis using seurat clusters and bulk rna-seq count matrix [[Code]](Code/Analysis/snRNAseq/)
- Plotting estimated immune fractions [[Code]](Code/Analysis/RNAseq/RNAseq-immune-deconvolution-bisque.R)
- Estimated immune fractions versus mutation status [[Code]](Code/Analysis/RNAseq/RNAseq-immune-deconvolution-mutation-correlation-summary-results.R)
