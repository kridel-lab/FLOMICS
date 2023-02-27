# CAPSEQ DNAseq Workflows

> Author: Victoria Shelton

This folder contains scripts for the analysis of Capture-Seq/Targeted DNAseq 

--------------------------
### Folder Descriptions

#### >[1_Calling_Variants_Pipeline](1_Calling_Variants_Pipeline/)
Scripts for executing the CAPSEQ pipeline of variant calling and gathers coverage data.

#### >[2_Variant_Filtering_Workflow](2_Variant_Filtering_Workflow/)
Scripts for filtering variant calls.

#### >[3_SNV_Clustering_Workflow](3_SNV_Clustering_Workfloww/)
Scripts for executing the clustering of variant calls.

--------------------------
## 1_Calling_Variants_Pipeline Files:
[All files described in this section are stored in the 1_Calling_Variants_Pipeline folder](1_Calling_Variants_Pipeline/)

#### Snakemake file
- CAPSEQ_Pipeline_BC.smk 
  - Produces sample coverage summary, SNV calls via Mutect2, and SNV annotations via Annovar. 
  - Requires inputs given in the config.yaml file

#### Configuration file
- [*config.yaml*](1_Calling_Variants_Pipeline/config.yaml)
  - Holds the paths to input files required for running the CAPSEQ Pipeline for variant calling.

#### Shell excution file
- snakemake_CAPSEQ_analysis.sh
  - File used to execute the CAPSEQ Pipeline 
  - Requires information given in the cluster.yaml file

#### Accompanying execution file
- cluster.yaml
  - File holding execution information to accompany snakemake_CAPSEQ_analysis.sh

#### Preparing BAM files for the CAPSEQ pipeline
- Highly recommend to execute the CAPSEQ pipeline upon symlinks and not on the raw/original BAM files themseleves. 
- Ensure the BAM files to be input do not contain any "-"/hyphens/dashes in their filenames. Otherwise the CAPSEQ pipeline will not execute.

## Running the CAPSEQ Pipeline: 
>Preparation:
1. Prepare the necessary CAPSEQ Pipeline input files:
    - BAM files (*.bam*) need to be stored together in a single directory.
    - It is highly recommended to create a directories of symbolic links pointing to your BAM files, and to not execute the CAPSEQ Pipeline directly upon your raw/original BAM files.
    - Ensure that the names of your sample BAM files do not contain any "-"/hyphens/dashes. Otherwise the CAPSEQ Pipeline will not be able to read your sample BAM files.
    - Determine which human reference genome was utilized to generate your sample BAM files and supply the paths to all other Pipeline inputs listed in [*config.yaml*](1_Calling_Variants_Pipeline/config.yaml). This includes:
      - The compressed and decompressed human reference genome fasta file
      - The gtf file of gene structure for the human reference genome
      - The raw sites file which contains population allele frequencies for human reference genome (*af-only-gnomad.raw.sites.b37.vcf.gz*)
      - The gene annotation directory utilized by Annovar. Ex.  `/../../../annovar/humandb`
      - The target and/or amplicon probe coordinate bed files
2. Create and supply an output directory path within *config.yaml* (after `outputDIR`). This is where all CAPSEQ Pipeline output files will be stored and organized for you.
3. Specify the path to an existing directory for the CAPSEQ Pipeline slurm files to be stored within *cluster.yaml* after `out:`, and within shell excution script after `#SBATCH -o`.
4. Within the shell execution script, also supply the path to the desired execution directory after `cd` on line 18, the path to where the CAPSEQ Pipeline snakemake file is stored after `snakemake` on line 20, and the path to where the CAPSEQ *cluster.yaml*" file is stored after `  --cluster-config` on line 23. 
  - Be sure to accurately name the job here as well. Example:
  `#SBATCH -J testrun-CAPSEQ`
5. Within the CAPSEQ Pipeline snakemake file (CAPSEQ_Pipeline_BC.smk) supply the path to where the *cluster.yaml* and *config.yaml* files are stored after each respective `configfile:` (line 10 and 11).

>Execution:
1. To execute the CAPSEQ Pipeline, modify the following command to indicate the shell execution script you will be using, and execute:
      `sbatch snakemake_CAPSEQ_analysis.sh`
2. To follow the progress of CAPSEQ Pipeline execute the `squeue` command and monitor the snakemake slurm files.

--------------------------
## 2_Variant_Filtering_Workflow:
[All files described in this section are stored in the 2_Variant_Filtering_Workflow folder](2_Variant_Filtering_Workflow/)

#### Step 1: Gathering variant calls and first stage filtering
- 001_Mutect2_prepare_matrix.R
  - Gathers variants calls into a single matrix
  - Excutes the first stage of filtering upon gathered variant calls
  - This occurs after variant calls have been annotated with Annovar
  - Requires "_annovar.vcf.gz.hg19_multianno.vcf" files

#### Step 2: Second stage filtering
- 002_filtering_Mutect2.R
  - Excutes the second stage of filtering upon gathered variant calls
  - This occurs after Step 1 of the 2_Variant_Filtering_Workflow
  - Requires the "sample_based_coverage_summary.csv" output by 1_Calling_Variants_Pipeline, "variant_calls.txt" output by 001_Mutect2_prepare_matrix.R, and any sample annotation necessary for binding the coverage summary and variant calls together. 

## Running the Variant Filtering Workflow: 
>Preparation:
1. Prepare the necessary CAPSEQ Pipeline input files:
    - Ensure the "_annovar.vcf.gz.hg19_multianno.vcf" files of all sample variant calls you wish to be gathered into a single matrix are placed into a single folder. ( This is for 001_Mutect2_prepare_matrix.R )
    - Know where the "sample_based_coverage_summary.csv" output by the 1_Calling_Variants_Pipeline is located, and where the "variant_calls.txt" output by 001_Mutect2_prepare_matrix.R is deposited. If additional information is required to connect the information in these two files by sample of origin is reuiqred, please provide this as well. ( This is for 002_filtering_Mutect2.R )

>Execution:
1. Run Step 1 first via `Rscript 001_Mutect2_prepare_matrix.R`
2. Run Step 2 via `Rscript 002_filtering_Mutect2.R`

--------------------------
## 3_SNV_Clustering_Workflow:
[All files described in this section are stored in the 3_SNV_Clustering_Workflow folder](3_SNV_Clustering_Workflow/)

#### Step 1: Creating SNV matrix (clustering analysis input)
- 001_make_SNV_matrix_and_plot.R
  - Summaries the variant calls output by 2_Variant_Filtering_Workflow
  - Outputs a summary matrix of SNV calls for each sample as a binary encoded 1/0
    - a 1 is denoted for any gene with a positive SNV status (aka; an SNV located in the respective gene was called for this sample), "_gene_vs_sample_mat_SNV_matrix.txt. 
  - Requires a file(s) listing the genes within the targeted panel, "panel1.csv", and the output "_exonic_filtered_MUTECT2_calls.txt" file containing the filtered variant calls output by 2_Variant_Filtering_Workflow. 

#### Step 2: Consensus clustering of the SNVs
- 002_SNV_ConsensusClustering_via_bootstrap.R
  - Performs consensus clustering as outlined in the Supplemental Methods to determine the optimal number of clusters exhibited by the data. 
  - Outputs the clustering seeds utilized in the consensus clustering analysis "clustering_seeds.txt", the consensus matrices "_aic_consensus_mat.txt" and "_bic_consensus_mat.txt", heatmaps visualizing the concensus matrices "_aic_consensuscluster.pdf" and "_bic_consensuscluster.pdf", and plots summarizing the the occurence of different cluster numbers "_aic_cluster_barplot.pdf" and "_bic_cluster_barplot.pdf".
  - Requries the "gene_vs_sample_SNV_mat_matrix.txt" matrix generated by 001_make_SNV_matrix_and_plot.R, as well the plotting functions scripted in "script-plot-all-mutations-all-cohorts_to_be_loaded.R" need to be loaded into the work environment. 

#### Step 3: Clustering and cluster stability analysis 
- 003_SNV_flexmix_GMM_Clustering_and_Stability_Analysis.R
  - Performs SNV clustering with cluster number determined to be of highest proportion from consensus clustering analysis(002_SNV_ConsensusClustering_via_bootstrap.R), as well as computing metrics for assessing the cluster stability and confidence of the clustering. 
  - Outputs the the clustering seeds utilized in the analysis ("_cluster_clustering_seeds.txt"), the components plot of each clustering run "_components-plot_5-clusters.pdf", the clustermap of the clustering output "_clustermap_aic.png", the significance scores of the genes input ("_sig_long_out_"), the cluster labels attributed to each sample by each clustering run "_GMM_Cluster_Labels_flexmix_clusters.csv", the number of samples attributed to each cluster ("_5_AIC_Cluster_Patient_Number_Matrix.csv"), the stability scores of each cluster ("_5_AIC_Cluster_Stability_Matrix.csv"), as well as the clustering stability score sheet ("aic_cluster_stability_score_sheet.csv"). 
  - Requries the "gene_vs_sample_SNV_mat_matrix.txt" matrix generated by 001_make_SNV_matrix_and_plot.R, file(s) listing the genes within the targeted panel(s) ("panel1.csv"), files listings which samples pertain to which cohort (if applicable), the aic and bic consensus matrices output by 002_SNV_ConsensusClustering_via_bootstrap.R ( "_aic_consensus_mat.txt" and "_bic_consensus_mat.txt"), and the plotting functions scripted in "script-plot-all-mutations-all-cohorts_to_be_loaded.R" need to be loaded into the work environment. 

#### Plotting functions file
-  script-plot-all-mutations-all-cohorts_to_be_loaded.R
  - This file contains the scripts for plotting functions which generate the cluster heatmap figure. Originally downloaded from https://github.com/ecsg-uoy/DLBCLGenomicSubtyping, only slightly modified to produce outputs.
  - Outputs the significance scores of genes included in the clustering ("_sig_long_out_"), the samples listed in the order they appear in the heatmap plotted for each cluster run ("_pat_order_"), and the dataframe used to plot the heatmap ("_df_long_"). 
  - Requires the variables 'i' and 'criteria' to be defined, alongside others. These variables are definied in 002_SNV_ConsensusClustering_via_bootstrap.R and 003_SNV_flexmix_GMM_Clustering_and_Stability_Analysis.R for their respective analyses. 

## Running the SNV Clustering Workflow: 
>Preparation:
1. Prepare the necessary required files listed above for each workflow component.
2. This workflow can be time and resource intensive. Please ensure enough time and resources are scheduled to run this analysis. The execution examples are listed in their simplest form below, but it is highly recommended to run this workflow in the background, for example by utilizing the SLURM job scheduler or similar schedulers at your disposal. 

>Execution:
1. Run Step 1 via `Rscript 001_make_SNV_matrix_and_plot.R`
2. Run Step 2 via `Rscript 002_SNV_ConsensusClustering_via_bootstrap.R` 
3. Run Step 3 via `Rscript 003_SNV_flexmix_GMM_Clustering_and_Stability_Analysis.R`


