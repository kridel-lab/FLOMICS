# CAPSEQDNAseq Pipeline

> Author: Victoria Shelton

This folder contains scripts for the analysis of Capture-Seq/Targeted DNAseq 

--------------------------
### Folder Descriptions

#### >[1_Calling_Variants_Pipeline](1_Calling_Variants_Pipeline/)


#### >[Scripts](Scripts/)
This folder contains the individual scripts/code on which the CAPSEQ Pipeline is based. Supplementary data analysis scripts are also stored here.

#### >[TargetedDNAseq_pipeline](TargetedDNAseq_pipeline/)
This folder contains the relevant scripts for executing the CAPSEQ pipeline.

#### >[TargetedDNAseq_pipeline/snakemake-tests](TargetedDNAseq_pipeline/snakemake-tests/)
Where all previous versions of the Snakemake pipeline are stored. These snakemake files were used when building and testing the CAPSEQ piepline. File *test003-pipeline.smk* is the final test version which spawned *CAPSEQ-PpielineBC.smk* and *CAPSEQ-Pipelilne_PLOS-E4402.smk*.

--------------------------
## Pipeline Files:
[All files described in this section are stored in the TargetedDNAseq_pipeline folder](TargetedDNAseq_pipeline/)

#### Snakemake file selection
> Tumor-Only analysis (only requries Tumor sample BAM files)
- CAPSEQ_Pipeline_BC.smk 
  - Intended for the BC/UHN and E2408 cohorts.
  - Produces sample coverage summary, SV calls via Manta, CNV calls via CNVkit, SNV calls via Mutect2, and SNV annotations via Annovar. 
  - Requires both target AND amplicon probe coordinate bed files.
- CAPSEQ_Pipeline_PLOSMED-E4402.smk 
  -  Intended for the PLOSMED and E4402 cohorts.
  -  Produces SNV calls via Mutect2 and SNV annotations via Annovar[^1]. 
  -  Requires either target OR amplicon probe coorindate bed files.
> Tumor-Normal analysis (requries both Tumor and Normal sample BAM files)
- CAPSEQ_Pipeline_Tumor-Normal_rev2.smk[^2]
  - Intended for the E2408 cohort.
  - Produces sample coverage summary, SV calls via Manta, SV calls via Gridss, CNV calls via CNVkit, SNV calls via Mutect2, and SNV annotations via Annovar.
  - Requires both target AND amplicon probe coordinate bed files.
- CAPSEQ_Pipieline_PLOSMED_Tumor-Normal.cmk[^3]
  - Intended for the PLOSMED cohort.
  - Produces SNV calls via Mutect2 and SNV annotations via Annovar[^1].
  - Requires either target OR amplicon probe coorindate bed files.

[^1]: This pipeline will be updated to include CNVkit copy number variant calling of BAM files.
[^2]: The name of this file will be changed to CAPSEQ_Pipeline_E2408_Tumor-Normal.smk
[^3]: This pipeline has yet to be run and tested (2022-01-03)

#### Config file selection
> Tumor-Only analysis
- config_2.yaml[^4]
  - Intended for all Tumor-Only analyses
- config.json[^4]
  - Is the same as config_2.yaml, but in json format
> Tumor-Normal analysis
- config_Tumor_Normal.yaml
  - Intended for E2408 Tumor-Normal analysis 
- config_PLOSMED_Tumor_Normal.yaml
  - Intended for PLOSMED Tumor-Normal analysis 

[^4]: These config files were modifed before conducting variant calling upon a different/new cohort.

#### Shell excution file selection
- snakemake_CAPSEQ_analysis.sh
  - Used when conducting Tumor-Only analysis
- snakemake_CAPSEQ_analysis_Tumor-Normal.sh 
  - Used when conducting Tumor-Normal analysis

#### Preparing BAM files for the CAPSEQ pipeline (To be moved to a new folder)
- Creating_symlinks.sh  
  - For creating symlinks that point to the BAM files. Highly recommend to execute the CAPSEQ pipeline upon symlinks and not on the BAM files themseleves. 
- renaming_BAMs.R 
  - For renaming the BAM files to ensure they do not contain any "-"/hyphens/dashes in the filename. Otherwise the CAPSEQ pipeline will not execute.

###### Testing and Obselite files & scripts (To be moved into a new folder)
- CAPSEW_Pipeline_PLOSMED_Tumor_Normal.smk
- CAPSEQ_Pipeline_Tumor-Normal_rev1.smk
- config.yaml
- TargetedDNAseq-pipeline.smk


## Running the CAPSEQ Pipeline: 
>Preparation:
1. Determine whether you are conducting Tumor-Only or Tumonr-Normal analysis, and whether target probe coordinates and/or amplicon probe coordinates are available for your cohort of samples. [Select the appropriate snakemake pipeline file.](https://github.com/kridel-lab/CapSeq_pipeline/blob/vic-test/README.md#snakemake-file-selection) 
2. Prepare the necessary CAPSEQ Pipeline input files:
    - BAM files (*.bam*) need to be stored together in a single directory. Tumor BAM files and Normal BAM files should be stored in separate directories.
    - It is highly recommended to create a directories of symbolic links pointing to your BAM files, and to not execute the CAPSEQ Pipeline directly upon your BAM files.
    - **!!!VERY IMPORTANT!!!** Ensure that the names of your sample BAM files do not contain any "-"/hyphens/dashes. Otherwise the CAPSEQ Pipeline will not be able to read your sample BAM files. See [Preparing BAM files for the CAPSEQ pipeline.](https://github.com/kridel-lab/CapSeq_pipeline/blob/vic-test/README.md#preparing-bam-files-for-the-capseq-pipeline-to-be-moved-to-a-new-folder)
    - Determine which human reference genome was utilized to generate your sample BAM files and supply the paths to all other Pipeline inputs listed in [*config_2.yaml*](TargetedDNAseq_pipeline/config_2.yaml). This includes:
      - The compressed and decompressed human reference genome fasta file
      - The gtf file of gene structure for the human reference genome
      - The raw sites file which contains population allele frequencies for human reference genome (*af-only-gnomad.raw.sites.b37.vcf.gz*)
      - The Manta software directory and *configManta.py* script
      - The CNVkit software directory, cnvkit.py and the gene annotations file used by CNVkit
      - The gene annotation directory utilized by Annovar. Ex.  `/../../../annovar/humandb`
      - The target and/or amplicon probe coordinate bed files
3. Create and supply an output directory path within *config_2.yaml* (after `outputDIR`). This is where all CAPSEQ Pipeline output files will be stored and organized for you.
- Example:
   
`outputDIR: /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_OUTPUT/E2408-Tumor-Only/run001/`

4. Specify the path to an existing directory for the CAPSEQ Pipeline slurm files to be stored within *cluster.yaml* after `out:`, and within [the appropriate shell excution script](https://github.com/kridel-lab/CapSeq_pipeline/tree/vic-test#shell-excution-file-selection) after `#SBATCH -o`.
- *cluster.yaml* example:

` out:  "/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_OUTPUT/E2408-Tumor-Normal/run001/slurm-outs/{rule}.{wildcards}-%j.out"`
- Shell execution script example:

`#SBATCH -o /cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/CAPSEQ_Pipeline_OUTPUT/E2408-Tumor-Only/run001/slurm-outs/%x-%j.out`

5. Within the shell execution script, also supply the path to the desired execution directory after `cd` on line 18, the path to where the CAPSEQ Pipeline snakemake file is stored after `snakemake` on line 20, and the path to where the CAPSEQ *cluster.yaml*" file is stored after `  --cluster-config` on line 23. 
- Example:

`cd /cluster/home/vshelton/CapSeq_pipeline/TargetedDNAseq_pipeline`

`snakemake -s /cluster/home/vshelton/CapSeq_pipeline/TargetedDNAseq_pipeline/CAPSEQ_Pipeline_BC.smk \`

`    --cluster-config /cluster/home/vshelton/CapSeq_pipeline/TargetedDNAseq_pipeline/cluster.yaml \`
  - Be sure to accurately name the job here as well. Example:
  `#SBATCH -J ppE2408-CAPSEQ`
  
8. Within the CAPSEQ Pipeline snakemake file supply the path to where the *cluster.yaml* and *config.yaml* files are stored after each respective `configfile:` (either line 10 and 11, OR, line 11 and 12).
- Example:

`configfile: "/cluster/home/vshelton/CapSeq_pipeline/TargetedDNAseq_pipeline/cluster.yaml"`

`configfile: "/cluster/home/vshelton/CapSeq_pipeline/TargetedDNAseq_pipeline/config_2.yaml"`


>Execution:
1. To execute the CAPSEQ Pipeline, modify the following command to indicate the shell execution script you will be using, and execute:
      `sbatch snakemake_CAPSEQ_analysis.sh`
2. To follow the progress of CAPSEQ Pipeline execute the `squeue` command and monitor the snakemake slurm files.

------------------------------------------
## Additional SV caller execution:
[**Gridss**](https://github.com/PapenfussLab/gridss) is an additional SV caller used. High confidence SVs are determined by comparison of Gridss SV calls with Manta SV calls
### Files
[These files can be found in the Scripts folder](Scripts/)[^5]

***Note***: Files without a cohort indicated in the filename were utlized to process BCAugust2020 and BCJuly2021 cohorts under Tumor-Only analysis. Files with E2408 indicated in the filename were utilized to process the E2408 cohort under Tumor-Only analysis. All other files have cohort and analysis-type indicated in the filename.
[^5]: All scripts used to run Gridss SV calling will be moved to a new folder within the respository.
#### Script Execution Order 
- 011_GRIDSS_order_commands_2021.sh
#### Indivudal command scripts
- 011_A_GRIDSS_image_file.sh
- 011_B_GRIDSS_setupreference.sh
- 011_C_GRIDSS_preprocess.sh
- 011_D_GRIDSS_assebmle.sh
- 011_E_GRIDSS_call.sh
#### Grouped command scripts (convenient)
- 011_C-E_GRIDSS_commands.sh
- 011_C-E_GRIDSS_commands_E2408.sh
- 011_C-E_GRIDSS_commands_E2408_Tumor-Normal.sh
#### Additional Commands that must be run for Tumor-Normal analysis, but not Tumor-Only analysis
- 011_GRIDSS_RepeatMasker_E2408_Tumor-Normal.sh
- 011_GRIDSS_SomaticFiltering_E2408_Tumor-Normal.sh
- 011_GRIDSS_SomaticFiltering_E2408_Tumor-Normal_1samp.sh


### Running Gridss SV calling
>Tumor-Only analysis
1. Read the [order of commands](Scripts/011_GRIDSS_order_commands_2021.sh), execute [011_A](Scripts/011_A_GRIDSS_image_file.sh) and [011_B](Scripts/011_B_GRIDSS_setupreference.sh) if the mandatory Gridss files are not present. (This can take ~3 hours)
   - You do not need to execute the order of commands script, instead use the individual batch scripts to run Gridss after reading the order of commands.
2. Create a directory to store the Gridss output. Don't forget to specifiy this directory in each Gridss script.
3. Modify and execute the appropriate [grouped command script](https://github.com/kridel-lab/CapSeq_pipeline/blob/vic-test/README.md#grouped-command-scripts-convenient).

>Tumor-Normal analysis
1. Follow the same steps as above. **Note:** CAPSEQ_Pipeline_Tumor-Normal_rev2.smk conducts these steps, you do not need to execute these scripts if you are using CAPSEQ_Pipeline_Tumor-Normal_rev2.smk
2. After the grouped command script has completed generating output, execute the RepeatMasker script. 
3. Once the RepeatMasker script has completed running, execute the appropriate Somatic Filtering script. (If conducting Gridss calling upon only 1 sample, use 011_GRIDSS_SomaticFiltering_E2408_Tumor-Normal_1samp.sh. If not, use 011_GRIDSS_SomaticFiltering_E2408_Tumor-Normal.sh).

----------------------------------
## Variant Analysis Files:
[These files can be found in the Scripts folder](Scripts/)
>COVERAGE
- 003A_summarize_probe_coverage_2021.R
- 003B_summarize_probe_coverage_2021.R
>Mutect2
- 008_Mutect2_prepare_matrix_July2021.R
- 010_A_filtering_Mutect2_SNVs.R
- 010_B_combining_SNV_matrices.R  (used to combine the BCJuly2021 and BCAugust2020 cohorts)
- 010_C_gene_vs_sample_SNV_matrix.R
>Manta
- 004_C_processing_manta_job_2021.R (executes 004_B_processing_manta_2021.R)
- 004_D_make_matrix_manta_2021.R
- counting_SVs.R
>Gridss
#### Processing Gridss output
- 011_F_processing_GRIDSS_2021.R
- 011_F_processing_GRIDSS_E2408_2021.R
- 011_F_processing_GRIDSS_Tumor-Normal_2021.R
- 011_G_processing_GRIDSS_job_2021.sh
- 011_G_processing_GRIDSS_job_E2408_2021.sh
- 011_G_processing_GRIDSS_job_E2408_Tumor-Normal_2021.sh
#### Wrangling the results into a utilizable matrix format
- 011_H_make_matrix_GRIDSS_2021.R
- 011_H_make_matrix_GRIDSS_E2408_2021.R
- 011_H_make_matrix_GRIDSS_E2408_Tumor-Normal.R
#### Filtering for the SV calls we are interested in, output as tables
- 011_I_GRIDSS_sv_breakapart_predictions_BCAugust2020.R
- 011_I_GRIDSS_sv_breakapart_predictions_BCJuly2021.R
- 011_I_GRIDSS_sv_breakapart_predictions_E2408_2021.R
- 011_I_GRIDSS_sv_breakapart_predictions_E2408_Tumor-Normal_2021.R

>SV comparion

>SNV Clustering

## Variant Analysis Steps:
>COVERAGE

>Mutect2

>Manta

>Gridss
1. Modify both of the appropriate [011_F and 011_G files](https://github.com/kridel-lab/CapSeq_pipeline/tree/vic-test#processing-gridss-output) (two files in total for one cohort), and execute only the 011_G file. (011_G will execute 011_F).
2. Modify and execute the appropriate [011_H file](https://github.com/kridel-lab/CapSeq_pipeline/tree/vic-test#wrangling-the-results-into-a-utilizable-matrix-format) to wrangle SV results into a matrix
3. Use the correct [011_I file](https://github.com/kridel-lab/CapSeq_pipeline/tree/vic-test#filtering-for-the-sv-calls-we-are-interested-in-output-as-tables) to generate sv breakapart prediction tables. 
>SV comparion

>SNV Clustering
--------------------------------
--------------------------------
#### Incoming additions:
- [x] Write README description for CAPSEQ_Pipeline_Tumor-Normal_rev2.smk
- [x] Add steps for Gridss execution.
- [x] Add links to README sections.
- [ ] Upload SV comparison and merging scripts.
- [ ] Upload SNV Clustering scripts.
- [ ] Change filename of CAPSEQ_Pipeline_Tumor-Normal.smk to CAPSEQ_Pipeline_E2408_Tumor-Normal.smk
- [ ] Change filename of snakemake_CAPSEQ_analysis.sh to snakemake_CAPSEQ_Tumor-Only_analysis.sh
- [ ] Move Gridss scripts to a new folder within the repository.
- [ ] Write README descriptions and instructions for the remainder of Variant Analysis scripts.

