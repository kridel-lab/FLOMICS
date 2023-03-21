FLOMICS
├─ DNAseq
│  ├─ 1_Calling_Variants_Pipeline
│  │  ├─ CAPSEQ_Pipeline.smk
│  │  ├─ cluster.yaml
│  │  ├─ config.yaml
│  │  └─ snakemake_CAPSEQ_analysis.sh
│  ├─ 2_Variant_Filtering_Workflow
│  │  ├─ 001_Mutect2_prepare_matrix.R
│  │  └─ 002_filtering_Mutect2.R
│  ├─ 3_SNV_Clustering_Workflow
│  │  ├─ 001_make_SNV_matrix_and_plot.R
│  │  ├─ 002_SNV_ConsensusClustering_via_bootstrap.R
│  │  ├─ 003_SNV_flexmix_GMM_Clustering_and_Stability_Analysis.R
│  │  └─ script-plot-all-mutations-all-cohorts_to_be_loaded.R
│  ├─ README.md
│  └─ archive
│     ├─ 003_SNV_flexmix_GMM_Clustering_and_Stability_Analysis copy.R
│     └─ DNAseqREADME.md