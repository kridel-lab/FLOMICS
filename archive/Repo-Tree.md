
```
Generated 2023-01-25
FLOMICS
├─ AllSamples_22July2020.png
├─ Code
│  ├─ Analysis
│  │  ├─ DNAseq
│  │  │  ├─ 000_RNAseq_QC_summary.R
│  │  │  ├─ 001_local_exploration_mutations.R
│  │  │  ├─ 002_mutational_patterns.R
│  │  │  ├─ 003_maftools_analysis.R
│  │  │  ├─ get_library_size_reads.R
│  │  │  ├─ library_size_vs_mean_coverage.R
│  │  │  ├─ plos_med_comparing_new_old_calls.R
│  │  │  ├─ script-Bayesian-prop-test.R
│  │  │  ├─ script_analysis_breakapart_predictions.R
│  │  │  ├─ script_analysis_mutations.R
│  │  │  └─ script_analysis_mutations_KI_new_PLOS.R
│  │  ├─ DNAseq_Methylation_Intersect
│  │  │  └─ SomeScript.R
│  │  ├─ IntegrativeAnalysis
│  │  │  ├─ 33_SNFClustering.R
│  │  │  └─ SNF_to_Cytoscape.R
│  │  ├─ Methylation
│  │  │  ├─ 0_PackageCheck.R
│  │  │  ├─ 10_BoxPlotsMethylationRNAseq
│  │  │  ├─ 11_HeatPlotMethylation.R
│  │  │  ├─ 12_EstimateCellCountsMethylation.R
│  │  │  ├─ 13_tSNEPlot.R
│  │  │  ├─ 14_Clustering_Kmeans_Medoids_Hierarchical_RPMM_InfiniumClust.R
│  │  │  ├─ 15_SurvivalAnalysis.R
│  │  │  ├─ 16_Huet23GeneModel.R
│  │  │  ├─ 17_LinePlots.R
│  │  │  ├─ 18_StandardDeviation.R
│  │  │  ├─ 19_IsolateEntries.R
│  │  │  ├─ 1_QCRemoveSamples.R
│  │  │  ├─ 20_GlmnetFeatureSelection.R
│  │  │  ├─ 21_ProportionVisualization.R
│  │  │  ├─ 22_MDSPlots.R
│  │  │  ├─ 23_BarplotMethylation.R
│  │  │  ├─ 24_Tumor_purity_check.R
│  │  │  ├─ 25_Variance.R
│  │  │  ├─ 26_Gprofiler.R
│  │  │  ├─ 27_CopyNumberAnalysis.R
│  │  │  ├─ 28_VennDiagramCreater.R
│  │  │  ├─ 29_TOASTtesting.R
│  │  │  ├─ 2_LoadingMethylData.R
│  │  │  ├─ 30_PieChart.R
│  │  │  ├─ 31_MeanSDPlot.R
│  │  │  ├─ 32_ClusterConfidencePlot.R
│  │  │  ├─ 33_SNFClustering.R
│  │  │  ├─ 34_MedianAbsoluteDeviation.R
│  │  │  ├─ 35_QCRNAseq.R
│  │  │  ├─ 36_DifferentialExpressionRNAseq.R
│  │  │  ├─ 37_ViolinPlotsRNAseq.R
│  │  │  ├─ 38_RNAseqToMutationCalls01.R
│  │  │  ├─ 39_RNAseqVsMethFC.R
│  │  │  ├─ 3_DensityPlotMethylation.R
│  │  │  ├─ 40_RESETRNAseqVsMethFCscript.R
│  │  │  ├─ 41_MethylMixRNAseqVsMethFCscript.R
│  │  │  ├─ 43_epiCMIT.R
│  │  │  ├─ 44_AlluvialPlots.R
│  │  │  ├─ 45_BioCircos.R
│  │  │  ├─ 46_scRNAseqBisque.R
│  │  │  ├─ 4_SummarizationMvalues.R
│  │  │  ├─ 5_SummarizationBeta.R
│  │  │  ├─ 6_SummaryStatistics.R
│  │  │  ├─ 7_DifferentialMethylation.R
│  │  │  ├─ 8_DifferentialVariability.R
│  │  │  ├─ 9_DifferentiallyMethylatedRegions.R
│  │  │  ├─ MainAnjaliMethylation.R
│  │  │  ├─ RD_minfi_Final.R
│  │  │  ├─ README.md
│  │  │  ├─ RK-plot-purity.R
│  │  │  ├─ RK-script-DNAmethyl-vs-EZH2.R
│  │  │  ├─ RK-script-FL-vs-norm.R
│  │  │  ├─ RK-script_replotGSEA.R
│  │  │  └─ RK_script-Bayesian-prop-test.R
│  │  ├─ MiXCR
│  │  │  ├─ 001_coxuni_mixcr.R
│  │  │  ├─ 002_mixcr_correlations-order-edits.R
│  │  │  ├─ 002_mixcr_correlations.R
│  │  │  ├─ correlation_source-order-edits.R
│  │  │  ├─ correlation_source.R
│  │  │  └─ run_correlations.sh
│  │  ├─ README.md
│  │  ├─ RK-outcome-analysis.R
│  │  ├─ RK-script-methylmix-enrichr.R
│  │  ├─ RNAseq
│  │  │  ├─ 002_local_exploration_gene_expression.R
│  │  │  ├─ 33_SNFClustering.R
│  │  │  ├─ 35_QCRNAseq.R
│  │  │  ├─ 36_DifferentialExpressionRNAseq.R
│  │  │  ├─ 37_ViolinPlotsRNAseq.R
│  │  │  ├─ 38_RNAseqToMutationCalls01.R
│  │  │  ├─ RNAseq-gene-mut-correlation.R
│  │  │  ├─ RNAseq-immune-deconvolution-bisque-table.R
│  │  │  ├─ RNAseq-immune-deconvolution-bisque-tier-based.R
│  │  │  ├─ RNAseq-immune-deconvolution-bisque.R
│  │  │  ├─ RNAseq-immune-deconvolution-gene-exp-correlation.R
│  │  │  ├─ RNAseq-immune-deconvolution-mutation-correlation-bisque-only.R
│  │  │  ├─ RNAseq-immune-deconvolution-mutation-correlation-summary-results.R
│  │  │  ├─ RNAseq-immune-deconvolution-mutation-correlation.R
│  │  │  ├─ RNAseq-immune-deconvolution.R
│  │  │  ├─ RNAseq-methylation-RESET-summary.R
│  │  │  ├─ RNAseq-methylation-RESET.R
│  │  │  ├─ generate_TMM_matrix_from_star_counts.R
│  │  │  └─ script_replotGSEA.R
│  │  ├─ RNAseq_Methylation_Intersect
│  │  │  ├─ 39_RNAseqVsMethFC.R
│  │  │  └─ 40_RESET_RNAseqVsMethFC_script.R
│  │  ├─ load_scripts_data_KI.R
│  │  └─ snRNAseq
│  │     ├─ bisque_001_sarah_final_code.R
│  │     ├─ singleR_annotation_AndyZhang.R
│  │     └─ singleR_annotation_LZ_DZ_affydata.R
│  ├─ BioinformaticsProcessing
│  │  ├─ 42_xlsxToGMTformat.R
│  │  ├─ DNAseq
│  │  │  ├─ BC_Manta_results
│  │  │  │  ├─ variants_005_summary_manta.R
│  │  │  │  ├─ variants_006_summary_manta.sh
│  │  │  │  ├─ variants_007_summary_manta.R
│  │  │  │  └─ variants_008_summary_manta.R
│  │  │  ├─ PLOS_MED
│  │  │  │  ├─ 001_Prepare_Amplicons_Targets.R
│  │  │  │  ├─ 002_CollectTargetedPcrMetrics.sh
│  │  │  │  ├─ 003_Run_Mutect2.sh
│  │  │  │  ├─ 004_Mutect2_filter_variants.sh
│  │  │  │  ├─ 005_Mutect2_Annotate_via_Annovar.sh
│  │  │  │  ├─ 006_Mutect2_prepare_matrix.R
│  │  │  │  └─ 007_Mutect2_add_sample_info.R
│  │  │  ├─ TargetSeq_call_mutations
│  │  │  │  ├─ 000_total_coverage_vs_mean_coverage.sh
│  │  │  │  ├─ 001_Prepare_Amplicons_Targets.R
│  │  │  │  ├─ 002_CollectTargetedPcrMetrics.sh
│  │  │  │  ├─ 003A_summarize_probe_coverage.R
│  │  │  │  ├─ 003B_summarize_probe_coverage.R
│  │  │  │  ├─ 004_A_platypus.sh
│  │  │  │  ├─ 004_B_platypus_cleanup.R
│  │  │  │  ├─ 004_B_platypus_cleanup.sh
│  │  │  │  ├─ 004_C_platypus_cleanup.sh
│  │  │  │  ├─ 005_A_manta.sh
│  │  │  │  ├─ 005_B_processing_manta.R
│  │  │  │  ├─ 005_C_processing_manta_job.sh
│  │  │  │  ├─ 005_D_make_matrix_manta.R
│  │  │  │  ├─ 006_Mutect2.sh
│  │  │  │  ├─ 007_Mutect2_filter_variants.sh
│  │  │  │  ├─ 008_Mutect2_Annotate_via_Annovar.sh
│  │  │  │  ├─ 009_Mutect2_prepare_matrix.R
│  │  │  │  ├─ 010_Mutect2_merged_wPlatypus_old_data.R
│  │  │  │  └─ order-scripts.sh
│  │  │  └─ summarizing_coverage_libraries_sequencing.R
│  │  ├─ RNAseq
│  │  │  ├─ 2022_Uniform_QC
│  │  │  │  ├─ 2022_Sep_RNA-seq_sample_QC.R
│  │  │  │  ├─ STAR_QC
│  │  │  │  │  ├─ extract_samples.R
│  │  │  │  │  ├─ get_QC_STAR_stats.R
│  │  │  │  │  ├─ read_STAR.R
│  │  │  │  │  ├─ read_sample_files.R
│  │  │  │  │  └─ rename_commd.bk
│  │  │  │  ├─ bk
│  │  │  │  │  ├─ RSeQC_check_rrna_hg19_forloop_QC.sh
│  │  │  │  │  ├─ RSeQC_check_rrna_hg19_sbatch_QC.sh
│  │  │  │  │  ├─ deal_qualimap_generate_exonic_summary.sh
│  │  │  │  │  ├─ generate_RSeQC_rrna_summary.sh
│  │  │  │  │  ├─ qualimap_bamQC.sh
│  │  │  │  │  ├─ qualimap_sbatch_BAM_QC.sh
│  │  │  │  │  └─ qualimap_sbatch_RNA_QC.sh
│  │  │  │  ├─ cal_insert_size
│  │  │  │  │  ├─ deal_qualimap_generate_insertsize_summary.sh
│  │  │  │  │  └─ qualimap_bamQC.sh
│  │  │  │  ├─ coding_bases_collectRnaSeqMetrics
│  │  │  │  │  ├─ deal_picard_RnaMetrics_summary_v2.sh
│  │  │  │  │  └─ picard_collectRnaSeqMetrics.sh
│  │  │  │  ├─ comBat-seq_batch_adj.R
│  │  │  │  ├─ commd.bk
│  │  │  │  └─ rRNA_cont_cal
│  │  │  │     ├─ deal_rRNA_contam_output_summary.sh
│  │  │  │     ├─ mapping2rRNA_E4402_v2.sh
│  │  │  │     └─ mapping2rRNA_TGL_n_OICR2022_v2.sh
│  │  │  ├─ 2022_combined_mapping_counting
│  │  │  │  ├─ E4402
│  │  │  │  │  ├─ STAR_parallel_sbatch_v37.sh
│  │  │  │  │  ├─ dataprep.readme
│  │  │  │  │  ├─ fastp_parellel_sbatch.sh
│  │  │  │  │  ├─ htseq_parallel_sbatch_v3_grch37.sh
│  │  │  │  │  ├─ rename_rawsymlink_sampleid_R1.sh
│  │  │  │  │  ├─ rename_rawsymlink_sampleid_R2.sh
│  │  │  │  │  └─ trimmomatic-0.39-2_conda_QC_parallel.sh
│  │  │  │  └─ TGL_n_OICR2022
│  │  │  │     ├─ STAR_parallel_sbatch_v37.sh
│  │  │  │     ├─ check_unfinished_bamfiles.sh
│  │  │  │     ├─ collect_counts_2022.R
│  │  │  │     ├─ commd.bk
│  │  │  │     ├─ htseq_forloop_sbatch_v3_grch37_pending16.sh
│  │  │  │     ├─ htseq_parallel_sbatch_v3_grch37.sh
│  │  │  │     ├─ rename_symlink.sh
│  │  │  │     └─ trimmomatic-0.39-2_conda_QC_parallel.sh
│  │  │  ├─ AlignmentGeneCounts
│  │  │  │  ├─ 001_generate_genome.sh
│  │  │  │  ├─ 002_run_STAR_job.sh
│  │  │  │  ├─ 002_run_STAR_script.sh
│  │  │  │  ├─ 003_generate_count_matrix.R
│  │  │  │  └─ 004_get_QC_from_STAR.R
│  │  │  ├─ ERVsDetection
│  │  │  │  ├─ 001_B_running_TELESCOPE_on_BAM_files_veryhimem.sh
│  │  │  │  ├─ 001_running_TELESCOPE_on_BAM_files.sh
│  │  │  │  ├─ 002_concatenating_TELESCOPE_results.sh
│  │  │  │  ├─ 003_A_get_ERV_counts_into_matrix.R
│  │  │  │  ├─ 003_A_get_ERV_counts_into_matrix_job.sh
│  │  │  │  ├─ 003_B_Clusters_DE_source.R
│  │  │  │  ├─ 003_B_Clusters_differential_expression_analysis.R
│  │  │  │  ├─ 003_B_Stages_DE_source.R
│  │  │  │  ├─ 003_B_Stages_differential_expression_analysis.R
│  │  │  │  ├─ 003_run_DE_analysis.sh
│  │  │  │  ├─ ERVs_diffMethylation
│  │  │  │  ├─ edgeR_GeneExpression.R
│  │  │  │  ├─ findOverlaps_ERV_CpG.R
│  │  │  │  ├─ format_DM_CpGs.R
│  │  │  │  ├─ plotting_heatmap_telescope_results.R
│  │  │  │  ├─ rGREAT_gprofiler_enrichment.R
│  │  │  │  └─ telescope_erv_visualization.R
│  │  │  ├─ FASTQ_prep
│  │  │  │  ├─ 002_A_merging_fastq_files_tgl.sh
│  │  │  │  ├─ 002_B_merging_fastq_files_tgl_second_upload.sh
│  │  │  │  ├─ 002_C_merging_fastq_files_tgl_nontopups.sh
│  │  │  │  ├─ 002_D_merging_fastq_files_tgl_second_upload_notopups.sh
│  │  │  │  └─ 002_E_merging_fastq_files_tgl_missing_files.sh
│  │  │  ├─ Kallisto
│  │  │  │  ├─ kallisto_001_generate_index.sh
│  │  │  │  ├─ kallisto_002_quantification_job.sh
│  │  │  │  ├─ kallisto_002_quantification_script.sh
│  │  │  │  ├─ tximport_001_import_kallisto_data.R
│  │  │  │  └─ tximport_001_import_kallisto_job.sh
│  │  │  ├─ MiXCR
│  │  │  │  ├─ 000_run_mixcr.sh
│  │  │  │  ├─ 001_run_mixcr_merge.sh
│  │  │  │  ├─ 002_mixcr_concatenate_all_clonotypes.R
│  │  │  │  ├─ 003_count_mixcr_processing.R
│  │  │  │  ├─ 004_clonal_fractions_mixcr_processing.R
│  │  │  │  ├─ 005_reads_unique_CDR3_mixcr_processing.R
│  │  │  │  ├─ 006_immunarch_diversity_indices.R
│  │  │  │  ├─ 007_vegan_diversity_indices.R
│  │  │  │  ├─ RK-script-IGH-isotype-analysis.R
│  │  │  │  ├─ RK-script-immunarch.R
│  │  │  │  ├─ RK-script-plot-diversity-from-Arash-script.R
│  │  │  │  ├─ script-Arash-RK-diversity-indices.R
│  │  │  │  └─ script-Arash-RK-diversity-indices_AS.R
│  │  │  ├─ VariantCalling
│  │  │  │  ├─ 001_call_variants_opossum_main_job.sh
│  │  │  │  ├─ 001_call_variants_opossum_start_snake.sh
│  │  │  │  ├─ 001_call_variants_snakemake.snakefile
│  │  │  │  ├─ 002_variant_annotation_cleanup.sh
│  │  │  │  ├─ 003_soft_variant_filtering_process.R
│  │  │  │  ├─ 003_soft_variant_filtering_process.sh
│  │  │  │  ├─ 003_soft_variant_filtering_process_job.sh
│  │  │  │  └─ 004_prepare_matrix.R
│  │  │  ├─ sort_FASTQ_files.sh
│  │  │  └─ sort_FASTQ_files_job.sh
│  │  └─ snRNAseq
│  │     ├─ 001_Seurat_sarah_final_code_get_clusters.R
│  │     ├─ 001_Seurat_submit_cluster.sh
│  │     ├─ 002_Seurat_sarah_final_code_visualize_genes_across_clusters.R
│  │     ├─ 002_Seurat_sarah_final_code_visualize_genes_across_clusters_local.R
│  │     ├─ 002_Seurat_submit_cluster.sh
│  │     ├─ 003_Seurat_rename_clusters_to_cell_types.R
│  │     ├─ 004_prepare_seurat_object_for_bisque.R
│  │     ├─ 005_Seurat_differential_expression_analysis_clusters.R
│  │     ├─ 006_Seurat_sarah_final_code_cell_cycle_analysis.R
│  │     ├─ 007_Seurat_B_cell_Heatmap.R
│  │     ├─ 007_Seurat_B_cell_Heatmap_local_version.R
│  │     ├─ 007_Seurat_B_cells_pathways_enrichment_local.R
│  │     ├─ 007_Seurat_T_cell_Heatmap.R
│  │     ├─ 007_Seurat_T_cell_Heatmap_local.R
│  │     ├─ 008_Seurat_pseudotime_analysis_local.R
│  │     ├─ doSeuratProc.R
│  │     ├─ seurat_cell_clusters_visualize_marker.R
│  │     └─ seurat_script_order.sh
│  ├─ Figures
│  │  └─ ImmuneDeconvolution_001.R
│  ├─ README.md
│  ├─ Tools
│  │  └─ RESET
│  │     ├─ DESCRIPTION
│  │     ├─ NAMESPACE
│  │     ├─ R
│  │     │  ├─ .Rhistory
│  │     │  ├─ FDRcal.R
│  │     │  ├─ eventScore.R
│  │     │  ├─ methNorSel.R
│  │     │  ├─ methStatus.R
│  │     │  ├─ reset.R
│  │     │  └─ shared_functions.R
│  │     ├─ README.md
│  │     └─ man
│  │        ├─ FDRcal.Rd
│  │        ├─ eventScore.Rd
│  │        ├─ methStatus.Rd
│  │        └─ reset.Rd
│  └─ archive
│     ├─ archive_001_Seurat.R
│     ├─ archive_bisque_001.R
│     ├─ variants_001_compiling_VCFs.sh
│     ├─ variants_001_compiling_VCFs_lofreq_only.sh
│     ├─ variants_001_compiling_VCFs_lofreq_only_Indels.sh
│     ├─ variants_001_compiling_VCFs_master.sh
│     ├─ variants_001_compiling_VCFs_master_platypus_only.sh
│     ├─ variants_001_compiling_VCFs_mutect2_only.sh
│     ├─ variants_001_compiling_VCFs_mutect2_only_Indels.sh
│     ├─ variants_002_helper_script.sh
│     ├─ variants_002_helper_script_lofreq_Indels_only.sh
│     ├─ variants_002_helper_script_lofreq_SNVonly.sh
│     ├─ variants_002_helper_script_master.sh
│     ├─ variants_002_helper_script_master_platypus_only.sh
│     ├─ variants_002_helper_script_mutect2_Indels_only.sh
│     ├─ variants_002_helper_script_mutect2_SNVonly.sh
│     ├─ variants_003_read_in_VCFs_into_matrix.R
│     ├─ variants_003_read_in_VCFs_into_matrix_mutect2_only_get_VAFs.R
│     ├─ variants_004_summary_visualize.R
│     └─ variants_004_summary_visualize_Mutect2_only.R
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
│  │  ├─ 003_SNV_flexmix_GMM_Clustering_and_Stability_Analysis copy.R
│  │  ├─ 003_SNV_flexmix_GMM_Clustering_and_Stability_Analysis.R
│  │  └─ script-plot-all-mutations-all-cohorts_to_be_loaded.R
│  ├─ DNAseqREADME.md
│  ├─ README.md
│  └─ archive
├─ Data
│  ├─ .Rhistory
│  └─ readme
├─ Figures
│  ├─ RK-Fig1.R
│  ├─ RK-Fig2B_E.R
│  ├─ RK-Fig3.R
│  ├─ RK-Fig5.R
│  ├─ RK-Fig6.R
│  └─ RK-snRNAseq-processing.R
├─ README.md
├─ RNAseq
│  ├─ 001_trimmomatic-0.39-2_conda_QC_parallel.sh
│  ├─ 002_STAR_parallel_sbatch_v37.sh
│  ├─ 003_htseq_parallel_sbatch_v3_grch37.sh
│  ├─ 004_collect_counts_2022.R
│  ├─ RNAseqREADME.md
│  ├─ STAR_QC
│  │  ├─ extract_samples.R
│  │  ├─ get_QC_STAR_stats.R
│  │  ├─ read_STAR.R
│  │  └─ read_sample_files.R
│  └─ archive
├─ RNAseq_2022_QC.png
├─ Repo-Tree.md
├─ UsedSamples_22July2020.png
├─ archive
├─ methylation
│  ├─ README.md
│  └─ rajesh
│     ├─ 00_minfi_analysis
│     │  ├─ RD_minfi_PQ.R
│     │  ├─ RD_minfi_PQ_All_Probes.R
│     │  └─ RD_minfi_SWAN.R
│     └─ README.md
└─ test_sarah_code.R

```