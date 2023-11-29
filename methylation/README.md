# FLOMICS Methylation Analysis

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
