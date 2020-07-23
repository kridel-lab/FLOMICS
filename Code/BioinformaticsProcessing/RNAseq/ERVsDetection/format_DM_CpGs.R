#----------------------------------------------------------------------------------
#format_DM_CpGs.R
#----------------------------------------------------------------------------------

#Sarah Russell
#Date:July 20th, 2020

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###LIBRARIES###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


date = Sys.Date()

# Loading needed packages
library(limma)
library(stringi)
library(data.table)
library(readxl)
library(dplyr)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###DATA###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/")

DE_ervs_stage=read.csv("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope/_2020-07-16_all_DE_ervs_stage.csv")
DE_ervs_cluster=read.csv("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope/_2020-07-16_all_DE_ervs_cluster.csv")
erv_annotations = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/hg19_clean_repeatmasker_annotations.bed")

sample_info = as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.xlsx"))

cpg_annotations=read.csv("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/Ann_IlluminaHumanMethylationEPICanno.ilm10b2.hg19.csv")
DE_cpgs_stage=read.csv("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/_2020-07-22_DMR_CpG_clusters.csv")
DE_cpgs_stage=read.csv("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/_2020-07-22_DMR_CpG_clusters.csv")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###ANALYSIS###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#format and merge DE ERVs with ERV annotations
####for stage comparisons

#include only Tier 3 QC samples
DE_ervs_stage_all=as.data.frame(filter(DE_ervs_stage, tier == "Tier 3"))
rownames(DE_ervs_stage_all)=DE_ervs_stage_all$transcript

#keep only following columns for DE ERVs
keeps=c("logFC", "FDR", "transcript")
telescope_data_keep=DE_ervs_stage_all[keeps]

#keep only chr, chromStart, chromEnd, strand for ERV annotation file
erv_annotations_keep=data.frame(erv_annotations$V1, erv_annotations$V2, erv_annotations$V3, erv_annotations$V6)
rownames(erv_annotations_keep)=erv_annotations$V10
colnames(erv_annotations_keep)=c("chr", "chromStart", "chromEnd", "strand")

merged_stage_telescope = merge(as.data.frame(erv_annotations_keep), as.data.frame(telescope_data_keep), by='row.names', all=TRUE)
#write.csv(merged_stage_telescope, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "merged_stage_telescope.csv", sep="_"), quote=F, row.names=F)


####for cluster comparisons

#include only Tier 3 QC samples
DE_ervs_cluster_all=as.data.frame(filter(DE_ervs_cluster, tier == "Tier 3"))
rownames(DE_ervs_cluster_all)=DE_ervs_cluster_all$transcript

#keep only following columns for DE ERVs
keeps=c("logFC", "FDR", "transcript")
telescope_data_keep=DE_ervs_cluster_all[keeps]

#keep only chr, chromStart, chromEnd, strand for ERV annotation file
erv_annotations_keep=data.frame(erv_annotations$V1, erv_annotations$V2, erv_annotations$V3, erv_annotations$V6)
rownames(erv_annotations_keep)=erv_annotations$V10
colnames(erv_annotations_keep)=c("chr", "chromStart", "chromEnd", "strand")

merged_cluster_telescope = merge(as.data.frame(erv_annotations_keep), as.data.frame(telescope_data_keep), by='row.names', all=TRUE)
#write.csv(merged_cluster_telescope, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "merged_cluster_telescope.csv", sep="_"), quote=F, row.names=F)

###############################################################################################################
#format and merge DM CpGs with CpG annotations
####for both comparisons

keeps=c("X", "logFC", "adj.P.Val")
cpg_data_keep_stage=DE_cpgs_stage[keeps]
cpg_data_keep_cluster=DE_cpgs_cluster[keeps]

rownames(cpg_data_keep_stage)=cpg_data_keep_stage$X
rownames(cpg_data_keep_cluster)=cpg_data_keep_cluster$X

cpg_annotations_keep=data.frame(cpg_annotations$V1, cpg_annotations$chr, cpg_annotations$pos, cpg_annotations$strand, cpg_annotations$UCSC_RefGene_Name)
rownames(cpg_annotations_keep)=cpg_annotations$V1
colnames(cpg_annotations_keep)=c("cpg", "chr", "chromStart", "strand", "UCSC_RefGene_Name)")

cpg_annotations_stage=filter(cpg_annotations_keep, cpg %in% cpg_data_keep_stage$X)
cpg_annotations_cluster=filter(cpg_annotations_keep, cpg %in% cpg_data_keep_cluster$X)

merged_stage_cpg = merge(as.data.frame(cpg_annotations_stage), as.data.frame(cpg_data_keep_stage), by='row.names', all=TRUE)
merged_cluster_cpg = merge(as.data.frame(cpg_annotations_cluster), as.data.frame(cpg_data_keep_cluster), by='row.names', all=TRUE)


#write.csv(merged_stage_cpg, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "merged_stage_cpg.csv", sep="_"), quote=F)
#write.csv(merged_cluster_cpg, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "merged_cluster_cpg.csv", sep="_"), quote=F)

###############################################################################################################
####format merged files to make GRanges object
de_erv_stage=merged_stage_telescope
dm_cpg_stage=merged_stage_cpg

de_erv_cluster=merged_cluster_telescope
dm_cpg_cluster=merged_cluster_cpg


#remove NAs, filter FDR for DE ERVs
de_erv_stage$Row.names=NULL
de_erv_stage=na.omit(de_erv_stage)
de_erv_stage=filter(de_erv_stage, FDR < 0.05)
colnames(de_erv_stage)[5]="ERV_logFC"

de_erv_cluster$Row.names=NULL
de_erv_cluster=na.omit(de_erv_cluster)
de_erv_cluster=filter(de_erv_cluster, FDR < 0.05)
colnames(de_erv_cluster)[5]="ERV_logFC"

write.csv(de_erv_stage, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "de_erv_stage", sep="_"), quote=F, row.names=F)
write.csv(de_erv_cluster, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "de_erv_cluster", sep="_"), quote=F, row.names=F)


#format de cpgs, end chromEnd
dm_cpg_stage$X.1=NULL
dm_cpg_stage$Row.names=NULL
dm_cpg_stage$X=NULL
#dm_cpg_stage=filter(dm_cpg_stage, adj.P.Val < 0.05)
  #not sure if I shoukd further filter these since already
    #filtered for p.val < 0.01 in DM analysis
colnames(dm_cpg_stage)[6]="CpG_logFC"

dm_cpg_stage$chromEnd=""
for(i in 1:nrow(dm_cpg_stage)){
  dm_cpg_stage$chromEnd[i] = dm_cpg_stage$chromStart[i]
}

dm_cpg_cluster$X.1=NULL
dm_cpg_cluster$Row.names=NULL
dm_cpg_cluster$X=NULL
#dm_cpg_cluster=filter(dm_cpg_cluster, adj.P.Val < 0.05)
#not sure if I shoukd further filter these since already
  #filtered for p.val < 0.01 in DM analysis
colnames(dm_cpg_cluster)[6]="ERV_logFC"

dm_cpg_cluster$chromEnd=""
for(i in 1:nrow(dm_cpg_cluster)){
  dm_cpg_cluster$chromEnd[i] = dm_cpg_cluster$chromStart[i]
}

write.csv(dm_cpg_stage, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "dm_cpg_stage", sep="_"), quote=F, row.names=F)
write.csv(dm_cpg_cluster, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "dm_cpg_cluster", sep="_"), quote=F, row.names=F)
