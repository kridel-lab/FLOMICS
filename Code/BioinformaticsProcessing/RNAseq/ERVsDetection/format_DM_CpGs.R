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

DE_ervs=fread("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope/_2020-08-12_de_ervs_STAGE_CLUSTER.csv")
erv_annotations = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/hg19_clean_repeatmasker_annotations.bed")

sample_info = as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.xlsx"))

cpg_annotations=fread("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/Ann_IlluminaHumanMethylationEPICanno.ilm10b2.hg19.csv")
DE_cpgs_stage=read.csv("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/_2020-07-22_DMR_CpG_stages.csv")
DE_cpgs_cluster=read.csv("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/_2020-07-22_DMR_CpG_clusters.csv")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###ANALYSIS###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#format and merge DE ERVs with ERV annotations
####for stage comparisons

#include all tiers and contrasts

#keep only following columns for DE ERVs, filter for sig FDR
de_keep=DE_ervs %>% select(logFC, FDR, transcript,contrast,tier)
dim(de_keep)
#[1] 29050     5
de_keep=filter(de_keep, FDR < 0.05)
dim(de_keep)
#[1] 5195    5

#keep only chr, chromStart, chromEnd, transcript for ERV annotation file
erv_annotations_keep=data.frame(erv_annotations$V1, erv_annotations$V2, erv_annotations$V3, erv_annotations$V10)
#rownames(erv_annotations_keep)=erv_annotations$V10
colnames(erv_annotations_keep)=c("chr", "chromStart", "chromEnd", "transcript")

erv_annotations_keep_filt=filter(erv_annotations_keep, transcript %in% de_keep$transcript)
dim(erv_annotations_keep_filt)
#[1] 2334    5

de_merged = merge(de_keep, erv_annotations_keep_filt, by='transcript', all=TRUE)
dim(de_merged)
#[1] 5195    9
write.csv(de_merged, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "de_ervs_with_annotations.csv", sep="_"), quote=F, row.names=F)


###############################################################################################################
#format and merge DM CpGs with CpG annotations
####for both comparisons

#include all tiers
DE_cpgs_stage
#dim(DE_cpgs_stage)
#[1] 131959      7

DE_cpgs_cluster
#> dim(DE_cpgs_cluster)
#[1] 58699     7

DE_cpgs_stage$contrast="ADVANCED_LIMITED"
DE_cpgs_cluster$contrast="Cluster1_Cluster2"

DE_cpgs=as.data.frame(rbind(DE_cpgs_stage,DE_cpgs_cluster))
dim(DE_cpgs)
#[1] 190658      8
DE_cpgs=filter(DE_cpgs, adj.P.Val < 0.05)
#dim(DE_cpgs)
#[1] 190658      8
DE_cpgs=DE_cpgs %>% select(X,logFC,adj.P.Val, contrast)
colnames(DE_cpgs)=c("cpg", "DM_logFC", "adj.P.Val", "contrast")

cpg_annotations_keep=cpg_annotations %>% select(V1, chr, pos, UCSC_RefGene_Name)
#dim(cpg_annotations_keep)
#[1] 866836      4
colnames(cpg_annotations_keep)=c("cpg", "chr", "chromStart", "UCSC_RefGene_Name")

cpg_annotations_keep_filt=filter(cpg_annotations_keep, cpg %in% DE_cpgs$cpg)
dim(cpg_annotations_keep_filt)
#[1] 164479      4
    
cpg_merged = merge(DE_cpgs, cpg_annotations_keep_filt, by='cpg', all=TRUE)
dim(cpg_merged)
#[1] 190658      7

write.csv(cpg_merged, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "dm_cpgs_with_annotations.csv", sep="_"), quote=F, row.names=F)
###############################################################################################################
####format merged files to make GRanges object

#de_merged
#cpg_merged

#remove NAs, filter FDR for DE ERVs
colnames(de_merged)[2]="ERV_logFC"

#format de cpgs, add end chromEnd
cpg_merged$chromEnd=cpg_merged$chromStart

write.csv(de_merged, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "DE_ERVS_FOR_HITMATRIX.csv", sep="_"), quote=F, row.names=F)
write.csv(cpg_merged,paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "DM_CPGS_FOR_HITMATRIX.csv", sep="_"), quote=F, row.names=F)
