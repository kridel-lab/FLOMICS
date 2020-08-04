#----------------------------------------------------------------------------------
#findOverlaps_ERV_CpG.R
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
library(GenomicRanges)
library(ggpubr)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###DATA###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/")

de_erv_stage=as.data.table(read.csv("_2020-07-21_de_erv_stage.csv"))
de_cpg_stage=as.data.table(read.csv("_2020-07-22_dm_cpg_stage"))
#de_erv_cluster=as.data.table(read.csv("_2020-07-22_de_erv_cluster.csv"))
#de_cpg_cluster=as.data.table(read.csv("_2020-07-22_mdm_cpg_cluster.csv"))
#only performed for stage comparisons

erv_coords = makeGRangesFromDataFrame(de_erv_stage, keep.extra.columns=TRUE)
cpg_coords = makeGRangesFromDataFrame(de_cpg_stage, keep.extra.columns=TRUE)

hits <- findOverlaps(vcf_dat_coords, all_genes_coords, ignore.strand=TRUE)
hits_overlap = cbind(as.data.table(vcf_dat[queryHits(hits),]), as.data.table(all_genes_coords)[subjectHits(hits),])
print(head(hits_overlap))

hits <- findOverlaps(erv_coords, cpg_coords, ignore.strand=TRUE)
hits_overlap = cbind(as.data.table(erv_coords[queryHits(hits),]), as.data.table(cpg_coords)[subjectHits(hits),])

write.csv(hits_overlap, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "hits_overlap_stage.csv",), quote=F, row.names=F)

pdf("/cluster/home/srussell/test_J22.pdf", width=5, height=4)
g=ggscatter(hits_overlap, x = "ERV_logFC", y = "CpG_logFC",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
          )+
          geom_vline(xintercept = 0) +
          geom_hline(yintercept = 0) +
  stat_cor(method = "pearson", label.x = 1, label.y = 1.5)  # Add correlation coefficient
  ggpar(g, y.lim=c(-2,2),x.lim=c(-6,6))
dev.off()
