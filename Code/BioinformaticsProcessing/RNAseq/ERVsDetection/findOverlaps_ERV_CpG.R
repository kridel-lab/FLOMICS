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

ervs=fread("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/_2020-08-12_DE_ERVS_FOR_HITMATRIX.csv")
cpgs=fread("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/_2020-08-12_DM_CPGS_FOR_HITMATRIX.csv")
#> dim(ervs)
#[1] 5195    8

#> dim(cpgs)
#[1] 190658      8

ervs_t3_FL=ervs[ervs$tier =="tier_3" & ervs$contrast =="ADVANCED_LIMITED",]
#> dim(ervs_t3_FL)
#[1] 381   8
ervs_t3_CL=ervs[ervs$tier =="tier_3" & ervs$contrast =="Cluster1_Cluster2",]
#dim(ervs_t3_CL)
#[1] 1700    8

cpgs_FL=cpgs[cpgs$contrast=="ADVANCED_LIMITED",]
#> dim(cpgs_FL)
#[1] 131959      8
cpgs_CL=cpgs[cpgs$contrast=="Cluster1_Cluster2",]
#> dim(cpgs_CL)
#[1] 58699     8

FL_erv_coords = makeGRangesFromDataFrame(ervs_t3_FL, keep.extra.columns=TRUE)
FL_cpg_coords = makeGRangesFromDataFrame(cpgs_FL, keep.extra.columns=TRUE)

CL_erv_coords = makeGRangesFromDataFrame(ervs_t3_CL, keep.extra.columns=TRUE)
CL_cpg_coords = makeGRangesFromDataFrame(cpgs_CL, keep.extra.columns=TRUE)

FL_hits <- findOverlaps(FL_erv_coords, FL_cpg_coords, ignore.strand=TRUE)
CL_hits <- findOverlaps(CL_erv_coords, CL_cpg_coords, ignore.strand=TRUE)

FL_hits_overlap = cbind(as.data.table(FL_erv_coords[queryHits(FL_hits),]), as.data.table(FL_cpg_coords)[subjectHits(FL_hits),])
print(head(FL_hits_overlap))
#> dim(FL_hits_overlap)
#[1]  4 20

CL_hits_overlap = cbind(as.data.table(CL_erv_coords[queryHits(CL_hits),]), as.data.table(CL_cpg_coords)[subjectHits(CL_hits),])
print(head(CL_hits_overlap))

#> dim(CL_hits_overlap)
#[1] 11 20

write.csv(FL_hits_overlap, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "FL_tier3_hits_overlap.csv"), quote=F, row.names=F)
write.csv(CL_hits_overlap, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "CL_hits_overlap.csv"), quote=F, row.names=F)

combined_hits=as.data.frame(rbind(FL_hits_overlap,CL_hits_overlap))

pdf("/cluster/home/srussell/test.pdf", width=5, height=4)
g=ggscatter(combined_hits, x = "logFC", y = "DM_logFC",
          color = "contrast", palette = "jco",
          shape = "contrast"
          )+
          geom_vline(xintercept = 0) +
          geom_hline(yintercept = 0) +
  ggpar(g, y.lim=c(-1,1),x.lim=c(-6,6))
dev.off()
