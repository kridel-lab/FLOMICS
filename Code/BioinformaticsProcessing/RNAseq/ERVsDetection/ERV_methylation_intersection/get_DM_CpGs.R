#----------------------------------------------------------------------------------
#get_DM_CpGs.R
#----------------------------------------------------------------------------------

#Sarah Russell
#Date:July 20th, 2020

date = Sys.Date()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###LIBRARIES###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
mvaluematrix=read.csv("2_MvalueMatrix_updSamples_Ordered_T1_FilteredProbes.csv", row.names=1)
sample_info = as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.xlsx"))
  #contains information on samples used in RNA-seq, Telescope analysis
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###ANALYSIS###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Differential methylation between ADVANCED-LIMITED

# Select advanced and limited stage cases only from the FL cases (all Teir 3 QC)
telescope_stage_samples=as.data.frame(filter(sample_info, STAGE %in% c("ADVANCED", "LIMITED")))
stage_status = factor(telescope_stage_samples$STAGE, levels=c("ADVANCED", "LIMITED"))

#construct design matrix
design_epic2 <- model.matrix (~ 0 + stage_status)
colnames(design_epic2) <- levels(stage_status)

#filter and double check only FL ADV_LIM samples included in mvaluematrix
samples=telescope_stage_samples$SAMPLE_ID
filtered_mval=mvaluematrix[ , samples]
print(colnames(filtered_mval)==telescope_stage_samples$SAMPLE_ID)

#perform DM CpG analysis
fit.reduced <- lmFit(filtered_mval,design_epic2)
contrasts <- makeContrasts("ADVANCED-LIMITED", levels = design_epic2)
fit.reduced <- eBayes(contrasts.fit(fit.reduced, contrasts))

#fit.reduced_MultipleTesting <- decideTests(fit.reduced, adjust.method = "BH", p.value = 0.01)
#Significant_Probes <- rownames(fit.reduced_MultipleTesting)[which(fit.reduced_MultipleTesting != 0)]
  #didn't include above two lines in analysis since printing summary results below
  #with adjust.method="BH"

print(Summary_Results <- summary(decideTests(fit.reduced, adjust.method = "BH", p.value = 0.01)))
  #should I keep p.value at 0.01? seems to be used in Anjali's analysis
DMR_CpG_stages = data.frame(topTable(fit.reduced, number = nrow(filtered_mval), adjust.method = "BH", p.value = 0.01))

#save results and upload to OneDrive
write.csv(DMR_CpG_stages, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/",date, "DMR_CpG_stages.csv", sep="_"), quote=F)

#################################################################################################################

# Differential methylation between Cluster1-Cluster2
# Select cases from sample_info file, no need to filter (all Tier 3 QC)

#contruct design matrix
cluster_status = factor(sample_info$Cluster, levels=c("1","2"))
design_epic2 <- model.matrix (~ 0 + cluster_status)
colnames(design_epic2) <- c("C1","C2")

#filter and double check only FL C1_C2 samples included in mvaluematrix
samples=sample_info$SAMPLE_ID
filtered_mval=mvaluematrix[ , samples]
print(colnames(filtered_mval)==sample_info$SAMPLE_ID)

#perform DM CpG analysis
fit.reduced <- lmFit(filtered_mval,design_epic2)
contrasts <- makeContrasts("C1-C2", levels = design_epic2)
fit.reduced <- eBayes(contrasts.fit(fit.reduced, contrasts))

#fit.reduced_MultipleTesting <- decideTests(fit.reduced, adjust.method = "BH", p.value = 0.01)
#Significant_Probes <- rownames(fit.reduced_MultipleTesting)[which(fit.reduced_MultipleTesting != 0)]
  #didn't include above two lines in analysis since printing summary results below
  #with adjust.method="BH"

print(Summary_Results <- summary(decideTests(fit.reduced, adjust.method = "BH", p.value = 0.01)))
DMR_CpG_clusters = data.frame(topTable(fit.reduced, number=nrow(filtered_mval), adjust.method = "BH", p.value = 0.01))

#save results and upload to OneDrive
write.csv(DMR_CpG_clusters, paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/diffmethylation/", date, "DMR_CpG_clusters.csv", sep="_"), quote=F)
