setwd("/Users/kisaev/FLOMICS_Anjali/FLOMICS_Anajli/Data/qc")

library(data.table)
library(plyr)
library(dplyr)
library(ggpubr)
library(EnvStats)
library(magrittr)

files=list.files(pattern=".captureKit.unique.tsv")
get_cov = function(file){
  print(file)
  f = fread(file)
  f=melt(f)
  f$pats = ""
  f$pats = sapply(as.character(f$variable), function(x){unlist(strsplit(x, "_"))[3]})
  f$type_measure = ""
  f$type_measure = sapply(as.character(f$variable), function(x){paste(unlist(strsplit(x, "_"))[1:2], collapse="_")})
  f$variable = NULL
  return(f)
}

all_covs = as.data.table(ldply(llply(files, get_cov)))

#are there regions that are never covered in anyone?
all_covs$region_pres = ""
all_covs$region_pres[all_covs$value == 0] = "NoCov"
all_covs$region_pres[all_covs$value >= 100] = "HighCov"
all_covs$region_pres[all_covs$region_pres == ""] = "SomeCov"
mean_cov = filter(all_covs, type_measure == "Mean_Coverage")
write.table(mean_cov, file="BC_TargetedSeq_summary_coverage_regions.csv", quote=F, row.names=F, sep=";")

#---------------------
#PATIENT BASED SUMMARY#
#---------------------

#summarize type of coverage per sample
pats_sum = as.data.table(table(all_covs$pats, all_covs$region_pres, all_covs$type_measure))
colnames(pats_sum) = c("Patient", "Type_Coverage", "Measure_Cov", "Number_regions")
pats_sum$Type_Coverage = factor(pats_sum$Type_Coverage, levels = c("NoCov", "SomeCov", "HighCov"))
pats_sum = pats_sum[order(-Number_regions)]
pats_sum$Patient = factor(pats_sum$Patient, levels=unique(pats_sum$Patient))

plots = pats_sum %>% group_by(Measure_Cov) %>% do(plots=ggplot(data=.) +
                                    aes(x=Patient, y=Type_Coverage) +
                                      geom_tile(aes(fill = Number_regions), colour = "grey50") + 
                                      rotate_x_text(90) + scale_fill_gradient2(low = "black", mid = "grey", 
                                                                               midpoint = 200, high = "orange", 
                                                                               na.value="transparent")
                                    + ggtitle(unique(.$Measure_Cov)))

pdf("all_libs_coverage_summary.pdf", width=16, height=7)
for(i in 1:4){
  print(plots$plots[i])
}
dev.off()

#---------------------
#REGION BASED SUMMARY#
#---------------------

#summarize type of coverage per sample
reg_sum = as.data.table(table(all_covs$Exon_Genomic_Coordinates, all_covs$region_pres, all_covs$type_measure))
colnames(reg_sum) = c("Region", "Type_Coverage", "Measure_Cov", "Number_regions")
reg_sum$Type_Coverage = factor(reg_sum$Type_Coverage, levels = c("NoCov", "SomeCov", "HighCov"))
reg_sum = reg_sum[order(-Number_regions)]
reg_sum$Region = factor(reg_sum$Region, levels=unique(reg_sum$Region))

#lowly covered regions = more than 20% of patients have poor coverage in regions 
filter(reg_sum, Type_Coverage == "NoCov", Number_regions >= 26)









