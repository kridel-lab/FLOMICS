#----------------------------------------------------------------------
#library_size_vs_mean_coverage.R
#----------------------------------------------------------------------

#Kisaev

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/Users/kisaev/UHN/kridel-lab - Documents (1)/FLOMICS/DNAseq/Targeted_Seq_Coverage/08-07-2020")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
"plyr", "data.table",
"ggrepel", "stringr", "ggpubr")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#get data
all_lib_sizes = fread("FL_Aug2020_131_samples_library_sizes.txt")

#coverage file
coverage_all = fread("2020-08-07_picard_tools_coverage_summary_targets_DNA_sequencing.csv")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#get mean on target coverage
# +++++++++++++++++++++++++++
target_cov = as.data.table(coverage_all %>%
  group_by(gene, Library) %>%
  dplyr::summarize(mean_cov = mean(mean_coverage)))
target_cov = target_cov[order(-mean_cov)]
target_cov = merge(target_cov, all_lib_sizes, by="Library")
target_cov$Library_Reads_Total_Million =  target_cov$Library_Reads_Total/1000000

# Basic plot
# +++++++++++++++++++++++++++
g1 = ggscatter(target_cov, x = "Library_Reads_Total_Million", y = "mean_cov",
   color = "black", shape = 21, size = 3)+theme_bw()+ylab("Mean coverage across genes per sample")

g2 = ggscatter(filter(target_cov, Library_Reads_Total_Million < 150), x = "Library_Reads_Total_Million", y = "mean_cov",
  color = "black", shape = 21, size = 3)+theme_bw()+ylab("Mean coverage across genes per sample")

pdf("Library_sizes_vs_mean_gene_coverage_per_sample.pdf")
print(g1)
print(g2)
dev.off()

#get mean on target coverage
# +++++++++++++++++++++++++++
sample_cov = as.data.table(coverage_all %>%
  group_by(External_identifier, TYPE, STAGE, Library) %>%
  dplyr::summarize(mean_cov = mean(mean_coverage)))
sample_cov = sample_cov[order(-mean_cov)]
sample_cov = merge(sample_cov, all_lib_sizes, by="Library")
sample_cov$Library_Reads_Total_Million =  sample_cov$Library_Reads_Total/1000000

# Basic plot
# +++++++++++++++++++++++++++
g1 = ggscatter(target_cov, x = "Library_Reads_Total_Million", y = "mean_cov",
   color = "black", shape = 21, size = 3)+theme_bw()+ylab("Mean coverage per sample")

g2 = ggscatter(filter(target_cov, Library_Reads_Total_Million < 150), x = "Library_Reads_Total_Million", y = "mean_cov",
  color = "black", shape = 21, size = 3)+theme_bw()+ylab("Mean coverage per sample")

pdf("Library_sizes_vs_mean_sample_coverage_per_sample.pdf")
print(g1)
print(g2)
dev.off()

#Summary of Library sizes
# +++++++++++++++++++++++++++
pdf("Library_sizes_distribution.pdf")
all_lib_sizes$Library_Reads_Total_Million =  all_lib_sizes$Library_Reads_Total/1000000
all_lib_sizes$All_Samples = "AllSamples"
print(ggboxplot(all_lib_sizes, x="All_Samples", y="Library_Reads_Total_Million")+theme_bw())
dev.off()

#add library size column to main file
coverage_all = merge(coverage_all, all_lib_sizes, by="Library")
write.csv(coverage_all, file="2021-06-24_picard_tools_coverage_summary_targets_DNA_sequencing_w_Library_Size.csv")
