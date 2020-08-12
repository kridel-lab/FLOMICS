#----------------------------------------------------------------------------------
#telescope_erv_visualization.R
#----------------------------------------------------------------------------------

#Sarah Russell
#Date:July 16th, 2020

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###SOURCE###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

options(stringsAsFactors=F)
date = Sys.Date()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###LIBRARIES###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "tidyr", "mclust", "data.table", "plyr",
              "ggrepel", "stringr")
lapply(packages, require, character.only = TRUE)

## loading packages
library(data.table)
library(dplyr)
library(plyr)
library(edgeR)
library(tidyr)
library(ggpubr)
library(readxl)
library(tidyverse)
library(rstatix)
library(datarium)


#WD
setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope/")

#----------------------------------------------------------------------------------
#DATA
#----------------------------------------------------------------------------------
ervs_de=fread("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope/_2020-08-12_de_ervs_STAGE_CLUSTER.csv")
dim(ervs_de)
#[1] 29050      8

ervs_de_sig = filter(ervs_de, FDR < 0.05)
dim(ervs_de_sig)
#[1] 5195    8

#label neg/pos logFC (FDR < 0.05)
ervs_de_sig$diff_exp = ""
ervs_de_sig$diff_exp[ervs_de_sig$logFC >0 ] ="Upregulated"
ervs_de_sig$diff_exp[ervs_de_sig$logFC <0 ] ="Downregulated"
#----------------------------------------------------------------------------------
#ANALYSIS
#----------------------------------------------------------------------------------
#Generate density plot comparing FC and FDR of all ERVs that were evaluated

# :::::::::::::::::::::::::::::::::::::::::::::::::::

#Only DE ERVs
pdf("/cluster/home/srussell/DE_ERVs_totaldensities_FL.pdf", width=7, height=4)
ggdensity(ervs_de_sig[ervs_de_sig$contrast=="ADVANCED_LIMITED",], x = "logFC",
          add = "mean", rug = TRUE, alpha=0.1,
          color = "diff_exp", fill = "diff_exp") + ggtitle("Fold changes all QC tiers for DE ERVs in ADV-LIM (FDR < 0.05), n=798")+theme_bw()
dev.off()

pdf("/cluster/home/srussell/DE_ERVs_totaldensities_CLUS.pdf", width=7, height=4)
ggdensity(ervs_de_sig[ervs_de_sig$contrast=="Cluster1_Cluster2",], x = "logFC",
          add = "mean", rug = TRUE, alpha=0.1,
          color = "diff_exp", fill = "diff_exp") + ggtitle("Fold changes all QC tiers for DE ERVs in C1-C2 (FDR < 0.05), n=4397")+theme_bw()
dev.off()

#make comparison matrix
ervs_de_dat = as.data.table(ervs_de_sig)
barplot_sum = as.data.table(table(ervs_de_dat$contrast, ervs_de_dat$tier, ervs_de_dat$diff_exp))
colnames(barplot_sum) = c("contrast", "tier", "Expression", "num_ervs")
barplot_sum$Expression = factor(barplot_sum$Expression, levels=c("Upregulated", "Downregulated"))
barplot_sum$num_ervs[barplot_sum$Expression == "Downregulated"] = -(barplot_sum$num_ervs[barplot_sum$Expression == "Downregulated"])
barplot_sum$N = abs(barplot_sum$num_ervs)

z = which(str_detect(barplot_sum$contrast, "ADVANCED_LIMITED"))
barplot_sum$contrast[z] = "Adv_vs_Lim"
z = which(str_detect(barplot_sum$contrast, "Cluster1_Cluster2"))
barplot_sum$contrast[z] = "C1_vs_C2"


pdf("/cluster/home/srussell/Down_Up_regulated_ERVs.pdf", width=5, height=4)
g=ggbarplot(barplot_sum, x="contrast", y="num_ervs", label = barplot_sum$N , lab.size=2, fill="Expression", facet.by=c("tier"),
  palette=c("red", "blue"))+
theme_bw()+
theme(text = element_text(size=10)) + ylab("Number of ERVs") + xlab("Contrast")
ggpar(g, legend="bottom", main="ERV Expression in ADV v. LIM/ C1 v. C2") + geom_hline(yintercept=0)
dev.off()

# :::::::::::::::::::::::::::::::::::::::::::::::::::

#for all sig ERVs, make histo plot of FCs
ervs_de_histo = ervs_de_dat
ervs_de_histo$logFC = abs(ervs_de_histo$logFC)

plot_cons_stage = c("ADVANCED_LIMITED")
plot_cons_cluster= c("Cluster1_Cluster2")

plot_cons_stage = as.data.table(filter(ervs_de_histo, contrast %in% plot_cons_stage))
plot_cons_cluster = as.data.table(filter(ervs_de_histo, contrast %in% plot_cons_cluster))

plot_cons_stage$contrast = factor(plot_cons_stage$contrast, levels=c("ADVANCED_LIMITED"))
plot_cons_cluster$contrast = factor(plot_cons_cluster$contrast, levels=c("Cluster1_Cluster2"))

plot_cons_stage$diff_exp = factor(plot_cons_stage$diff_exp, levels=c("Upregulated", "Downregulated"))
plot_cons_cluster$diff_exp = factor(plot_cons_cluster$diff_exp, levels=c("Upregulated", "Downregulated"))

plot_cons_stage$tier = factor(plot_cons_stage$tier, levels=c("tier_1", "tier_2", "tier_3"))
plot_cons_cluster$tier = factor(plot_cons_cluster$tier, levels=c("tier_1", "tier_2", "tier_3"))
#could probably put above in function


pdf("/cluster/home/srussell/stage_histo.pdf", width=7, height=4)
g=gghistogram(plot_cons_stage, x = "logFC",
   add = "median", alpha=0.5,
   color = "diff_exp",
   fill = "diff_exp", palette = c("red", "blue"),
   add_density = TRUE, ggtheme=theme_bw(),
   facet.by = "tier")
ggpar(g, legend="bottom", legend.title = "Expression", main="Distribution of ERV logFC in ADV v. LIM")+theme(text = element_text(size=12))+ylab("Count")+xlab("Absolute logFC")
dev.off()

pdf("/cluster/home/srussell/cluster_histo.pdf", width=7, height=4)
g=gghistogram(plot_cons_cluster, x = "logFC",
   add = "median", alpha=0.5,
   color = "diff_exp",
   fill = "diff_exp", palette = c("red", "blue"),
   add_density = TRUE, ggtheme=theme_bw(),
   facet.by = "tier")
ggpar(g, legend="bottom", legend.title = "Expression", main="Distribution of ERV logFC in C1 v. C2" )+theme(text = element_text(size=12))+ylab("Count")+xlab("Absolute logFC")
dev.off()
# :::::::::::::::::::::::::::::::::::::::::::::::::::
