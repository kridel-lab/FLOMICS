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
all_de_stages=fread("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope/_2020-07-16_all_DE_ervs_stage.csv")
all_de_clusters=fread("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope/_2020-07-16_all_DE_ervs_cluster.csv")
ervs_de=rbind(all_de_stages, all_de_clusters)

ervs_de_sig = unique(filter(ervs_de, FDR < 0.05)$transcript)
#get all significant ERVs (FDR < 0.05)
ervs_de$diff_exp = ""
ervs_de$diff_exp[ervs_de$logFC >0 ] ="Upregulated"
ervs_de$diff_exp[ervs_de$logFC <0 ] ="Downregulated"
#----------------------------------------------------------------------------------
#ANALYSIS
#----------------------------------------------------------------------------------
#Generate density plot comparing FC and FDR of all ERVs that were evaluated

# :::::::::::::::::::::::::::::::::::::::::::::::::::

#Only DE ERVs
ggdensity(filter(ervs_de, transcript_id %in% ervs_de_sig), x = "logFC",
          add = "mean", rug = TRUE, alpha=0.1,
          color = "diff_exp", fill = "diff_exp") + ggtitle("Fold changes across QC tiers for DE ERVs (FDR < 0.05), n=X")+theme_bw()
dev.off()

#make comparison matrix
ervs_de_dat = as.data.table(filter(ervs_de, FDR < 0.05))
barplot_sum = as.data.table(table(ervs_de_dat$contrast, ervs_de_dat$tier, ervs_de_dat$diff_exp))
colnames(barplot_sum) = c("contrast", "tier", "Expression", "num_ervs")
barplot_sum$Expression = factor(barplot_sum$Expression, levels=c("Upregulated", "Downregulated"))
barplot_sum$num_ervs[barplot_sum$Expression == "Downregulated"] = -(barplot_sum$num_ervs[barplot_sum$Expression == "Downregulated"])
barplot_sum$N = abs(barplot_sum$num_ervs)

z = which(str_detect(barplot_sum$contrast, "ADVANCED_LIMITED"))
barplot_sum$contrast[z] = "Adv vs Lim"
z = which(str_detect(barplot_sum$contrast, "Cluster1_Cluster2"))
barplot_sum$contrast[z] = "C1 vs C2"


pdf("Down_Up_regulated_ERVs.pdf", width=5, height=4)
g=ggbarplot(barplot_sum, x="contrast", y="num_ervs", label = barplot_sum$N , lab.size=2, fill="Expression", facet.by=c("tier"),
  palette=c("red", "blue"))+
theme_bw()+
theme(text = element_text(size=10)) + ylab("Number of ERVs") + xlab("Contrast")
ggpar(g, legend="bottom") + geom_hline(yintercept=0)
dev.off()

#how did you get the N label for downregulated below the x axis?
# :::::::::::::::::::::::::::::::::::::::::::::::::::
#having an issue getting the distribution of this graph right.. any suggestions?

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

plot_cons_stage$tier = factor(plot_cons_stage$tier, levels=c("Tier 1", "Tier 2", "Tier 3"))
plot_cons_cluster$tier = factor(plot_cons_cluster$tier, levels=c("Tier 1", "Tier 2", "Tier 3"))
#could probably put above in function


pdf("/cluster/home/srussell/stage_histo.pdf", width=7, height=4)
g=gghistogram(plot_cons_stage, x = "logFC",
   add = "median", alpha=0.5,
   color = "diff_exp",
   fill = "diff_exp", palette = c("red", "blue"),
   add_density = TRUE, ggtheme=theme_bw(),
   facet.by = "tier")+xlim(0,10)
ggpar(g, legend="bottom")+theme(text = element_text(size=12))+ylab("Count")+xlab("logFC")
dev.off()

pdf("/cluster/home/srussell/cluster_histo.pdf", width=7, height=4)
g=gghistogram(plot_cons_cluster, x = "logFC",
   add = "median", alpha=0.5,
   color = "diff_exp",
   fill = "diff_exp", palette = c("red", "blue"),
   add_density = TRUE, ggtheme=theme_bw(),
   facet.by = "tier")+xlim(0,10)
ggpar(g, legend="bottom")+theme(text = element_text(size=12))+ylab("Count")+xlab("logFC")
dev.off()
# :::::::::::::::::::::::::::::::::::::::::::::::::::
