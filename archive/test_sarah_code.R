#----------------------------------------------------------------------
#Sarah Rusell
#Use MiXCR metrics and to plot survival statistics
#done locally
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------
date = Sys.Date()

library(data.table)
library(dplyr)
library(plyr)
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)
library(purrr)

#----------------------------------------------------------------------
#load data & filtering
#----------------------------------------------------------------------
setwd("/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/Analysis-Files/Mixcr")

# All FLOMICS samples included - load sample information
all.samples.RNAseq.FLOMICS <- fread("/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/metadata/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")

#survival data
clin_dat=as.data.frame(fread("/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/metadata/clinical_data_rcd11Aug2020.csv"))

tcr_metrics=as.data.frame(fread("/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/Analysis-Files/Mixcr/_2020-10-02_tcr_metrics.csv"))
bcr_metrics=as.data.frame(fread("/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/Analysis-Files/Mixcr/_2020-10-02_bcr_metrics.csv"))
div_met=as.data.frame(fread("/Users/kisaev/UHN/kridel-lab - Documents/FLOMICS/Analysis-Files/Mixcr/_2020-10-08_repetoire_diversity.csv"))
colnames(div_met)[1]="rna_seq_file_sample_ID"

div_met=merge(div_met, all.samples.RNAseq.FLOMICS %>%
                   select(LY_FL_ID,rna_seq_file_sample_ID),
                 by="rna_seq_file_sample_ID")

surv_dat=filter(clin_dat %>% select(LY_FL_ID,CODE_TTP,TTP),LY_FL_ID %in% div_met$LY_FL_ID)

cov_bcr=merge(bcr_metrics %>% select(sample,stage,total_clones,clone_abundance,unique_CDR3),
               div_met %>% select(sample=rna_seq_file_sample_ID,LY_FL_ID,contains("BCR")),
               by="sample") %>% filter(., LY_FL_ID %in% surv_dat$LY_FL_ID) %>%
               merge(., surv_dat, by="LY_FL_ID")

cov_bcr=merge(bcr_metrics %>% select(sample,stage,total_clones,clone_abundance,unique_CDR3),
              div_met %>% select(sample=rna_seq_file_sample_ID,LY_FL_ID,contains("BCR")),
              by="sample") %>% merge(surv_dat, ., by="LY_FL_ID")
cov_bcr_sep=list(advanced=filter(cov_bcr, stage == "ADVANCED"),
                 limited=filter(cov_bcr, stage == "LIMITED"))

cov_tcr=merge(tcr_metrics %>% select(sample,stage,total_clones,clone_abundance,unique_CDR3),
              div_met %>% select(sample=rna_seq_file_sample_ID,LY_FL_ID,contains("TCR")),
              by="sample") %>% merge(surv_dat, ., by="LY_FL_ID")
cov_tcr_sep=list(advanced=filter(cov_tcr, stage == "ADVANCED"),
                 limited=filter(cov_tcr, stage == "LIMITED"))

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------
tcr.covariate_names <- c('total_clones' = "number of unique clones",
                     'clone_abundance'= "total number of clones",
                     'unique_CDR3' = "unique CDR3 sequences",
                     'TCR.chao' = "chao diversity" ,
                     'TCR.true.diversity' = "true diversity",
                     'TCR.gini' = "gini coefficient" ,
                     'TCR.gini.simp' = "gini-simpson diversity" ,
                     'TCR.d50' = "d50 diversity",
                     'TCR.shannon' =  "shannon diversity",
                     'TCR.invsimpson'= "inverse simpson diversity",
                     'TCR.pielou' = "pielous eveness index")
bcr.covariate_names <- c('total_clones' = "number of unique clones",
                         'clone_abundance'= "total number of clones",
                         'unique_CDR3' = "unique CDR3 sequences",
                         'BCR.chao' = "chao diversity" ,
                         'BCR.true.diversity' = "true diversity",
                         'BCR.gini' = "gini coefficient" ,
                         'BCR.gini.simp' = "gini-simpson diversity" ,
                         'BCR.d50' = "d50 diversity",
                         'BCR.shannon' =  "shannon diversity",
                         'BCR.invsimpson'= "inverse simpson diversity",
                         'BCR.pielou' = "pielous eveness index")

get_bcr_forest = function(dat, comparison){
map(vars(total_clones, clone_abundance, unique_CDR3,BCR.chao,BCR.true.diversity,
         BCR.gini, BCR.gini.simp, BCR.d50,BCR.shannon, BCR.invsimpson,
         BCR.pielou), function(by)
         {
           analyse_multivariate(dat,
                                vars(TTP, CODE_TTP),
                                covariates = list(by), # covariates expects a list
                                covariate_name_dict = bcr.covariate_names)
         }) %>%
  forest_plot(factor_labeller = bcr.covariate_names,
              endpoint_labeller = c(time = "TTP"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor", "n"),
              ggtheme = ggplot2::theme_bw(base_size = 10),
              HR_x_breaks = c(0, 0.25, 0.5, 1, 2, 5, 10, 20),
              title=paste("Univariate Cox regression analysis BCR metrics and TPP",comparison,sep=" "))
  ggsave(filename = paste(date,comparison,"bcr_TPP_coxuni.pdf",sep="_"), width=10, height=4)
}

get_bcr_forest(cov_bcr_sep$advanced,"advanced")
get_bcr_forest(cov_bcr_sep$limited,"limited")


get_tcr_forest = function(dat, comparison){
  map(vars(total_clones, clone_abundance, unique_CDR3,TCR.chao,TCR.true.diversity,
           TCR.gini, TCR.gini.simp, TCR.d50,TCR.shannon, TCR.invsimpson,
           TCR.pielou), function(by)
           {
             analyse_multivariate(dat,
                                  vars(TTP, CODE_TTP),
                                  covariates = list(by), # covariates expects a list
                                  covariate_name_dict = tcr.covariate_names)
           }) %>%
    forest_plot(factor_labeller = tcr.covariate_names,
                endpoint_labeller = c(time = "TTP"),
                orderer = ~order(HR),
                labels_displayed = c("endpoint", "factor", "n"),
                ggtheme = ggplot2::theme_bw(base_size = 10),
                HR_x_breaks = c(0, 0.25, 0.5, 1, 2, 5, 10, 20),
                title=paste("Univariate Cox regression analysis TCR metrics and TPP",comparison,sep=" "))
  ggsave(filename = paste(date,comparison,"tcr_TPP_coxuni.pdf",sep="_"), width=10, height=4)
}

get_tcr_forest(cov_tcr_sep$advanced,"advanced")
get_tcr_forest(cov_tcr_sep$limited,"limited")
