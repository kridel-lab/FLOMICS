# Updated 3 Aug 2021
# 18 February 2021
# Edited by Anjali Silva, Developed by Robert Kridel 


packages <- c("dplyr", "tidyr", "BayesianFirstAid", "cowplot")
lapply(packages, require, character.only = TRUE)



date <- Sys.Date()

setwd("~/github/FLOMICS/")

SNF.clust <- read.csv("Cluster Labels/InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv") %>%
  select(ID, SNFClust = SNFClust10Feb2021)



table(SNF.clust$SNFClust) # C1, n = 56; C2, n = 45

# Here, doing C1 mutation / (C1 total patients) 
# NOT C1 mutation / (total mutations)

mutations <- read.csv("DNAseq/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl.csv") %>%
  pivot_longer(!Hugo_Symbol, names_to = "ID", values_to = "mutated") %>%
  right_join(SNF.clust) %>%
  filter(!is.na(SNFClust) & mutated == "1") %>%
  group_by(Hugo_Symbol, SNFClust) %>%
  dplyr::summarize(n = n()) %>%
  pivot_wider(names_from = SNFClust, values_from = n) %>%
  dplyr::rename(C1.mut = '1', C2.mut = '2') %>%
  replace(is.na(.), 0)

set.seed(123)
Bayesian.prop.test.mutations <- mutations %>%
  rowwise() %>%
  mutate(group.diff = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5,1],
         sd = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5,2],
         prob.C1.greater = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5,8],
         prob.C2.greater = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5,7],
         HDIlo = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5,5],
         HDIup = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5,6])

write.csv(Bayesian.prop.test.mutations, "DNAseq/Bayesian.prop.test.mutations.csv", row.names = FALSE)

plot1 <- mutations %>%
  pivot_longer(!Hugo_Symbol, names_to = "SNFClust", values_to = "nb_mutations") %>%
  mutate(SNFClust = ifelse(SNFClust == "C1.mut", "C1", "C2")) %>%
  mutate(perc.mutated = ifelse(SNFClust == "C1", nb_mutations/56, nb_mutations/45)) %>%
  left_join(Bayesian.prop.test.mutations[,c("Hugo_Symbol", "group.diff")]) %>%
  ggplot(aes(x = reorder(Hugo_Symbol, -group.diff), y = perc.mutated, fill = SNFClust)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  labs(y = "Percentage Mutated") +
  scale_fill_manual(values = c("#4363D8", "#F58231")) +
  theme(axis.text.x = element_text(face = "italic", 
                                   angle = 90, vjust = 0.5, 
                                   hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.box = "horizontal") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

plot2 <- Bayesian.prop.test.mutations %>%
  select(Hugo_Symbol, group.diff, HDIlo, HDIup) %>%
  ggplot(aes(x = reorder(Hugo_Symbol, -group.diff), 
             y = group.diff, ymin = HDIlo, ymax = HDIup)) +
  geom_linerange(position = position_dodge(0.75)) +
  geom_point(position = position_dodge(0.75)) +
  theme_bw() +
  labs(y = "Group Difference") +
  scale_fill_manual(values = c("#4363D8", "#F58231")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5)

plot_grid(plot2, plot1, ncol = 1, rel_heights = c(1, 1.5))

ggsave(paste0("img/", date, " Mutations_prob_group.pdf"), width = 8, height = 4)





require("devtools")
devtools::install_github("rasmusab/bayesian_first_aid")
library("BayesianFirstAid")

install.packages("rjags")
library(rjags)

install.packages("coda")
library("coda")

install.packages("MASS")
library(MASS)

install.packages("mnormt")
library(mnormt)

install.packages("stringr")
library(stringr)

install.packages("dplyr")
library(dplyr)

install.packages("tidyr")
library(tidyr)

install.packages("cowplot")
library(cowplot)

library(ggplot2)


# ERROR: lazy loading failed for package ‘BayesianFirstAid’
# * removing ‘/Library/Frameworks/R.framework/Versions/4.0/Resources/library/BayesianFirstAid’
# Error: Failed to install 'BayesianFirstAid' from GitHub:
#   (converted from warning) installation of package 
# ‘/var/folders/vs/wc7d5v9x7x3g40tqs3yd7j3c0000gp/T//RtmpRxgml0/file52a0146a7340/BayesianFirstAid_0.1.tar.gz’ had non-zero exit status

# change directory to bayesian_first_aid
# source("bayes_binom_test.R")
# source("bayes_cor_test.R")
# source("bayes_poisson_test.R")
# source("bayes_prop_test.R")
# source("bayes_t_test.R")
# source("BayesianFirstAid-package.R")
# source("generic_functions.R")
# source("utility_functions.R")

SNF.clust <- read.csv("InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv") %>%
  dplyr::select(ID, SNFClust = SNFClust10Feb2021)
table(SNF.clust$SNFClust) # C1, n = 56; C2, n = 45

TargetedDNAseqDirPath <- "~/Desktop/UHN/FLOMICS/TargetedSequencing"

mutations <- read.csv(paste0(TargetedDNAseqDirPath, 
                             "/Mutation_and_BA_matrices/mut.merged.df.T1.poor.cov.excl_10Feb2021.csv")) %>%
  pivot_longer(!Hugo_Symbol, names_to = "ID", values_to = "mutated") %>%
  right_join(SNF.clust) %>%
  filter(!is.na(SNFClust) & mutated == "1") %>%
  group_by(Hugo_Symbol, SNFClust) %>%
  dplyr::summarize(n = n()) %>%
  pivot_wider(names_from = SNFClust, values_from = n) %>%
  dplyr::rename(C1.mut = '1', C2.mut = '2') %>%
  replace(is.na(.), 0)


# Here, doing C1 mutation / (total mutations)
# NOT C1 mutation / (C1 total patients) 

set.seed(123)
Bayesian.prop.test.mutations <- mutations %>%
  rowwise() %>%
  mutate(group.diff = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5,1],
         sd = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5,2],
         prob.C1.greater = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5, 8],
         prob.C2.greater = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5, 7],
         HDIlo = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5,5],
         HDIup = bayes.prop.test(c(C1.mut, C2.mut), c(56, 45))$stats[5,6])

#write.csv(Bayesian.prop.test.mutations, "DNAseq/Bayesian.prop.test.mutations.csv", row.names = FALSE)


# barplot
plot1 <- mutations %>%
  mutate(sumRow = sum(C1.mut,C2.mut)) %>%
  pivot_longer(!c(Hugo_Symbol, sumRow), 
               names_to = "SNFClust", values_to = "nb_mutations") %>%
  mutate(SNFClust = ifelse(SNFClust == "C1.mut", "C1", "C2")) %>%
  mutate(perc.mutated = ifelse(SNFClust == "C1", 
                               nb_mutations/sumRow, nb_mutations/sumRow)) %>%
  left_join(Bayesian.prop.test.mutations[,c("Hugo_Symbol", "group.diff")]) %>%
  ggplot2::ggplot(aes(x = reorder(Hugo_Symbol, -group.diff), 
                      y = perc.mutated, fill = SNFClust)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  labs(y = "Percentage Mutated") +
  scale_fill_manual(values = c("#4363D8", "#F58231")) +
  theme(axis.text.x = element_text(face = "italic", 
                                   angle = 90, vjust = 0.5, 
                                   hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.box = "horizontal") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

plot2 <- Bayesian.prop.test.mutations %>%
         dplyr::select(Hugo_Symbol, group.diff, HDIlo, HDIup) %>%
         ggplot(aes(x = reorder(Hugo_Symbol, -group.diff), 
                    y = group.diff, ymin = HDIlo, ymax = HDIup)) +
         geom_linerange(position = position_dodge(0.75)) +
         geom_point(position = position_dodge(0.75)) +
         theme_bw() +
         labs(y = "Group Difference") +
         scale_fill_manual(values = c("#4363D8", "#F58231")) +
         theme(axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(),
               axis.title.x = element_blank(),
               axis.text.y = element_text(color = "black"),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) +
         geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5)

cowplot::plot_grid(plot2, plot1, ncol = 1, rel_heights = c(1, 1.5))

ggplot2::ggsave(paste0("img/", date, " Mutations_prob_group2.pdf"), width = 8, height = 4)
# [END]

