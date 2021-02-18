
packages <- c("dplyr", "tidyr", "BayesianFirstAid", "cowplot")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

setwd("~/github/FLOMICS/")

SNF.clust <- read.csv("Cluster Labels/InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv") %>%
  select(ID, SNFClust = SNFClust10Feb2021)

table(SNF.clust$SNFClust) # C1, n = 56; C2, n = 45

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
  scale_fill_manual(values = c("#4363D8", "#F58231")) +
  theme(axis.text.x = element_text(face = "italic", angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.box = "horizontal") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
  
plot2 <- Bayesian.prop.test.mutations %>%
  select(Hugo_Symbol, group.diff, HDIlo, HDIup) %>%
  ggplot(aes(x = reorder(Hugo_Symbol, -group.diff), y = group.diff, ymin = HDIlo, ymax = HDIup)) +
  geom_linerange(position = position_dodge(0.75)) +
  geom_point(position = position_dodge(0.75)) +
  theme_bw() +
  scale_fill_manual(values = c("#4363D8", "#F58231")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.5)

plot_grid(plot2, plot1, ncol = 1, rel_heights = c(1, 1.5))

ggsave(paste0("img/", date, " Mutations_prob_group.pdf"), width = 8, height = 4)


