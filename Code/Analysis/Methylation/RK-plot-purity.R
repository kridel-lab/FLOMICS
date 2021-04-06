
library(dplyr)
library(ggplot2)
library(ggpubr)

date <- Sys.Date()

SNF.clust <- read.csv("Cluster Labels/InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv") %>%
  select(ID, SNFClust = SNFClust10Feb2021)

purity <- read.csv("methylation/Purity_281probes_10Jan2020.csv") %>%
  left_join(SNF.clust[,c("ID", "SNFClust")], by = c("X" = "ID")) %>%
  filter(!is.na(SNFClust)) %>%
  mutate(purity = 100*purity) %>%
  select(Sample_ID = X, purity = purity, SNFClust)

purity %>%
  group_by(SNFClust) %>%
  summarize(mean = mean(purity), median = median(purity))

p <- purity %>%
  select(Sample_ID, 'purity (%)' = purity, 'SNF cluster' = SNFClust) %>%
  ggboxplot(x = "SNF cluster", y = "purity (%)", color = "SNF cluster",
            add = "jitter", palette = c("1" = "#4363d8", "2" = "#f58231")) +
  ggtitle("Purity estimated from methylation") +
  theme(legend.position = "none")
p + stat_compare_means(method = "t.test")
ggsave("img/purity_SNFClust.pdf", width = 3.5, height = 3.5)
