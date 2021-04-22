
packages <- c("enrichR", "dplyr")

lapply(packages, require, character.only = TRUE)

setwd("~/github/FLOMICS/")

hypometh_overexpr <- read.csv("RNAseq-Methylation-Intersect/41_MethylMix_n=121_HypoGenes235.csv") %>% .$HypoGenes
hypermeth_overexpr <- read.csv("RNAseq-Methylation-Intersect/41_MethylMix_n=121_HyperGenes135.csv") %>% .$HyperGenes

dbs <- c("BioPlanet_2019")

enrichr_results_hypo <- enrichr(hypometh_overexpr, databases = dbs)
enrichr_results_hypo <- rbind(enrichr_results_hypo[[1]])
colnames(enrichr_results_hypo) <- c("Term",	"Overlap",	"P-value",	"Adjusted P-value",	"Old P-value",	"Old Adjusted P-value",	"Odds Ratio",	"Combined Score", "Genes")
enrichr_results_hypo <- enrichr_results_hypo %>%
  mutate(nb.genes.overlap = stringr::str_extract(enrichr_results_hypo$Overlap, "[^/]+")) %>%
  mutate(nb.genes.overlap = as.numeric(nb.genes.overlap)) %>%
  filter(nb.genes.overlap > 2) %>%
  select(-nb.genes.overlap)

write.table(enrichr_results_hypo, "RNAseq-Methylation-Intersect/enrichr_results_hypo.txt", sep = "\t", row.names = FALSE)

enrichr_results_hyper <- enrichr(hypermeth_overexpr, databases = dbs)
enrichr_results_hyper <- rbind(enrichr_results_hyper[[1]])
colnames(enrichr_results_hyper) <- c("Term",	"Overlap",	"P-value",	"Adjusted P-value",	"Old P-value",	"Old Adjusted P-value",	"Odds Ratio",	"Combined Score", "Genes")
enrichr_results_hyper <- enrichr_results_hyper %>%
  mutate(nb.genes.overlap = stringr::str_extract(enrichr_results_hyper$Overlap, "[^/]+")) %>%
  mutate(nb.genes.overlap = as.numeric(nb.genes.overlap)) %>%
  filter(nb.genes.overlap > 2) %>%
  select(-nb.genes.overlap)
write.table(enrichr_results_hyper, "RNAseq-Methylation-Intersect/enrichr_results_hyper.txt", sep = "\t", row.names = FALSE)