
###
# This script explores MANTA predictions based on targeted DNA sequencing done at BC Cancer
# And compiles list of cases with BCL2 and/or BCL6 translocations
# (non-coding regions of BCL2 and BCL6 were targeted with specific purpose to detect translocations)
###

# load libraries 
packages <- c("dplyr", "ggplot2")
lapply(packages, require, character.only = TRUE)              

date <- Sys.Date()

# setwd("~/FLOMICS")

# read in meta data
meta.data <- read.table("Git/Methylation/Pipeline/sample_annotations_rcd6Nov2019.txt", head = TRUE, sep = "\t")

# load library identifiers
library_identifiers <- read.table("Targeted_sequencing_BCGSC/DNAseq_library_identifiers.txt", header = TRUE, sep = "\t")
# B48248 and B48278 failed library construction
pts_capseq <- library_identifiers %>%
  filter(!Library %in% c("B48248", "B48278")) %>%
  .$External_ID %>% as.character() # n = 131

# load coverage statistics
coverage <- read.table("Targeted_sequencing_BCGSC/Results/BC_TargetSeq_summary_mean_coverage_regions.txt", header = T, sep = "\t") %>%
  group_by(pats) %>%
  summarize(mean = mean(value)) %>%
  data.frame()

# load MANTA rearrangements predictions
df <- read.csv("Targeted_sequencing_BCGSC/Results/Manta/2019-12-03_FLOMICS_MANTA_BC_SVs_with_annotations.csv")

# Remember: IGK on chr2, IGL on chr 22 .. but in the end we don't see BCL2 rearrangements to chr2, 22.

# list genes with most predictions
sort(table(df$hgnc_symbol), decreasing = TRUE)[1:20] # Only BCL2 and BCL6 translocations seem real

# Most predictions have very little coverage
# coverage MATE_BND_DEPTH for BCL2 and BCL6 rearrangements are higher than coverage for all predictions
# hence, reasonable to set filter for MATE_BND_DEPTH to search for putative novel rearrangements
# although important to bear in mind that most predictions are in regions that were not actually targeted with capture panel
# hence method largely imperfect to detect novel rearrangements

# BND_DEPTH = Read depth at local translocation breakend
# MATE_BND_DEPTH = Read depth at remote translocation mate breakend

df %>%
  ggplot(aes(x = MATE_BND_DEPTH)) +
  geom_histogram()

df %>%
  filter(hgnc_symbol == "BCL2") %>%
  filter(grepl("14:", ALT)) %>%
  ggplot(aes(x = MATE_BND_DEPTH)) +
  geom_histogram()

df %>%
  filter(hgnc_symbol == "BCL6") %>%
  ggplot(aes(x = MATE_BND_DEPTH)) +
  geom_histogram()

# putative true rearrangements:
df %>%
  filter(BND_DEPTH > 20 & MATE_BND_DEPTH > 20) %>%
  group_by(hgnc_symbol) %>%
  summarize(count = n()) %>%
  arrange(desc(count))

# can further filter by identifying those rearrangements that have systematically identical start/end coordinates
# these predictions are very likely artifacts
# in the end, candidate novel rearrangements:
# FUT8, PRDM6, CIITA, ANKUB1, ERCC6L2, LNX2, RNF13, TACR1
# FUT8 and PRDM6 likely artifacts as many filtered out predictions mapping to same coordinates
df %>%
  filter(BND_DEPTH > 20 & MATE_BND_DEPTH > 20) %>%
  group_by(hgnc_symbol, SV_start) %>%
  summarize(count = n()) %>% 
  filter(count < 2) %>%
  ungroup() %>%
  group_by(hgnc_symbol) %>%
  summarize(count = n()) %>% 
  arrange(desc(count)) %>%
  print(n = 20)

# Identify BCL2 rearrangements to chr14, i.e. t(14;18)
df_BCL2_BA <- df %>%
  filter(hgnc_symbol == "BCL2") %>%
  filter(grepl("14:", ALT)) %>%
  mutate(BCL2_BA_MANTA = 1) %>%
  select(External_ID, BCL2_BA_MANTA) %>%
  unique()

merged <-
  df_BCL2_BA %>%
  right_join(library_identifiers) %>%
  filter(!Library %in% c("B48248", "B48278")) %>%
  mutate(BCL2_BA_MANTA = ifelse(is.na(BCL2_BA_MANTA), 0, BCL2_BA_MANTA)) %>%
  right_join(meta.data[,c("SAMPLE_ID", "TYPE", "STAGE", "TRANSLOC_14_18")], by = c("External_ID" = "SAMPLE_ID")) %>%
  filter(TYPE == "FL") %>%
  mutate(BCL2_BA_consensus = TRANSLOC_14_18) %>%
  # mutate(BCL2_BA_consensus = ifelse(is.na(TRANSLOC_14_18) & BCL2_BA_MANTA == 1, 1, BCL2_BA_consensus))
  mutate(BCL2_BA_consensus = ifelse(is.na(TRANSLOC_14_18), BCL2_BA_MANTA, BCL2_BA_consensus))

table(merged$BCL2_BA_MANTA, merged$TRANSLOC_14_18)
# 6 false negatives by Manta

BCL2_BA_MANTA_FALSE_NEG <- merged %>%
  filter(TRANSLOC_14_18 == 1 & BCL2_BA_MANTA == 0) %>%
  left_join(coverage, by = c("Library" = "pats"))
# coverage overall rather high except LY_FL_292_T1  

table(merged$BCL2_BA_MANTA, merged$STAGE)
table(merged$TRANSLOC_14_18, merged$STAGE)
table(merged$BCL2_BA_consensus, merged$STAGE)

# Identify BCL6 rearrangements
BCL6.calls.pts <- df %>%
  filter(BND_DEPTH > 50 & MATE_BND_DEPTH > 50) %>%
  filter(hgnc_symbol == "BCL6", SVTYPE == "BND") %>%
  .$External_ID %>% as.character %>% unique()

BC.BCL6 <- read.table("~/tfl_capseq/meta-data/clinical-pathology-data/clinical-FF-FFPET.txt", sep = "\t", header = T) %>%
  filter(patient_id %in% meta.data$OTHER_ID) %>%
  select(patient_id, T1.BCL6.BA) %>%
  left_join(meta.data[,c("SAMPLE_ID", "OTHER_ID")], by = c("patient_id" = "OTHER_ID"))
# only 3 BC cases out of 31 had BCL6 translocations:
# FL1190 = 	LY_FL_235_T1, FL2008 = 	LY_FL_237_T1, FL2121 = LY_FL_247_T1

merged <- merged %>%
  mutate(BCL6_BA_MANTA = ifelse(External_ID %in% pts_capseq, 0, NA)) %>%
  mutate(BCL6_BA_MANTA = ifelse(External_ID %in% BCL6.calls.pts, 1, BCL6_BA_MANTA)) %>%
  left_join(BC.BCL6[,c("T1.BCL6.BA", "SAMPLE_ID")], by = c("External_ID" = "SAMPLE_ID")) %>%
  rename(BCL6_BA_BC = T1.BCL6.BA) %>%
  mutate(BCL6_BA_BC = ifelse(BCL6_BA_BC == "NO", 0,
                      ifelse(BCL6_BA_BC == "YES", 1, NA))) %>%
  mutate(BCL6_BA_consensus = BCL6_BA_MANTA) %>%
  mutate(BCL6_BA_consensus = ifelse(!is.na(BCL6_BA_BC), BCL6_BA_BC, BCL6_BA_consensus)) %>%
  mutate(External_ID = as.character(External_ID)) %>%
  arrange(External_ID)
  
table(merged$BCL6_BA_consensus) # 13 BCL6 translocations out of 148 cases, i.e. 8%                                

merged <- merged %>%
  select(External_ID, BCL2_BA_consensus, BCL6_BA_consensus)

# Write BCL2 and BCL6 break apart results out
write.table(merged, file = paste0("Targeted_sequencing_BCGSC/Results/Manta/", date, "_BCL2_BCL6_rearrangements.txt"),
                                  sep = "\t", row.names = FALSE)
  
            