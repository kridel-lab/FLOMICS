
###
# This script analyses MANTA predictions based on targeted DNA sequencing done at BC Cancer
# And compiles list of cases with BCL2 and/or BCL6 translocations
# (non-coding regions of BCL2 and BCL6 were targeted with specific purpose to detect translocations)
###

# load libraries
packages <- c("dplyr", "ggplot2", "data.table")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

# setwd("~/github/FLOMICS/")

# Read in sample annotation
sample.annotation <- read.csv("metadata/sample_annotations_rcd6Nov2019.csv", header = T)

sample.annotation.FLOMICS.131 <- sample.annotation %>%
  filter(CAPSEQ_DATA == TRUE & CAPSEQ_DATA_YEAR == "2019" & CAPSEQ_INCLUDE == "YES") %>%
  select(SAMPLE_ID, STAGE, TRANSLOC_14_18, INSTITUTION)

# load MANTA rearrangements predictions
predictions <- fread("DNAseq/BC_SV_MANTA/08-17-2020/2020-08-17_Manta_SVs_pass.txt")

# Identify BCL2 rearrangements to chr14, i.e. t(14;18)
BCL2.BA.predictions <- predictions %>%
  filter(gene == "BCL2") %>%
  filter(grepl("14:", ALT)) %>%
  mutate(BCL2_BA_MANTA = 1) %>%
  select(External_identifier, BCL2_BA_MANTA) %>%
  unique()

# Merge BCL2 BA predictions with sample annotation to look for concordance
BCL2.BA.merged <- BCL2.BA.predictions %>%
  right_join(sample.annotation.FLOMICS.131, by = c("External_identifier" = "SAMPLE_ID")) %>%
  mutate(BCL2_BA_MANTA = ifelse(is.na(BCL2_BA_MANTA), 0, BCL2_BA_MANTA)) %>%
  mutate(BCL2_BA_consensus = TRANSLOC_14_18) %>%
  mutate(BCL2_BA_consensus = ifelse(is.na(TRANSLOC_14_18), BCL2_BA_MANTA, BCL2_BA_consensus))

table(BCL2.BA.merged$BCL2_BA_MANTA, BCL2.BA.merged$TRANSLOC_14_18)
# 5 false negatives by Manta
# 1 with MANTA prediction but no t(14:18) (LY_FL_294_T1, JGH)

# Recode LY_FL_294_T1 consensus as 1 because MANTA prediction seems real
# and generate consensus results for BCL2
BCL2.BA.merged <- BCL2.BA.merged %>%
  mutate(BCL2_BA_consensus = ifelse(External_identifier == "LY_FL_294_T1", 1, BCL2_BA_consensus))

# Identify BCL6 rearrangements
BCL6.BA.predictions <- predictions %>%
  filter(BND_DEPTH > 50 | MATE_BND_DEPTH > 50) %>%
  filter(gene == "BCL6", SVTYPE == "BND") %>%
  mutate(BCL6_BA_MANTA = 1) %>%
  select(External_identifier, BCL6_BA_MANTA) %>%
  unique()

# only 3 BC cases out of 31 had BCL6 translocations:
# FL1190 = 	LY_FL_235_T1, FL2008 = 	LY_FL_237_T1, FL2121 = LY_FL_247_T1

final.BA.results <- sample.annotation %>%
  filter(TYPE == "FL") %>%
  select(SAMPLE_ID, STAGE, TIME_POINT, TRANSLOC_14_18) %>%
  left_join(BCL2.BA.merged[,c("External_identifier", "BCL2_BA_consensus")], by = c("SAMPLE_ID" = "External_identifier")) %>%
  left_join(BCL6.BA.predictions, by = c("SAMPLE_ID" = "External_identifier")) %>%
  mutate(BCL2_BA_consensus = ifelse(is.na(BCL2_BA_consensus) & !is.na(TRANSLOC_14_18), TRANSLOC_14_18, BCL2_BA_consensus)) %>%
  mutate(BCL6_BA_consensus = ifelse(SAMPLE_ID %in% sample.annotation.FLOMICS.131$SAMPLE_ID & is.na(BCL6_BA_MANTA), 0, BCL6_BA_MANTA)) %>%
  mutate(BCL6_BA_consensus = ifelse(SAMPLE_ID %in% c("LY_FL_235_T1", "LY_FL_237_T1", "LY_FL_247_T1"), 1, BCL6_BA_consensus)) %>%
  mutate(FLOMICS.131 = ifelse(SAMPLE_ID %in% sample.annotation.FLOMICS.131$SAMPLE_ID, "YES", "NO")) %>%
  select(SAMPLE_ID, TIME_POINT, FLOMICS.131, STAGE, BCL2_BA_consensus, BCL6_BA_consensus)

# # file for Daniel Xia:
# final.BA.results <- sample.annotation %>%
#   filter(TYPE == "FL") %>%
#   select(SAMPLE_ID, STAGE, TIME_POINT, TRANSLOC_14_18) %>%
#   left_join(BCL2.BA.merged[,c("External_identifier", "BCL2_BA_consensus")], by = c("SAMPLE_ID" = "External_identifier")) %>%
#   left_join(BCL6.BA.predictions, by = c("SAMPLE_ID" = "External_identifier")) %>%
#   mutate(BCL2_BA_consensus = ifelse(is.na(BCL2_BA_consensus) & !is.na(TRANSLOC_14_18), TRANSLOC_14_18, BCL2_BA_consensus)) %>%
#   mutate(BCL6_BA_consensus = ifelse(SAMPLE_ID %in% sample.annotation.FLOMICS.131$SAMPLE_ID & is.na(BCL6_BA_MANTA), 0, BCL6_BA_MANTA)) %>%
#   mutate(BCL6_BA_consensus = ifelse(SAMPLE_ID %in% c("LY_FL_235_T1", "LY_FL_237_T1", "LY_FL_247_T1"), 1, BCL6_BA_consensus)) %>%
#   mutate(FLOMICS.131 = ifelse(SAMPLE_ID %in% sample.annotation.FLOMICS.131$SAMPLE_ID, "YES", "NO")) %>%
#   filter(TIME_POINT == "T1") %>%
#   select(SAMPLE_ID, STAGE, TRANSLOC_14_18, BCL2_BA_consensus)
# write.table(final.BA.results, file = "DNAseq/Mutation_and_BA_matrices/BCL2_transloc_for_Daniel.txt",
#             sep = "\t", row.names = FALSE)

# 8 cases in FLOMICS have T2 and T1 sequenced
# for LY_FL_449, there is BCL2 prediction for T2 but not for T1
# we know that t(14;18) is always ancestral, hence can use T2 data to encode BCL2 BA result for T1
final.BA.results <- final.BA.results %>%
  mutate(BCL2_BA_consensus = ifelse(SAMPLE_ID == "LY_FL_449_T1", 1, BCL2_BA_consensus))

nrow(final.BA.results) # n = 183

# Correct these results with FISH data from Adam Smith
final.BA.results.corr <- final.BA.results %>%
  mutate(BCL2_BA_consensus = ifelse(SAMPLE_ID %in% c("LY_FL_050_T1", "LY_FL_219_T1"), "1", BCL2_BA_consensus)) %>% # false negatives by MANTA
  mutate(BCL2_BA_consensus = ifelse(SAMPLE_ID == "LY_FL_119_T1", "0", BCL2_BA_consensus)) %>% # there is a translocation to chr14 based on Manta but not into IGH locus
  mutate(BCL6_BA_consensus = ifelse(SAMPLE_ID == "LY_FL_076_T1", "1", BCL6_BA_consensus)) %>% 
  # another discrepancy is LY_FL_449_T2 where MANTA predicts BA and FISH does not show this, but looks like real translocation when inspecting MANTA prediction -> keep prediction of rearrangement
  mutate(BCL6_BA_consensus = ifelse(SAMPLE_ID %in% c("LY_FL_134_T1", "LY_FL_442_T2", "LY_FL_449_T2"), "0", BCL6_BA_consensus)) # look like false predictions in MANTA as translocation to some repetitive sequence or to same chromosome
  
final.BA.results.T1 <- final.BA.results %>% filter(TIME_POINT == "T1") # n = 173
final.BA.results.corr.T1 <- final.BA.results.corr %>% filter(TIME_POINT == "T1") # n = 173

table(final.BA.results.T1$BCL2_BA_consensus)
table(final.BA.results.corr.T1$BCL2_BA_consensus)
table(final.BA.results.T1$BCL6_BA_consensus)
table(final.BA.results.corr.T1$BCL6_BA_consensus)

table(final.BA.results.T1$STAGE, final.BA.results.T1$BCL2_BA_consensus)
table(final.BA.results.corr.T1$STAGE, final.BA.results.corr.T1$BCL2_BA_consensus)
table(final.BA.results.T1$STAGE, final.BA.results.T1$BCL6_BA_consensus)
table(final.BA.results.corr.T1$STAGE, final.BA.results.corr.T1$BCL6_BA_consensus)

write.csv(final.BA.results, file = "DNAseq/Mutation_and_BA_matrices/BA.results.csv", row.names = FALSE)
write.csv(final.BA.results.corr, file = "DNAseq/Mutation_and_BA_matrices/BA.results.corr.csv", row.names = FALSE)
write.csv(final.BA.results.T1, file = "DNAseq/Mutation_and_BA_matrices/BA.results.T1.csv", row.names = FALSE)
write.csv(final.BA.results.corr.T1, file = "DNAseq/Mutation_and_BA_matrices/BA.results.corr.T1.csv", row.names = FALSE)
