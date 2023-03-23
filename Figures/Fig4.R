# Load packages
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(dplyr)
library(ggfortify)
library(tidyr)
library(reshape)
library(ggpubr)
library(GenomicRanges)
library(data.table)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(stringr)

# System settings
options(stringsAsFactors = F, error = NULL)
date <- Sys.Date()
set.seed(1234)

# Change directory
setwd("<working directory path>")

# Get the 850k annotation data
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
EPIC_Annotations <- data.frame(ann850k)
remove(ann850k)

# Define CpGs of interestBetaMatrix
CpG_islands <- EPIC_Annotations %>% filter(Relation_to_Island == "Island") %>% data.frame() %>% .$Name
OpenSea <- EPIC_Annotations %>% filter(Relation_to_Island == "OpenSea") %>% data.frame() %>% .$Name

# Load epiCMIT Data
load("Estimate.epiCMIT.RData")

# Read in methylation data
bval <- readRDS("beta_combat.rds")
mval <- readRDS("mval_combat.rds")

# Read in sample annotations
Sample_Annotations <- read.table(file = "20221228_sample_annotations.txt",
                                 sep = "\t", header = TRUE, na.strings=c("","NA")) %>% 
subset(TYPE!="TFL") %>%
subset(TIME_POINT!="T2") %>%
filter(SAMPLE_ID %in% colnames(bval))

# Subset of methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
mval <- mval[, Sample_Annotations$SAMPLE_ID]

# Read in cluster information
flexmix_clust <- read.table(file = "2023-01-19_14_GMM_Cluster_Labels_flexmix_clusters.txt",
                            sep = "\t",
                            header = TRUE, na.strings=c("","NA")) %>%
                            dplyr::select(SAMPLE_ID, ClusterAIC)

# Read in sample sheet with batch information
targets <- read.csv(file = "SampleSheet.csv",
                    skip = 1, header = TRUE, na.strings=c("","NA"))

targets$ID <- gsub(" ", "", paste(targets$Sample_Name,".",targets$Sentrix_Position))

# Remove repeated samples to obtain correct batch information
rep_rm_samples <- c("LY_FL_158_T1.R01C01", "LY_FL_479_T1.R02C01", "LY_FL_488_T1.R06C01",
                    "LY_FL_498_T1.R06C01", "LY_FL_523_T1.R01C01", "LY_FL_524_T1.R02C01",
                    "LY_FL_525_T1.R03C01", "LY_FL_527_T1.R04C01", "LY_FL_529_T1.R01C01",
                    "LY_FL_1156_T1.R06C01", "LY_FL_571_T1.R01C01", "HCT116_DKO_methylated.R08C01",
                    "HCT116_DKO_methylated.R04C01", "LY_FL_311_T1.R04C01", "LY_FL_159_T1.R08C01",
                    "LY_FL_535_T1.R08C01", "LY_FL_536_T1.R01C01")

targets <- targets[ ! targets$ID %in% rep_rm_samples, ]
targets$Sample_Name <- str_replace(targets$Sample_Name, "LY_FL_159_T1_rep", "LY_FL_159_T1")
targets <- targets %>% dplyr::select(Sample_Name, Batch)
colnames(targets)[1]  <- "SAMPLE_ID"

# Read in clinical data
clinical_data <- read.csv(file = "20220525_clinical_data.csv",
                          header = TRUE, na.strings=c("","NA")) %>%
dplyr::select(LY_FL_ID, ANN_ARBOR_STAGE)

# Merge sample annotation, cluster information, sample sheet
df_merged <- merge(x = targets,
             y = flexmix_clust,
             by.x = "SAMPLE_ID",
             by.y = "SAMPLE_ID",
             all.x = TRUE,
             all.y = TRUE)

df_merged_1 <- merge(x = Sample_Annotations,
               y = df_merged,
               by.x = "SAMPLE_ID",
               by.y = "SAMPLE_ID",
               all.x = TRUE)

Sample_Annotations <- merge(x = df_merged_1,
                      y = clinical_data,
                      by.x = "LY_FL_ID",
                      by.y = "LY_FL_ID",
                      all.x = TRUE)

# Create a additional column with ClusterAIC, RLN, DLBCL
Sample_Annotations$ClusterAIC_TYPE <- Sample_Annotations$ClusterAIC
Sample_Annotations <- Sample_Annotations %>% mutate(ClusterAIC_TYPE = coalesce(ClusterAIC_TYPE,TYPE))
Sample_Annotations$ClusterAIC_TYPE[Sample_Annotations$ClusterAIC_TYPE == "FL"] <- NA

# Create a additional column with Ann Arbor Stage, RLN, DLBCL
Sample_Annotations$ANN_ARBOR_STAGE_TYPE <- Sample_Annotations$ANN_ARBOR_STAGE
Sample_Annotations <- Sample_Annotations %>% mutate(ANN_ARBOR_STAGE_TYPE = coalesce(as.character(ANN_ARBOR_STAGE_TYPE), TYPE))
Sample_Annotations$ANN_ARBOR_STAGE_TYPE[Sample_Annotations$ANN_ARBOR_STAGE_TYPE == "FL"] <- NA

row.names(Sample_Annotations) <- Sample_Annotations$SAMPLE_ID

## bval - All CpGs
# bval_sub <- bval[sample(nrow(bval), 10000), ] # Change For Subset
bval_sub <- bval
bval_sub_t <- t(bval_sub)
bval_sub_t <- merge(bval_sub_t, Sample_Annotations, by=0)
names(bval_sub_t)[1] <- ""
bval_sub_t <- data.frame(bval_sub_t[,-1], row.names=bval_sub_t[,1])
bval_sub_t_df <- bval_sub_t[1:765769] # Change For Subset

## mval - All CpGs
# mval_sub <- mval[sample(nrow(mval), 10000), ] # Change For Subset
mval_sub <- mval
mval_sub_t <- t(mval_sub)
mval_sub_t <- merge(mval_sub_t, Sample_Annotations, by=0)
names(mval_sub_t)[1] <- ""
mval_sub_t <- data.frame(mval_sub_t[,-1], row.names=mval_sub_t[,1])
mval_sub_t_df <- mval_sub_t[1:765769] # Change For Subset

# epiCMIT Analysis
DNAm.epiCMIT <- DNAm.to.epiCMIT(DNAm = data.frame(bval),
                                DNAm.genome.assembly = "hg19",
                                map.DNAm.to = "Illumina.450K.epiCMIT",
                                min.epiCMIT.CpGs = 800)

epiCMIT.Illumina <- epiCMIT(DNAm.epiCMIT = DNAm.epiCMIT,
                            return.epiCMIT.annot = FALSE,
                            export.results = TRUE,
                            export.results.dir = "",
                            export.results.name = "<working directory path>")

head(epiCMIT.Illumina$epiCMIT.scores)
epiCMIT.Illumina$epiCMIT.run.info

## PLOT - A
## ========
# M-Value Clustering - All CpGs
t <- mval_sub_t %>% drop_na(ClusterAIC_TYPE)
df <- mval_sub_t_df[intersect(row.names(t), row.names(mval_sub_t_df)), ]
pca_res <- prcomp(df, scale. = TRUE)

# pdf("mval-pca-plot-ClusterAIC-All-CpGs.pdf")
P1 <- autoplot(pca_res,
           data = t,
           colour = 'ClusterAIC_TYPE',
           frame = TRUE
      ) +
      theme_bw() +
      theme(title =element_text(size=20, face='bold'),
            axis.title.x = element_text(size = 20),
            axis.text.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            legend.text=element_text(size=20),
            legend.position="bottom",
            legend.title = element_text(size=20)) + 
      ggtitle("Methylation - All CpGs") +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_point(aes(colour = factor(ClusterAIC_TYPE)), size = 1.5) +
      scale_color_manual(values = c(RLN="orange", C1="#F8766D", C2="#A3A500", C3="#00BF7D", C4="#00B0F6", C5="#E76BF3", DLBCL="grey20")) +
      scale_fill_manual(values = c(RLN="orange", C1="#F8766D", C2="#A3A500", C3="#00BF7D", C4="#00B0F6", C5="#E76BF3", DLBCL="grey20")) +
      coord_cartesian(clip = 'off') +
      theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
      labs(fill='Cluster') +
      labs(color='Cluster')
# dev.off()

## PLOT - B
## ========
# Mean methylation flexmix - All probes
my_comparisons <- list( c("RLN", "DLBCL"), 
                        c("RLN", "C5"), 
                        c("RLN", "C4"),
                        c("RLN", "C3"),
                        c("RLN", "C2"),
                        c("RLN", "C1"))

# pdf("Mean-Methylation-All-CpGs-ClusterAIC.pdf")
df <- bval_sub %>%
      data.frame() %>%
      mutate(probe = row.names(.)) %>%
      reshape2::melt() %>%
      left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
      group_by(variable) %>%
      summarize(mean = mean(value)) %>%
      left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
      left_join(Sample_Annotations, by = c("variable" = "SAMPLE_ID")) %>%
      filter(!is.na(ClusterAIC_TYPE)) %>%
      mutate(ClusterAIC_TYPE = factor(ClusterAIC_TYPE, levels = c("RLN", "C1", "C2", "C3", "C4", "C5", "DLBCL")))

      P2 <- ggboxplot(df, x = "ClusterAIC_TYPE", y = "mean", color = "ClusterAIC_TYPE", add = "jitter", size = 0.5) +
        # geom_hline(yintercept = mean(df$mean), linetype="dashed", color="blue2", linewidth=0.8) +
        scale_color_manual(values = c(RLN="orange", C1="#F8766D", C2="#A3A500", C3="#00BF7D", C4="#00B0F6", C5="#E76BF3", DLBCL="grey20")) +
        scale_fill_manual(values = c(RLN="orange", C1="#F8766D", C2="#A3A500", C3="#00BF7D", C4="#00B0F6", C5="#E76BF3", DLBCL="grey20")) +
        theme_bw() +
        theme(title = element_text(size=20, face='bold'),
              axis.title.x = element_text(size = 20),
              axis.text.x = element_text(size = 20),
              axis.title.y = element_text(size = 20),
              axis.text.y = element_text(size = 20),
              legend.position = "none") +
        ylab("Mean Methylation / Beta Value") +
        ggtitle("All CpGs") +
        theme(plot.title = element_text(hjust = 0.5)) +
        stat_compare_means(size = 6, comparisons = my_comparisons,
                           label.y = c(max(df$mean),
                                       max(df$mean) + 0.002,
                                       max(df$mean) + 0.004,
                                       max(df$mean) + 0.006,
                                       max(df$mean) + 0.0085,
                                       max(df$mean) + 0.0105)) +
        stat_compare_means(size = 7, label.x = 3.8, label.y = max(df$mean) + 0.011) +
        coord_cartesian(clip = 'off') +
        theme(plot.margin = margin(0.4, 1, 0.4, 0.4, "cm")) +
        scale_x_discrete(name ="Cluster") +
        scale_y_continuous(limits=c(min(df$mean) - 0.002, max(df$mean) + 0.012))
# dev.off()

## PLOT - C
## ========
# Mean methylation flexmix - CpG probes
my_comparisons <- list( c("RLN", "DLBCL"), 
                        c("RLN", "C5"), 
                        c("RLN", "C4"),
                        c("RLN", "C3"),
                        c("RLN", "C2"),
                        c("RLN", "C1"))

# pdf("Mean-Methylation-CpG-Islands-ClusterAIC.pdf")
df <- bval_sub[intersect(CpG_islands, row.names(bval_sub)),] %>%
        data.frame() %>%
        mutate(probe = row.names(.)) %>%
        reshape2::melt() %>%
        left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
        group_by(variable) %>%
        summarize(mean = mean(value)) %>%
        left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
        left_join(Sample_Annotations, by = c("variable" = "SAMPLE_ID")) %>%
        filter(!is.na(ClusterAIC_TYPE)) %>%
        mutate(ClusterAIC_TYPE = factor(ClusterAIC_TYPE, levels = c("RLN", "C1", "C2", "C3", "C4", "C5", "DLBCL")))
  
        P3 <- ggboxplot(df, x = "ClusterAIC_TYPE", y = "mean", color = "ClusterAIC_TYPE", add = "jitter", size = 0.5) +
          # geom_hline(yintercept = mean(df$mean), linetype="dashed", color="blue2", linewidth=0.8) +
          scale_color_manual(values = c(RLN="orange", C1="#F8766D", C2="#A3A500", C3="#00BF7D", C4="#00B0F6", C5="#E76BF3", DLBCL="grey20")) +
          scale_fill_manual(values = c(RLN="orange", C1="#F8766D", C2="#A3A500", C3="#00BF7D", C4="#00B0F6", C5="#E76BF3", DLBCL="grey20")) +
          theme_bw() +
          theme(title = element_text(size=20, face='bold'),
                axis.title.x = element_text(size = 20),
                axis.text.x = element_text(size = 20),
                axis.title.y = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.position = "none") +
          ylab("Mean Methylation / Beta Value") +
          theme(plot.title = element_text(hjust = 0.5))+
          stat_compare_means(size = 6, comparisons = my_comparisons,
                               label.y = c(max(df$mean),
                                           max(df$mean) + 0.004,
                                           max(df$mean) + 0.008,
                                           max(df$mean) + 0.013,
                                           max(df$mean) + 0.018,
                                           max(df$mean) + 0.023)) +
          stat_compare_means(size = 7, label.x = 3.8, label.y = max(df$mean) + 0.024) +
          ggtitle("CpG Islands") +
          coord_cartesian(clip = 'off') +
          theme(plot.margin = margin(0.4, 1, 0.4, 0.4, "cm")) +
          scale_x_discrete(name ="Cluster") +
          scale_y_continuous(limits=c(min(df$mean) - 0.002, max(df$mean) + 0.026))
# dev.off()

## PLOT - D
## ========
my_comparisons <- list( c("RLN", "DLBCL"), 
                        c("RLN", "FL"),
                        c("FL", "DLBCL"))

# pdf("epiCMIT_Type.pdf")
P4 <- epiCMIT.Illumina$epiCMIT.scores %>%
          left_join(Sample_Annotations, by = c("Samples" = "SAMPLE_ID")) %>%
          mutate(TYPE = factor(TYPE, levels = c("RLN", "FL", "DLBCL"))) %>%
          ggpubr::ggboxplot(x = "TYPE", y = "epiCMIT", add = "jitter", palette = "jco", color = "TYPE", size = 0.5) +
          scale_color_manual(values = c(RLN="orange", FL="blue1", DLBCL="grey20")) +
          scale_fill_manual(values = c(RLN="orange", FL="blue1", DLBCL="grey20")) +
          theme_bw() +
          theme(title = element_text(size=20, face='bold'),
                axis.title.x = element_text(size = 20),
                axis.text.x = element_text(size = 20),
                axis.title.y = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.position = "none") +
          ylab("epiCMIT Score") +
          theme(plot.title = element_text(hjust = 0.5)) +
          ggtitle("Proliferative History") +
          # geom_hline(yintercept = mean(epiCMIT.Illumina$epiCMIT.scores$epiCMIT), linetype="dashed", color="blue2", linewidth=0.8) +
          stat_compare_means(size = 6, comparisons = my_comparisons,
                             label.y = c(max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT),
                                         max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT) + 0.06,
                                         max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT) + 0.10)) +
          stat_compare_means(size = 7, label.x = 1, label.y = max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT) + 0.16) +
          coord_cartesian(clip = 'off') +
          theme(plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm")) +
          scale_y_continuous(limits=c(min(epiCMIT.Illumina$epiCMIT.scores$epiCMIT)-0.02,
                                      max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT)+0.18)) +
          scale_x_discrete(name ="Type")
# dev.off()


## PLOT - E
## ========
my_comparisons <- list( c("RLN", "C5"),
                        c("RLN", "C4"), 
                        c("RLN", "C3"),
                        c("RLN", "C2"),
                        c("RLN", "C1"))

# pdf("epiCMIT_ClusterAIC.pdf")
P5 <- epiCMIT.Illumina$epiCMIT.scores %>%
          inner_join(Sample_Annotations, by = c("Samples" = "SAMPLE_ID")) %>%
          filter(!is.na(ClusterAIC_TYPE)) %>%
          filter(!grepl("DLBCL", ClusterAIC_TYPE)) %>%
          mutate(ClusterAIC_TYPE = factor(ClusterAIC_TYPE, levels = c("RLN", "C1", "C2", "C3", "C4", "C5"))) %>%
          ggpubr::ggboxplot(x = "ClusterAIC_TYPE", y = "epiCMIT", add = "jitter", palette = "jco", color = "ClusterAIC_TYPE", size = 0.5) +
          scale_color_manual(values = c(RLN="orange", C1="#F8766D", C2="#A3A500", C3="#00BF7D", C4="#00B0F6", C5="#E76BF3")) +
          scale_fill_manual(values = c(RLN="orange", C1="#F8766D", C2="#A3A500", C3="#00BF7D", C4="#00B0F6", C5="#E76BF3")) +
          theme_bw() +
          theme(title = element_text(size=20, face='bold'),
                axis.title.x = element_text(size = 20),
                axis.text.x = element_text(size = 20),
                axis.title.y = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.position = "none") +
          ylab("epiCMIT Score") +
          theme(plot.title = element_text(hjust = 0.5)) +
          ggtitle("Proliferative History") +
          # geom_hline(yintercept = mean(epiCMIT.Illumina$epiCMIT.scores$epiCMIT), linetype="dashed", color="blue2", linewidth=0.8) +
          stat_compare_means(size = 6, comparisons = my_comparisons,
                             label.y = c(max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT),
                                         max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT) + 0.04,
                                         max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT) + 0.08,
                                         max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT) + 0.12,
                                         max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT) + 0.16)) +
          stat_compare_means(size = 7, label.x = 3.6, label.y = max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT)+0.16) +
          coord_cartesian(clip = 'off') +
          theme(plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm")) +
          scale_y_continuous(limits=c(min(epiCMIT.Illumina$epiCMIT.scores$epiCMIT)-0.02,
                                      max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT)+0.18)) +
          scale_x_discrete(name ="Cluster")
# dev.off()

## PLOT - F
## ========
my_comparisons <- list( c("RLN", "4"),
                        c("RLN", "3"), 
                        c("RLN", "2"),
                        c("RLN", "1"))

# pdf("epiCMIT_Ann_Arbor_Stage.pdf")
P6 <- epiCMIT.Illumina$epiCMIT.scores %>%
        left_join(Sample_Annotations, by = c("Samples" = "SAMPLE_ID")) %>%
        filter(!is.na(ANN_ARBOR_STAGE_TYPE)) %>%
        filter(!grepl("DLBCL", ANN_ARBOR_STAGE_TYPE)) %>%
        mutate(ANN_ARBOR_STAGE_TYPE = factor(ANN_ARBOR_STAGE_TYPE, levels = c("RLN", "1", "2", "3", "4"))) %>%
        ggpubr::ggboxplot(x = "ANN_ARBOR_STAGE_TYPE", y = "epiCMIT", add = "jitter", palette = "jco", color = "ANN_ARBOR_STAGE_TYPE", size = 0.5) +
        scale_color_manual(values = c(RLN="orange", "1"="yellow4", "2"="green4", "3"="cyan4", "4"="tan4")) +
        scale_fill_manual(values = c(RLN="orange", "1"="yellow4", "2"="green4", "3"="cyan4", "4"="tan4")) +
        theme_bw() +
        theme(title = element_text(size=20, face='bold'),
              axis.title.x = element_text(size = 20),
              axis.text.x = element_text(size = 20),
              axis.title.y = element_text(size = 20),
              axis.text.y = element_text(size = 20),
              legend.position = "none") +
        ylab("epiCMIT Score") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle("Proliferative History") +
        # geom_hline(yintercept = mean(epiCMIT.Illumina$epiCMIT.scores$epiCMIT), linetype="dashed", color="blue2", linewidth=0.8) +
        stat_compare_means(size = 6, comparisons = my_comparisons,
                           label.y = c(max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT),
                                       max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT) + 0.045,
                                       max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT) + 0.085,
                                       max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT) + 0.125)) +
        stat_compare_means(size = 7, label.x = 2.9, label.y = max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT)+0.14) +
        coord_cartesian(clip = 'off') +
        theme(plot.margin = margin(0.4, 1, 0.4, 0.4, "cm")) +
        scale_y_continuous(limits=c(min(epiCMIT.Illumina$epiCMIT.scores$epiCMIT)-0.02,
                                    max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT)+0.16)) +
        scale_x_discrete(name ="Ann Arbor Stage")
# dev.off()

## Arrange Plots
## =============
pdf("Fig4.pdf", width = 21.02, height = 14.01)
ggarrange(P1, P2, P3, P4, P5, P6,
          labels = c("A", "B", "C", "D", "E", "F"),
          font.label = list(size = 32, face = "bold"),
          nrow = 2, ncol = 3)
dev.off()
