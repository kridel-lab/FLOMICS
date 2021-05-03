#---
# This script analyzes Blueprint data, intersects diff methylated CpGs with out data,
# explores chromatin states
# This script needs to be run from /FLOMICS folder
#---

packages <- c("dplyr", "ggplot2", "minfi", "missMethyl", "factoextra", "limma", "reshape2",
              "GenomicRanges", "ggpubr", "gridExtra", "DMRcate", "enrichR", "forcats", "regioneR")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

# setwd("~/github/FLOMICS/")

AnnotationFile <- data.table::fread(file = "methylation/Ann_IlluminaHumanMethylationEPICanno.ilm10b2.hg19.csv")

# Load Blueprint sample sheet
SampleSheet <- read.table(file = "methylation/Blueprint/EGAS00001001196.txt", sep = "\t", header = T)
SampleSheet$Raw_data_File <- substr(SampleSheet$Raw_data_File, 1, 63)
SampleSheet$Raw_data_File <- gsub("/", "_", SampleSheet$Raw_data_File)
SampleSheet <- SampleSheet %>%
  mutate(TYPE = case_when(grepl("^naiBC", Sample.Name) ~ "naiBC",
                          grepl("^t-naiBC", Sample.Name) ~"t-naiBC",
                          grepl("^memBC", Sample.Name) ~"memBC",
                          grepl("^gcBC", Sample.Name) ~"gcBC",
                          grepl("^t-PC", Sample.Name) ~"t-PC",
                          grepl("^bm-PC", Sample.Name) ~"bm-PC"))

# Load Blueprint data
Blueprint_RGset <- read.metharray.exp(base = "methylation/Blueprint", force = TRUE)
colnames(Blueprint_RGset) <- SampleSheet$Sample.Name[match(colnames(Blueprint_RGset), SampleSheet$Raw_data_File)]

# Load our data
ClinicalFile <- data.table::fread(file = "metadata/sample_annotations_rcd6Nov2019.csv") %>%
  filter(EPIC_QC != "Bad" & EPIC_INCLUDE == "YES" & TIME_POINT == "T1" & !is.na(SAMPLE_ID_TGL))

# Make combined phenodata file
pData_1 <- ClinicalFile[,c("SAMPLE_ID", "TYPE")]
pData_2 <- SampleSheet %>%
  select(SAMPLE_ID = Sample.Name, TYPE) 
pData <- rbind(pData_1, pData_2)

load("methylation/Output_Data_2_28Jan2020.RData")
UHN_RGset <- Output_Data_2$RGChannelSet
# UHN_RGset <- UHN_RGset[,ClinicalFile$SAMPLE_ID_TGL]
UHN_RGset <- UHN_RGset[,ClinicalFile$SAMPLE_ID_TGL[1:176]]
colnames(UHN_RGset) <- ClinicalFile$SAMPLE_ID[match(colnames(UHN_RGset), ClinicalFile$SAMPLE_ID_TGL)]

# Combine arrays
Combined_RGset <- combineArrays(UHN_RGset, Blueprint_RGset, outType = "IlluminaHumanMethylation450k")

# Normalize
# Combined_mSet_Sw <- preprocessNoob(Combined_RGset)
Combined_mSet_Sw <- preprocessSWAN(Combined_RGset)

# Order Combined_M_valueMatrix_sub based on pData file
Combined_mSet_Sw <- Combined_mSet_Sw[,pData$SAMPLE_ID]

# Get M value matrix
Combined_M_valueMatrix <- getM(Combined_mSet_Sw)
Combined_M_valueMatrix_sub <- Combined_M_valueMatrix[sample(nrow(Combined_M_valueMatrix), 100000), ]

# Perform PCA - only on Blueprint samples
BLueprint_M_valueMatrix_sub <- Combined_M_valueMatrix_sub[,as.character(pData_2$SAMPLE_ID)]
Blueprint.M.matrix.pca <- prcomp(t(BLueprint_M_valueMatrix_sub), scale = TRUE)
fviz_eig(Blueprint.M.matrix.pca)
fviz_pca_ind(Blueprint.M.matrix.pca)
Blueprint.M.matrix.pca.ind <- get_pca_ind(Blueprint.M.matrix.pca)
Blueprint.M.matrix.pca.ind <- Blueprint.M.matrix.pca.ind$coord
Blueprint.M.matrix.pca.ind <- data.frame(Blueprint.M.matrix.pca.ind)
Blueprint.M.matrix.pca.ind$SAMPLE_ID <- row.names(Blueprint.M.matrix.pca.ind)
Blueprint.M.matrix.pca.ind <- Blueprint.M.matrix.pca.ind %>%
  select(SAMPLE_ID, M.matrix.Dim.1 = Dim.1, M.matrix.Dim.2 = Dim.2, M.matrix.Dim.3 = Dim.3) %>%
  left_join(pData_2)

Blueprint.M.matrix.pca.ind %>%
  mutate(TYPE = ifelse(TYPE == "naiBC", "peripheral blood B cells",
                       ifelse(TYPE == "t-naiBC", "tonsil naive B cells",
                              ifelse(TYPE == "gcBC", "germinal centre B cells",
                                     ifelse(TYPE == "t-PC", "tonsil plasma cells",
                                            ifelse(TYPE == "bm-PC", "bone marrow plasma cells",
                                                   ifelse(TYPE == "memBC", "memory B cells", NA))))))) %>%
  mutate(TYPE = factor(TYPE, levels = c("peripheral blood B cells",
                                        "tonsil naive B cells",
                                        "germinal centre B cells",
                                        "tonsil plasma cells",
                                        "bone marrow plasma cells",
                                        "memory B cells"))) %>%
  ggplot(aes(M.matrix.Dim.1, M.matrix.Dim.2, col = TYPE),
         main = "PCA plot") +
  geom_point(size = 3) +
  theme_bw() +
  theme_bw() +
  theme(axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_color_brewer(palette = "Spectral") +
  xlab("Dimension 1 (38.12%)") +
  ylab("Dimension 2 (12.28%)")
ggsave(paste0("img/",  date, " RK-PCA-Blueprint.pdf"), width = 14, height = 8, units = "cm")

# Perform PCA
M.matrix.pca <- prcomp(t(Combined_M_valueMatrix), scale = TRUE)
fviz_eig(M.matrix.pca)
M.matrix.pca.ind <- get_pca_ind(M.matrix.pca)
M.matrix.pca.ind <- M.matrix.pca.ind$coord
M.matrix.pca.ind <- data.frame(M.matrix.pca.ind)
M.matrix.pca.ind$SAMPLE_ID <- row.names(M.matrix.pca.ind)
M.matrix.pca.ind <- M.matrix.pca.ind %>%
  select(SAMPLE_ID, M.matrix.Dim.1 = Dim.1, M.matrix.Dim.2 = Dim.2, M.matrix.Dim.3 = Dim.3) %>%
  left_join(pData)

plots_to_save <- list()
plots_to_save[[1]] <- fviz_eig(M.matrix.pca)
plots_to_save[[2]] <- M.matrix.pca.ind %>%
  ggplot(aes(M.matrix.Dim.1, M.matrix.Dim.2, col = TYPE), main = "PCA plot") +
  geom_point() + theme_bw()
plots_to_save[[3]] <- M.matrix.pca.ind %>%
  ggplot(aes(M.matrix.Dim.1, M.matrix.Dim.3, col = TYPE), main = "PCA plot") +
  geom_point() + theme_bw()
plots_to_save[[4]] <- M.matrix.pca.ind %>%
  ggplot(aes(M.matrix.Dim.2, M.matrix.Dim.3, col = TYPE), main = "PCA plot") +
  geom_point() + theme_bw()

ggsave(filename = paste0("img/", date, " RK-PCA-BLueprint-ourFLs.png"),
       marrangeGrob(grobs = plots_to_save, ncol = 2, nrow = 2),
       device = "png", width = 25, height = 18, units = "cm")

#---
# Differentially methylated probes - model NOT taking purity into account
#---

design <- model.matrix(~0 + TYPE, data = pData[,c("SAMPLE_ID", "TYPE")])
colnames(design) <- c("bmPC", "DLBCL", "FL", "gcBC", "memBC", "naiBC", "RLN", "tnaiBC", "tPC")
fit <- lmFit(Combined_M_valueMatrix, design)
contMatrix <- makeContrasts(FL-naiBC,
                            FL-tnaiBC,
                            FL-gcBC,
                            FL-tPC,
                            FL-bmPC,
                            FL-memBC,
                            FL-DLBCL,
                            FL-RLN,
                            levels = design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2)) # DM CpGs at FDR < 0.05

log.FC = 2
adj.P.Value = 0.05

# FL vs. naiBC
diff_meth_FL_naiBC <- topTable(fit2, coef = "FL - naiBC", adjust.method = "BH", number = Inf)
diff_meth_FL_naiBC$Probe <- row.names(diff_meth_FL_naiBC)
diff_meth_FL_naiBC$type <- "naiBC"

# FL vs. tnaiBC
diff_meth_FL_tnaiBC <- topTable(fit2, coef = "FL - tnaiBC", adjust.method = "BH", number = Inf)
diff_meth_FL_tnaiBC$Probe <- row.names(diff_meth_FL_tnaiBC)
diff_meth_FL_tnaiBC$type <- "tnaiBC"

# FL vs. gcBC
diff_meth_FL_gcBC <- topTable(fit2, coef = "FL - gcBC", adjust.method = "BH", number = Inf)
diff_meth_FL_gcBC$Probe <- row.names(diff_meth_FL_gcBC)
diff_meth_FL_gcBC$type <- "gcBC"

# FL vs. tPC
diff_meth_FL_tPC <- topTable(fit2, coef = "FL - tPC", adjust.method = "BH", number = Inf)
diff_meth_FL_tPC$Probe <- row.names(diff_meth_FL_tPC)
diff_meth_FL_tPC$type <- "tPC"

# FL vs. bmPC
diff_meth_FL_bmPC <- topTable(fit2, coef = "FL - bmPC", adjust.method = "BH", number = Inf)
diff_meth_FL_bmPC$Probe <- row.names(diff_meth_FL_bmPC)
diff_meth_FL_bmPC$type <- "bmPC"

# FL vs. memBC
diff_meth_FL_memBC <- topTable(fit2, coef = "FL - memBC", adjust.method = "BH", number = Inf)
diff_meth_FL_memBC$Probe <- row.names(diff_meth_FL_memBC)
diff_meth_FL_memBC$type <- "memBC"

# FL vs. RLN
diff_meth_FL_RLN <- topTable(fit2, coef = "FL - RLN", adjust.method = "BH", number = Inf)
diff_meth_FL_RLN$Probe <- row.names(diff_meth_FL_RLN)
diff_meth_FL_RLN$type <- "RLN"

# FL vs. DLBCL
diff_meth_FL_DLBCL <- topTable(fit2, coef = "FL - DLBCL", adjust.method = "BH", number = Inf)
diff_meth_FL_DLBCL$Probe <- row.names(diff_meth_FL_DLBCL)
diff_meth_FL_DLBCL$type <- "DLBCL"

# Merged
rbind(diff_meth_FL_naiBC, diff_meth_FL_tnaiBC, diff_meth_FL_gcBC,
      diff_meth_FL_tPC, diff_meth_FL_bmPC, diff_meth_FL_memBC,
      diff_meth_FL_RLN, diff_meth_FL_DLBCL) %>%
  mutate(sig = ifelse(logFC > log.FC & adj.P.Val < adj.P.Value, "sig.up",
               ifelse(logFC < -log.FC & adj.P.Val < adj.P.Value, "sig.down", "not.sig"))) %>%
  mutate(type = factor(type, levels = c("tnaiBC", "naiBC", "gcBC", "tPC", "memBC", "bmPC", "RLN", "DLBCL"))) %>%
  mutate(sig = factor(sig, levels = c("sig.up", "not.sig", "sig.down"))) %>%
  group_by(type, sig) %>%
  summarize(count = n()) %>%
  mutate(type.cat = ifelse(type %in% c("RLN", "DLBCL"), "RLN, DLBCL", "BLUEPRINT")) %>%
  filter(!sig == "not.sig") %>%
  ggplot(aes(x = type, y = count, fill = sig)) +   
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(~ type.cat, scales = "free", space = "free") +
  theme_bw() +
  scale_fill_manual(values = c("#d6604d", "#4393c3"))

rbind(diff_meth_FL_naiBC, diff_meth_FL_tnaiBC, diff_meth_FL_gcBC,
      diff_meth_FL_tPC, diff_meth_FL_bmPC, diff_meth_FL_memBC) %>%
  mutate(sig = ifelse(logFC > log.FC & adj.P.Val < adj.P.Value, "sig.up",
                      ifelse(logFC < -log.FC & adj.P.Val < adj.P.Value, "sig.down", "not.sig"))) %>%
  mutate(type = factor(type, levels = c("tnaiBC", "naiBC", "gcBC", "tPC", "memBC", "bmPC"))) %>%
  mutate(sig = factor(sig, levels = c("sig.up", "not.sig", "sig.down"))) %>%
  group_by(type, sig) %>%
  summarize(count = n()) %>%
  filter(!sig == "not.sig") %>%
  ggplot(aes(x = type, y = count, fill = sig)) +   
  geom_bar(position = "stack", stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#d6604d", "#4393c3")) +
  scale_y_continuous(name = "Nb. differentially methylated probes .. FL vs.")

ggsave(paste0("img/",  date, " RK-Nb-diff-probes-Blueprint.pdf"), width = 12, height = 10, units = "cm")

# Hypermethylated in FL
hypermeth_FL_gcBC <- diff_meth_FL_gcBC %>%
  filter(logFC > log.FC & adj.P.Val < adj.P.Value) %>% .$Probe
hypermeth_FL_RLN <- diff_meth_FL_RLN %>%
  filter(logFC > log.FC & adj.P.Val < adj.P.Value) %>% .$Probe
hypermeth_FL_intersect <- intersect(hypermeth_FL_gcBC, hypermeth_FL_RLN)
hypermeth_FL_intersect_genes <- AnnotationFile %>%
  filter(V1 %in% hypermeth_FL_intersect)
hypermeth_FL_intersect_genes$UCSC_RefGene_Name <- gsub(";.*", "\\1", hypermeth_FL_intersect_genes$UCSC_RefGene_Name)
hypermeth_FL_intersect_genes <- hypermeth_FL_intersect_genes %>%
  filter(hypermeth_FL_intersect_genes$UCSC_RefGene_Name != "") %>%
  .$UCSC_RefGene_Name %>% unique()

# Hypomethylated in FL
hypometh_FL_gcBC <- diff_meth_FL_gcBC %>%
  filter(logFC < -log.FC & adj.P.Val < adj.P.Value) %>% .$Probe
hypometh_FL_RLN <- diff_meth_FL_RLN %>%
  filter(logFC < -log.FC & adj.P.Val < adj.P.Value) %>% .$Probe
hypometh_FL_intersect <- intersect(hypometh_FL_gcBC, hypometh_FL_RLN)
hypometh_FL_intersect_genes <- AnnotationFile %>%
  filter(V1 %in% hypometh_FL_intersect)
hypometh_FL_intersect_genes$UCSC_RefGene_Name <- gsub(";.*", "\\1", hypometh_FL_intersect_genes$UCSC_RefGene_Name)
hypometh_FL_intersect_genes <- hypometh_FL_intersect_genes %>%
  filter(hypometh_FL_intersect_genes$UCSC_RefGene_Name != "") %>%
  .$UCSC_RefGene_Name %>% unique()

not.diff.methylated <- diff_meth_FL_gcBC %>%
  filter(!Probe %in% c(hypermeth_FL_intersect_genes, hypometh_FL_intersect_genes)) %>%
  .$Probe

length(hypermeth_FL_intersect)
length(hypometh_FL_intersect)

length(hypermeth_FL_intersect_genes)
length(hypometh_FL_intersect_genes)

AnnotationFile %>%
  mutate(methyl.FL.status = ifelse(V1 %in% hypermeth_FL_intersect, "hypermethylated",
                            ifelse(V1 %in% hypometh_FL_intersect, "hypomethylated",
                            ifelse(V1 %in% not.diff.methylated, "other", NA)))) %>%
  filter(!is.na(methyl.FL.status)) %>%
  select(V1, Relation_to_Island, methyl.FL.status) %>%
  group_by(Relation_to_Island, methyl.FL.status) %>%
  summarize(proportion = n()) %>%
  mutate(methyl.FL.status = factor(methyl.FL.status, levels = c("hypomethylated", "other", "hypermethylated"))) %>%
  ggplot(aes(x = methyl.FL.status, y = proportion, fill = Relation_to_Island)) + 
  geom_bar(position = "fill", stat = "identity") +
  # geom_bar(position="stack", stat="identity") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values=c(
    "Island" = "#a6cee3", 
    "N_Shelf" = "#1f78b4", 
    "N_Shore" = "#fb9a99", 
    "OpenSea" = "#fdbf6f", 
    "S_Shelf" = "#ff7f00", 
    "S_Shore" = "#cab2d6"))

ggsave(paste0("img/",  date, " RK-diff-meth-probes-by-region.pdf"), width = 14, height = 10, units = "cm")

#---
# Testing for differential variability (DiffVar)
#---
pData.sub <- pData %>% filter(TYPE %in% c("FL", "RLN"))
FL.RLN <- pData.sub$SAMPLE_ID

design_sub <- design <- model.matrix(~0 + TYPE,
                                     data = pData.sub[,c("SAMPLE_ID", "TYPE")])
  
colnames(design_sub) <- c("FL", "RLN")

Combined_M_valueMatrix_sub_sub <- Combined_M_valueMatrix[,FL.RLN]

fitvar.contr <- varFit(Combined_M_valueMatrix_sub_sub, design = design_sub)

contMatrix_sub <- makeContrasts(FL-RLN, levels = design_sub)

fitvar.contr <- contrasts.varFit(fitvar.contr, contrasts = contMatrix_sub)
fitvar.contr <- eBayes(fitvar.contr)
summary(decideTests(fitvar.contr))

diff_meth_varfit <- topTable(fitvar.contr, adjust.method = "BH", number = Inf)

#--
# Enrichment of diff methyl CpGs within publicly available databases
#--

# Make GRanges
hypermeth_FL_intersect_df <- AnnotationFile %>%
  filter(V1 %in% hypermeth_FL_intersect) %>%
  select(chr = chr, start = pos, end = pos)
hypermeth_FL_intersect_df <- makeGRangesFromDataFrame(hypermeth_FL_intersect_df)

hypometh_FL_intersect_df <- AnnotationFile %>%
  filter(V1 %in% hypometh_FL_intersect) %>%
  select(chr = chr, start = pos, end = pos)
hypometh_FL_intersect_df <- makeGRangesFromDataFrame(hypometh_FL_intersect_df)

# Intersection with chromatin states
chrom_states_df <- read.table(file = "methylation/wgEncodeBroadHmmGm12878HMM_hg19.bed", header = F)
# downloaded from: http://compbio.mit.edu/ENCODE_chromatin_states/
chrom_states_vector <- chrom_states_df %>% .$V4 %>% unique() %>% as.character()

chrom_states_list <- list()
for (i in 1:length(chrom_states_vector)) {
  chrom_states_df_sub <- chrom_states_df %>%
    filter(V4 == chrom_states_vector[i]) %>%
    select(chr = V1, start = V2, end = V3)
  chrom_states_df_sub <- makeGRangesFromDataFrame(chrom_states_df_sub)
  chrom_states_list[i] <- chrom_states_df_sub
  print(chrom_states_vector[i])
}  

chrom_states_output_hypo <- list()
for (i in 1:length(chrom_states_vector)) {
  x <- countOverlaps(hypometh_FL_intersect_df, chrom_states_list[[i]])
  if (sum(x) == 0) {
    x <- table(x)
    x[2] <- 0
    names(x) <- c(0,1)
    x$class <- chrom_states_vector[i]
    chrom_states_output_hypo[[i]] <- x
    print(chrom_states_vector[i])
  } else {
    x <- table(x)
    x$class <- chrom_states_vector[i]
    chrom_states_output_hypo[[i]] <- x
    print(chrom_states_vector[i])
  }
}
chrom_states_output_hypo = do.call(rbind, chrom_states_output_hypo)
chrom_states_output_hypo <- data.frame(chrom_states_output_hypo)
chrom_states_output_hypo$direction <- "hypo"

chrom_states_output_hyper <- list()
for (i in 1:length(chrom_states_vector)) {
  x <- countOverlaps(hypermeth_FL_intersect_df, chrom_states_list[[i]])
  if (sum(x) == 0) {
    x <- table(x)
    x[2] <- 0
    names(x) <- c(0,1)
    x$class <- chrom_states_vector[i]
    chrom_states_output_hyper[[i]] <- x
    print(chrom_states_vector[i])
  } else {
    x <- table(x)
    x$class <- chrom_states_vector[i]
    chrom_states_output_hyper[[i]] <- x
    print(chrom_states_vector[i])
  }
}
chrom_states_output_hyper = do.call(rbind, chrom_states_output_hyper)
chrom_states_output_hyper <- data.frame(chrom_states_output_hyper)
chrom_states_output_hyper$direction <- "hyper"

chrom_states_output <- rbind(chrom_states_output_hyper, chrom_states_output_hypo)
colnames(chrom_states_output) <- c("not_overlap", "overlap", "class", "direction")
chrom_states_output <- chrom_states_output[,c(3,4,1,2)]
chrom_states_output$class <- as.character(chrom_states_output$class)
chrom_states_output$direction <- as.character(chrom_states_output$direction)
chrom_states_output$not_overlap <- as.numeric(chrom_states_output$not_overlap)
chrom_states_output$overlap <- as.numeric(chrom_states_output$overlap)
chrom_states_output$overlap_prop <- (chrom_states_output$overlap/(chrom_states_output$not_overlap+chrom_states_output$overlap))*100

# Plot percentage CpGs falling into given chromatin states
chrom_states_output %>%
  mutate(direction = factor(direction, levels = c("hypo", "hyper"))) %>%
  mutate(class = factor(class, levels = c("15_Repetitive/CNV",
                                          "14_Repetitive/CNV",
                                          "13_Heterochrom/lo",
                                          "12_Repressed",
                                          "11_Weak_Txn",
                                          "10_Txn_Elongation",
                                          "9_Txn_Transition",
                                          "8_Insulator",
                                          "7_Weak_Enhancer","6_Weak_Enhancer",
                                          "5_Strong_Enhancer",
                                          "4_Strong_Enhancer",
                                          "3_Poised_Promoter",
                                          "2_Weak_Promoter",
                                          "1_Active_Promoter"))) %>%
  ggplot(aes(direction, class, fill = overlap_prop)) +
  geom_raster() +
  geom_text(aes(label = round(overlap_prop, 1)), color = "white") +
  scale_fill_gradient(low = "#4393c3", high = "#b2182b") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(paste0("img/",  date, " RK-chrom-states.pdf"), width = 12, height = 8, units = "cm")

# Perform fisher test
chrom_states_output_fisher <- list()

for (i in 1:length(chrom_states_vector)) {
  a <- chrom_states_output %>%
    filter(class == chrom_states_vector[i])
  a <- cbind(a[,3],a[,4])
  fisher <- fisher.test(a)
  chrom_states_output_fisher[[i]] <- c("class" = chrom_states_vector[i],
                                       fisher$estimate,
                                       "P" = fisher$p.value)
} 

chrom_states_output_fisher <- do.call("rbind", chrom_states_output_fisher)
chrom_states_output_fisher <- data.frame(chrom_states_output_fisher)
chrom_states_output_fisher$P <- as.numeric(as.character(chrom_states_output_fisher$P))
chrom_states_output_fisher$fdr <- p.adjust(chrom_states_output_fisher$P, method = "fdr")

#--
# create bed files to use in GREAT
#--

hyper.CpGs.bed <- AnnotationFile %>%
  filter(V1 %in% hypermeth_FL_intersect) %>%
  select(chrom = chr, chromStart = pos, chromEnd = pos)
write.table(hyper.CpGs.bed, "hyper.CpGs.file.bed", row.names = F, col.names = F, quote = FALSE)

hypo.CpGs.bed <- AnnotationFile %>%
  filter(V1 %in% hypometh_FL_intersect) %>%
  select(chrom = chr, chromStart = pos, chromEnd = pos)
write.table(hypo.CpGs.bed, "hypo.CpGs.file.bed", row.names = F, col.names = F, quote = FALSE)

#--
# Enrichment of diff methyl CpGs within H3K27me3
#--
H3K27me3 <- read.table("methylation/merged_peaks_h3k27me3_karpas_wsudlcl2_ocily.bed", sep = "\t", header = FALSE) %>%
  select(chr = V1, start = V2, end = V3)
H3K27me3_GRanges <- makeGRangesFromDataFrame(H3K27me3)

hypermeth_FL_intersect_GRanges <- AnnotationFile %>%
  filter(V1 %in% hypermeth_FL_intersect) %>%
  select(chr = chr, start = pos, end = pos)
hypermeth_FL_intersect_GRanges <- makeGRangesFromDataFrame(hypermeth_FL_intersect_GRanges)

hypometh_FL_intersect_GRanges <- AnnotationFile %>%
  filter(V1 %in% hypometh_FL_intersect) %>%
  select(chr = chr, start = pos, end = pos)
hypometh_FL_intersect_GRanges <- makeGRangesFromDataFrame(hypometh_FL_intersect_GRanges)

pt <- overlapPermTest(A = H3K27me3_GRanges, B = hypermeth_FL_intersect_GRanges, ntimes = 100)
numOverlaps(A = H3K27me3_GRanges, B = hypermeth_FL_intersect_GRanges, count.once = TRUE)
lz <- localZScore(pt = pt, A = H3K27me3_GRanges, B = hypermeth_FL_intersect_GRanges)
pdf(paste0("img/",  date, " RK-overlapPermTest-H3K27me3-hypermethFL_1.pdf"), width = 6, height = 4)
plot(pt)
dev.off()
pdf(paste0("img/",  date, " RK-overlapPermTest-H3K27me3-hypermethFL_2.pdf"), width = 6, height = 4)
plot(lz)
dev.off()

# library(VennDiagram) 
# venn.diagram(list(B = 1:10375, A = 8920:30699), fill = c("lightblue", "green"), 
#              alpha = c(0.5, 0.5), lwd =0, "venn_diagram.tiff")

pt <- overlapPermTest(A = H3K27me3_GRanges, B = hypometh_FL_intersect_GRanges, ntimes = 100)
lz <- localZScore(pt = pt, A = H3K27me3_GRanges, B = hypometh_FL_intersect_GRanges)
pdf(paste0("img/",  date, " RK-overlapPermTest-H3K27me3-hypomethFL_1.pdf"), width = 6, height = 4)
plot(pt)
dev.off()
pdf(paste0("img/",  date, " RK-overlapPermTest-H3K27me3-hypomethFL_2.pdf"), width = 6, height = 4)
plot(lz)
dev.off()

library("venneuler")
plot(venneuler(c(A = 30699, B = 10375, "A&B" = 1456)))

#---
# Differentially methylated regions 
#---
myAnnotation_FL_gcBC <- cpg.annotate(object = Combined_M_valueMatrix, datatype = "array", what = "M", 
                                     analysis.type = "differential", design = design, 
                                     contrasts = TRUE, cont.matrix = contMatrix, 
                                     coef = "FL - gcBC", arraytype = "EPIC")
DMRs_FL_gcBC <- dmrcate(myAnnotation_FL_gcBC, lambda = 1000, C = 2)
DMRs_FL_gcBC_ranges <- extractRanges(DMRs_FL_gcBC, genome = "hg19")
DMRs_FL_gcBC_results <- data.frame(DMRs_FL_gcBC_ranges)

myAnnotation_FL_RLN <- cpg.annotate(object = Combined_M_valueMatrix, datatype = "array", what = "M", 
                                    analysis.type = "differential", design = design, 
                                    contrasts = TRUE, cont.matrix = contMatrix, 
                                    coef = "FL - RLN", arraytype = "EPIC")
DMRs_FL_RLN <- dmrcate(myAnnotation_FL_RLN, lambda = 1000, C = 2)
DMRs_FL_RLN_ranges <- extractRanges(DMRs_FL_RLN, genome = "hg19")
DMRs_FL_RLN_results <- data.frame(DMRs_FL_RLN_ranges)

# Hypermethylated regions in FL
meandiff.cutoff = 0.2
Stouffer.P.Value = 0.05 # 1.0e-20

protein.coding.genes <- read.table("annotables_grch37.txt", sep = "\t", header = TRUE) %>%
  filter(biotype == "protein_coding") %>% .$symbol %>% as.character() %>% unique()

hypermeth_FL_gcBC_regions <- DMRs_FL_gcBC_results %>%
  mutate(overlapping.genes = gsub(",.*$", "", overlapping.genes)) %>%
  filter(!is.na(overlapping.genes)) %>%
  filter(meandiff > meandiff.cutoff & Stouffer < Stouffer.P.Value) %>% .$overlapping.genes
hypermeth_FL_RLN_regions <- DMRs_FL_RLN_results %>%
  mutate(overlapping.genes = gsub(",.*$", "", overlapping.genes)) %>%
  filter(!is.na(overlapping.genes)) %>%
  filter(meandiff > meandiff.cutoff & Stouffer < Stouffer.P.Value) %>% .$overlapping.genes
hypermeth_FL_intersect_regions_genes <- sort(unique(intersect(hypermeth_FL_gcBC_regions, hypermeth_FL_RLN_regions)))
hypermeth_FL_intersect_regions_genes <- hypermeth_FL_intersect_regions_genes[hypermeth_FL_intersect_regions_genes %in% protein.coding.genes]
length(hypermeth_FL_intersect_regions_genes)  

hypometh_FL_gcBC_regions <- DMRs_FL_gcBC_results %>%
  mutate(overlapping.genes = gsub(",.*$", "", overlapping.genes)) %>%
  filter(!is.na(overlapping.genes)) %>%
  filter(meandiff < -meandiff.cutoff & Stouffer < Stouffer.P.Value) %>% .$overlapping.genes
hypometh_FL_RLN_regions <- DMRs_FL_RLN_results %>%
  mutate(overlapping.genes = gsub(",.*$", "", overlapping.genes)) %>%
  filter(!is.na(overlapping.genes)) %>%
  filter(meandiff < -meandiff.cutoff & Stouffer < Stouffer.P.Value) %>% .$overlapping.genes
hypometh_FL_intersect_regions_genes <- sort(unique(intersect(hypometh_FL_gcBC_regions, hypometh_FL_RLN_regions)))
hypometh_FL_intersect_regions_genes <- hypometh_FL_intersect_regions_genes[hypometh_FL_intersect_regions_genes %in% protein.coding.genes]
length(hypometh_FL_intersect_regions_genes)  

intersect(hypometh_FL_intersect_regions_genes, hypermeth_FL_intersect_regions_genes)

# Using EnrichR to explore DMRcate results
dbs <- listEnrichrDbs()
dbs <- c("ChEA_2016", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "ENCODE_Histone_Modifications_2015", "ENCODE_TF_ChIP-seq_2015",
         "Epigenomics_Roadmap_HM_ChIP-seq")
enrichr.results <- enrichr(hypermeth_FL_intersect_regions_genes, databases = dbs)

enrichr1 <- data.frame(enrichr.results[1]) %>%
  select(term = ChEA_2016.Term, adj.P.val = ChEA_2016.Adjusted.P.value) %>%
  mutate(minus.log10.P.val = -log10(adj.P.val)) %>%
  mutate(database = "ChEA_2016") %>%
  filter(!grepl("Mouse", term)) %>%
  head(10)

enrichr2 <- data.frame(enrichr.results[2])  %>%
  select(term = ENCODE_and_ChEA_Consensus_TFs_from_ChIP.X.Term, adj.P.val = ENCODE_and_ChEA_Consensus_TFs_from_ChIP.X.Adjusted.P.value) %>%
  mutate(minus.log10.P.val = -log10(adj.P.val)) %>%
  mutate(database = "ENCODE_and_ChEA_Consensus_TFs") %>%
  filter(!grepl("Mouse", term)) %>%
  head(10)

enrichr3 <- data.frame(enrichr.results[3])  %>%
  select(term = ENCODE_Histone_Modifications_2015.Term, adj.P.val = ENCODE_Histone_Modifications_2015.Adjusted.P.value) %>%
  mutate(minus.log10.P.val = -log10(adj.P.val)) %>%
  mutate(database = "ENCODE_Histone_Modifications") %>%
  filter(!grepl("Mouse", term)) %>%
  head(10)

enrichr4 <- data.frame(enrichr.results[4])  %>%
  select(term = ENCODE_TF_ChIP.seq_2015.Term, adj.P.val = ENCODE_TF_ChIP.seq_2015.Adjusted.P.value) %>%
  mutate(minus.log10.P.val = -log10(adj.P.val)) %>%
  mutate(database = "ENCODE_TF_ChIP-seq") %>%
  filter(!grepl("Mouse", term)) %>%
  head(10)

enrichr5 <- data.frame(enrichr.results[5])  %>%
  select(term = Epigenomics_Roadmap_HM_ChIP.seq.Term, adj.P.val = Epigenomics_Roadmap_HM_ChIP.seq.Adjusted.P.value) %>%
  mutate(minus.log10.P.val = -log10(adj.P.val)) %>%
  mutate(database = "Epigenomics_Roadmap_HM_ChIP.seq") %>%
  filter(!grepl("Mouse", term)) %>%
  head(10)

# Plot EnrichR results
rbind(enrichr1, enrichr2, enrichr3, enrichr4, enrichr5) %>%
  mutate(term = fct_reorder(term, minus.log10.P.val)) %>%
  ggplot(aes(x = term, y = minus.log10.P.val)) +
  geom_point(stat = 'identity', aes(col = database), size = 3)  +
  coord_flip() +
  theme_bw()

