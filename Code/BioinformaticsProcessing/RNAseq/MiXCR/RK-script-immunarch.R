
library(dplyr)
library(data.table)
library(ggpubr)
library(immunarch)

date <- Sys.Date()

setwd("~/github/FLOMICS/")

immdata <- repLoad("Analysis-Files/Mixcr/mixcr_data_immunarch/")
immdata$meta$CLUSTER <- NULL

SNF.clust <- read.csv("Cluster Labels/InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv") %>%
  select(ID, SNFClust = SNFClust10Feb2021)

immdata$meta <- immdata$meta %>%
  left_join(SNF.clust[,c("ID", "SNFClust")], by = c("Sample" = "ID")) %>%
  mutate(SNF_incl = ifelse(is.na(SNFClust), "NO", "YES"))

clin <- read.csv("metadata/clinical_data_rcd11Aug2020.csv") %>%
  mutate(LY_FL_ID = paste0(LY_FL_ID, "_T1")) %>%
  mutate(POD24 = ifelse(CODE_TTP == 1 & TTP < 2, "YES", "NO")) %>%
  mutate(POD24_stage = paste0("POD24_", POD24, "_", TYPE))

immdata$meta <- immdata$meta %>%
  left_join(clin[,c("LY_FL_ID", "ANN_ARBOR_STAGE", "POD24_stage")], by = c("Sample" = "LY_FL_ID"))

#--
# Quality control metrics
#--

tiers <- read.csv("RNAseq/qc/Tiers1to3RNAseqFLOMICS_15Jan2021ASilva.csv")
tier1 <- tiers$T1[!is.na(tiers$T1)] # n = 81
tier2 <- tiers$T2[!is.na(tiers$T2)] # n = 104
tier3 <- tiers$T3[!is.na(tiers$T3)] # n = 132

qc <- read.csv("RNAseq/qc/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")

#--
# Filter dataset based on quality control metrics
#--

# immdata$data <- immdata$data[tier1]
# immdata$meta <- immdata$meta %>%
#   filter(Sample %in% tier1)

#--
# Filter for TRA/TRB V chains only, just looking at TCR set of data
#--

tr_imm = immdata
x = ''
for(i in 1:length(tr_imm$data)){
  # tr_imm$data[[i]] <- tr_imm$data[[i]][grepl("TRA|TRB|TRG|TRD",tr_imm$data[[i]]$V.name),]
  tr_imm$data[[i]] <- tr_imm$data[[i]][grepl("TRB",tr_imm$data[[i]]$V.name),]
  x[i] = nrow(tr_imm$data[[i]]) > 0
  z = names(tr_imm$data[as.logical(x)])
  dat = tr_imm$data[as.logical(x)]
  met = filter(tr_imm$meta, Sample %in% z)
  repertoire = list(data = dat, meta = met)
}

# tr_imm$data <- tr_imm$data[tier1]
# tr_imm$meta <- tr_imm$meta %>%
#   filter(Sample %in% tier1)

#--
# Descriptive
#---

# Number of unique clonotypes by sample
exp_vol <- repExplore(tr_imm$data, .method = "volume", .coding = TRUE)
# vis(exp_vol)
# p1 <- vis(exp_vol, .by = c("TYPE"), .meta = tr_imm$meta)
# p2 <- vis(exp_vol, .by = c("STAGE"), .meta = tr_imm$meta)
# p3 <- vis(exp_vol, .by = c("SNFClust"), .meta = tr_imm$meta)
# p1 + p2 + p3

# Number of clones by sample
clones.nb.sample <- list()
for (i in (1:length(tr_imm$data))) {
  x = names(tr_imm$data[i])
  y <- sum(data.frame(tr_imm$data[i])[,1])
  clones.nb.sample[[i]] <- c(Sample = x, nb = y)
  }
clones.nb.sample <- do.call(rbind, clones.nb.sample)
clones.nb.sample <- data.frame(clones.nb.sample)
clones.nb.sample$nb <- as.numeric(clones.nb.sample$nb)

# Repertoire overlap
# imm_ov1 <- repOverlap(tr_imm$data, .method = "public", .verbose = F)
# imm_ov2 <- repOverlap(tr_imm$data, .method = "morisita", .verbose = F)

# Gene usage computation
# imm_gu <- geneUsage(tr_imm$data, "hs.trbv")

#--
# Nb of clones vs. RNAseq nb uniquely mapped reads
#---

clones.nb.sample_qc <- clones.nb.sample %>%
  left_join(qc[,c("SAMPLE_ID", "Uniquely.mapped")], by = c("Sample" = "SAMPLE_ID")) %>%
  select(Sample, nb_clones = nb, nb_uniquely_mapped_reads = Uniquely.mapped)

clones.nb.sample_qc %>%
  ggscatter(x = "nb_uniquely_mapped_reads", y = "nb_clones", add = "reg.line", conf.int = TRUE,                                
            add.params = list(color = "blue", fill = "lightgray")) +
  stat_cor(method = "spearman")
# R = 0.53, i.e. correlation between nb_uniquely_mapped_reads and nb_unique_clonotypes

# -> normalize nb_clones by nb_uniquely_mapped_reads

clones.nb.sample_qc <- clones.nb.sample_qc %>%
  mutate(norm_factor = nb_uniquely_mapped_reads/sum(.$nb_uniquely_mapped_read)*nrow(.)) %>%
  mutate(nb_clones_norm = nb_clones/norm_factor) %>%
  left_join(immdata$meta)

p <- clones.nb.sample_qc %>%
  filter(Sample %in% tier3) %>%
  filter(TYPE != "RLN") %>%
  ggboxplot(x = "TYPE", y = "nb_clones_norm", color = "TYPE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

clones.nb.sample_qc %>%
  filter(Sample %in% tier3) %>%
  filter(TYPE != "RLN") %>%
  group_by(TYPE) %>%
  summarize(mean = mean(nb_clones_norm))

p <- clones.nb.sample_qc %>%
  filter(STAGE %in% c("ADVANCED", "LIMITED")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "STAGE", y = "nb_clones_norm", color = "STAGE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- clones.nb.sample_qc %>%
  filter(SNFClust %in% c("1", "2")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "SNFClust", y = "nb_clones_norm", color = "SNFClust", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- clones.nb.sample_qc %>%
  filter(POD24_stage %in% c("POD24_YES_ADVANCED", "POD24_NO_ADVANCED")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "POD24_stage", y = "nb_clones_norm", color = "POD24_stage", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

#--
# Unique clonotypes vs. RNAseq nb uniquely mapped reads
#---

exp_vol_qc <- exp_vol %>%
  left_join(qc[,c("SAMPLE_ID", "Uniquely.mapped")], by = c("Sample" = "SAMPLE_ID")) %>%
  select(Sample, nb_unique_clonotypes = Volume, nb_uniquely_mapped_reads = Uniquely.mapped)

exp_vol_qc %>%
  ggscatter(x = "nb_uniquely_mapped_reads", y = "nb_unique_clonotypes", add = "reg.line", conf.int = TRUE,                                
            add.params = list(color = "blue", fill = "lightgray")) +
  stat_cor(method = "spearman")
# R = 0.55, i.e. correlation between nb_uniquely_mapped_reads and nb_unique_clonotypes

# -> normalize nb_unique_clonotypes by nb_uniquely_mapped_reads

exp_vol_qc <- exp_vol_qc %>%
  mutate(norm_factor = nb_uniquely_mapped_reads/sum(.$nb_uniquely_mapped_read)*nrow(.)) %>%
  mutate(nb_unique_clonotypes_norm = nb_unique_clonotypes/norm_factor) %>%
  left_join(immdata$meta)

p <- exp_vol_qc %>%
  filter(Sample %in% tier3) %>%
  filter(TYPE != "RLN") %>%
  ggboxplot(x = "TYPE", y = "nb_unique_clonotypes_norm", color = "TYPE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

exp_vol_qc %>%
  filter(Sample %in% tier3) %>%
  filter(TYPE != "RLN") %>%
  group_by(TYPE) %>%
  summarize(mean = mean(nb_unique_clonotypes_norm))

p <- exp_vol_qc %>%
  filter(STAGE %in% c("ADVANCED", "LIMITED")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "STAGE", y = "nb_unique_clonotypes_norm", color = "STAGE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- exp_vol_qc %>%
  filter(SNFClust %in% c("1", "2")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "SNFClust", y = "nb_unique_clonotypes_norm", color = "SNFClust", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- exp_vol_qc %>%
  filter(POD24_stage %in% c("POD24_YES_ADVANCED", "POD24_NO_ADVANCED")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "POD24_stage", y = "nb_unique_clonotypes_norm", color = "POD24_stage", add = "jitter")
p + stat_compare_means(method = "wilcox.test")
  
#--
# Calculate diversity metrics
#---

div_chao <- as.data.frame(repDiversity(repertoire$data, "chao1"))
div_chao$Sample = rownames(div_chao)
div_chao = merge(div_chao,met,by="Sample")

div_hill <- repDiversity(repertoire$data, "hill")

div_div <- as.data.frame(repDiversity(repertoire$data, .col="aa", "div"))
div_div=merge(div_div,met,by="Sample")

div_gini.simp <- as.data.frame(repDiversity(repertoire$data, "gini.simp"))
div_gini.simp=merge(div_gini.simp,met,by="Sample")

div_gini <- as.data.frame(repDiversity(repertoire$data, .col="aa", "gini"))
div_gini$Sample=rownames(div_gini)
div_gini=merge(div_gini,met,by="Sample")

div_inv.simp <- as.data.frame(repDiversity(repertoire$data, .col="aa", "inv.simp"))
div_inv.simp=merge(div_inv.simp,met,by="Sample")

div_d50 <- as.data.frame(repDiversity(repertoire$data, .col="aa", "d50"))
div_d50$Sample=rownames(div_d50)
div_d50=merge(div_d50,met,by="Sample")

all_div = list(div_chao=div_chao,
               div_div=div_div,
               div_gini.simp=div_gini.simp,
               div_gini=div_gini,
               div_inv.simp=div_inv.simp,
               div_d50=div_d50)

for(i in 1:length(all_div)){
  g1=ggboxplot(all_div[[i]], x = "TYPE", y = names(all_div[[i]][2]),
               title = "sample diversity estimation", ylab = names(all_div[i]), xlab = "Sample Type",
               color = "TYPE", palette = "jco")
  my_comparisons1 <- list(c("DLBCL", "FL"))
  
  #######
  
  g2=ggboxplot(all_div[[i]][all_div[[i]]$TYPE=="FL",], x = "STAGE", y = names(all_div[[i]][2]),
               title = "sample diversity estimation", ylab = names(all_div[i]), xlab = "Sample Stage",
               color = "STAGE", palette = "jco")
  my_comparisons2 <- list( c("ADVANCED", "LIMITED"))
  
  ######
  
  g3=ggboxplot(all_div[[i]][all_div[[i]]$SNF_incl=="YES",], x = "SNFClust", y = names(all_div[[i]][2]),
               title = "sample diversity estimation", ylab = names(all_div[i]), xlab = "SNF Cluster",
               color = "SNFClust", palette = "jco")
  my_comparisons3 <- list( c(1, 2))
  g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed")
  
  ######
  
  (g1 + stat_compare_means(comparisons = my_comparisons1) + grids(linetype = "dashed")) +
  (g2 + stat_compare_means(comparisons = my_comparisons2) + grids(linetype = "dashed")) +
  (g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed"))
  
  ggsave(filename = paste0("img/TCR_diversity_", names(all_div[i]), ".pdf"), width=10, height=5)
}


################################################################################################################################################################################
#try sampling clones

#downsample would be similar to the sampling I had done with the vegan package
##basically expanding clone counts and choosing from a population of all clones, not
##all clonotypes - unclear if this is with or without replacement.. might sample with replacement

#sample chooses from clonotypes not clones, according to size if prob = T.

#should try downsampling, see what results look like, could also try sampling population similar to vegan package
#try n = 50,100,150

# Downsampling to 1000 clones (not clonotypes!)


x=''
y=''
z=''
for(i in 1:length(repertoire$data)){
  x[i]=sum(repertoire$data[[i]]$Clones) >= 50
  y[i]=sum(repertoire$data[[i]]$Clones) >= 100
  z[i]=sum(repertoire$data[[i]]$Clones) >= 150
}
pop=list(x,y,z)
get_samp_pops = function(pops){
  z=names(repertoire$data[as.logical(pops)])
  dat=repertoire$data[as.logical(pops)]
  met=filter(repertoire$meta, Sample %in% z)
  sample.pop=list(data=dat,meta=met)
  return(sample.pop)
}
sample_pops=plyr::llply(pop,get_samp_pops)

samp50 <- repSample(sample_pops[[1]]$data[c(1:length(sample_pops[[1]]$data))], .method = "downsample", .n=50)
#> length(sample_pops[[1]]$data) == length(samp50)  ##94
#TRUE
samp100 <- repSample(sample_pops[[2]]$data[c(1:length(sample_pops[[2]]$data))], .method = "downsample", .n=100)
#> length(sample_pops[[2]]$data) == length(samp100)   ##80
#[1] TRUE
samp150 <- repSample(sample_pops[[3]]$data[c(1:length(sample_pops[[3]]$data))], .method = "downsample", .n=150)
#> length(sample_pops[[3]]$data) == length(samp150)   ##60
#[1] TRUE
#is as expected based off previous sampling with vegan package

#calculate diversity metrics
get_sample_metrics = function(sampling,
                              metadata){
  
  div_chao <- as.data.frame(repDiversity(sampling, "chao1"))
  div_chao$Sample=rownames(div_chao)
  div_chao=merge(div_chao,metadata,by="Sample")
  
  #div_hill <- repDiversity(sampling, "hill")
  
  div_div <- as.data.frame(repDiversity(sampling, .col="aa", "div"))
  div_div=merge(div_div,metadata,by="Sample")
  
  div_gini.simp <- as.data.frame(repDiversity(sampling, "gini.simp"))
  div_gini.simp=merge(div_gini.simp,metadata,by="Sample")
  
  div_gini <- as.data.frame(repDiversity(sampling, .col="aa", "gini"))
  div_gini$Sample=rownames(div_gini)
  div_gini=merge(div_gini,metadata,by="Sample")
  
  div_inv.simp <- as.data.frame(repDiversity(sampling, .col="aa", "inv.simp"))
  div_inv.simp=merge(div_inv.simp,metadata,by="Sample")
  
  div_d50 <- as.data.frame(repDiversity(sampling, .col="aa", "d50"))
  div_d50$Sample=rownames(div_d50)
  div_d50=merge(div_d50,metadata,by="Sample")
  
  
  all_div = list(div_chao=div_chao,
                 div_div=div_div,
                 div_gini.simp=div_gini.simp,
                 div_gini=div_gini,
                 div_inv.simp=div_inv.simp,
                 div_d50=div_d50)
  return(all_div)
}
diversity.50 = get_sample_metrics(samp50,sample_pops[[1]]$meta)
diversity.100 = get_sample_metrics(samp100,sample_pops[[2]]$meta)
diversity.150 = get_sample_metrics(samp150,sample_pops[[3]]$meta)

#visualize sampling

get_sample_plots = function(samp_mets){
  
  for(i in 1:length(samp_mets)){
    g1=ggboxplot(samp_mets[[i]], x = "TYPE", y = names(samp_mets[[i]][2]),
                 title = "sample diversity estimation", ylab = names(samp_mets[i]), xlab = "Sample Type",
                 color = "TYPE", palette = "jco")
    my_comparisons1 <- list( c("DLBCL", "FL"))
    
    #######
    g2=ggboxplot(samp_mets[[i]][samp_mets[[i]]$TYPE=="FL",], x = "STAGE", y = names(samp_mets[[i]][2]),
                 title = "sample diversity estimation", ylab = names(samp_mets[i]), xlab = "Sample Stage",
                 color = "STAGE", palette = "jco")
    my_comparisons2 <- list( c("ADVANCED", "LIMITED"))
    
    ######
    g3=ggboxplot(samp_mets[[i]][samp_mets[[i]]$SNF_incl=="YES",], x = "SNFClust", y = names(samp_mets[[i]][2]),
                 title = "sample diversity estimation", ylab = names(samp_mets[i]), xlab = "Sample Cluster",
                 color = "SNFClust", palette = "jco")
    my_comparisons3 <- list( c(1, 2))
    ######
    (g1 + stat_compare_means(comparisons = my_comparisons1) + grids(linetype = "dashed")) +
      (g2 + stat_compare_means(comparisons = my_comparisons2) + grids(linetype = "dashed")) +
      (g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed"))
    

    ggsave(filename = paste0("img/TCR_diversity_downsampled_100", names(samp_mets[i]), ".pdf"), width=10, height=5)
  }
}
get_sample_plots(diversity.50)
get_sample_plots(diversity.100)
get_sample_plots(diversity.150)

####
# BCR
####

#--
# Filter for IG V chains only, just looking at IG set of data
#--

br_imm = immdata
x = ''
for(i in 1:length(br_imm$data)){
  br_imm$data[[i]] <- br_imm$data[[i]][grepl("IG",br_imm$data[[i]]$V.name),]
  x[i] = nrow(br_imm$data[[i]]) > 0
  z = names(br_imm$data[as.logical(x)])
  dat = br_imm$data[as.logical(x)]
  met = filter(br_imm$meta, Sample %in% z)
  repertoire = list(data = dat, meta = met)
}

# br_imm$data <- br_imm$data[tier1]
# br_imm$meta <- br_imm$meta %>%
#   filter(Sample %in% tier1)

#--
# Descriptive
#---

# Number of unique clonotypes by sample
exp_vol <- repExplore(br_imm$data, .method = "volume", .coding = TRUE)
vis(exp_vol)
p1 <- vis(exp_vol, .by = c("TYPE"), .meta = br_imm$meta)
p2 <- vis(exp_vol, .by = c("STAGE"), .meta = br_imm$meta)
p3 <- vis(exp_vol, .by = c("SNFClust"), .meta = br_imm$meta)
p1 + p2 + p3

# Number of clones by sample
clones.nb.sample <- list()
for (i in (1:length(br_imm$data))) {
  x = names(br_imm$data[i])
  y <- sum(data.frame(br_imm$data[i])[,1])
  clones.nb.sample[[i]] <- c(Sample = x, nb = y)
}
clones.nb.sample <- do.call(rbind, clones.nb.sample)
clones.nb.sample <- data.frame(clones.nb.sample)
clones.nb.sample$nb <- as.numeric(clones.nb.sample$nb)

# Repertoire overlap
# imm_ov1 <- repOverlap(br_imm$data, .method = "public", .verbose = F)
# imm_ov2 <- repOverlap(br_imm$data, .method = "morisita", .verbose = F)

# Gene usage computation
# imm_gu <- geneUsage(br_imm$data, "macmul.IGHV")

#--
# Nb of clones vs. RNAseq nb uniquely mapped reads
#---

clones.nb.sample_qc <- clones.nb.sample %>%
  left_join(qc[,c("SAMPLE_ID", "Uniquely.mapped")], by = c("Sample" = "SAMPLE_ID")) %>%
  select(Sample, nb_clones = nb, nb_uniquely_mapped_reads = Uniquely.mapped)

clones.nb.sample_qc %>%
  ggscatter(x = "nb_uniquely_mapped_reads", y = "nb_clones", add = "reg.line", conf.int = TRUE,                                
            add.params = list(color = "blue", fill = "lightgray")) +
  stat_cor(method = "spearman")
# R = 0.6, i.e. correlation between nb_uniquely_mapped_reads and nb_unique_clonotypes

# -> normalize nb_clones by nb_uniquely_mapped_reads

clones.nb.sample_qc <- clones.nb.sample_qc %>%
  mutate(norm_factor = nb_uniquely_mapped_reads/sum(.$nb_uniquely_mapped_read)*nrow(.)) %>%
  mutate(nb_clones_norm = nb_clones/norm_factor) %>%
  left_join(immdata$meta)

p <- clones.nb.sample_qc %>%
  filter(Sample %in% tier3) %>%
  filter(TYPE != "RLN") %>%
  ggboxplot(x = "TYPE", y = "nb_clones_norm", color = "TYPE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

clones.nb.sample_qc %>%
  filter(Sample %in% tier3) %>%
  filter(TYPE != "RLN") %>%
  group_by(TYPE) %>%
  summarize(mean = mean(nb_clones_norm))

p <- clones.nb.sample_qc %>%
  filter(STAGE %in% c("ADVANCED", "LIMITED")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "STAGE", y = "nb_clones_norm", color = "STAGE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- clones.nb.sample_qc %>%
  filter(SNFClust %in% c("1", "2")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "SNFClust", y = "nb_clones_norm", color = "SNFClust", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- clones.nb.sample_qc %>%
  filter(POD24_stage %in% c("POD24_YES_ADVANCED", "POD24_NO_ADVANCED")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "POD24_stage", y = "nb_clones_norm", color = "POD24_stage", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

#--
# Unique clonotypes vs. RNAseq nb uniquely mapped reads
#---

exp_vol_qc <- exp_vol %>%
  left_join(qc[,c("SAMPLE_ID", "Uniquely.mapped")], by = c("Sample" = "SAMPLE_ID")) %>%
  select(Sample, nb_unique_clonotypes = Volume, nb_uniquely_mapped_reads = Uniquely.mapped)

exp_vol_qc %>%
  ggscatter(x = "nb_uniquely_mapped_reads", y = "nb_unique_clonotypes", add = "reg.line", conf.int = TRUE,                                
            add.params = list(color = "blue", fill = "lightgray")) +
  stat_cor(method = "spearman")
# R = 0.6, i.e. correlation between nb_uniquely_mapped_reads and nb_unique_clonotypes

# -> normalize nb_unique_clonotypes by nb_uniquely_mapped_reads

exp_vol_qc <- exp_vol_qc %>%
  mutate(norm_factor = nb_uniquely_mapped_reads/sum(.$nb_uniquely_mapped_read)*nrow(.)) %>%
  mutate(nb_unique_clonotypes_norm = nb_unique_clonotypes/norm_factor) %>%
  left_join(immdata$meta)

p <- exp_vol_qc %>%
  filter(Sample %in% tier3) %>%
  filter(TYPE != "RLN") %>%
  ggboxplot(x = "TYPE", y = "nb_unique_clonotypes_norm", color = "TYPE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

exp_vol_qc %>%
  filter(Sample %in% tier3) %>%
  filter(TYPE != "RLN") %>%
  group_by(TYPE) %>%
  summarize(mean = mean(nb_unique_clonotypes_norm))

p <- exp_vol_qc %>%
  filter(STAGE %in% c("ADVANCED", "LIMITED")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "STAGE", y = "nb_unique_clonotypes_norm", color = "STAGE", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- exp_vol_qc %>%
  filter(SNFClust %in% c("1", "2")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "SNFClust", y = "nb_unique_clonotypes_norm", color = "SNFClust", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

p <- exp_vol_qc %>%
  filter(POD24_stage %in% c("POD24_YES_ADVANCED", "POD24_NO_ADVANCED")) %>%
  filter(Sample %in% tier3) %>%
  ggboxplot(x = "POD24_stage", y = "nb_unique_clonotypes_norm", color = "POD24_stage", add = "jitter")
p + stat_compare_means(method = "wilcox.test")

#--
# Calculate diversity metrics
#---

div_chao <- as.data.frame(repDiversity(repertoire$data, "chao1"))
div_chao$Sample = rownames(div_chao)
div_chao = merge(div_chao,met,by="Sample")

div_hill <- repDiversity(repertoire$data, "hill")

div_div <- as.data.frame(repDiversity(repertoire$data, .col="aa", "div"))
div_div=merge(div_div,met,by="Sample")

div_gini.simp <- as.data.frame(repDiversity(repertoire$data, "gini.simp"))
div_gini.simp=merge(div_gini.simp,met,by="Sample")

div_gini <- as.data.frame(repDiversity(repertoire$data, .col="aa", "gini"))
div_gini$Sample=rownames(div_gini)
div_gini=merge(div_gini,met,by="Sample")

div_inv.simp <- as.data.frame(repDiversity(repertoire$data, .col="aa", "inv.simp"))
div_inv.simp=merge(div_inv.simp,met,by="Sample")

div_d50 <- as.data.frame(repDiversity(repertoire$data, .col="aa", "d50"))
div_d50$Sample=rownames(div_d50)
div_d50=merge(div_d50,met,by="Sample")

all_div = list(div_chao=div_chao,
               div_div=div_div,
               div_gini.simp=div_gini.simp,
               div_gini=div_gini,
               div_inv.simp=div_inv.simp,
               div_d50=div_d50)

for(i in 1:length(all_div)){
  g1=ggboxplot(all_div[[i]], x = "TYPE", y = names(all_div[[i]][2]),
               title = "sample diversity estimation", ylab = names(all_div[i]), xlab = "Sample Type",
               color = "TYPE", palette = "jco")
  my_comparisons1 <- list( c("DLBCL", "FL"), c("FL", "RLN"), c("DLBCL", "RLN") )
  
  #######
  
  g2=ggboxplot(all_div[[i]][all_div[[i]]$TYPE=="FL",], x = "STAGE", y = names(all_div[[i]][2]),
               title = "sample diversity estimation", ylab = names(all_div[i]), xlab = "Sample Stage",
               color = "STAGE", palette = "jco")
  my_comparisons2 <- list( c("ADVANCED", "LIMITED"))
  
  ######
  
  g3=ggboxplot(all_div[[i]][all_div[[i]]$SNF_incl=="YES",], x = "SNFClust", y = names(all_div[[i]][2]),
               title = "sample diversity estimation", ylab = names(all_div[i]), xlab = "SNF Cluster",
               color = "SNFClust", palette = "jco")
  my_comparisons3 <- list( c(1, 2))
  g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed")
  
  ######
  
  (g1 + stat_compare_means(comparisons = my_comparisons1) + grids(linetype = "dashed")) +
    (g2 + stat_compare_means(comparisons = my_comparisons2) + grids(linetype = "dashed")) +
    (g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed"))
  
  ggsave(filename = paste0("img/BCR_diversity_", names(all_div[i]), ".pdf"), width=10, height=5)
}


################################################################################################################################################################################
#try sampling clones

#downsample would be similar to the sampling I had done with the vegan package
##basically expanding clone counts and choosing from a population of all clones, not
##all clonotypes - unclear if this is with or without replacement.. might sample with replacement

#sample chooses from clonotypes not clones, according to size if prob = T.

#should try downsampling, see what results look like, could also try sampling population similar to vegan package
#try n = 50,100,150

# Downsampling to 1000 clones (not clonotypes!)


x=''
y=''
z=''
for(i in 1:length(repertoire$data)){
  x[i]=sum(repertoire$data[[i]]$Clones) >= 50
  y[i]=sum(repertoire$data[[i]]$Clones) >= 100
  z[i]=sum(repertoire$data[[i]]$Clones) >= 150
}
pop=list(x,y,z)
get_samp_pops = function(pops){
  z=names(repertoire$data[as.logical(pops)])
  dat=repertoire$data[as.logical(pops)]
  met=filter(repertoire$meta, Sample %in% z)
  sample.pop=list(data=dat,meta=met)
  return(sample.pop)
}
sample_pops=plyr::llply(pop,get_samp_pops)

samp50 <- repSample(sample_pops[[1]]$data[c(1:length(sample_pops[[1]]$data))], .method = "downsample", .n=50)
#> length(sample_pops[[1]]$data) == length(samp50)  ##123
#TRUE
samp100 <- repSample(sample_pops[[2]]$data[c(1:length(sample_pops[[2]]$data))], .method = "downsample", .n=100)
#> length(sample_pops[[2]]$data) == length(samp100)   ##117
#[1] TRUE
samp150 <- repSample(sample_pops[[3]]$data[c(1:length(sample_pops[[3]]$data))], .method = "downsample", .n=150)
#> length(sample_pops[[3]]$data) == length(samp150)   ##114
#[1] TRUE
#is as expected based off previous sampling with vegan package

#calculate diversity metrics
get_sample_metrics = function(sampling,
                              metadata){
  
  div_chao <- as.data.frame(repDiversity(sampling, "chao1"))
  div_chao$Sample=rownames(div_chao)
  div_chao=merge(div_chao,metadata,by="Sample")
  
  #div_hill <- repDiversity(sampling, "hill")
  
  div_div <- as.data.frame(repDiversity(sampling, .col="aa", "div"))
  div_div=merge(div_div,metadata,by="Sample")
  
  div_gini.simp <- as.data.frame(repDiversity(sampling, "gini.simp"))
  div_gini.simp=merge(div_gini.simp,metadata,by="Sample")
  
  div_gini <- as.data.frame(repDiversity(sampling, .col="aa", "gini"))
  div_gini$Sample=rownames(div_gini)
  div_gini=merge(div_gini,metadata,by="Sample")
  
  div_inv.simp <- as.data.frame(repDiversity(sampling, .col="aa", "inv.simp"))
  div_inv.simp=merge(div_inv.simp,metadata,by="Sample")
  
  div_d50 <- as.data.frame(repDiversity(sampling, .col="aa", "d50"))
  div_d50$Sample=rownames(div_d50)
  div_d50=merge(div_d50,metadata,by="Sample")
  
  
  all_div = list(div_chao=div_chao,
                 div_div=div_div,
                 div_gini.simp=div_gini.simp,
                 div_gini=div_gini,
                 div_inv.simp=div_inv.simp,
                 div_d50=div_d50)
  return(all_div)
}
diversity.50 = get_sample_metrics(samp50,sample_pops[[1]]$meta)
diversity.100 = get_sample_metrics(samp100,sample_pops[[2]]$meta)
diversity.150 = get_sample_metrics(samp150,sample_pops[[3]]$meta)

#visualize sampling

get_sample_plots = function(samp_mets){
  
  for(i in 1:length(samp_mets)){
    g1=ggboxplot(samp_mets[[i]], x = "TYPE", y = names(samp_mets[[i]][2]),
                 title = "sample diversity estimation", ylab = names(samp_mets[i]), xlab = "Sample Type",
                 color = "TYPE", palette = "jco")
    my_comparisons1 <- list( c("DLBCL", "FL"))
    
    #######
    g2=ggboxplot(samp_mets[[i]][samp_mets[[i]]$TYPE=="FL",], x = "STAGE", y = names(samp_mets[[i]][2]),
                 title = "sample diversity estimation", ylab = names(samp_mets[i]), xlab = "Sample Stage",
                 color = "STAGE", palette = "jco")
    my_comparisons2 <- list( c("ADVANCED", "LIMITED"))
    
    ######
    g3=ggboxplot(samp_mets[[i]][samp_mets[[i]]$SNF_incl=="YES",], x = "SNFClust", y = names(samp_mets[[i]][2]),
                 title = "sample diversity estimation", ylab = names(samp_mets[i]), xlab = "Sample Cluster",
                 color = "SNFClust", palette = "jco")
    my_comparisons3 <- list( c(1, 2))
    ######
    (g1 + stat_compare_means(comparisons = my_comparisons1) + grids(linetype = "dashed")) +
      (g2 + stat_compare_means(comparisons = my_comparisons2) + grids(linetype = "dashed")) +
      (g3 + stat_compare_means(comparisons = my_comparisons3) + grids(linetype = "dashed"))
    
    
    ggsave(filename = paste0("img/BCR_diversity_downsampled_100", names(samp_mets[i]), ".pdf"), width=10, height=5)
  }
}
get_sample_plots(diversity.50)
get_sample_plots(diversity.100)
get_sample_plots(diversity.150)

