#----------------------------------------------------------------------
#Strelka_009_processing_manta_results.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
              "ggrepel", "stringr", "maftools", "VariantAnnotation", "ggpubr")
lapply(packages, require, character.only = TRUE)
library(GenomicRanges)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
library(biomaRt)
library(openxlsx)
library(plotly)

date = Sys.Date()

print(date)

setwd("~/FLOMICS_Anjali/FLOMICS_Anajli/Data")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?

#----------------------------------------------------------------------
#data 
#----------------------------------------------------------------------

svs = fread(list.files(pattern="all_SVs_samples.txt")[1])

#id mapping
part1 = read.xlsx("part1_mapping_bc_dna.xlsx")
part2 = read.xlsx("part2_mapping_bc_dna.xlsx")
all_parts= rbind(part1, part2)
colnames(all_parts)[8] = "pat"

#sample information
sample_annotations = read.xlsx("sample_annotations_rcd6Nov2019.xlsx")
colnames(sample_annotations)[1] = "External_ID"
#clean up column
sample_annotations$External_ID = sapply(sample_annotations$External_ID, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
sample_annotations$External_ID = paste(sample_annotations$External_ID, sample_annotations$TIME_POINT, sep="_")

#merge datasets first by library ID then actual patient ID 
svs = merge(svs, all_parts, by="pat") ; svs = merge(svs, sample_annotations, by = "External_ID")

write.csv(svs, paste(date,"FLOMICS_MANTA_BC_SVs_with_annotations.csv", sep="_"), quote=F, row.names=F)

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#clean up translocations or any SVs with mates 
svs$id = as.character(svs$id)
snvs_wmates = svs[which(!(is.na(svs$MATE_BND_DEPTH))),]
snvs_orig = svs[which(svs$id %in% snvs_wmates$MATEID)]
all_mate_ids = as.character(unique(unique(snvs_wmates$id, snvs_orig$id)))
z = which(svs$id %in% all_mate_ids)
all_others = svs[-z,]

#distribution of types of structural variation
pats_svs = as.data.table(table(svs$SVTYPE, svs$pat)) ; pats_svs = pats_svs[order(-N)]

ggbarplot(pats_svs, x = "V2", y="N", fill="V1") +theme_bw() + 
  rotate_x_text(75) + ylab("Types of SVs") + xlab("Sample")

#genes most impacted by SVs 
genes_svs = as.data.table(table(svs$hgnc_symbol, svs$pat, svs$SVTYPE)) ; genes_svs = genes_svs[order(-N)]
genes_svs = as.data.table(filter(genes_svs, N >0, !(V1 == "")))
theme_update(text = element_text(size=7))

genes_sum = as.data.table(table(genes_svs$V1))
genes_sum = genes_sum[order(-N)]
genes_sum = as.data.table(filter(genes_sum, N > 6))
genes_svs$V1 = factor(genes_svs$V1, levels = genes_sum$V1)

p = ggplot(filter(genes_svs, V1 %in% genes_sum$V1), aes(V2, V1)) +
  geom_tile(aes(fill = V3), colour = "grey50") + 
  rotate_x_text(90)

#ggplotly(p)
pdf(paste(date,"FLOMICS_MANTA_BC_SVs_with_annotations.pdf", sep="_"), width=12, height=15) 
print(p)
dev.off()

