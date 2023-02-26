
#--
# This script executes second-pass filtering upon variant calls
# post Mutect2 and Annovar analysis and prepares the SNV-by-sample matrix
# Author: Victoria Shelton
# Last modified: February 24th, 2023
# Tested on R/4.1.0
#--

#--
# Load in packages
#--

packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust",
 "data.table", "plyr", "ggrepel", "stringr", "maftools", "VariantAnnotation",
 "EnvStats", "ggpubr", "GenomicRanges", "Homo.sapiens",
 "TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts", "biomaRt")
lapply(packages, require, character.only = TRUE)

#--
# Setup working space
#--

date <- Sys.Date()

# set directory to where the "variant_calls.txt" matrix output by
## 001_Mutect2_prepare_matrix.R was deposited
setwd("/your working directory/")
dir <- "/your working directory/"

#--
# Read in data
#--

coverage_info <- read.delim("sample_based_coverage_summary.csv",
 header = TRUE, sep = ",")
snv_calls <- read.delim("variant_calls.txt", header = TRUE, sep = ";")
samp_info <- read.delim("sample_annotations.csv", header = TRUE, sep = ",")


#--
# Filtering out samples with low mean coverage,
## variants with SNP annotation and non-exonic variants/exonic variants
#--

#renaming samples/columns:
names(snv_calls)[names(snv_calls) == "Tumor_Sample_Barcode"] <- "Library"
names(samp_info)[1] <- "Library"

#merging datasets
snv_samp <- merge(samp_info, snv_calls, by = "Library", all = TRUE)
snv_samp_cov <- merge(snv_samp, coverage_info, by = "Library", all = TRUE)
dim(snv_samp_cov)

#filtering by sample coverage:
snv_samp_cov$mean_cov <- as.numeric(snv_samp_cov$mean_cov)
snv_samp_cov_filtered <- as.data.table(filter(snv_samp_cov, mean_cov >= 50))
dim(snv_samp_cov_filtered)
samples <- unique(snv_samp_cov_filtered$Library)
print("This is the number of samples with >50 coverage")
length(samples)

#filtering by SNP annotation
snv_samp_snp_filtered <- filter(snv_samp_cov_filtered,
 avsnp142 == "." | (avsnp142 != "." & cosmic68 != "."))

#re-labelingbelieved to be mis-annotated genes
snv_samp_snp_filtered$Hugo_Symbol[snv_samp_snp_filtered$Hugo_Symbol ==
 "CTD-2369P2.2"] <- "S1PR2"
snv_samp_snp_filtered$Hugo_Symbol[snv_samp_snp_filtered$Hugo_Symbol ==
 "RP11-211G3.2"] <- "BCL6"
snv_samp_snp_filtered$Hugo_Symbol[snv_samp_snp_filtered$Hugo_Symbol ==
 "RP11-132N15.3"] <- "BCL6"
snv_samp_snp_filtered$Hugo_Symbol[snv_samp_snp_filtered$Hugo_Symbol ==
 "RP11-132N15.1"] <- "BCL6"
snv_samp_snp_filtered$Hugo_Symbol[snv_samp_snp_filtered$Hugo_Symbol ==
 "RP4-655J12.5"] <- "CD58"

dim(snv_samp_snp_filtered)
print("This is the number of variants after removing SNPs")
length(snv_samp_snp_filtered$Library)

#------------------------------------------------------------------------------

#filtering out variants that ARE NOT exonic (non-exonic),
  #ARE "exonic synonymous SNVs",
  #& retaining variants that ARE "splicing ." => "EXONIC SNV DATA"
mut_flomics <- snv_samp_snp_filtered %>%
 filter(grepl("^exonic", Variant_Classification))

#selecting for exonic splice site variants
mut_flomics_splice <- filter(snv_samp_snp_filtered,
 Variant_Classification == "splicing .")
mut_flomics <- rbind(mut_flomics, mut_flomics_splice)
mut_flomics <- filter(mut_flomics,
 Variant_Classification != "exonic synonymous_SNV")
dim(mut_flomics)
print("This is the number of variants in exonic data matrix")
length(mut_flomics$Library)

#number of samples after filtering
  ##exonic SNVs (excluding synonymous SNVs)
samples <- unique(mut_flomics$Library)
print("This is the number of samples in EXONIC SNV DATA")
length(samples)

#number of samples LIMITED and ADVANCED
  ##exonic SNVs (excluding synonymous SNVs)
samples_mat <- dplyr::distinct(mut_flomics, Library, STAGE.y)
print("This is the number of LIMITED & ADVANCED samples in EXONIC SNV DATA")
dplyr::count(samples_mat, STAGE.y)

#write out tables
table_path_exonic <- paste0(dir, date, "_exonic_filtered_MUTECT2_calls.txt")
write.table(mut_flomics,
 file = table_path_exonic, quote = FALSE, row.names = FALSE, sep = ";")
print("EXONIC SNV DATA table created here:")
print(table_path_exonic)


#------------------------------------------------------------------------------

#--
#Analyze EXONIC SNV DATA in more detail
#--
count_table <- mut_flomics %>%
  group_by(Library) %>%
  dplyr::summarize(count = n())

mut_flomics_counts <- merge(mut_flomics,
 count_table, by = "Library", all = TRUE)

mut_flomics_summary <- unique(dplyr::select(mut_flomics_counts,
 Library, RES_ID, mean_cov, count))

# Plot number of variants per sample
a_plot <- ggplot(mut_flomics_summary,
 aes(Library, count)) +
  geom_col() +
   theme(axis.text.x = element_text(angle = 90)) +
    scale_y_continuous(breaks = seq(0, 35, by = 1)) +
     ggtitle("# of SNVs per Sample")
print(a_plot)

#Plot number of variants vs. coverage
b_plot <- ggplot(mut_flomics_summary,
 aes(mean_cov, count)) +
  geom_point() +
   geom_smooth(method = lm) +
    theme(axis.text.x = element_text(angle = 90)) +
     scale_y_continuous(breaks = seq(0, 35, by = 1)) +
      scale_x_continuous(breaks = seq(0, 1500, by = 200)) +
       ggtitle("# of SNVs vs. Sample Coverage")
print(b_plot)

mean_coverage <- summarize(mut_flomics_counts, mean = mean(mean_cov))
print("this is the mean coverage of samples within EXONIC SNV DATA")
print(mean_coverage)

#Plot number of mutations per gene
c_plot <- mut_flomics_counts %>%
  group_by(Hugo_Symbol, Library) %>%
  dplyr::summarize(count = n()) %>%
  group_by(Hugo_Symbol) %>%
  dplyr::summarize(count = n()) %>%
  arrange(dplyr::desc(count)) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, Hugo_Symbol)) %>%
  ggplot(aes(x = Hugo_Symbol, y = count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle("# of SNVs per Gene")
print(c_plot)

#write out exonic snv data plots
save_path <- paste0(dir, date, "_exonic_variants_plot.pdf")
pdf(save_path, width = 14, height = 8)
print(b_plot)
print(c_plot)
print(mean_coverage)
print(d_plot)
dev.off()

message("exonic snv data plots created")
