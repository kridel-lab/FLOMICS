#----------------------------------------------------------------------
#karin isaev
#post mutect2 and annovar soft filtering and matrix prep
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)
date=Sys.Date()

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)
library(readxl)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#mutect2 run on tumour only mode
#annotated by annovar

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#mutations
muts=fread(list.files(pattern="Mutect2")[length(list.files(pattern="Mutect2"))])

#samples
samp_info = as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/sample_annotations_rcd6Nov2019.xlsx"))

#bc mutation data
bc_mut_data = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/BC_mutation_data/BC_Cancer_capseq_data.csv")

#old mutation matrix
old_muts = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq/mutect2_full_data_w_AF_DP_2020-02-13_.txt")

#platypus matrix 

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------
