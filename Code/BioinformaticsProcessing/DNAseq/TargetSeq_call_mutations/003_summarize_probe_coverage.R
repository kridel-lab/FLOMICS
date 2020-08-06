#----------------------------------------------------------------------
#prepare intervals for picard collect metrics tool
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)
date=Sys.Date()

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TargetedDNAseq")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
"plyr",
"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#Interval files provided by IDT and processed by Robert
#----------------------------------------------------------------------

#Sample info
samp_info = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/library_mapping_BC.csv")

#detailed sample info
more_samp_info = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.csv")
more_samp_info = more_samp_info[,c("SAMPLE_ID", "rna_seq_file_sample_ID",
"RES_ID", "LY_FL_ID", "TIME_POINT", "PILOT_30_cases", "SEX", "INSTITUTION",
"TYPE", "STAGE" ,"COO", "Cluster", "RNAseq_DATA", "EPIC_INCLUDE")]

#coding
c_amp = fread("coding_genes_probe_coordinates_n_1917.txt")
c_amp$V1 = sapply(c_amp$V1, function(x){unlist(strsplit(x, "chr"))[2]})
c_amp = c_amp[,c(1:4)]
c_amp$type = "coding"

c_target = fread("coding_genes_target_coordinates_n_838.txt")
c_target$V1 = sapply(c_target$V1, function(x){unlist(strsplit(x, "chr"))[2]})
c_target = c_target[,c(1:4)]
c_target$type = "coding"

#noncoding
nc_amp = fread("non_coding_probe_coordinates_n_1481.txt")
nc_amp = nc_amp[,c(3:5,1)]
nc_amp$type = "noncoding"
colnames(nc_amp) = colnames(c_amp)
#nc_amp$V1 = paste("chr", nc_amp$V1, sep="")
nc_target = fread("non_coding_target_coordinates_n_35.txt")
nc_target$V1 = sapply(nc_target$V1, function(x){unlist(strsplit(x, "chr"))[2]})
nc_target = nc_target[,c(1:4)]
nc_target$type = "noncoding"

#collect all amplicon intervals
all_amps = rbind(c_amp, nc_amp)
colnames(all_amps) = c("chr", "start", "stop", "region", "type")

#collect all target intervals
all_targets = rbind(c_target, nc_target)
colnames(all_targets) = c("chr", "start", "stop", "region", "type")

#all files with coverage per target
all_res = list.files("/cluster/projects/kridelgroup/FLOMICS/DATA/BC_TargetSeq_Calls", pattern="coverage")

read_file = function(file_name){

  f = paste("/cluster/projects/kridelgroup/FLOMICS/DATA/BC_TargetSeq_Calls/", file_name, sep="")
  f = fread(f)
  sample = unlist(strsplit(file_name, "_"))[1]
  sample_check = unlist(strsplit(file_name, "-"))[1]

    if(!(sample_check == file_name)){
      sample = sample_check
    }

  f$Library = sample
  return(f)
}

all_res = as.data.table(ldply(llply(all_res, read_file)))
all_res$id = paste(all_res$chrom, all_res$end, sep="_")
all_targets$id = paste(all_targets$chr, all_targets$stop, sep="_")

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

all_res = merge(all_res, all_targets, by = "id")
colnames(all_res)[3] = "start_picard"
colnames(all_res)[18] = "start_target"

#merge with sample name
all_res = merge(all_res, samp_info, by="Library")

#save all_res
write.csv(all_res, file=paste(date,
  "picard_tools_coverage_summary_targets_DNA_sequencing.csv", sep="_"), quote=F, row.names=F)

#make boxplot comparing coding vs noncoding
pdf("/cluster/projects/kridelgroup/FLOMICS/DATA/001_coverage_mutation_summary.pdf")

# Change box plot colors by groups
p<-ggplot(all_res, aes(x=type, y=mean_coverage, fill=type)) +
  geom_boxplot()
p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

p<-ggplot(all_res, aes(x=type, y=normalized_coverage, fill=type)) +
  geom_boxplot()
p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

p<-ggplot(all_res, aes(x=type, y=min_coverage, fill=type)) +
  geom_boxplot()
p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

p<-ggplot(all_res, aes(x=type, y=min_normalized_coverage, fill=type)) +
  geom_boxplot()
p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

p<-ggplot(all_res, aes(x=type, y=max_coverage, fill=type)) +
  geom_boxplot()
p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

p<-ggplot(all_res, aes(x=type, y=max_normalized_coverage, fill=type)) +
  geom_boxplot()
p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

p<-ggplot(all_res, aes(x=type, y=read_count, fill=type)) +
  geom_boxplot()
p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dev.off()
