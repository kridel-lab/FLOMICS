#----------------------------------------------------------------------------------
#003_B_Stages_differential_expression_analysis.R - Stages
#----------------------------------------------------------------------------------

#Karin Isaev
#Date: January 16th, 2020
#This script takes in individual files obtained by telescope
#for each individual BAM file from RNA-Seq FL TGL13 samples
#and conducts EdgeR differential expression analysis of these trancripts

setwd("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/concatenated_results") #or where ever the 136 tsv files are stored
source("/cluster/home/srussell/FLOMICS/Code/BioinformaticsProcessing/RNAseq/ERVsDetection/003_B_Stages_DE_source.R")
#----------------------------------------------------------------------------------
#DATA
#----------------------------------------------------------------------------------

#The file below was created by script 003_A_get_ERV_counts_into_matrix.R
#It read in all 136 tsv files and put it together into one dataframe
#which we can use further
all_telescope = fread("/cluster/projects/kridelgroup/FLOMICS/DATA/_2020-07-08_TELESCOPE_OUTPUT_WITH_SAMPLE_ANNOTATION.csv")

#information regarding each sample and which stage of disease and cluster they are part of
sample_info = as.data.table(read_excel("/cluster/projects/kridelgroup/FLOMICS/DATA/Sample_Info/FL_TGL_STAR_logQC_2020-06-18_summary_KI_ClusterContamAdded.xlsx"))

#subset telescope results based on QC tiers
tiers=c("tier1", "tier2", "tier3")
#apply downstream code to each tier and then compare results
get_telescope_tier_summary = function(tier){
  if(tier == "tier3"){
    dat=sample_info
    print(length(unique(dat$SAMPLE_ID)))}
  if(tier == "tier2"){
    dat=sample_info
    z1 = which(dat$rRNAcontam > 40)
    z2 = which(dat$Uniquely.mapped < 10000000)
    z3 = which(dat$'Uniquely.mapped.reads..' < 50)
    z4 = which(dat$'X..of.reads.mapped.to.multiple.loci' > 20)
    total_lost = unique(c(z1, z2, z3, z4))
    dat=sample_info[-total_lost]
    print(length(unique(dat$SAMPLE_ID)))}
  if(tier == "tier1"){
    dat=sample_info
    z1 = which(dat$rRNAcontam > 40)
    z2 = which(dat$Uniquely.mapped < 10000000)
    z3 = which(dat$'Uniquely.mapped.reads..' < 70)
    z4 = which(dat$'X..of.reads.mapped.to.multiple.loci' > 20)
    total_lost = unique(c(z1, z2, z3, z4))
    dat=sample_info[-total_lost]
    print(length(unique(dat$SAMPLE_ID)))}
    tier_telescope_filter=as.data.table(filter(all_telescope, sample %in% dat$rna_seq_file_sample_ID))
    return(tier_telescope_filter)
    print("done tier analysis")
}
alltiers_telescope = llply(tiers, get_telescope_tier_summary)

vert_alltiers_telescope = llply(alltiers_telescope, vert_format)
#----------------------------------------------------------------------------------
#Differential expression analysis using EdgeR adapted from RNA-Seq code
#----------------------------------------------------------------------------------
#format for EdgeR
xtiers=llply(vert_alltiers_telescope,filter_removeNA)

print(ncol(xtiers[[1]]))
print(ncol(xtiers[[2]]))
print(ncol(xtiers[[3]]))

group_tiers=llply(vert_alltiers_telescope,groups_on_tiers)

#create DGEList object, normalize
tier_numbers=c(1,2,3)
df_tiers=llply(tier_numbers,make_DGEs)

#convert_CPM=function(dftiers){
#  exprs.df <- data.frame(cpm(dftiers, normalized.lib.sizes = FALSE, log = FALSE))
#  return(exprs.df)
#  ###should I write this out and to what directory?
#}
#exprs.df_tiers=llply(df_tiers,convert_CPM)

#construct design matrix
design_tiers=llply(tier_numbers,design_matrix)

#estimate dispersion
fit_tiers=llply(tier_numbers,get_disp_fit)

#make contrasts
my.contrasts_tier=llply(tier_numbers,get_myconstrasts)

#get differential expression results summary from all contrasts
all_contrasts_tier = llply(tier_numbers,get_allcontrasts)

# Set cut-offs for logFC and FDR
x <- 1 # logFCx
y <- 1 #keep all p-values for now can filter significant hits later

all_de_ervs_stage = llply(tier_numbers, get_res)

#colnames(all_de_ervs)[6] = "transcript"
#could just replace "mutate(ensembl_gene_id" with "mutate(transcript"
 #and exclude the above line?

#all_de_ervs = join(all_de_ervs,telescope_annotations)
#not sure what is meant by "telescope_annotations" above <-- this can be done
#to add back all the ERV coordinates and transcript info to main result file
#but not necessary!

#save results and upload to OneDrive
write.csv(as.data.table(all_de_ervs_stage[[1]]), paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope/", date, "tier1_de_ervs_stage.csv", sep="_"), quote=F, row.names=F)
write.csv(as.data.table(all_de_ervs_stage[[2]]), paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope/", date, "tier2_de_ervs_stage.csv", sep="_"), quote=F, row.names=F)
write.csv(as.data.table(all_de_ervs_stage[[3]]), paste("/cluster/projects/kridelgroup/FLOMICS/ANALYSIS/TELESCOPE_ANALYSIS/edgeR_telescope/", date, "tier3_de_ervs_stage.csv", sep="_"), quote=F, row.names=F)
