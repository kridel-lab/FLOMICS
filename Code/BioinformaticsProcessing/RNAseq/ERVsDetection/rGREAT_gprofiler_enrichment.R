#----------------------------------------------------------------------------------
#rGREAT_gprofiler_enrichment.R
#----------------------------------------------------------------------------------

#Sarah Russell
#Date: Sep 3, 2020
#This script takes differentially expressed ERVs from previous analysis and performs
#rGREAT analysis to obtain cis associated genes. Overlap between these genes and DEGs
#is then performed and used in gprofiler anaysis for pathway enrichment.

#----------------------------------------------------------------------------------
#PRE-PROCESSING
#----------------------------------------------------------------------------------
#Differential Expression for Genes: FLOMICS/Code/BioinformaticsProcessing/RNAseq/ERVsDetection/edgeR_GeneExpression.R

#To obtain Gene names for Ensembl IDs (Ensembl IDs used in DEG analysis):

#library(biomaRt)
#RNAseqCountMatrix=fread("STAR_quantmode_counts_matrix_FL_136_patients.txt")
#en_genes=data.frame(gene_id=RNAseqCountMatrix$gene)

#mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#genes <-  en_genes$gene_id
#gene_IDs <- getBM(attributes= c("ensembl_gene_id","external_gene_name","gene_biotype"),
#          mart= mart,values=genes)

#----------------------------------------------------------------------------------
#PACKAGES
#----------------------------------------------------------------------------------

# Loading needed packages
library(dplyr)
library(plyr)
library(data.table)
library(ggpubr)
library(readxl)
library(rGREAT)
library(gprofiler2)
#----------------------------------------------------------------------------------
#DATA
#----------------------------------------------------------------------------------

setwd("/Users/sarahrussell/R")

mat=fread("rgreat_coords.csv")
#concatenated file of DE ERVs for Tier 3, Adv-Lim
#> dim(mat)
#[1] 26754  10

#for DEGs
degs_stage=fread("_2020-08-26_AllDEIDs_stage.csv")
colnames(degs_stage)=c("logFC","logCPM","F","Pvalue","FDR","ensembl_gene_id")
degs_stage$contrast="ADVANCED_LIMITED"
degs_cluster=fread("_2020-08-26_AllDEIDs_cluster.csv")
colnames(degs_cluster)=c("logFC","logCPM","F","Pvalue","FDR","ensembl_gene_id")
degs_cluster$contrast="Cluster1_Cluster2"
degs=rbind(degs_stage,degs_cluster)
degs$diff_exp = ""
degs$diff_exp[degs$logFC >0 ] ="Upregulated"
degs$diff_exp[degs$logFC <0 ] ="Downregulated"
#> dim(degs)
#[1] 3920    8

#Format Gene Annotations
gene_annos=fread("gene_IDs_all.csv", header=T)
gene_annos$V1=NULL
#> dim(gene_annos)
#[1] 67130     3
gene_annos_filt=filter(gene_annos, ensembl_gene_id %in% degs$ensembl_gene_id & gene_biotype=="protein_coding")#
#> dim(gene_annos_filt)
#[1] 680   3

degs=filter(merge(degs,gene_annos_filt,by="ensembl_gene_id",all=T), gene_biotype == "protein_coding")
colnames(degs)[9]="gene"

#----------------------------------------------------------------------------------
#ANALYSIS
#----------------------------------------------------------------------------------
#how many protein coding DEGs
con=c("ADVANCED_LIMITED","Cluster1_Cluster2")
reg=c("Upregulated","Downregulated")
length(con)  
for(i in 1:length(con)){
cat("\n Contrast:", con[i], "has total of", nrow(degs[degs$contrast==con[i]]),
    "Protein-coding ENSEMBL IDs; # of Upregulated genes=", 
    nrow(degs[degs$contrast==con[i] & degs$diff_exp=="Upregulated"]), 
    "and Downregulated genes=", nrow(degs[degs$contrast==con[i] & degs$diff_exp=="Downregulated"]),".")
}
#Contrast: ADVANCED_LIMITED has total of 1960 Protein-coding ENSEMBL IDs; # of Upregulated genes= 562 and Downregulated genes= 1398 .
#Contrast: Cluster1_Cluster2 has total of 1960 Protein-coding ENSEMBL IDs; # of Upregulated genes= 912 and Downregulated genes= 1048 .

#########################
#Get rGREAT jobs for ERV contrasts
job_stage_up = submitGreatJob(mat[mat$contrast=="ADVANCED_LIMITED" & mat$diff_exp=="UP", ], request_interval=5)
job_cluster_up = submitGreatJob(mat[mat$contrast=="Cluster1_Cluster2" & mat$diff_exp=="UP", ], request_interval=5)
job_stage_down = submitGreatJob(mat[mat$contrast=="ADVANCED_LIMITED" & mat$diff_exp=="DOWN", ], request_interval=5)
job_cluster_down = submitGreatJob(mat[mat$contrast=="Cluster1_Cluster2" & mat$diff_exp=="DOWN", ], request_interval=5)

  
jobs=list(stage_up=job_stage_up,cluster_up=job_cluster_up,stage_down=job_stage_down,cluster_down=job_cluster_down)
jn=names(jobs)
get_rGREAT_genes = function(jobname){
  #peaks input is chr, start, end
	#pdf(paste("/Users/sarahrussell/",jobname,"_rGREAT.pdf"), width=9, height=4)
	j=jobs[[jobname]]
  res = plotRegionGeneAssociationGraphs(j)
	#dev.off()
	#tb = getEnrichmentTables(job, category = c("Genes"))
	#genes = as.data.table(tb[[1]])
	#genes = genes[order(-Binom_Fold_Enrichment)]

	genes = as.data.table(res)
	z = which(is.na(genes$gene))
	if(!(length(z)==0)){genes=genes[-z,]}
	#genes = as.data.table(filter(genes, abs(distTSS) < 5000))
	tot_genes = length(unique(genes$gene))
	print(tot_genes)
	return(genes)
}
all_genes = llply(jn,get_rGREAT_genes)
names(all_genes)=names(jobs)

get_rGREAT_pathways = function(jobname){
  j=jobs[[jobname]]
  
  #peaks input is chr, start, end
  tb = getEnrichmentTables(j, request_interval=5)

  #molecular function
  mf = tb[[1]] ; mf = as.data.table(mf)
  mf = as.data.table(filter(mf, Binom_Adjp_BH < 0.05))

  #biological process
  bp = tb[[2]] ; bp = as.data.table(bp)
  bp = as.data.table(filter(bp, Binom_Adjp_BH < 0.05))

  #pathways
  pathways = rbind(mf, bp) ;
  pathways = pathways[order(-Binom_Fold_Enrichment)]

  return(pathways)
}
all_pathways = llply(jn,get_rGREAT_pathways)
names(all_pathways)=names(jobs)
#########################
#Get overlaps for ERV associated Genes and DEGs
get_overlaps = function(con = NA,
                        expr = NA,
                        deg_con = NA){
  dat = degs[degs$contrast==con & degs$diff_exp==expr,]
  overlap = merge(all_genes[[deg_con]], dat, by = "gene")
  cat("\n Analysis Contrast:", con,"in", expr, "ERVs has total of", nrow(overlap),
      "Overlapping ERV associated genes (rGREAT) with DEGs in the following contrast:",
      deg_con,".", "unique overlapping genes=",length(unique(overlap$gene)),".")
  return(overlap)
}
Uervs_Udegs_stage=get_overlaps(con = "ADVANCED_LIMITED",expr = "Upregulated", deg_con = "stage_up" )
# Analysis Contrast: ADVANCED_LIMITED in Upregulated ERVs has total of 116 Overlapping ERV associated genes (rGREAT) with DEGs 
#in the following contrast: stage_up . unique overlapping genes= 44 .

Uervs_Ddegs_stage=get_overlaps(con = "ADVANCED_LIMITED",expr = "Upregulated", deg_con = "stage_down" )
# Analysis Contrast: ADVANCED_LIMITED in Upregulated ERVs has total of 29 Overlapping ERV associated genes (rGREAT) with DEGs 
#in the following contrast: stage_down . unique overlapping genes= 17 .

Dervs_Udegs_stage=get_overlaps(con = "ADVANCED_LIMITED",expr = "Downregulated", deg_con = "stage_up" )
# Analysis Contrast: ADVANCED_LIMITED in Downregulated ERVs has total of 184 Overlapping ERV associated genes (rGREAT) with DEGs 
#in the following contrast: stage_up . unique overlapping genes= 68 .

Dervs_Ddegs_stage=get_overlaps(con = "ADVANCED_LIMITED",expr = "Downregulated", deg_con = "stage_down" )
# Analysis Contrast: ADVANCED_LIMITED in Downregulated ERVs has total of 106 Overlapping ERV associated genes (rGREAT) with DEGs 
#in the following contrast: stage_down . unique overlapping genes= 50 .

Uervs_Udegs_cluster=get_overlaps(con = "Cluster1_Cluster2",expr = "Upregulated", deg_con = "cluster_up" )
# Analysis Contrast: Cluster1_Cluster2 in Upregulated ERVs has total of 417 Overlapping ERV associated genes (rGREAT) with DEGs 
#in the following contrast: cluster_up . unique overlapping genes= 71 .

Uervs_Ddegs_cluster=get_overlaps(con = "Cluster1_Cluster2",expr = "Upregulated", deg_con = "cluster_down" )
# Analysis Contrast: Cluster1_Cluster2 in Upregulated ERVs has total of 35 Overlapping ERV associated genes (rGREAT) with DEGs 
#in the following contrast: cluster_down . unique overlapping genes= 21 .

Dervs_Udegs_cluster=get_overlaps(con = "Cluster1_Cluster2",expr = "Downregulated", deg_con = "cluster_up" )
# Analysis Contrast: Cluster1_Cluster2 in Downregulated ERVs has total of 131 Overlapping ERV associated genes (rGREAT) with DEGs 
#in the following contrast: cluster_up . unique overlapping genes= 47 .

Dervs_Ddegs_cluster=get_overlaps(con = "Cluster1_Cluster2",expr = "Downregulated", deg_con = "cluster_down" )
# Analysis Contrast: Cluster1_Cluster2 in Downregulated ERVs has total of 59 Overlapping ERV associated genes (rGREAT) with DEGs 
#in the following contrast: cluster_down . unique overlapping genes= 27 .

all_overlaps=list(Uervs_Udegs_stage=Uervs_Udegs_stage,Uervs_Ddegs_stage=Uervs_Ddegs_stage,
                  Dervs_Udegs_stage=Dervs_Udegs_stage,Dervs_Ddegs_stage=Dervs_Ddegs_stage,
                  Uervs_Udegs_cluster=Uervs_Udegs_cluster,Uervs_Ddegs_cluster=Uervs_Ddegs_cluster,
                  Dervs_Udegs_cluster=Dervs_Udegs_cluster,Dervs_Ddegs_cluster=Dervs_Ddegs_cluster)

#########################
#functional enrichment analysis of overlap gene lists
overlap_con=names(all_overlaps)
get_enr_analysis = function(overlapdat){
  dat=all_overlaps[[overlapdat]]
  gostres = gprofiler2::gost(unique(dat$gene), organism = "hsapiens",
                             significant = TRUE, user_threshold = 0.05,
                             correction_method = "fdr",
                             sources=c("GO:BP", "REAC","TF"))
  gostres$result$fdr = p.adjust(gostres$result$p_value, method="fdr")
  return(gostres)
}
enriched_pathways = llply(overlap_con,get_enr_analysis)
names(enriched_pathways)=overlap_con


#overlap_con=names(all_overlaps)
get_plots_analysis = function(gostdat){
  #Manhattan-like-plot
  dat=enriched_pathways[[gostdat]]
  publish_gostplot(gostplot(dat, capped = T, interactive = F), width = NA, height = NA, filename = paste("/Users/sarahrussell/R/",gostdat,".pdf",sep = "_") )
  
  paths=dat$result
  paths$log10 = -log10(paths$p_value)
  paths = as.data.table(paths)
  paths_full = paths
  paths = paths[order(-log10)][1:20,]
  paths = paths[!is.na(paths$term_name), ]
  
  
  if(!(dim(paths)[1]==0)){
    g=ggdotchart(paths, x = "term_name", y = "log10",
                 color = "source",                                # Color by groups
                 palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
                 sorting = "descending",                       # Sort value in descending order
                 rotate = TRUE,                                # Rotate vertically
                 dot.size = 2,                                 # Large dot size
                 y.text.col = TRUE,                            # Color y text by groups
                 ggtheme = theme_pubr()                        # ggplot2 theme
    )+
      theme_cleveland()
    g=ggpar(g, font.y=6, font.tickslab=2)  + theme_bw() + ggtitle("")
    ggsave(filename = paste("/Users/sarahrussell/R/", "gprofiler_matched_genes_plot",gostdat, ".pdf",sep="_"), width=8, height=4)
  }
  return(paths)
}
path_top_20=llply(overlap_con,get_plots_analysis)

#plots uploaded to Teams: /FLOMICS/Analysis-Files/Telescope_ERVs/_2020_09_03_pathenrichment.pdf 
