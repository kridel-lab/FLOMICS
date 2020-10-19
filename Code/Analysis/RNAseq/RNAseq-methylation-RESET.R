#-------------------------------------------------------------------------------
#RNA-seq-methylation-RESET.R
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#load functions and libraries
#-------------------------------------------------------------------------------

#getwd() --> FLOMICS teams folder
#cd /Users/kisaev/UHN/kridel-lab - Documents/FLOMICS

#-------------------------------------------------------------------------------
#data filtering
#-------------------------------------------------------------------------------

#methylation probes+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
meth = readRDS("methylation/2_BetaMatrix_updSamples_Ordered_T1_FilteredProbes.rds")
#only keep RLN and FL, remove dlbcl
z = which(str_detect(colnames(meth), "DLC"))
meth = meth[,-z] #155 FL + 5RN

#load the probe list++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load("Analysis-Files/RESET/promoter-probes-list.rdata")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#prepare input for resest

#1. first function is to prepare TSS probes that are consistently silenced or enhanced
#across normal samples

probes=promoter.probes.list
normal_matrix = meth[,c(which(colnames(meth) %in% c("LY_RLN_001","LY_RLN_002","LY_RLN_003",
"LY_RLN_004","LY_RLN_005")))]
dim(normal_matrix)
#[1] 595564      5

methNorSel_res = methNorSel(normal_matrix, probes) #output is a list of two (normalmeth.sil.probes and normalmeth.enh.probes)

#2. define methylation matrix for tumour samples
tum_matrix = meth[,c(which(!(colnames(meth) %in% c("LY_RLN_001","LY_RLN_002","LY_RLN_003",
"LY_RLN_004","LY_RLN_005"))))]
dim(tum_matrix)
#[1] 595564    155

#3. run main analysis using the reset function to get silencing/enhancing scores

#do everything using just the reset function
# INPUTS:
  # 1. selected normal-sample datasets (the probe should have either sil or enh probes condition as specified in 'methNorSel' function)
    ## rows: probe index
    ## columns: normal samples
  # 2. tumor methylome datasets
    ## columns: tumor samples
    ## probeIDs
  # 3. Transcriptome dataset of the tumor samples
    ## columns: tumor samples
    ## Gene name (as specified in the probe index)
  # 4. no.permutation for FDR calculation
  # 5. Methylation event type: either "sil" or "enh"

#reset=function(normal.db, meth.tumor, transcriptome,
#    methylation.event=c('sil','enh'),
#    FDR.permutation.no=100, seed=100)

reset_res_sil = reset(methNorSel_res[[1]], tum_matrix, tpm, FDR.permutation.no=100, "sil")
reset_res_enh = reset(methNorSel_res[[2]], tum_matrix, tpm, FDR.permutation.no=100, "enh")

saveRDS(reset_res_sil, "Analysis-Files/RESET/silening_events_results.rds")
saveRDS(reset_res_enh, "Analysis-Files/RESET/enhancing_events_results.rds")
