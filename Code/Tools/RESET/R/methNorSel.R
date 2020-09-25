###################################################################################################################################################################

# This function selects the hyper/hypomethylated probes (Sil/Enh probes) from the normal-sample dataset to study the DNA methylation:

# INPUTS:
  # 1. normal-samples data-set: matrix of beta-values
    ## rows: probeIDs
    ## columns: sample names
  # 2. list of selected probes: a data-frame consist of two colomns: (This can have a defult dataset attached to that)
    ## first column: list of probeID
    ## second column: list of the associated probe.index
# OUTPUT:
  # this function first extract the data regarding the selected probes, and later based on the beta-values prepare two seperate normal-sample data-sets:
    ## 1. hypermethylated probes in normal condition: Enh.probes
    ## 2. hypomethylated probes in normal condition: Sil.probes
###################################################################################################################################################################
#' Probe selection
#'
#' This function selects the hyper/hypomethylated probes (Sil/Enh probes) from the normal-sample dataset to study the DNA methylation:
#' @param normal.mtx matrix of beta-values (rows: probeIDs, columns: sample names)
#' @param probe.list a data-frame consist of at least two colomns: first column: list of probeID, second column: list of the associated probe.index
#' @return Two datasets: 1. hypermethylated probes in normal condition: Enh.probes, 2. hypomethylated probes in normal condition: Sil.probes
#' @export

methNorSel=function(normal.mtx,probe.list){
  # Selecting out the selected probes
    no=which(as.character(probe.list[,1])%in%rownames(normal.mtx))
  # generating the normal matrix of the selected probes
    row.n=as.character(probe.list[,1])[no]
    mtx.selected=normal.mtx[row.n,]
    rownames(mtx.selected)=as.character(probe.list[,2])[no]
  # Selecting the hypermethylated and hypomethylated probes in normal condition
    mean.probes=rowMeans(mtx.selected,na.rm = TRUE)
    var.probes=apply(mtx.selected, 1, function(x) {return(var(x,na.rm = TRUE))})
    
    sil.probes=which(mean.probes<=0.1&var.probes<=0.005)
    enh.probes=which(mean.probes>0.8&var.probes<=0.005)
    
  # generating the datasets
    normalmeth.sil.probes=mtx.selected[sil.probes,]
    normalmeth.enh.probes=mtx.selected[enh.probes,]
    
    normalmeth.list=list(normal.sil.probes=normalmeth.sil.probes,normal.enh.probes=normalmeth.enh.probes)
  return(normalmeth.list)
}


