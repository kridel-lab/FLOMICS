###################################################################################################################################################################

# This function evaluates the methylation status of tumor samples comparing to their normal condition.

# INPUTS:
  # 1. selected normal-sample datasets (the probe should have either sil or enh probes condition as specified in 'methNorSel' function)
    ## rows: probe index
    ## columns: normal samples
  # 2. tumor methylome datasets
    ## columns: tumor samples
    ## probeIDs
# OUTPUT: a list composed of 4 datasets and one vector: (all datasets has probeIndex in the rows and samples in column)
    # 1. Normal samples methylome data
    # 2. Tumor samples methylome data
    # 3. Tumor samples methylation status
    # 4. Specific Beta distribution (based on each specific probes value)
    # 5. Universal Beta estimation of normal-samples distribution (output of 'fitdistr' function)
###################################################################################################################################################################
# source('shared_functions.R')

#' Methylation Status Determination
#'
#' This function evaluates the methylation status of tumor samples comparing to their normal condition.
#' @param normal.db selected normal-sample datasets (the probe should have either sil or enh probes condition as specified in 'methNorSel' function. rows: probe index, columns: normal samples)
#' @param meth.tumor tumor methylome datasets (columns: tumor samples, rows: probeIDs)
#' @param methylation.event methylation event type. Silencing = 'sil',Enhancing = 'enh'
#' @return A list contains: 1. Normal samples methylome data, 2. Tumor samples methylome data, 3. Tumor samples methylation status, 4. Specific Beta distribution (based on each specific probes value), 5. Universal Beta estimation of normal-samples distribution (output of 'fitdistr' function)
#' @export
methStatus=function(normal.db,meth.tumor,methylation.event=c('sil','enh')){
  ## Generating the corespinding normal and tumor datasets
  probeId=sapply(rownames(normal.db), function(x){return(strsplit(x,'@')[[1]][1])})
  no=which(probeId%in%rownames(meth.tumor))
    indexname=rownames(normal.db)[no]
    probeid=probeId[no]
  
  normal=normal.db[indexname,]
  tumor=meth.tumor[probeid,]
    rownames(tumor)=rownames(normal)
  
  ##Beta distribution estimation
    #1)for each probe specific: from definition (Mean, Var => a,b)
      mean_normalsamples=rowMeans(normal,na.rm = TRUE)
      var_normalsamples=apply(normal,1,function(x){return(var(x,na.rm = TRUE))})
      alfa_normalsamples=(mean_normalsamples*(mean_normalsamples-(mean_normalsamples^2)-(var_normalsamples)))/var_normalsamples
      beta_normalsamples=(alfa_normalsamples/mean_normalsamples)-(alfa_normalsamples)
      df=data.frame(mean_normalsamples,var_normalsamples,alfa_normalsamples,beta_normalsamples)
      matrix_beta_dist_Probs_sep=as.matrix(df) #Specific Beta distribution (based on each specific probes value)
    #2) Universal distribution
      # to make the normal datasets without the repeated probeIDa
      normal.uniq=normal[which(!duplicated(probeid)),]
      meth.values=as.vector(normal.uniq)
      mean_values=mean(meth.values,na.rm = TRUE)
      var_values=var(meth.values,na.rm = TRUE)
      alfa_unmet=(mean_values*(mean_values-(mean_values^2)-(var_values)))/var_values
      beta_unmet=(alfa_unmet/mean_values)-(alfa_unmet)
      #fitdist function
      # library(MASS)
      beta_values = tryCatch(
      {
        beta_values=fitdistr(meth.values[which(!is.na(meth.values))],'beta',start=list(shape1=alfa_unmet,shape2=beta_unmet)) ## Universal estimation of normal-samples distribution 
      }, error = function(err) {
        # print(err)
        print('Reverting to empirical estimation of global beta density function using moments...')
        beta_values = list(c(alfa_unmet, beta_unmet))
      })
  
  
  
  ## extracting the methylation status of the tumor samples based on normal distributions
    #UNIVERSAL
      status.uni=t(apply(tumor,1, function(x){
        stat=sapply(x,function(y){
          if(methylation.event=='sil'){
            bin.call=ifelse(pbeta(y,beta_values[[1]][1],beta_values[[1]][2])>0.995,1,0)
          }else if(methylation.event=='enh'){
            bin.call=ifelse(pbeta(y,beta_values[[1]][1],beta_values[[1]][2])<0.005,1,0)
          }
          return(bin.call)
        })
        return(stat)
      }))
  
      
    #Specific
      status.spec=lapply(1:nrow(tumor), function(x){
        stat=sapply(tumor[x,],function(y){
          if(methylation.event=='sil'){
            bin.call=ifelse(pbeta(y,matrix_beta_dist_Probs_sep[x,3],matrix_beta_dist_Probs_sep[x,4])>0.995,1,0)
          }else if(methylation.event=='enh'){
            bin.call=ifelse(pbeta(y,matrix_beta_dist_Probs_sep[x,3],matrix_beta_dist_Probs_sep[x,4])<0.005,1,0)
          }
          return(bin.call)})
        return(stat)
      })
      status.spec=as.matrix(do.call(rbind,status.spec))
      rownames(status.spec)=rownames(status.uni)
  
    # generate the final methylation satatus: 1 in both univeral and specific.
      meth.status=ifelse((status.spec+status.uni)==2,1,0)

  ## RESULTS
    final.result=list(normal.meth=normal,tumor.meth=tumor,tumor.meth.status=meth.status,beta.dist.normal.specific=matrix_beta_dist_Probs_sep,beta.dist.normal.universal=beta_values)
    
    return(final.result)
}


