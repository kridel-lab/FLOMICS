###################################################################################################################################################################

# This function estimates the false discovery rate (FDR) of the calculated scores, we adopted the strategy previously proposed in Significance Analysis of Microarrays (SAM) by Tusher et al.It permutates (n=100) the transcriptome data. 

# INPUTS:
  # 1. score.values : the score calculated for each probe (output of "eventScore" function)
  # 2. meth.tumor : methylome data used in score calculation
  # 3. meth.tumor.status : methylation status data used in score calculation
  # 4. transcriptome : the transcriptome data used for score calculation (ranked normally distributed transcriptome data)
  # 5. number of desired permutation
# OUTPUT: 
  # 1. FDR.res : table containing the FDR calculated for each observed score
  # 2. permutation.res : the results of (n=100) permutations + mean permutated score + observed scores
  # 3. score.cutoff : the suggested score cutoff: score > 1.5 & FDR < 0.1 & false score < 0.5
###################################################################################################################################################################
# source('shared_functions.R')

#' FDR Calculation
#'
#' This function estimates the false discovery rate (FDR) of the calculated scores, we adopted the strategy previously proposed in Significance Analysis of Microarrays (SAM) by Tusher et al.It permutates (n=100) the transcriptome data.
#' @param score.values The score calculated for each probe (output of "eventScore" function)
#' @param meth.tumor methylome data used in score calculation
#' @param meth.tumor.status methylation status data used in score calculation
#' @param transcriptome the transcriptome data used for score calculation (ranked normally distributed transcriptome data)
#' @param no.permutation number of desired permutation(default=100)
#' @param methylation.event methylation event type. Silencing = 'sil',Enhancing = 'enh'
#' @return A list contains: 1. FDR.res : table containing the FDR calculated for each observed score, 2. permutation.res : the results of (n=100) permutations + mean permutated score + observed scores, 3. score.cutoff : the suggested score cutoff: score > 1.5 & FDR < 0.1 & false score < 0.5
#' @export
FDRcal=function(score.values,meth.tumor,meth.tumor.status,transcriptome,no.permutation=100,methylation.event=c('sil','enh'),seed=100){
    
  ###################################################################################################################
  #Define the datasets for the calculation: transcriptome//methylome//meth.status
    Samplesize=ncol(meth.tumor)
    EXP_T=transcriptome
    METH_S=meth.tumor.status
    METH_T=meth.tumor
  ## make 4 transcriptome datasets pooling the transcriptome values of Bin1-4
    B1=lapply(1:nrow(METH_T), function(x){
      Meth_status=METH_T[x,]<=0.25
      data=EXP_T[x,which(Meth_status)]
      return(data)
    })
    B1.exp.val=na.omit(unlist(B1))
  
  
  
    B2=lapply(1:nrow(METH_T), function(x){
      Meth_status=(METH_T[x,]>0.25)&(METH_T[x,]<=0.5)
      data=EXP_T[x,which(Meth_status)]
      return(data)
    })
    B2.exp.val=na.omit(unlist(B2))
  
  
    B3=lapply(1:nrow(METH_T), function(x){
      Meth_status=(METH_T[x,]>0.5)&(METH_T[x,]<=0.75)
      data=EXP_T[x,which(Meth_status)]
      return(data)
    })
    B3.exp.val=na.omit(unlist(B3))
  
  
    B4=lapply(1:nrow(METH_T), function(x){
      Meth_status=(METH_T[x,]>0.75)
      data=EXP_T[x,which(Meth_status)]
      return(data)
    })
    B4.exp.val=na.omit(unlist(B4))

  ##########################################################################################
  #Permutation
    set.seed(seed)
    prob.no=nrow(METH_T)
    y.rep=rep(NA,times=ncol(METH_T))
    
    exp.score=lapply(1:no.permutation, function(U){
      print(paste('Permutation --> ',U,sep = ''))
      scores=sapply(1:nrow(METH_S), function(x){
        #Generating the permutated transcriptome data
          y=y.rep
          first=which(METH_T[x,]<=0.25) 
          y.1=sample(B1.exp.val,length(first),replace = TRUE)
          second=which((METH_T[x,]>0.25)&(METH_T[x,]<=0.5))
          y.2=sample(B2.exp.val,length(second),replace = TRUE)
          third=which((METH_T[x,]>0.5)&(METH_T[x,]<=0.75))
          y.3=sample(B3.exp.val,length(third),replace = TRUE)
          forth=which(METH_T[x,]>0.75)
          y.4=sample(B4.exp.val,length(forth),replace = TRUE)
          
          y[first]=y.1
          y[second]=y.2
          y[third]=y.3
          y[forth]=y.4
          
        # constructing the z= methylome data and meth_status
          Meth_status=(METH_S[x,]==1)
          z=as.numeric(METH_T[x,])
          
          
          
          
        #exclude NAs
          omit.y=is.na(EXP_T[x,])
          omit.z=is.na(z)
          omit=omit.y+omit.z
          omit=omit>0
          
          y=y[which(!omit)]
          z=z[which(!omit)]
          Meth_status=Meth_status[which(!omit)]
          
        ##### continue if and only if still any value left by NA omition
        if(length(y)<5){
            return(NA)
        }else{
          #Normalize the data  
          miny=min(y)
          maxy=max(y)
          if((maxy-miny)!=0){
            y=(y-miny)/(maxy-miny)}
          
          minz=min(z)
          maxz=max(z)
          if((maxz-minz)!=0){
            z.cor=(z-minz)/(maxz-minz)}else{
              z.cor=z
            }
          
          m=as.matrix(data.frame(z,z.cor,y))
          ##### continue if and only if methylation-different samples are more than 5
          if(sum(Meth_status,na.rm = TRUE)<5|sum(!Meth_status,na.rm = TRUE)<5){Score=NA}else{
            #C1: average of no event samples // average of event samples
              c1=c(mean(z[which(!Meth_status)],na.rm = TRUE),mean(y[which(!Meth_status)],na.rm = TRUE))
              c2=c(mean(z[which(Meth_status)],na.rm = TRUE),mean(y[which(Meth_status)],na.rm = TRUE))
            
            # ranking the samples by their expression values
              rank.y=length(y)+1-rank(y)
            # info samples in bin.1
              B1=z<=0.25
              rank.y.B1=rank.y[which(B1)]
              B1.No=length(rank.y.B1)
              B1.BetaVal=z[which(B1)]
              OutL.B1=rank.y.B1>B1.No
              
            
            # info samples in bin.2
              B2=(z>0.25)&(z<=0.5)
              rank.y.B2=rank.y[which(B2)]
              B2.No=length(rank.y.B2)
              B2.BetaVal=z[which(B2)]
              OutL.B2=(rank.y.B2>(B2.No+B1.No))|(rank.y.B2<=B1.No)
              
            
            # info samples in bin.3
              B3=(z>0.5)&(z<=0.75)
              rank.y.B3=rank.y[which(B3)]
              B3.No=length(rank.y.B3)
              B3.BetaVal=z[which(B3)]
              OutL.B3=(rank.y.B3>(B3.No+B2.No+B1.No))|(rank.y.B3<=(B2.No+B1.No))
            
            # info samples in bin.4 
              B4=(z>0.75)
              rank.y.B4=rank.y[which(B4)]
              B4.No=length(rank.y.B4)
              B4.BetaVal=z[which(B4)]
              OutL.B4=rank.y.B4<=(B3.No+B2.No+B1.No)
            
              
            ## Score calculations components based on the the bin info 
              if(length(which(OutL.B1))==0){W.out.B1=0
            } else {
                W.out.B1=sapply(which(OutL.B1), function(q) {
                  a=Bin.No(x=rank.y.B1[q],a1=B1.No,a2=(B2.No+B1.No),a3=(B1.No+B2.No+B3.No),a4=(B1.No+B2.No+B3.No+B4.No))
                  w=abs(1-a)/4
                  return(w)
                })
            }
              if(length(which(OutL.B2))==0){W.out.B2=0}else{
                W.out.B2=sapply(which(OutL.B2), function(q){
                  a=Bin.No(x=rank.y.B2[q],a1=B1.No,a2=(B2.No+B1.No),a3=(B1.No+B2.No+B3.No),a4=(B1.No+B2.No+B3.No+B4.No))
                  w=abs(2-a)/4
                  return(w)
                })
              }
              if(length(which(OutL.B3))==0){W.out.B3=0}else{
                W.out.B3=sapply(which(OutL.B3), function(q){
                  a=Bin.No(x=rank.y.B3[q],a1=B1.No,a2=(B2.No+B1.No),a3=(B1.No+B2.No+B3.No),a4=(B1.No+B2.No+B3.No+B4.No))
                  w=abs(3-a)/4
                  return(w)
                })
              }
              if(length(which(OutL.B4))==0){W.out.B4=0}else{
                W.out.B4=sapply(which(OutL.B4), function(q){
                  a=Bin.No(x=rank.y.B4[q],a1=B1.No,a2=(B2.No+B1.No),a3=(B1.No+B2.No+B3.No),a4=(B1.No+B2.No+B3.No+B4.No))
                  w=abs(4-a)/4
                  return(w)
                })
              }
              
              
              if(length(which(B1))==0){score.B1=0}else{score.B1=(1+(length(which(!OutL.B1))/(B1.No+4))-(sum(W.out.B1)/B1.No))}
              if(length(which(B2))==0){score.B2=0}else{score.B2=(1+(length(which(!OutL.B2))/(B2.No+4))-(sum(W.out.B2)/B2.No))}
              if(length(which(B3))==0){score.B3=0}else{score.B3=(1+(length(which(!OutL.B3))/(B3.No+4))-(sum(W.out.B3)/B3.No))}
              if(length(which(B4))==0){score.B4=0}else{score.B4=(1+(length(which(!OutL.B4))/(B4.No+4))-(sum(W.out.B4)/B4.No))}
              
              
              
              Score=sqrt(abs(c1[1]-c2[1])*abs(c1[2]-c2[2]))*(score.B1+score.B2+score.B3+score.B4)
              
            # if the correlation is in the wrong direction returns NA
              if((c1[2]<c2[2]&methylation.event=='sil')|(c1[2]>c2[2]&methylation.event=='enh')){Score=NA}
          }
          return(Score)
        }
      })
      return(scores)
    })
  ###########################################################################################  
  #FDR calculation
    exp.score.all=do.call(cbind,exp.score)
    exp.score.mean=rowMeans(exp.score.all,na.rm = TRUE)
    obs.score.pool=na.omit(score.values)
      obs.score.pool=sort(obs.score.pool,decreasing = TRUE)
    size.effect.FDR=sapply(1:length(obs.score.pool), function(x){
      delta=obs.score.pool[x]
      delta.exp=(exp.score.all>=delta)
        delta.exp=colSums(delta.exp,na.rm = TRUE)
      no.false=mean(delta.exp)
      no.obs=sum(obs.score.pool>=delta,na.rm = TRUE)
      y=no.false/no.obs
      return(y)
    })
    
    FDR.res=data.frame(observed.scores=obs.score.pool,FDR=size.effect.FDR)
    permutation.res=data.frame(exp.score.all,exp.score.mean,observed.scores=score.values)
      rownames(permutation.res)=rownames(meth.tumor)
  
  # FDR cutoff
      no.false=1:nrow(FDR.res)*FDR.res$FDR
      cutoff=which((FDR.res$observed.scores>1.5)&(FDR.res$FDR<=0.1)&(no.false<0.5))
      if(length(cutoff)>0){
        score.cutoff=FDR.res$observed.scores[max(cutoff)]}else{
          score.cutoff=NA}
  ##Results
    final.res=list(FDR.res=FDR.res,permutation.res=permutation.res,score.cutoff=score.cutoff)
    return(final.res)
}




