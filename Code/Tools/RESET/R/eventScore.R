###################################################################################################################################################################

# This function calculates the silencing/enhancing methylation score; which evaluates the level of correlation between hyper/hypomethylation and gene expression

# INPUTS:
  # 1. Processed tumor methylome data; beta values (the probe should have either sil or enh probes condition as specified in 'methNorSel' function)
    ## rows: probe index
    ## columns: normal samples
  # 2. Processed tumor methylation status (this dataset should have the coresponding methylation status of the tumor samples)
    ## rows: probe index
    ## columns: normal samples
  # 3. Transcriptome dataset of the tumor samples
    ## columns: tumor samples
    ## Gene name (as specified in the probe index)
# OUTPUT: 
    # 1. Score.report : the score reported on each event (probe methylation event)
    # 2. transcriptome : the transcriptome data used in the analysis (all the data are match with methylome and methylome staus)
    # 3. transcriptome.nor.dis : the ranked normally distributed transcriptome data
    # 4. meth.tumor : methylome data used in score calculation
    # 5. meth.tumor.status : methylation status data used in score calculation
###################################################################################################################################################################

# source('shared_functions.R')

#' Score Calculation
#'
#' This function calculates the silencing/enhancing methylation score; which evaluates the level of correlation between hyper/hypomethylation and gene expression
#' @param meth.tumor Processed tumor methylome data; beta values (the probe should have either sil or enh probes condition as specified in 'methNorSel' function. rows: probe index, columns: normal samples)
#' @param meth.tumor.status Processed tumor methylation status (this dataset should have the coresponding methylation status of the tumor samples. rows: probe index, columns: normal samples)
#' @param transcriptome Transcriptome dataset of the tumor samples (columns: tumor samples, Gene name (as specified in the probe index))
#' @param methylation.event methylation event type. Silencing = 'sil',Enhancing = 'enh'
#' @return A list contains: 1. Score.report : the score reported on each event (probe methylation event), 2. transcriptome : the transcriptome data used in the analysis (all the data are match with methylome and methylome staus), 3. transcriptome.nor.dis : the ranked normally distributed transcriptome data, 4. meth.tumor : methylome data used in score calculation, 5. meth.tumor.status : methylation status data used in score calculation
#' @export
eventScore=function(meth.tumor,meth.tumor.status,transcriptome,methylation.event=c('sil','enh')){
  # Transcriptome and methylome sample selection
    #Samples sorting
      s=sortection(list(colnames(meth.tumor),colnames(transcriptome)))
      # print(s[1:20,])
      meth.tumor=meth.tumor[,s[,1]]
      meth.tumor.status=meth.tumor.status[,s[,1]]
      transcriptome=transcriptome[,s[,2]]
    #transcriptome rows sorting
      #gene name methylome
        gene.meth=toupper(sapply(rownames(meth.tumor), function(x){
          return(strsplit(x,'@')[[1]][3])
        }))
      #sorting and excluding duplicate rows of transcriptome
        overal.exp=rowSums(transcriptome,na.rm = TRUE)
        transcriptome=transcriptome[order(overal.exp,decreasing = TRUE),]
        row.transc=toupper(rownames(transcriptome))
        no.exc=which(duplicated(row.transc))
          if(length(no.exc)>0){
            transcriptome=transcriptome[-no.exc,]
            row.transc=row.transc[-no.exc]
          }
        rownames(transcriptome)=row.transc
      # generate final transcriptome dataset each row containing the expression level of coresponding methylome dataset's gene
          tr=as.data.frame(transcriptome)
        transcriptome.final=as.matrix(tr[gene.meth,])
          rownames(transcriptome.final)=rownames(meth.tumor)
          
        ### TRANSCRIPTOME NEW dataset
          transcriptome.new=t(apply(transcriptome.final, 1, function(x){
            if(all(is.na(x))){return(x)}else{
              r = rank(x)
              qq = qqnorm(r,plot.it=F)
              qq = qq$x
              qq[is.na(x)] = NA # useless, will leave here just to keep the original code version
              return(qq)
            }
          }))
          colnames(transcriptome.new)=colnames(transcriptome.final)
  #Score calculation
    Score=sapply(1:nrow(meth.tumor.status), function(x){
      # probe data extraction : meth_status, y= expression, z= methylome
        Meth_status=meth.tumor.status[x,]==1
        y=as.numeric(transcriptome.new[x,])
        z=as.numeric(meth.tumor[x,])
      #exclude NAs
        omit.y=is.na(y)
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
            if(length(which(OutL.B1))==0){W.out.B1=0}else{
              W.out.B1=sapply(which(OutL.B1), function(q){
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
  ##Results
    Gene=toupper(sapply(rownames(meth.tumor), function(x){
      return(strsplit(x,'@')[[1]][3])
    }))
    Promoter.no=toupper(sapply(rownames(meth.tumor), function(x){
      return(strsplit(x,'@')[[1]][2])
    }))
    ProbID=toupper(sapply(rownames(meth.tumor), function(x){
      return(strsplit(x,'@')[[1]][1])
    }))
    No.Methylation.Events=rowSums(meth.tumor.status,na.rm = TRUE)
    Score.report=data.frame(ProbID,Gene,Promoter.no,Score,No.Methylation.Events)
  
  final.result=list(Score.report=Score.report,transcriptome=transcriptome.final,transcriptome.norm.dis=transcriptome.new,meth.tumor=meth.tumor,meth.tumor.status=meth.tumor.status)
  return(final.result)
}


