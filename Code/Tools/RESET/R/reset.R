###################################################################################################################################################################

# This function is the sequence of methStatus, eventScore and FDRcal

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

# OUTPUT:
  # 1. Normal samples methylome data ("methStatus" output)
  # 2. Tumor samples methylome data ; all samples ("methStatus" output)
  # 3. Tumor samples methylation status; all samples ("methStatus" output)
  # 2. Specific Beta distribution (based on each specific probes value) ("methStatus" output)
  # 3. Universal Beta estimation of normal-samples distribution (output of 'fitdistr' function) ("methStatus" output)
  # 4. Score.report : the score reported on each event (probe methylation event) ('eventScore' output)
  # 5. transcriptome : the transcriptome data used in the analysis (all the data are match with methylome and methylome staus) ('eventScore' output)
  # 6. transcriptome.nor.dis : the ranked normally distributed transcriptome data ('eventScore' output)
  # 7. meth.tumor : methylome data used in score calculation ; samples matched with transcriptome data ('eventScore' output)
  # 8. meth.tumor.status : methylation status data used in score calculation; samples matched with transcriptome data ('eventScore' output)
  # 9. FDR.res : table containing the FDR calculated for each observed score ("FDRcal" output)
  # 10. permutation.res : the results of (n=100) permutations + mean permutated score + observed scores ("FDRcal" output)
  # 11. score.cutoff : the suggested score cutoff: score > 1.5 & FDR < 0.1 & false score < 0.5 ("FDRcal" output)

###################################################################################################################################################################
# source('shared_functions.R')
# source('methStatus.R')
# source('eventScore.R')
# source('FDRcal.R')

#' RESET
#'
#' This function runs the comprehensive RESET analysis on the provided datasets
#'
#' @import MASS
#'
#' @param normal.db normal methylation
#' @param meth.tumor tumor methylation
#' @param transcriptome gene expression data
#' @param methylation.event =c('sil','enh')
#' @param FDR.permutation.no =100
#'
#'
#' @return A list contains all the info regarding the analysis
#'
#' #@examples
#' #pipelineCore=function(normal.db,meth.tumor,transcriptome,methylation.event=c('sil','enh'),FDR.permutation.no=100)
#'
#' @export

reset=function(normal.db, meth.tumor, transcriptome,
  methylation.event=c('sil','enh'),
  FDR.permutation.no=100, seed=100){

  ###################################################################################################################
  ###Functions defintion


  ###################################################################################################################
  # Methylation status
    print('Tumor samples Methylation status determination')
    meth.stat=methStatus(normal.db,meth.tumor,methylation.event)
      # Output:
        # 1. Normal samples methylome data (FINAL)
        # 2. Tumor samples methylome data
        # 3. Tumor samples methylation status
        # 4. Specific Beta distribution (based on each specific probes value) (FINAL)
        # 5. Universal Beta estimation of normal-samples distribution (output of 'fitdistr' function) (FINAL)

  # Score calculation
    print('Evaluation of the events Scores')
    score=eventScore(meth.tumor=meth.stat$tumor.meth,meth.tumor.status=meth.stat$tumor.meth.status,transcriptome,methylation.event)
      # Output:
        # 1. Score.report : the score reported on each event (probe methylation event) (FINAL)
        # 2. transcriptome : the transcriptome data used in the analysis (all the data are match with methylome and methylome staus) (FINAL)
        # 3. transcriptome.nor.dis : the ranked normally distributed transcriptome data (FINAL)
        # 4. meth.tumor : methylome data used in score calculation (FINAL)
        # 5. meth.tumor.status : methylation status data used in score calculation (FINAL)

  # FDR Calculation
    print('FDR estimation')
    fdr=FDRcal(score.values=score$Score.report[,4],meth.tumor=score$meth.tumor,meth.tumor.status=score$meth.tumor.status,transcriptome=score$transcriptome.norm.dis,no.permutation=FDR.permutation.no,methylation.event, seed=seed)
      # Output
        # 1. FDR.res : table containing the FDR calculated for each observed score
        # 2. permutation.res : the results of (n=100) permutations + mean permutated score + observed scores
        # 3. score.cutoff : the suggested score cutoff: score > 1.5 & FDR < 0.1 & false score < 0.5


  ##Results
    final.res=list(
      #1 methStatus
        normal.meth=meth.stat$normal.meth,
        meth.tumor.all=meth.stat$tumor.meth,
        meth.tumor.status.all=meth.stat$tumor.meth.status,
        beta.dist.normal.specific=meth.stat$beta.dist.normal.specific,
        beta.dist.normal.universal=meth.stat$beta.dist.normal.universal,
      #2 eventScore
        Score.report=score$Score.report,
        transcriptome=score$transcriptome,
        transcriptome.norm.dis=score$transcriptome.norm.dis,
        meth.tumor.matched=score$meth.tumor,
        meth.tumor.status.matched=score$meth.tumor.status,
      #3 FDRcal
        FDR.res=fdr$FDR.res,
        permutation.res=fdr$permutation.res,
        score.cutoff=fdr$score.cutoff
                   )
    print('DONE')
    return(final.res)
}
