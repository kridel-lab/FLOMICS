# Updated 3 Aug 2021
# 26 Feb 2020
# Function: Return probes that has a Median Absolute Deviation value above a given threshold across samples, 
#           using Beta values OR return X probes with the highest MAD value.
#           (Similar work has been done in 
#           https://link-springer-com.myaccess.library.utoronto.ca/article/10.1007%2Fs00432-019-02895-2 
#           and https://clinicalsarcomaresearch.biomedcentral.com/articles/10.1186/s13569-019-0113-6 )
# Author: Anjali Silva

# Input:
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns.
# MvalueMatrix: A matrix of probes x patients, with probes in rows and patients as columns; if not proivded
#               set to NA (default)
# CutOff: A numeric indicating the cutoff for MAD (median absolute deviation).  Probes with MAD greater than
#         this cutoff will be selected. 
# NumberofProbes: An integer specifying the number of probes to be returned with the highest MAD value. Example,
#                NumberofProbes = 1000. 

# Output:
# BetaMatrix_MAD_Filtered: Matrix of beta values for probes over the specified median absolute deviation x patients. 
# MvalueMatrix_MAD_Filtered Matrix of M values for probes over the specified median absolute deviation x patients. 
# All_MAD_Values: A vector specifying MAD value for each probe in the Beta Matrix.

MedianAbsoluteDeviation34 <- function(BetaMatrix, 
                                      MvalueMatrix = NA, 
                                      CutOff = NA, 
                                      NumberofProbes = NA) {
    
  # Checking input values
  if(! is.na(CutOff) && ! is.na(NumberofProbes)) {
    stop("\n Cannot provide values for both CutOff and NumberofProbes. Leave one argument as NA.")
  }
  
  # Calculate MAD for each probe
  set.seed(1)
  MADEachProbe <- apply(BetaMatrix, 1, mad)
  
    # If using cut off
    if (! is.na(CutOff)) {
      BetaMatrix_MAD_Filtered <- BetaMatrix[which(MADEachProbe > CutOff), ]
      if(! is.na(MvalueMatrix)[1]) {
        MvalueMatrix_MAD_Filtered <- MvalueMatrix[which(MADEachProbe > CutOff), ]
      } else {
        MvalueMatrix_MAD_Filtered <- MvalueMatrix
      }
    }
  
    # If using number of probes
    if (! is.na(NumberofProbes)) {
      BetaMatrix_MAD_Filtered <- BetaMatrix[match(names(sort(MADEachProbe, 
                                 decreasing = TRUE)[1:NumberofProbes]), rownames(BetaMatrix)), ]
      if(! is.na(MvalueMatrix)[1]) {
        MvalueMatrix_MAD_Filtered <- MvalueMatrix[match(rownames(BetaMatrix_MAD_Filtered), 
                                                        rownames(MvalueMatrix)), ]
      } else {
        MvalueMatrix_MAD_Filtered <- MvalueMatrix
      }
    }

  
  RESULTS <- list(BetaMatrix_MAD_Filtered = BetaMatrix_MAD_Filtered,
                  MvalueMatrix_MAD_Filtered = MvalueMatrix_MAD_Filtered,
                  All_MAD_Values = MADEachProbe)
  
  class(RESULTS) <- "MedianAbsoluteDeviation_ASilva"
  return(RESULTS)
}
# [END]
