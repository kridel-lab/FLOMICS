# Updated 26 April 2019
# Function: Obtain probes that has a standard deviation (SD) value across samples, above a given cutoff. 
# Values calculated using Beta values! (Similar work has been done by Hinoue, T et al., 2012: 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3266034/)
# Author: Anjali Silva

# Input:
# BetaMatrix: A matrix of beta values for probes x patients, with probes in rows and patients as columns.
# MvalueMatrix: A matrix of probes x patients, with probes in rows and patients as columns; if not proivded
#               set to NA (default)
# CutOff: A numeric value for cutoff. E.g., 0.2. Probes with a SD value greater this cutoff will be selected. 
#         Default NA. 
# NumberofProbes: An integer specifying the number of probes to be returned with the highest SD value. Example,
#                NumberofProbes = 1000. 


# Output:
# BetaMatrix_SD_Filtered: Matrix of beta values for probes over the specified standard deviation x patients. 
# MvalueMatrix_SD_Filtered: Matrix of M values for probes over the specified standard deviation x patients, 
#                           if provided as an argument MvalueMatrix. 
# All_SD_Values: A vector specifying standard deviation value for each probe in the Beta Matrix.

SDeviation <- function(BetaMatrix, 
                       MvalueMatrix = NA, 
                       CutOff = NA, 
                       NumberofProbes = NA) {
  
  # Checking input values
  if(! is.na(CutOff) && ! is.na(NumberofProbes)) {
    stop("\n Cannot provide values for both CutOff and NumberofProbes. Leave one argument as NA.")
  }
  
  # Calculate SD for each probe
  set.seed(1)
  SDEachProbe <- apply(BetaMatrix, 1, sd)
  
  # If using cut off
  if (! is.na(CutOff)) {
    BetaMatrix_SD_Filtered <- BetaMatrix[which(SDEachProbe > CutOff), ]
    if(! is.na(MvalueMatrix)[1]) {
      MvalueMatrix_SD_Filtered <- MvalueMatrix[which(SDEachProbe > CutOff), ]
    } else {
      MvalueMatrix_SD_Filtered <- MvalueMatrix
    }
  }
  
  # If using number of probes
  if (! is.na(NumberofProbes)) {
    BetaMatrix_SD_Filtered <- BetaMatrix[match(names(sort(SDEachProbe, decreasing = TRUE)[1:NumberofProbes]), 
                                               rownames(BetaMatrix)), ]
    if(! is.na(MvalueMatrix)[1]) {
      MvalueMatrix_SD_Filtered <- MvalueMatrix[match(rownames(BetaMatrix_SD_Filtered), rownames(MvalueMatrix)), ]
    } else {
      MvalueMatrix_SD_Filtered <- MvalueMatrix
    }
  }
  
  RESULTS <- list(BetaMatrix_SD_Filtered = BetaMatrix_SD_Filtered,
                  MvalueMatrix_SD_Filtered = MvalueMatrix_SD_Filtered,
                  All_SD_Values = SDEachProbe)
  
  class(RESULTS) <- "SDeviation_ASilva"
  return(RESULTS)
}
