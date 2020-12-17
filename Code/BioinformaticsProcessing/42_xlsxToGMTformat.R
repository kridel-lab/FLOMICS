# Date: 5 October 2020
# Function: A function that read in Staudt_SignatureDB_version.xlsx and 
#           produce the corresponding .gmt file. The GMT file format is 
#           a tab delimited file format used for gene sets. An example is
#           provided below. 
# Author: Anjali Silva 


# Input:
# StaudtSignatureDB: Data frame or table resulting from reading Staudt_SignatureDB_version.xlsx
#                     via readxl::read_excel(). Must contain fields "Short.signature.name"
#                     and "Gene symbol concensus" as column titles.  
# writeFile: Logical TRUE or FALSE. If TRUE, then a file in .gmt format will be saved in current 
#            working directory.
# nameFile: A unique character string indicating the name to be used for file, if writeFile == TRUE. 

# Output:
# DataFrameSignatureGenes: A dataframe



xlsxToGMT <- function(StaudtSignatureDB, 
                      writeFile = TRUE, 
                      nameFile = NA) {
  
  # Create a list to save entries
  listSignatureGenes <- vector(mode = "list", 
         length = length(unique(StaudtSignatureDB$`____Short signature name`)))
  names(listSignatureGenes) <- unique(StaudtSignatureDB$`____Short signature name`)
  
  # For each unique signature name, do the following
  for (i in 1:length(listSignatureGenes)) { 
    
    listSignatureGenes[[i]] <-
      c(NA, base::unname(unlist(StaudtSignatureDB_Excel %>% # NA introduced to provide extra space
        data.frame() %>%
        dplyr::filter(X____Short.signature.name == names(listSignatureGenes)[i]) %>%
        dplyr::select("X____Gene.symbol.concensus"))))
  }
    
  # save into a dataframe
  DataFrameSignatureGenes <- data.frame(lapply(listSignatureGenes, 
                                               "length<-", max(lengths(listSignatureGenes))),
                                        check.names = FALSE)
  # This function was obtained from 
  # https://intellipaat.com/community/10208/the-simplest-way-to-convert-a-list-with-various-length-vectors-to-a-data-frame-in-r
  # October 2020 
  
  # option check.names = FALSE ensure '-' in name are not converted to '.'
  
  
  # transpose data frame
  transpStaudtSignatureDBGMT <- t(StaudtSignatureDBGMT)
  
  if(writeFile == TRUE) {
    if(is.na(nameFile) == FALSE) {
      utils::write.table(testing, 
                         file = paste0("StaudtSig_",nameFile, "_", Sys.Date(), ".gmt"), 
                         sep = "\t",
                         na = "", 
                         quote = FALSE,
                         row.names = TRUE, 
                         col.names = FALSE)
    } else if(is.na(nameFile) == TRUE) {
      utils::write.table(testing, 
                         file = paste0("StaudtSig030920_", Sys.Date(),".gmt"), 
                         sep = "\t",
                         na = "", 
                         quote = FALSE,
                         row.names = TRUE, 
                         col.names = FALSE)
    }
  }
  
  return(transpStaudtSignatureDBGMT)
}





# Examples
# https://lymphochip.nih.gov/signaturedb/ 
# library(readxl)
# RNAseqDirPath <- "~/Desktop/UHN/FLOMICS/RNAseq/ExtendedStudy2020/"
# Staudt_SignatureDB <- readxl::read_excel(path = paste0(RNAseqDirPath, "GSEA/Staudt_SignatureDB_030920.xlsx"))
# head(Staudt_SignatureDB)
# class(Staudt_SignatureDB)

# dataFrameStaudtSigDB <- xlsxToGMT(StaudtSignatureDB = Staudt_SignatureDB, 
#                                  writeFile = TRUE,
#                                  nameFile = "030920")




