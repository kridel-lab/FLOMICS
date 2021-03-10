
###
# This script is based on original script from Arash (Pugh lab).
# It reads in output from MiXCR
# and then runs iNext to calculate diversity indices.
# https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html
###

setwd("~/github/FLOMICS/")

date <- Sys.Date()

library(iNEXT)
library(ggplot2)

#' List clone counts for diversity measures
#' @description This funtion removes non-productive aaCDR3s and lists read counts (column "cloneCount" on mixcr output)
#' for dowstream diversity calculations. Note that this function does not include files with only one clone.
#' This is just to remove the shallow samples which happens quite often if using RNAseq data. This should not cause any issues
#' when using the function with captured data.
#'
#' @param datapath path to mixcr files
#' @param chain any of: TRA, TRB, TRD, TRG, IGH, IGL, IGK
#'
#' @examples immunelistfx("~/git/PLTK/PLTK/data-raw/", "TRB")

immunelistfx <- function(datapath, chain){
# File pattern is CLONES_TRB*.txt. Ensure you output each chain seperately.
  # file_list<- list.files(datapath, pattern = paste0("CLONES_", chain))
  file_list<- list.files(datapath, pattern = "TRB")
  readlist = list()
  i <- 1
  for(f in file_list){
    mixcrfle <- read.table(paste0(datapath, f),header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, na.strings = c("", "NA"))
    message("number of clones in file:", f)
    print(nrow(mixcrfle))
    # f <- substr(f, 11, nchar(f)-4)
    f <- paste0(gsub("\\..*", "", f), ".", chain)
    print(f)
    if(nrow(mixcrfle) <= 1){next()}
    message("number of non prodcutive CDR3s:")
    print(length(mixcrfle$aaSeqCDR3[grepl("[*]", mixcrfle$aaSeqCDR3) | grepl("_", mixcrfle$aaSeqCDR3)]))
    mixcrfle <- mixcrfle[!grepl("_", mixcrfle$aaSeqCDR3) & !grepl("[*]", mixcrfle$aaSeqCDR3),]
    message("nonproductive aaCDR3 removed")
    readlist[[i]] <- mixcrfle$cloneCount
    names(readlist)[i] <- f
    i <- i + 1
  }
  return(readlist)
}

lst <- immunelistfx("Analysis-Files/Mixcr/specific_clonotypes/", "TRB")

#' Matrix for diversity measures and related descriptives
#' @description This function takes the list generated from immunelistfx and generates a dataframe for descriptives and diversity estimators
#' @param lst list output from immunelistfx
#' @param chain any of: TRA, TRB, TRD, TRG, IGH, IGL, IGK
#' @param batchname output csv filename will be in this form: "divstats_chain_batchname.csv"
#' @param outpath path to save output file
#'
#' @examples data(celllines_list)
#' Divstats.fx(celllines_list, "TRB", "thisismytest", "~/")
#'
Divstats.fx <- function(lst, chain, batchname, outpath){

  div_stats <- matrix(ncol = 19, nrow = length(lst))
  
  colnames(div_stats) <- c(chain,"Reads", "CPKR","Average_reads", "VMR","Max_reads","Singletons", "Doubletons",
                           "qD", "Sample_Coverage",
                           "observed_Richness", "estimated_Richness", "SE_Richness",
                           "observed_Shannon", "estimated_Shannon", "SE_Shannon",
                           "observed_Simpson", "estimated_Simpson", "SE_Simpson")
  rownames(div_stats) <- names(lst)
  
  #Descriptives
  div_stats[,chain] <- unlist(lapply(lst, length))
  div_stats[,"Reads"] <- unlist(lapply(lst, sum))
  div_stats[,"CPKR"] <- unlist(lapply(lst, function(x){(length(x)/sum(x))*1000}))
  div_stats[,"Average_reads"] <- unlist(lapply(lst, mean))
  div_stats[,"VMR"] <- unlist(lapply(lst, function(x){sd(x)/mean(x)}))
  div_stats[,"Max_reads"] <- unlist(lapply(lst, max))
  div_stats[,"Singletons"] <- unlist(lapply(lst, function(x){length(which(x==1))}))
  div_stats[,"Doubletons"] <- unlist(lapply(lst, function(x){length(which(x==2))}))
  
  # Estimators
  out <- iNEXT(lst, 0, datatype="abundance")
  # Estimations based on Extrapolation
  est <- out$iNextEst
  qDlist <- lapply(est, "[[", "qD")
  div_stats[, "qD"] <- unlist(lapply(qDlist, max))
  
  SClist <- lapply(est, "[[", "SC")
  div_stats[, "Sample_Coverage"] <- unlist(lapply(SClist, max))
  
  # Asymptotic estimators
  
  AsyEst <- out$AsyEst
  div_stats[, c("observed_Richness", "estimated_Richness", "SE_Richness")] <- as.matrix(AsyEst[AsyEst$Diversity == "Species richness",
                                                                                                c("Observed", "Estimator", "s.e.")])
  div_stats[, c("observed_Shannon", "estimated_Shannon", "SE_Shannon")] <- as.matrix(AsyEst[AsyEst$Diversity == "Shannon diversity",
                                                                                            c("Observed", "Estimator", "s.e.")])
  div_stats[, c("observed_Simpson", "estimated_Simpson", "SE_Simpson")] <- as.matrix(AsyEst[AsyEst$Diversity == "Simpson diversity",
                                                                                            c("Observed", "Estimator", "s.e.")])
  write.csv(div_stats,
            file = paste0(outpath, "divstats_", chain, batchname, ".csv"),
            row.names = TRUE)
}

Divstats.fx(lst, "TRB", "", "Analysis-Files/Mixcr/Arash_script_output/")

out <- iNEXT(lst, 0, datatype="abundance")

ggiNEXT(out, type=2, facet.var="none", color.var="site") +
  ylim(0,1) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(paste0("img/", date, " iNext_Sample_coverage.pdf"), width = 8, height = 8)
