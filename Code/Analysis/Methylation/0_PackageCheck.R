# Updated 10 Feb 2020
# Function: check whether a package is available in the system, if not downloads
# Author: Marcelo Ponce and Anjali Silva


RegularPckgs <- c("ade4",
                  "BioCircos",
                  "BiocManager",
                  #"BootValidation",
                  "BiocGenerics",
                  "BisqueRNA",
                  "caret",
                  "cluster",
                  "corrplot",
                  "csSAM",
                  "cowplot",
                  "data.table",
                  "devtools",
                  "dplyr",
                  "enrichR",
                  "EnvStats",
                  "fpc",
                  "factoextra",
                  "ggpubr",
                  "gprofiler2",
                  "gplots",
                  "ggplot2",
                  "ggrepel",
                  "InfiniumPurify",
                  "magrittr",
                  "mclust",
                  "minfi",
                  "MASS",
                  "pastecs",
                  "pheatmap",
                  "plotly",
                  "plotmo",
                  "pROC",
                  "Rtsne",
                  "readr",
                  "reshape2",
                   "rlang",
                   "rjags",
                   "RPMM",
                   "RColorBrewer",
                   "SNFtool",
                   "scales",
                   "Seurat",
                   "stringi",
                   "stats",
                   "survival",
                   "survminer",
                   "tidyverse",
                   "tibble",
                   "data.table",
                   "tidyr")

BioManagerPckgs <- c("AnnotationDbi", 
                     "Biobase",
                     "bumphunter",
                     "DMRcate",
                     "EnhancedVolcano",
                     "FlowSorted.Blood.EPIC",
                     "eulerr",
                     "edgeR",
                     "Gviz",
                     "org.Hs.eg.db",
                     "IlluminaHumanMethylationEPICmanifest",
                     "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
                     "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                     "KEGGprofile",
                     "limma",
                     "missMethyl",
                     "shinyMethyl", 
                     "TOAST")
                     #"WGCNA")

GitHubPckgs <- c("brentp/celltypes450",
                 "rasmusab/bayesian_first_aid")

# A function that loads the packages if not found - not used 
LoadCheckPkg <- function(RegularPckgs = NA, BioManagerPckgs = NA, GitHubPckgs = NA) {
  
  # Code developed by A.Silva based on code provided by Dr. Marcelo Ponce, April 2019
  if (! is.na(RegularPckgs)[1]) {
    for (pckg in RegularPckgs) {
      cat("Checking package", pckg, "\n")
      # check whether the package is NOT loaded
      if (! paste('package:', pckg, sep = "") %in% search()) {
        # check whether the package is available in the system
        if (pckg %in%  .packages(all.available = TRUE)) {
          # load the package
          cat("Loading library", pckg,"... \n")
          library(pckg, character.only = TRUE)
        } else {
          cat("Installing package", pckg,"... \n")
          install.packages(pckg)
          library(pckg, character.only = TRUE)
        }
      }
    }
  }
  
  if (! is.na(BioManagerPckgs)[1]) {
    for (pckg in BioManagerPckgs) {
      cat("Checking package", pckg, "\n")
      # check whether the package is NOT loaded
      if (! paste('package:', pckg, sep = "") %in% search()) {
        # check whether the package is available in the system
        if (pckg %in%  .packages(all.available = TRUE)) {
          # load the package
          cat("Loading library", pckg,"... \n")
          library(pckg, character.only=TRUE)
        } else {
          if(getRversion() >= 3.5) {
            # If R version is equal to greater than 3.5 
            install.packages("BiocManager")
            library("BiocManager")
            BiocManager::install(pckg)
          } else {
            # If R version is less than 3.5
            source("https://bioconductor.org/biocLite.R")
            cat("Installing package",pckg,"... \n")
            biocLite(pckg)
            library(pckg, character.only = TRUE)
          }
        }
      }
    }
  }
  
  if (! is.na(GitHubPckgs)[1]) {
    for (pckg in GitHubPckgs) {
      pckg1 <- strsplit(pckg, split = "/" )[[1]][2]
      cat("Checking package", pckg1, "\n")
      # check whether the package is NOT loaded
      if (! paste('package:', pckg1, sep = "") %in% search()) {
        # check whether the package is available in the system
        if (pckg1 %in%  .packages(all.available = TRUE)) {
          # load the package
          cat("Loading library", pckg,"... \n")
          library(pckg1, character.only = TRUE)
        } else {
            install.packages("devtools")
            library("devtools")
            devtools::install_github(pckg)
        }
      }
    }
  }
  
  # Developed by Anjali Silva and Marcelo Ponce
}

CheckPackageOnly0 <- function(RegularPckgs = NA, BioManagerPckgs = NA) {
  fail = FALSE
  # Code developed by A.Silva based on code provided by Dr. Marcelo Ponce, April 2019
  if (! is.na((RegularPckgs)[1])) {
    for (pckg in RegularPckgs) {
      cat("Checking package", pckg, "\n")
      # check whether the package is NOT loaded
      if (! paste('package:', pckg, sep = "") %in% search()) {
        # check whether the package is available in the system
        if (pckg %in%  .packages(all.available = TRUE)) {
          # load the package
          cat("Loading library", pckg, "... \n")
          library(pckg, verbose = FALSE, character.only = TRUE)
        } else {
          cat("**************** Package: ", pckg, "is not found. Install using",
              paste("install.packages('", pckg, "') **************** \n", sep = "")) 
          fail = TRUE
         }
      }
    }
  }
  
  if (! is.na((BioManagerPckgs)[1])) {
    for (pckg in BioManagerPckgs) {
      cat("Checking package", pckg, "\n")
      # check whether the package is NOT loaded
      if (! paste('package:', pckg, sep = "") %in% search()) {
        # check whether the package is available in the system
        if (pckg %in%  .packages(all.available = TRUE)) {
          # load the package
          cat("Loading library", pckg, "... \n")
          library(pckg, character.only = TRUE, verbose = FALSE)
        } else {
          cat("**************** Package: ", pckg, "is not found. Install using", 
              paste("BiocManager::install('", pckg, "') ****************\n", sep = "")) 
          fail =TRUE
        }
      }
    }
  }
  
  if (! is.na((GitHubPckgs)[1])) {
    for (pckg in GitHubPckgs) {
      pckg1 <- strsplit(pckg, split = "/" )[[1]][2]
      cat("Checking package", pckg1, "\n")
      # check whether the package is NOT loaded
      if (! paste('package:', pckg1, sep = "") %in% search()) {
        # check whether the package is available in the system
        if (pckg1 %in%  .packages(all.available = TRUE)) {
          # load the package
          cat("Loading library", pckg1, "... \n")
          library(pckg1, character.only = TRUE, verbose = FALSE)
        } else {
          cat("**************** Package: ", pckg1, "is not found. Install using", 
              paste("devtools::install_github('", pckg, "') ****************\n", sep = "")) 
          fail = TRUE
        }
      }
    }
  }
  
  if (fail) stop()
  # Developed by Anjali Silva and Marcelo Ponce
}

# Run the function CheckPackageOnly0- for auto check
CheckPackageOnly0(RegularPckgs = RegularPckgs, 
                 BioManagerPckgs = BioManagerPckgs)

# [END]
