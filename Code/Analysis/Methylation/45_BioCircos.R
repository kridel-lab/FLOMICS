# Updated 3 Aug 2021
# 19 March 2021
# Function: Generating circular multi-track plots
# Author: Anjali Silva

# https://cran.r-project.org/web/packages/BioCircos/vignettes/BioCircos.html
# Input:
# PNGorPDF: Output format of the image, options = "png" or "pdf". Default: png.
# ProduceImages: Produce images or not, options = "Yes" or "No". Default: Yes

# Output:
# FLResultsDataFrame: A data frame containing patients x categories, categories containing "epiCMIT", "Type",
#                    "Stage", "Sex", "BCL2Translocation", "Clusters", "EZH2Mut", "BCL2Mut", "KMT2DMut", "CREBBPMut",
#                    "EP300Mut", "epiCMIThyper", "epiCMIThypo", "indicator"


library(BioCircos)

### Trying example ####

BioCircos(yChr = FALSE)
BioCircos(yChr = FALSE, xChr = FALSE, chrPad = 0, genomeFillColor = "Blues")



# Background track
tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.5, maxRadius = 0.8,
                                     borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#FFBBBB")  
BioCircos(tracklist, 
          genomeFillColor = "PuOr",
          chrPad = 0.05, 
          displayGenomeBorder = FALSE, 
          genomeTicksDisplay = FALSE,  
          genomeLabelTextSize = "9pt", 
          genomeLabelDy = 0)


# Heatmap tracks
# Define a custom genome
genomeChr = LETTERS[1:10]
lengthChr = 5*1:length(genomeChr)
names(lengthChr) <- genomeChr

# Define boxes positions
boxPositions = unlist(sapply(lengthChr, seq))
boxChromosomes = rep(genomeChr, lengthChr)

# Define values for two heatmap tracks
boxVal1 = boxPositions %% 13 / 13
boxVal2 = (7 + boxPositions) %% 17 / 17

tracks = BioCircosHeatmapTrack("heatmap1", boxChromosomes, boxPositions - 1, boxPositions,
                               boxVal1, minRadius = 0.6, maxRadius = 0.75)
tracks = tracks + BioCircosHeatmapTrack("heatmap1", boxChromosomes, boxPositions - 1, 
                                        boxPositions, boxVal2, minRadius = 0.75, maxRadius = 0.9, color = c("#FFAAAA", "#000000"))

BioCircos(tracks, genome = as.list(lengthChr), genomeTicksDisplay = F, genomeLabelDy = 0, 
          HEATMAPMouseOverColor = "#F3C73A")

# [END]
