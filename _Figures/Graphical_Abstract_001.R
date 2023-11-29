#---
# making files for input into the 3d contour plot
# which generates the graphical abstract
# Author: Victoria Shelton
#---

#--
# Load in packages
#--

packages <- c("dplyr", "ggplot2", "starsExtra")
lapply(packages, require, character.only = TRUE)

#--
# Setup working space
#--

date <- Sys.Date()
setwd("~/your working directory/")

#--
# Read in Data
#--

sig_file <- read.table("sig_long_out_AIC.txt", sep = ";", header = TRUE)
df_file_long <- read.table("df_long_AIC.txt", sep = ";", header = TRUE)
pat_order <- read.table("pat_order_AIC.txt", sep = ";", header = TRUE) %>%
 unlist()
z <- read.table("SNV_vs_sample_matrix.txt", sep = ";")
gene_order <- c("CREBBP", "STAT6", "TNFRSF14", "TNFAIP3", "TP53", "MYD88",
 "PIM1", "EP300", "HIST1H1E", "BTG2", "STAT3", "IL4R", "NOTCH2", "FAS",
 "SOCS1", "CD79B", "GNA13", "MEF2B", "HVCN1", "BCL7A", "BCL2", "KLHL6",
 "EEF1A1", "EBF1", "P2RY8", "HIST1H1C", "DTX1", "FOXO1", "HIST1H2AM",
 "EZH2", "S1PR2", "SGK1", "ARID1A", "B2M", "MEF2C", "IRF8", "MYC", "BTG1",
 "CD83", "CCND3", "POU2AF1", "RHOA", "HLA-DMB", "NOTCH1", "ATP6AP1",
 "RRAGC", "CTSS", "ATP6V1B2", "KMT2D", "BCL6")

#reorder matrix based on cluster map organization
z <- z[gene_order, pat_order]

#--
# defining the function which will create cluster hills
#--

make_hill = function(k, y_scaled) {
  # Create matrix
  m <- matrix(1, ncol = k, nrow = k)
  # Set focal cell to zero
  m[(nrow(m) + 1) / 2, (ncol(m) + 1) / 2] <- 0
  # Calculate azimuths
  s <- matrix_to_stars(m)
  names(s) <- "value"
  u <- st_as_sf(s, as_points = TRUE)
  ctr <- u[u$value == 0, ]
  #-- first, most outer, radius
  radius <- unname(diff(st_bbox(s)[c(1, 3)]) / 2)
  circle <- st_buffer(ctr, radius + 0.001)
  in_circle <- st_intersects(u, circle, sparse = FALSE)[, 1]
  u$value[in_circle] <- 1
  u$value[!in_circle] <- 0
  #-- all other radii
  for(i in elevparam){
    radius <- unname(i)
    circle <- st_buffer(ctr, radius + 0.001)
    in_circle <- st_intersects(u, circle, sparse = FALSE)[, 1]
    u$value[in_circle] <- y_scaled[(((k / 2) + 1) - i)]
  }
  s <- st_rasterize(u[, "value"], s)
  m <- layer_to_matrix(s)
  m[(nrow(m) + 1) / 2, (ncol(m) + 1) / 2] <- max(y_scaled)
  # Return
  return(m)
}

#--
#Generating traces for each cluster
#--

#
#~AR subtype
sampnum <- 94 #number of samples in the subtype
x <- seq(0, sampnum, by = 0.1)
x <- seq(0, sampnum, by = 0.4)
#calculating what the sd of the normal distribution will be
clus_enrich_sum <- sig_file %>%
 dplyr::filter(cluster == "C5") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(val) %>%
    sum()
clus_stability <- 1.062215
clus_prop <- sig_file %>%
 dplyr::filter(cluster == "C5") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(cluster_prop) %>%
    sum()
overall_prop <- sig_file %>%
 dplyr::filter(cluster == "C5") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(overall_prop) %>%
    sum()
sd <- clus_enrich_sum + clus_stability + clus_prop + overall_prop

y <- dnorm(x, mean = (sampnum / 2), sd = (20))
#lets make the max value of the normal distribution equal
# to the number samples in the subtype
ymax <- y[which.max(y)] # extracting the largest value
AR_y_scaled <- y * (sd / ymax) # this will scale the normal distribution
           #so the peak is equal to the number of samples in the subtype.

#checking the distribution
ARplot <- plot(x, AR_y_scaled, main = "Normal Distribution", col = "blue")

#making the hill matrix
AR_k <- as.numeric(length(AR_y_scaled))
elevparam <- (((AR_k - 2):1)[!(((AR_k - 2):1) %% 2 == 0)]) / 2
 #keeping only odd numbers

qAR <- make_hill(k = AR_k, y_scaled = AR_y_scaled)

#
#~TT subtype~
sampnum <- 101 #number of samples in the subtype
x <- seq(0, sampnum, by = 0.4)
#calculating what the sd of the normal distribution will be
clus_enrich_sum <- sig_file %>%
 dplyr::filter(cluster == "C2") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(val) %>%
    sum()
clus_stability <- 0.585277
clus_prop <- sig_file %>%
 dplyr::filter(cluster == "C2") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(cluster_prop) %>%
    sum()
overall_prop <- sig_file %>%
 dplyr::filter(cluster == "C2") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(overall_prop) %>%
    sum()
sd <- clus_enrich_sum + clus_stability + clus_prop + overall_prop

y <- dnorm(x, mean = (sampnum / 2), sd = (20))
#lets make the max value of the normal distribution equal
# to the number samples in the subtype
ymax <- y[which.max(y)] # extracting the largest value
TT_y_scaled <- y * (sd / ymax) # this will scale the normal distribution
           # so the peak is equal to the number of samples in the subtype.

#checking the distribution
TTplot <- plot(x, TT_y_scaled, main = "Normal Distribution", col = "blue")

#making the hill matrix
TT_k <- as.numeric(length(TT_y_scaled))
elevparam <- (((TT_k - 2):1)[!(((TT_k - 2):1) %% 2 == 0)]) / 2
 #keeping only odd numbers

qTT <- make_hill(k = TT_k, y_scaled = TT_y_scaled)

#
#~Q subtype~
sampnum <- 124 #number of samples in the subtype
x <- seq(0, sampnum, by = 0.4)
#calculating what the sd of the normal distribution will be
clus_enrich_sum <- sig_file %>%
 dplyr::filter(cluster == "C4") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(val) %>%
    sum()
clus_stability <- 1.005452
clus_prop <- sig_file %>%
 dplyr::filter(cluster == "C4") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(cluster_prop) %>%
    sum()
overall_prop <- sig_file %>%
 dplyr::filter(cluster == "C4") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(overall_prop) %>%
    sum()
sd <- clus_enrich_sum + clus_stability + clus_prop + overall_prop
sd <- (clus_enrich_sum + clus_stability + clus_prop + overall_prop) * 4

y <- dnorm(x, mean = (sampnum / 2), sd = (20))
#lets make the max value of the normal distribution equal
# to the number samples in the subtype
ymax <- y[which.max(y)] # extracting the largest value
Q_y_scaled <- y * (sd / ymax) # this will scale the normal distribution
           # so the peak is equal to the number of samples in the subtype.

#checking the distribution
Qplot <- plot(x, Q_y_scaled, main = "Normal Distribution", col = "blue")

#making the hill matrix
Q_k <- as.numeric(length(Q_y_scaled))
elevparam <- (((Q_k - 2):1)[!(((Q_k - 2):1) %% 2 == 0)]) / 2
 #keeping only odd numbers

qQ <- make_hill(k = Q_k, y_scaled = Q_y_scaled)
qQ_new <- replace(qQ, qQ == 1, 0) #getting rid of excess 1s

#
#~GM subtype~
sampnum <- 137 #number of samples in the subtype
x <- seq(0, sampnum, by = 0.4)
#calculating what the sd of the normal distribution will be
clus_enrich_sum <- sig_file %>%
 dplyr::filter(cluster == "C3") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(val) %>%
    sum()
clus_stability <- 0.865383
clus_prop <- sig_file %>%
 dplyr::filter(cluster == "C3") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(cluster_prop) %>%
    sum()
overall_prop <- sig_file %>%
 dplyr::filter(cluster == "C3") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(overall_prop) %>%
    sum()
sd <- clus_enrich_sum + clus_stability + clus_prop + overall_prop

y <- dnorm(x, mean = (sampnum / 2), sd = (20))
#lets make the max value of the normal distribution equal
# to the number samples in the subtype
ymax <- y[which.max(y)] # extracting the largest value
GM_y_scaled <- y * (sd / ymax) # this will scale the normal distribution
            #so the peak is equal to the number of samples in the subtype.

#checking the distribution
GMplot <- plot(x, GM_y_scaled, main = "Normal Distribution", col = "blue")

#making the hill matrix
GM_k <- as.numeric(length(GM_y_scaled))
elevparam <- (((GM_k - 2):1)[!(((GM_k - 2):1) %% 2 == 0)]) / 2
 #keeping only odd numbers

qGM <- make_hill(k = GM_k, y_scaled = GM_y_scaled)

#
#~CS subtype~
sampnum <- 257 #number of samples in the subtype
x <- seq(0, sampnum, by = 0.4)
#calculating what the sd of the normal distribution will be
clus_enrich_sum <- sig_file %>%
 dplyr::filter(cluster == "C1") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(val) %>%
    sum()
clus_stability <- 0.619218
clus_prop <- sig_file %>%
 dplyr::filter(cluster == "C1") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(cluster_prop) %>%
    sum()
overall_prop <- sig_file %>%
 dplyr::filter(cluster == "C1") %>%
  dplyr::filter(val > 10) %>%
   dplyr::select(overall_prop) %>%
    sum()
sd <- clus_enrich_sum + clus_stability + clus_prop + overall_prop

y <- dnorm(x, mean = (sampnum / 2), sd = (50))
#lets make the max value of the normal distribution equal
# to the number samples in the subtype
ymax <- y[which.max(y)] # extracting the largest value
CS_y_scaled <- y * (sd / ymax) # this will scale the normal distribution
            #so the peak is equal to the number of samples in the subtype.

#checking the distribution
CSplot <- plot(x, CS_y_scaled, main = "Normal Distribution", col = "blue")

#making the hill matrix
CS_k <- as.numeric(length(CS_y_scaled))
elevparam <- (((CS_k - 2):1)[!(((CS_k - 2):1) %% 2 == 0)]) / 2
 #keeping only odd numbers

qCS <- make_hill(k = CS_k, y_scaled = CS_y_scaled)

#--
# Scaling and positioning cluster traces according to axis parameters
#--

cs_padded <- rbind(matrix(0, ncol = 643, nrow = 643),
 qCS, matrix(0, ncol = 643, nrow = 500))
dim(cs_padded)
cs_padded_2 <- cbind(matrix(0, ncol = 589, nrow = 1786),
 cs_padded, matrix(0, ncol = 554, nrow = 1786))
dim(cs_padded_2)

tt_padded <- rbind(matrix(0, ncol = 253, nrow = 1049),
 qTT, matrix(0, ncol = 253, nrow = 484))
dim(tt_padded)
tt_padded_2 <- cbind(matrix(0, ncol = 811, nrow = 1786),
 tt_padded, matrix(0, ncol = 722, nrow = 1786))
dim(tt_padded_2)

gm_padded <- rbind(matrix(0, ncol = 343, nrow = 641),
 qGM, matrix(0, ncol = 343, nrow = 802))
dim(gm_padded)
gm_padded_2 <- cbind(matrix(0, ncol = 128, nrow = 1786),
 gm_padded, matrix(0, ncol = 1315, nrow = 1786))
dim(gm_padded_2)

q_padded <- rbind(matrix(0, ncol = 311, nrow = 568),
 qQ_new, matrix(0, ncol = 311, nrow = 907))
dim(q_padded)
q_padded_2 <- cbind(matrix(0, ncol = 1037, nrow = 1786),
 q_padded, matrix(0, ncol = 438, nrow = 1786))
dim(q_padded_2)
q_padded_3 <- (q_padded_2 + 5.33)
#increased all values by 5 to lift this cluster,
# so it won't get cut off by the margins of the z-axis in the 3D-plot

ar_padded <- rbind(matrix(0, ncol = 236, nrow = 1025),
 qAR, matrix(0, ncol = 236, nrow = 525))
dim(ar_padded)
ar_padded_2 <- cbind(matrix(0, ncol = 1190, nrow = 1786),
 ar_padded, matrix(0, ncol = 360, nrow = 1786))
dim(ar_padded_2)

write.table(cs_padded_2, file = "cs_padded.txt",
 sep = ";", row.names = FALSE, col.names = FALSE)
write.table(tt_padded_2, file = "tt_padded.txt",
 sep = ";", row.names = FALSE, col.names = FALSE)
write.table(gm_padded_2, file = "gm_padded.txt",
 sep = ";", row.names = FALSE, col.names = FALSE)
write.table(q_padded_3, file = "q_padded.txt",
 sep = ";", row.names = FALSE, col.names = FALSE)
write.table(ar_padded_2, file = "ar_padded.txt",
 sep = ";", row.names = FALSE, col.names = FALSE)
