#---
# This script plots the Gallium dataset analysis for Figure 4
# Author: Victoria Shelton
#---

packages <- c("dplyr", "readr", "ggpubr", "survival", "survminer")
lapply(packages, require, character.only = TRUE)

date <- Sys.Date()

setwd("~/your working directory/GALLIUM/")

# Read in R objects of the Gallium cohort
p1 <- readRDS("MutationsClustermapGALLIUM.rds")
p2 <- readRDS("KMcurve_TT_vs_Other.rds")
p3 <- readRDS("MultivariateCox_TT_vs_Other.rds")

#Saving plot 1 -  the clustermap
ggsave(file = paste0(date, "_GALLIUM_Clustermap_aic.png"), p1,
 width = 17, height = 14, units = "cm")

# Updating plot 2 labels
p2[["plot"]][["labels"]][["x"]] <- "Time (years)"
p2[["plot"]][["labels"]][["title"]] <- "GALLIUM - overall cohort"
p2[["table"]][["labels"]][["x"]] <- "Time (years)"
p2[["table"]][["labels"]][["title"]] <- "Number at risk"

# Updating plot 2 break points
p2[["plot"]][["plot_env"]][["times"]] <- c(0, 2, 4, 6)
p2[["plot"]][["plot_env"]][["xlim"]] <- c(0, 6)
p2[["plot"]][["plot_env"]][["xticklabels"]] <- c(0, 2, 4, 6)
p2[["plot"]][["plot_env"]][["break.time.by"]] <- 2
p2[["table"]][["plot_env"]][["times"]] <- c(0, 2, 4, 6)

# removing legend in plot 2
p2$plot <- p2$plot + theme_survminer(
  legend = "none"
)
p2$table <- p2$table + theme_survminer(
  legend = "none"
)

#Saving plot 2 -  progrssion free survival & plot 3 - hazard ratio
gg <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p2$plot, p2$table,
                                heights = c(0.9, 0.4), ncol = 1, nrow = 2),
  p3, nrow = 1, ncol = 2, widths = c(0.60, 1))

ggsave(file = paste0(date, "_GALLIUM_plots_TT_vs_Other.pdf"), gg,
 width = 25, height = 15, units = "cm")
