# Date: 26 June 2020
# Function: Using output from gprofiler2() function, create a plot with term names and p value. 
#           If more than 50 terms are present, only the first 50 terms with smallest p value are plotted. 
# Author: Anjali Silva and Robert Kridel

# Input:
# GeneIDs: A list of genes to be input as query. 
# Organism: Character string indicating the name of organism as defined by gprofiler2::gost().
# OrderedQuery: logical indicating TRUE or FALSE, indicating if the gene list is ranked. 
#                Default FALSE. Default: "hsapiens".
# PositiveorNegFC: A character string indicating if genes corresponding to "negative" or 
#                  "positive" FC. Default NA. 
# PvalAlphaLevel: Alpha level as threshold for significance. Default 0.01.
# ConditionName: Character string indicating the name of condition or contrast used to generate
#                the list of genes. This will be the used on x-axis of graph. Default: "Condition".
# FigureGenerate: Produce images or not, options = "Yes" or "No"; default "Yes" 
# PNGorPDF: Output format of the image, options = "png" or "pdf"; default "png"

# Output:
# shortLink: A link that can be pasted on web to see results. 

# Visuals saved to img folder
# 26_gProfiler_", ConditionName.*

GProfilerAnalysis <- function(GeneIDs,
                      Organism = "hsapiens",
                      OrderedQuery = FALSE,
                      PvalAlphaLevel = 0.01,
                      PositiveorNegFC = NA,
                      ConditionName = "condition",
                      ProduceImages = "Yes", 
                      PNGorPDF = "png") {
  
  library(gprofiler2)
  library(grDevices)
  library(ggplot2)
  
  set.seed(1234)
  TablePlot <- Table <- gprofiler2::gost(query = GeneIDs, 
                                         organism = Organism, 
                                         ordered_query = OrderedQuery,
                                         user_threshold = PvalAlphaLevel)  %>%
                                         .$result %>%
                                         arrange(p_value) %>%
                                         mutate(log.p = -log10(p_value)) %>%
                                         mutate(term.name = paste(term_id, term_name, sep = "\n")) %>%
                                         mutate(term.name = fct_reorder(term.name, p_value, .desc = TRUE)) %>%
                                         dplyr::select(term_name, term_size, 
                                                       query_size, intersection_size, log.p) %>%
                                         mutate(condition = ConditionName)
  
  shortLink <- gprofiler2::gost(query = GeneIDs, 
                                organism = Organism, 
                                ordered_query = OrderedQuery,
                                user_threshold = PvalAlphaLevel, 
                                as_short_link = TRUE)
  # If more than 50 terms present, select top 50 only to plot based on p value (smallest first)
  if(nrow(TablePlot) > 50) {
    TablePlot <- TablePlot[c(1:50), ]
  }
  
  # define colours
  if((is.na(PositiveorNegFC) == FALSE) && (PositiveorNegFC == "positive")) {
    colours = c("#67000d", "#a50f15", "#cb181d", "#ef3b2c", "white")
  } else if((is.na(PositiveorNegFC) == FALSE) && (PositiveorNegFC == "negative")) {
    colours = c("#023858", "#045a8d", "#0570b0", "#74a9cf", "white")
  } else {
    colours = c("#e5f5f9", "#99d8c9", "#66c2a4", "#41ae76", "white")
  }
  
  par(mfrow = c(1, 1))
  TablePlot %>%
    ggplot2::ggplot(aes(x = factor(condition), y = term_name, fill = log.p, colour = factor(condition))) + 
          geom_tile(colour = "white") +
          labs(x = "Category", y = "gProfiler Term Name", fill = "-log10(P value)") +
          # geom_text(aes(label = intersection_size), size = 2) +
          theme_classic() +
          scale_fill_gradientn(values = c(1.0, 0.75, 0.5, 0.25, 0), 
                               colours = colours, 
                               na.value = "white") +  #facet_wrap(~ cell_line, scales = "free_x") + theme_classic() +
          theme(axis.title.x = element_blank(), 
                axis.text.x = element_text(size = 9, angle = 90, hjust = 1),
                axis.title.y = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = 6),
                panel.background = element_rect(colour = "black", size = 0.5), 
                strip.text.x = element_text(size = 8), strip.background = element_blank())

   pathNow <- getwd() # getting the path
   if (PNGorPDF == "pdf") {
      ggplot2::ggsave(paste0(pathNow, "/img/26_gProfiler_", ConditionName, ".pdf"))
    } else {
      ggplot2::ggsave(paste0(pathNow, "/img/26_gProfiler_", ConditionName, ".png"))
    }
  
    grDevices::dev.off()
    
    

   RESULTS <- list(shortLink = shortLink,
                   TableAllValues = Table)
   
   class(RESULTS) <- "GProfilerAnalysis_ASilva"
   return(RESULTS)

  }
