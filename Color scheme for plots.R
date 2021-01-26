#------------------------------------
# Color scheme for plots
#------------------------------------

# Create a custom color scale
suppressPackageStartupMessages({
  library(RColorBrewer)
  library(ggplot2)
})

## Color set, from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
colors <- c(Cluster1 = "#FFCC00", Cluster2 = "#990033", Cluster3 = "#33CC33", 
            Cluster4 = "#FF0066", Cluster5 = "#0066FF",
            Cluster6 = "#000080")

manual.scale <- scale_colour_manual(name = "", values = colors)
