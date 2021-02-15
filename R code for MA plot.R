##R code for MA plot with gene name
##Data: C1 single cell time series Data
## Authors Upasna Srivastava
## Affiliation: Tsuda lab, University of Tokyo
library(ggpubr)
library(gplots)
library(RColorBrewer)

## load Deseq2 object of colon c1 single cell time series Data

## MA plot for res_day3_day0
ggmaplot(res_day3_day0_filt, main = expression("day3" %->% "day0"),fdr = 0.05, fc = 2, size = 0.4,
palette = c("#B31B21", "#1465AC", "darkgray"),genenames = as.vector(res_day6_day0_filt$timepoint),
legend = "top", top = 20,font.label = c("bold", 11),font.legend = "bold",font.main = "bold",
ggtheme = ggplot2::theme_minimal())

## MA plot for res_day6_day0

ggmaplot(res_day6_day0_filt, main = expression("day6" %->% "day0"),fdr = 0.05, fc = 2, size = 0.4,
palette = c("#B31B21", "#1465AC", "darkgray"),genenames = as.vector(res_day6_day0_filt$timepoint),
legend = "top", top = 20,font.label = c("bold", 11),font.legend = "bold",font.main = "bold",
ggtheme = ggplot2::theme_minimal())


## MA plot for res_day9_day0
ggmaplot(res_day9_day0_filt, main = expression("day9" %->% "day0"),
fdr = 0.05, fc = 2, size = 0.4,palette = c("#B31B21", "#1465AC", "darkgray"),
genenames = as.vector(res_day6_day0_filt$timepoint),legend = "top", top = 20,
font.label = c("bold", 11),font.legend = "bold",font.main = "bold",
ggtheme = ggplot2::theme_minimal())
