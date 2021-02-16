##Heatmap plots R codes
##Authors: Upasna srivastava
##Affiliations:tsuda lab, University of Tokyo

##Load Dseq2 object of C1 scTimeSeries data and then run these commands

##for day3_day0
topVarianceGenes_day3_day0 <- head(order(rowVars(assay(vsd.sub_day0_day3)), decreasing=T),35)
matrix_day3_day0 <- assay(vsd.sub_day0_day3)[ topVarianceGenes_day3_day0, ]
matrix_day3_day0 <- matrix_day3_day0 - rowMeans(matrix_day3_day0)
annotation_data_day3_day0 <- as.data.frame(colData(vsd.sub_day0_day3)[c("timepoint")])
pheatmap(matrix_day3_day0, annotation_col=annotation_data_day3_day0)

library(pheatmap)
library(gplots)
if (nrow(matrix_day3_day0) > 100) stop("Too many rows for heatmap, who can read?!")
fontsize_row = 8 - nrow(matrix_day3_day0) / 15
pheatmap(matrix_day3_day0,annotation_col=annotation_data_day3_day0, col=greenred(256), main="day3_day0_Top35_DEG", cluster_cols=T,
         fontsize_row=fontsize_row, border_color=T)

#### for day6_vs_day0

topVarianceGenes_day6_day0 <- head(order(rowVars(assay(vsd.sub_day0_day6)), decreasing=T),35)
matrix_day6_day0 <- assay(vsd.sub_day0_day6)[ topVarianceGenes_day6_day0, ]
matrix_day6_day0 <- matrix_day6_day0 - rowMeans(matrix_day6_day0)
annotation_data_day6_day0 <- as.data.frame(colData(vsd.sub_day0_day6)[c("timepoint")])
pheatmap(matrix_day6_day0, annotation_col=annotation_data_day6_day0)

library(pheatmap)
library(gplots)
if (nrow(matrix_day6_day0) > 100) stop("Too many rows for heatmap, who can read?!")
fontsize_row = 8 - nrow(matrix_day6_day0) / 15
pheatmap(matrix_day6_day0,annotation_col=annotation_data_day6_day0, col=greenred(256), main="day6_day0_Top35_DEG", cluster_cols=T,
         fontsize_row=fontsize_row, border_color=T)

  ## day9_vs_day0
  topVarianceGenes_day9_day0 <- head(order(rowVars(assay(vsd.sub_day0_day9)), decreasing=T),35)
  matrix_day9_day0 <- assay(vsd.sub_day0_day9)[ topVarianceGenes_day9_day0, ]
  matrix_day9_day0 <- matrix_day9_day0 - rowMeans(matrix_day9_day0)
  annotation_data_day9_day0 <- as.data.frame(colData(vsd.sub_day0_day9)[c("timepoint")])
  pheatmap(matrix_day9_day0, annotation_col=annotation_data_day9_day0)

  library(pheatmap)
  library(gplots)
  if (nrow(matrix_day9_day0) > 100) stop("Too many rows for heatmap, who can read?!")
  fontsize_row = 8 - nrow(matrix_day9_day0) / 15
  pheatmap(matrix_day9_day0,annotation_col=annotation_data_day9_day0, col=greenred(256), main="day9_day0_Top35_DEG", cluster_cols=T,
           fontsize_row=fontsize_row, border_color=T)


> ## day3_vs_day6
 topVarianceGenes_day3_day6 <- head(order(rowVars(assay(vsd.sub_day3_day6)), decreasing=T),50)
matrix_day3_day6 <- assay(vsd.sub_day3_day6)[ topVarianceGenes_day3_day6, ]
 matrix_day3_day6 <- matrix_day3_day6 - rowMeans(matrix_day3_day6)
 annotation_data_day3_day6 <- as.data.frame(colData(vsd.sub_day3_day6)[c("timepoint")])
 pheatmap(matrix_day3_day6, annotation_col=annotation_data_day3_day6)

 library(pheatmap)
 library(gplots)
 if (nrow(matrix_day3_day6) > 100) stop("Too many rows for heatmap, who can read?!")
fontsize_row = 8 - nrow(matrix_day3_day6) / 15
pheatmap(matrix_day3_day6, annotation_col=annotation_data_day3_day6, col=greenred(256), main="day3_day6_Top50_DEG", cluster_cols=T, fontsize_row=fontsize_row, border_color=T)

###
## day6_vs_day9
topVarianceGenes_day6_day9 <- head(order(rowVars(assay(vsd.sub_day6_day9)), decreasing=T),50)
matrix_day6_day9 <- assay(vsd.sub_day6_day9)[ topVarianceGenes_day6_day9, ]
matrix_day6_day9 <- matrix_day6_day9 - rowMeans(matrix_day6_day9)
annotation_data_day6_day9 <- as.data.frame(colData(vsd.sub_day6_day9)[c("timepoint")])
pheatmap(matrix_day6_day9, annotation_col=annotation_data_day6_day9)

library(pheatmap)
library(gplots)
if (nrow(matrix_day6_day9) > 100) stop("Too many rows for heatmap, who can read?!")
fontsize_row = 8 - nrow(matrix_day6_day9) / 15
pheatmap(matrix_day6_day9, annotation_col=annotation_data_day6_day9, col=greenred(256), main="day6_day9_Top50_DEG", cluster_cols=T,
         fontsize_row=fontsize_row, border_color=T)


## day6_vs_day9
topVarianceGenes_day6_day9 <- head(order(rowVars(assay(vsd.sub_day6_day9)), decreasing=T),60)
matrix_day6_day9 <- assay(vsd.sub_day6_day9)[ topVarianceGenes_day6_day9, ]
matrix_day6_day9 <- matrix_day6_day9 - rowMeans(matrix_day6_day9)
annotation_data_day6_day9 <- as.data.frame(colData(vsd.sub_day6_day9)[c("timepoint")])
pheatmap(matrix_day6_day9, annotation_col=annotation_data_day6_day9)

library(pheatmap)
library(gplots)
if (nrow(matrix_day6_day9) > 100) stop("Too many rows for heatmap, who can read?!")
fontsize_row = 8 - nrow(matrix_day6_day9) / 15
pheatmap(matrix_day6_day9, annotation_col=annotation_data_day6_day9, col=greenred(256), main="day6_day9_Top60_DEG", cluster_cols=T,
         fontsize_row=fontsize_row, border_color=T)


 ## all days
 topVarianceGenes <- head(order(rowVars(assay(vsd)), decreasing=T),60)
 matrix <- assay(vsd)[ topVarianceGenes, ]
 matrix <- matrix - rowMeans(matrix)
 annotation_data <- as.data.frame(colData(vsd)[c("timepoint")])
pheatmap(matrix, annotation_col=annotation_data)

 library(pheatmap)
 library(gplots)
 if (nrow(matrix) > 100) stop("Too many rows for heatmap, who can read?!")
 fontsize_row = 8 - nrow(matrix) / 15
 pheatmap(matrix, annotation_col=annotation_data, col=greenred(256), main="CancerTimeSeries_Top60_DEG", cluster_cols=T,
     fontsize_row=fontsize_row, border_color=T)

 ## all days
 topVarianceGenes <- head(order(rowVars(assay(vsd)), decreasing=T),50)
 matrix <- assay(vsd)[ topVarianceGenes, ]
matrix <- matrix - rowMeans(matrix)
annotation_data <- as.data.frame(colData(vsd)[c("timepoint")])
pheatmap(matrix, annotation_col=annotation_data)

library(pheatmap)
library(gplots)
 if (nrow(matrix) > 100) stop("Too many rows for heatmap, who can read?!")
 fontsize_row = 8 - nrow(matrix) / 15
pheatmap(matrix, annotation_col=annotation_data, col=greenred(256), main="CancerTimeSeries_Top50_DEG", cluster_cols=T,
  fontsize_row=fontsize_row, border_color=T)
