---
output:
  pdf_document: default
  html_document: default
---
##title: "Gene expression Analysis pipeline for colon cancer C1 time series Single cell RNAseq dataset"
##author: "upasna srivastava"
##date: "feb 2021"
## Tsuda lab, University of Tokyo

output:
  pdf_document:
    fig_height: 7
    fig_width: 7
    toc: yes
  html_document:
    code_folding: hide
    fig_height: 7
    fig_width: 7
    toc: yes
    toc_float: yes
params:
  eval_no_correction: no
  show_preproc: no

#if(dev.cur() == 1) dev.new()

`{r options, echo=FALSE, results="hide",message=FALSE, error=FALSE, include=FALSE, autodep=TRUE}
knitr::opts_chunk$set(fig.align="center", cache=TRUE, error=FALSE, message=FALSE, warning=TRUE,echo=FALSE)
library(DESeq2)
library(RColorBrewer)
library(BiocParallel)
library(doParallel)
library(Matrix)
library(ggplot2)
library(gplots)
library(matrixStats)
library(zinbwave)
library(SummarizedExperiment)
library(clusterExperiment)
library(dplyr)
library(biomaRt)
library(scran)
library(edgeR)
library(ggpubr)
library(zingeR)
source("https://bioconductor.org/biocLite.R")
biocLite("GSEABase")
biocLite("GSVA")
biocLite("gage")
library(GSEABase)
library(GSVA)
library(gage)
library(mvoutlier)
library(limma)
library(destiny)
library(scater)
library(gridExtra)
library(BiocStyle)
library(rmarkdown)
library(geneplotter)
library(ggplot2)
library(plyr)
library(LSD)
library(gplots)
library(RColorBrewer)
library(stringr)
library(topGO)
library(genefilter)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)
library(S4Vectors)
library(org.Mm.eg.db)
library( "org.Hs.eg.db" )
library( "reactome.db" )
library("GenomicFeatures")
library( "genefilter" )
library("GenomicAlignments")
library("BiocParallel")
library("magrittr")
library(datasets)
library(knitr)
library(ggplot2)
library(GGally)
library(DESeq2)
library(edgeR)
library(RColorBrewer)
library(devtools)
install_github("raivokolde/pheatmap")
BiocManager::install("PoiClaClu")
library("PoiClaClu")
library(RColorBrewer)
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
devtools::install_github("kevinblighe/EnhancedVolcano")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
#Note: to install development version:
devtools::install_github("kevinblighe/EnhancedVolcano")
#2.2 2. Load the package into R session
library(EnhancedVolcano)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("heatmaps", version = "3.8")
BiocManager::install("ComplexHeatmap")
source("http://bioconductor.org/biocLite.R")
biocLite("CountClust")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap", version = "3.8")
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(timepoint))])

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
set.seed(20) #set seed, just in case
countdata <- read.table("/Users/srivastava/Desktop/clean_Day0_3_6_9_September_counts.tsv", header=TRUE ,  row.name="Geneid")
countdata
dim(countdata)

countdata = countdata[rowSums(countdata) >10,]
dim(countdata)

countdata <- as.matrix(countdata)

head(countdata)
(timepoint <- factor(c(rep("day0", 96), rep("day3", 96), rep("day6", 96), rep("day9", 96))))
(coldata <- data.frame(row.names=colnames(countdata), timepoint))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~timepoint)

dds <- estimateSizeFactors(dds)
dds$timepoint <- factor(dds$timepoint)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
res_day3_day0<-results(dds, contrast= c("timepoint", "day3", "day0"))
 res_day6_day0<-results(dds, contrast= c("timepoint", "day6", "day0"))
 res_day9_day0<-results(dds, contrast= c("timepoint", "day9", "day0"))
 
 #res_day3_day6<-results(dds, contrast= c("timepoint", "day3", "day6"))
 #res_day9_day6<-results(dds, contrast= c("timepoint", "day9", "day6"))
 
 ####################
 
 # removing of NA containg lines at padj
 
 res_day3_day0_filt <- res_day3_day0[!is.na(res_day3_day0$padj),]
 res_day6_day0_filt <- res_day6_day0[!is.na(res_day6_day0$padj),]
 res_day9_day0_filt <- res_day9_day0[!is.na(res_day9_day0$padj),]
 res_day3_day0_filt[1:5]

# order by adjusted P-value
res_day3_day0_filt_ordered <- res_day3_day0_filt [order(res_day3_day0_filt$padj),]
res_day6_day0_filt_ordered <- res_day6_day0_filt [order(res_day6_day0_filt$padj),]
res_day9_day0_filt_ordered <- res_day9_day0_filt [order(res_day9_day0_filt$padj),]

#################################################

#Summary for up and down regulated genes
 summary(res_day3_day0_filt, alpha = 0.05)

 summary(res_day6_day0_filt, alpha = 0.05)

 summary(res_day9_day0_filt, alpha = 0.05)

#######################################################################
head(res_day3_day0_filt_ordered)

head(res_day6_day0_filt_ordered)

#heatmap to get top 50 DE Genes
de.genes_day3_day0 <- rownames(res_day3_day0_filt_ordered) [1:35]
de.genes_day6_day0 <- rownames(res_day6_day0_filt_ordered) [1:35]
de.genes_day9_day0 <- rownames(res_day9_day0_filt_ordered) [1:35]

###
vsd <- vst(dds)
norm.data <- assay(vsd)
colnames(norm.data) <- colData(vsd)$timepoint
norm.de_day3_day0 <- norm.data[row.names(norm.data) %in% de.genes_day3_day0,]
norm.de_day6_day0 <- norm.data[row.names(norm.data) %in% de.genes_day6_day0,]
norm.de_day9_day0 <- norm.data[row.names(norm.data) %in% de.genes_day9_day0,]

######################################################################################################################

##heatmap plot
library(gplots)
heatmap_day3_day0 <- heatmap.2(norm.de_day3_day0, trace ="none", cexCol = 0.9, main = "heatmap Plot top 50 DEG day3_vs_day0")
dev.copy(heatmap_day3_day0.png,'heatmap_day3_day0.png')
heatmap_day6_day0 <- heatmap.2(norm.de_day6_day0, trace ="none", cexCol = 0.9, main = "heatmap Plot top 50 DEG day6_vs_day0")
dev.copy(heatmap_day6_day0.png,'heatmap_day6_day0.png')
heatmap_day9_day0 <- heatmap.2(norm.de_day9_day0, trace ="none", cexCol = 0.9, main = "heatmap Plot top 50 DEG day9_vs_day0")
dev.copy(heatmap_day9_day0.png,'heatmap_day9_day0.png')
##########################################################################
##cluster dendrogram of top 50 DE genes
d_genes <- dist(norm.de_day3_day0)
plot(hclust(d_genes)) ## Do clustering and split the samples into 2 groups
rect.hclust(hclust(d_genes),k=5)


d_genes <- dist(norm.de_day6_day0)
plot(hclust(d_genes))
## Do clustering and split the samples into 2 groups
rect.hclust(hclust(d_genes),k=5)

d_genes <- dist(norm.de_day9_day0)
plot(hclust(d_genes))
## Do clustering and split the samples into 2 groups
rect.hclust(hclust(d_genes),k=5)
 dev.off()




sum( res_day3_day0_filt$padj < 0.1, na.rm=TRUE )

sum( res_day6_day0_filt$padj < 0.1, na.rm=TRUE )

sum( res_day9_day0_filt$padj < 0.1, na.rm=TRUE )

table(res$padj<=0.05)
table(res_day3_day0$padj<=0.05)
table(res_day6_day0$padj<=0.05)
table(res_day9_day0$padj<=0.05)

###############################################

cat(summary(res))
significance= rep(NA, nrow(res))
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

names(resdata)[1] <- "gene"
resdata = data.frame(resdata,significance)
print(head(resdata))
resdata[1:10,1:7]
resdata <- resdata[!is.na(resdata$padj),]
resdata[1:10,1:7]
expression_matrix <- resdata[,c(1,8:ncol(resdata))]
print(head(expression_matrix))
mcols(res, use.names=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsdMat <- assay(vsd)
write.table(as.data.frame(vsdMat), file="normalized_expression_values.txt", sep="\t", quote=F, col.names=T, row.names = T)
vsd.sub_day0_day3 <-vsd[,vsd@colData$timepoint %in% c("day0","day3")]
vsd.sub_day0_day3
row.names(vsd.sub_day0_day3)
vsd.sub_day0_day6 <-vsd[,vsd@colData$timepoint %in% c("day0","day6")]
vsd.sub_day0_day9 <-vsd[,vsd@colData$timepoint %in% c("day0","day9")]
vsd.sub_day6_day9 <-vsd[,vsd@colData$timepoint %in% c("day6","day9")]
vsd.sub_day3_day6 <-vsd[,vsd@colData$timepoint %in% c("day3","day6")]
####################################################################
# order results by padj value (most significant to least)
res= subset(res, padj<0.05)
res <- res[order(res$padj),]
# order results by padj value (most significant to least)
res_day3_day0_filt= subset(res_day3_day0_filt, padj<0.05)
res_day3_day0_filt <- res_day3_day0_filt[order(res_day3_day0_filt$padj),]
# order results by padj value (most significant to least)
res_day6_day0_filt= subset(res_day6_day0_filt, padj<0.05)
res_day6_day0_filt <- res_day6_day0_filt[order(res_day6_day0_filt$padj),]
# order results by padj value (most significant to least)
res_day9_day0_filt= subset(res_day9_day0_filt, padj<0.05)
res_day9_day0_filt <- res_day9_day0_filt[order(res_day9_day0_filt$padj),]

##
sum( res$padj < 0.1, na.rm=TRUE )
resSig <- subset(res, res$padj < 0.1 )
head( resSig[ order( resSig$log2FoldChange ), ], 4) ## downregulated genes
##the strongest upregulation :
head( resSig[ order( -resSig$log2FoldChange ), ], 4) ##upregulated genes
#MA plot(The MA-plot represents each gene with a dot. The x axis is the average expression over all samples,
#the y axis the log2 fold change of normalized counts (i.e the average of counts normalized by size factor) between treatment and
#control. Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.
#This plot demonstrates that only genes with a large average normalized count contain sufficient information to yield a significant call.)

plotMA( res, ylim = c(-3, 3) )
#PCA plots
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(vsd))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))

# set condition
timepoint <- days
scores <- data.frame(pc$x, timepoint)

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(timepoint))))
 geom_point(size = 5)+
 ggtitle("Principal Components")
+ scale_colour_brewer(name = " ", palette = "Set1")
+ theme(
  plot.title = element_text(face = 'bold'),
  legend.position = c(.9,.2),
  legend.key = element_rect(fill = 'NA'),
  legend.text = element_text(size = 10, face = "bold"),
  axis.text.y = element_text(colour = "Black"),
  axis.text.x = element_text(colour = "Black"),
  axis.title.x = element_text(face = "bold"),
  axis.title.y = element_text(face = 'bold'),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_rect(color = 'black',fill = NA)
))

#ggsave(pcaplot,file=paste0(outputPrefix, "-ggplot2.pdf"))
#######################################################
plot(assay(vsd)[,1:3],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")


#Gene clustering
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(vsd) ), decreasing=TRUE ), 35 )
heatmap.2( assay(vsd)[ topVarGenes, ], scale="row",
     trace="none", dendrogram="column",
     col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
     ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[
        colData(vsd)$timepoint ] )


###################################################
## note: heatmap plot code using dds and vsd
##Heatmap
biocLite("pheatmap")
library("pheatmap")
vsd <- normTransform(dds)
head(assay(vsd), 3)
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing =TRUE)[1:35]
df <- as.data.frame(colData(dds)[,c("timepoint)])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, showrownames=TRUE, cluster_cols=FALSE, annotation_col=df)

### OR if you want to plot heatmap based on log2foldchange the try this:Instead of plotting assay(ntd), you'll want to instead plot lfcShrink(dds)$log2FoldChange.
#Realistically, use something like:

res = lfcShrink(dds)
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing =TRUE)[1:35]
df <- as.data.frame(colData(dds)[,c("timepoint)])
pheatmap(res$log2FoldChange,[select,], cluster_rows=FALSE, showrownames=TRUE, cluster_cols=FALSE, annotation_col=df)

########################################################################################################
##heatmap plot using ComplexHeatmap

require(ComplexHeatmap)
require(circlize)
require(cluster)

pamClusters <- pam(heat, k=5)

hmap <- Heatmap(heat,
        name="Transcript Z-score",
        col=colorRamp2(myBreaks, myCol),
        heatmap_legend_param=list(
              color_bar="continuous",
              legend_direction="horizontal",
              legend_width=unit(5,"cm"),
              title_position="topcenter",
              title_gp=gpar(fontsize=15, fontface="bold")),
        split=paste0("", pamClusters$clustering),
        row_title="Transcripts",
        row_title_side="left",
        row_title_gp=gpar(fontsize=15, fontface="bold"),
        show_row_names=FALSE,
        column_title="",
        column_title_side="top",
        column_title_gp=gpar(fontsize=15, fontface="bold"),
        column_title_rot=0,
        show_column_names=FALSE,
        clustering_distance_columns=function(x) as.dist(1-cor(t(x))),
        clustering_method_columns="ward.D2",
        clustering_distance_rows="euclidean",
        clustering_method_rows="ward.D2",
        row_dend_width=unit(30,"mm"),
        column_dend_height=unit(30,"mm"),
        top_annotation=colAnn,
        top_annotation_height=unit(1.75,"cm"),
        bottom_annotation=sampleBoxplot,
        bottom_annotation_height=unit(4, "cm"))

draw(hmap, heatmap_legend_side="top", annotation_legend_side="right")
#####################################################################
##most significant genes
mostSig <- res_day3_day0_filt[res_day3_day0_filt$padj < 0.05 & abs(res_day3_day0_filt$log2FoldChange) >= 1, ], n=35))

 upgenes_day3_day0 <- rownames(head(res_day3_day0_filt[ order( res_day3_day0_filt$log2FoldChange ), ], n=35))
row.names(upgenes_day3_day0)

down_genes_day3_day0 <- rownames(head(res_day3_day0_filt[ order( -res_day3_day0_filt$log2FoldChange ), ], n=35))
down_genes_day3_day0

##############################################################################
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("timepoint"))
#####################################################################
res_day3_day0_filt_order <- res_day3_day0_filt[order(res_day3_day0_filt$padj),]
res_day3_day0_filt_order
res_day3_day0_Ordered_upgenes <- rownames(head(res_day3_day0_filt_order[ order( res_day3_day0_filt_order$log2FoldChange >=1 ), ], n=35))
res_day3_day0_Ordered_upgenes
res_day3_day0_Ordered_downgenes <- rownames(head(res_day3_day0_filt_order[ order( -res_day3_day0_filt_order$log2FoldChange), ], n=35))
res_day3_day0_Ordered_downgenes

#########################################################################################
## Date 28 march 2019
author: "upasna srivastava"
##significant gene list at log2Foldchange(We subset the results table to
#these genes and then sort it by the log2 fold change estimate to get the significant
#genes with the strongest down-regulation.)
##for day3_day0
resSig_day3_day0 <- subset(res_day3_day0_filt, padj < 0.05)
resSig_log2foldchange_day3_day0_downregulated <- (resSig_day3_day0[ order( resSig_day3_day0$log2FoldChange ), ])
write.csv(resSig_log2foldchange_day3_day0_downregulated, file= "/Users/srivastava/Desktop/2031/resSig_log2foldchange_day3_day0_downregulated.csv")

resSig_day6_day0 <- subset(res_day6_day0_filt, padj < 0.05)
resSig_log2foldchange_day6_day0_downregulated <- (resSig_day6_day0[ order( resSig_day6_day0$log2FoldChange ), ])
resSig_log2foldchange_day6_day0_downregulated
write.csv(resSig_log2foldchange_day6_day0_downregulated, file= "/Users/srivastava/Desktop/2031/resSig_log2foldchange_day6_day0_downregulated.csv")

resSig_day9_day0 <- subset(res_day9_day0_filt, padj < 0.05)
resSig_log2foldchange_day9_day0_downregulated <- (resSig_day9_day0[ order( resSig_day9_day0$log2FoldChange ), ])
resSig_log2foldchange_day9_day0_downregulated
write.csv(resSig_log2foldchange_day9_day0_downregulated, file= "/Users/srivastava/Desktop/2031/resSig_log2foldchange_day9_day0_downregulated.csv")

##upregulated genes
resSig_day3_day0 <- subset(res_day3_day0_filt, padj < 0.05)
resSig_log2foldchange_day3_day0_upregulated <- (resSig_day3_day0[ order( -resSig_day3_day0$log2FoldChange ), ])
resSig_log2foldchange_day3_day0_upregulated
write.csv(resSig_log2foldchange_day3_day0_upregulated, file= "/Users/srivastava/Desktop/2031/resSig_log2foldchange_day3_day0_upregulated.csv")

resSig_day6_day0 <- subset(res_day6_day0_filt, padj < 0.05)
resSig_log2foldchange_day6_day0_upregulated <- (resSig_day6_day0[ order( -resSig_day6_day0$log2FoldChange ), ])
resSig_log2foldchange_day6_day0_upregulated
write.csv(resSig_log2foldchange_day6_day0_upregulated, file= "/Users/srivastava/Desktop/2031/resSig_log2foldchange_day6_day0_upregulated.csv")

resSig_day9_day0 <- subset(res_day9_day0_filt, padj < 0.05)
resSig_log2foldchange_day9_day0_upregulated <- (resSig_day9_day0[ order( -resSig_day9_day0$log2FoldChange ), ])
resSig_log2foldchange_day9_day0_upregulated
write.csv(resSig_log2foldchange_day9_day0_upregulated, file= "/Users/srivastava/Desktop/2031/resSig_log2foldchange_day9_day0_upregulated.csv")

#################################################################################
vsd.sub_day3_day0 <-vsd[,vsd@colData$timepoint %in% c("day3","day0")]
 vsd.sub_day3_day0
 PCA_day3_day0 <- plotPCA(vsd.sub_day3_day0, intgroup=c("timepoint"), returnData=TRUE)
head(PCA_day3_day0)
row.names(vsd.sub_day3_day0)
mat_day3_day0_downregulated <- assay(vsd.sub_day3_day0) [ row.names(resSig_log2foldchange_day3_day0_downregulated),]
head(mat_day3_day0_downregulated)
library(AnnotationDbi)
library(pheatmap)
library(gplots)
mat_day3_day0_downregulated <- mat_day3_day0_downregulated -rowMeans(mat_day3_day0_downregulated)
mat_day3_day0_downregulated



############################################################################################
mat_day3_day0 <- assay(vsd.sub_day3_day0) [ row.names(resSig_log2foldchange_day3_day0),]
head(mat_day3_day0,5)
topgene_day3_day0 <- rbind[(head(mat_day3_day0,15),tail(mat_day3_day0,15))]
tail(topgene_day3_day0)
pheatmap(topgene_day3_day0, annotation_col = annotationdata)

###################################################################################################
 ## date: 30 march 2019
mat_day6_day0 <- assay(vsd.sub_day0_day6) [ row.names(resSig_log2foldchange_day6_day0),]
topgene_day6_day0 <- rbind (head(mat_day6_day0,15),tail(mat_day6_day0,15))
row.names (topgene_day6_day0)
write.csv(topgene_day6_day0, file= "/Users/srivastava/Desktop/2030/topgene_day6_day0_matrix.csv")
annotationdata <- as.data.frame(colData(vsd.sub_day0_day6) [c("timepoint")])
annotationdata
pheatmap(topgene_day6_day0, annotation_col = annotationdata)

mat_day9_day0 <- assay(vsd.sub_day0_day9) [ row.names(resSig_log2foldchange_day9_day0),]
topgene_day9_day0 <- rbind (head(mat_day9_day0,15),tail(mat_day9_day0,15))
row.names (topgene_day9_day0)
write.csv(topgene_day9_day0, file= "/Users/srivastava/Desktop/2030/topgene_day9_day0_matrix.csv")
pheatmap(topgene_day9_day0, annotation_col = annotationdata)

#######################################################################################################
## EnhancedVolcano plot

EnhancedVolcano(resSig_log2foldchange_day6_day9_downregulated,

                lab = rownames(resSig_log2foldchange_day6_day9_downregulated),

                x = "log2FoldChange",

                y = "pvalue", title = "day6_vs._day9")

################################################################################################
## Rcode for generation of 35 most variable genes using vsd

topVarianceGenes <- head(order(rowVars(assay(vsd.sub_day0_day3)), decreasing=T),35)
matrix <- assay(vsd.sub_day0_day3)[ topVarianceGenes, ]
 matrix <- matrix - rowMeans(matrix)
annotation_data <- as.data.frame(colData(vsd.sub_day0_day3)[c("timepoint")])
 pheatmap(matrix, annotation_col=annotation_data)


topVarianceGenes <- head(order(rowVars(assay(vsd.sub_day0_day6)), decreasing=T),35)
 matrix <- assay(vsd.sub_day0_day6)[ topVarianceGenes, ]
 matrix <- matrix - rowMeans(matrix)
 annotation_data <- as.data.frame(colData(vsd.sub_day0_day6)[c("timepoint")])
 pheatmap(matrix, annotation_col=annotation_data)


 topVarianceGenes <- head(order(rowVars(assay(vsd.sub_day0_day9)), decreasing=T),35)
 matrix <- assay(vsd.sub_day0_day9)[ topVarianceGenes, ]
 matrix <- matrix - rowMeans(matrix)
 annotation_data <- as.data.frame(colData(vsd.sub_day0_day9)[c("timepoint")])
pheatmap(matrix, annotation_col=annotation_data)


 topVarianceGenes <- head(order(rowVars(assay(vsd.sub_day3_day6)), decreasing=T),35)
matrix <- assay(vsd.sub_day3_day6)[ topVarianceGenes, ]
matrix <- matrix - rowMeans(matrix)
annotation_data <- as.data.frame(colData(vsd.sub_day3_day6)[c("timepoint")])
pheatmap(matrix, annotation_col=annotation_data)


topVarianceGenes <- head(order(rowVars(assay(vsd.sub_day6_day9)), decreasing=T),35)
 matrix <- assay(vsd.sub_day6_day9)[ topVarianceGenes, ]
 matrix <- matrix - rowMeans(matrix)
annotation_data <- as.data.frame(colData(vsd.sub_day6_day9)[c("timepoint")])
 pheatmap(matrix, annotation_col=annotation_data)


#############################################################################
 colnames(ann) <- c("timepoint")
 colours <- list("timepoint"=c("day0"="red2","day3"="royalblue","day6"="limegreen","day9"="gold"))
 colAnn <- HeatmapAnnotation(df=ann, which="col", col=colours, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
 boxplotCol <- HeatmapAnnotation(boxplot=anno_boxplot(data.matrix(mat_day9_day0), border=TRUE, gp=gpar(fill="#CCCCCC"), pch=".", size=unit(2, "mm"), axis=TRUE, axis_side="left", axis_gp=gpar(fontsize=12)), annotation_width=unit(c(1, 5.0), "cm"), which="col")
 hmap <- Heatmap(mat_day9_day0, name = "expression",  col = greenred(75),
 show_row_names = FALSE, show_column_names = FALSE, cluster_rows = TRUE,
cluster_columns = TRUE, show_column_dend = TRUE, show_row_dend = TRUE,
row_dend_reorder = TRUE, column_dend_reorder = TRUE, clustering_method_rows = "ward.D2",
clustering_method_columns = "ward.D2", width = unit(100, "mm"),
top_annotation_height=unit(1.0,"cm"), top_annotation=colAnn,
bottom_annotation_height=unit(3, "cm"), bottom_annotation=boxplotCol)
 draw(hmap, heatmap_legend_side="left", annotation_legend_side="right")

 #####
