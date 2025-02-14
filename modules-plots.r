# The MIT License (MIT)
# 
# Copyright (c) <2024> <DresdenConceptGenomeCenter>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# 
# contact: mathias.lesche(at)tu-dresden.de

# most of these plots are for expression analysis of DESeq2 and require a deseq object

# it's the same function as in the modules-deseq-report.r
# if you update one, update the other to keep it consistent
createfilename <- function(dir, firstpart, appendix, extension) {
  if (appendix != '') {
    filename <- file.path(dir, paste(firstpart, appendix, extension, sep = '.'))
  } else {
    filename <- file.path(dir, paste(firstpart, extension, sep = '.'))
  }
  return(filename)
}

# function removes various symbols from the deseq formula 
# and reduces it to the factors
# input
# ftemp ... deseqformula as string
getFormulavector <- function(ftemp) {
  require(magrittr)
  ftemp <- gsub('\\+', ',', ftemp) %>% gsub('~', '', .) %>% gsub(':', ',', .) %>%
    gsub('\\*', ',', .) %>% gsub(' ', '', .) %>% strsplit(',') %>% unlist()
  return(ftemp)
}


# add and change settings of a ggplot object for pca, mds ...
# ggtemp ... ggplot object
# title ... string for ggtitle
addStandards <- function(ggtemp, title, scaleaxis = T) {
  require(ggplot2)
  axislabelsize <- 0.9
  axistitlesize <- 0.9
  titlesize <- 0.95
  legendsize <- 0.85
  ggtemp <- ggtemp + ggtitle(title) + 
    theme(plot.title=element_text(vjust=1, size=rel(titlesize))) + # coord_fixed() +
    theme(legend.title=element_text(size=rel(legendsize)), legend.position = 'bottom', legend.key=element_blank(), legend.key.size=unit(0.7, "cm"), legend.text=element_text(size=rel(legendsize))) +
    theme(axis.text.x=element_text(size=rel(axislabelsize)), axis.text.y=element_text(size=rel(axislabelsize))) +
    theme(axis.title.x=element_text(size=rel(axistitlesize)), axis.title.y=element_text(size=rel(axistitlesize))) +
    theme(legend.box = 'vertical') +
    guides(col = guide_legend(ncol = 4, byrow = T))
  if (scaleaxis) {
    ggtemp <- ggtemp + scale_x_continuous(expand = c(.15, .15)) + scale_y_continuous(expand = c(.15, .15))
  }
  return(ggtemp)
}

addConditions <- function(ggtemp, conditions, shapes = 0) {
  require(ggplot2)
  if (length(conditions) == 1) {
    ggtemp <- ggtemp + geom_point(pch=1, size=2, aes_string(colour=conditions[1]))
  } else if (length(conditions)==2) {
    if (shapes < 7) {
      ggtemp <- ggtemp + geom_point(size=2, aes_string(colour=conditions[1], shape = conditions[2]))
    } else {
      ggtemp <- ggtemp + geom_point(size=2, aes_string(colour=conditions[1], shape = conditions[2])) + scale_shape_manual(values = 1:shapes)
    }
  } else if (length(conditions)==3) {
    if (shapes < 7) {
      ggtemp <- ggtemp + geom_point(size=2, aes_string(colour=conditions[1], shape = conditions[2], alpha = conditions[3])) + 
        scale_alpha_discrete(range=c(0.4, 1))
    } else {
      ggtemp <- ggtemp + geom_point(size=2, aes_string(colour=conditions[1], shape = conditions[2], alpha = conditions[3])) + 
        scale_alpha_discrete(range=c(0.4, 1)) + scale_shape_manual(values = 1:shapes)
    }
  } else {
    warning('Do not use more than 3 conditions!')
  }
  return(ggtemp)
}

addSampleLabel <- function(ggtemp, samplename) {
  require(ggplot2)
  if (samplename == "variable") {
    require(ggrepel)
    ggtemp <- ggtemp +  geom_text_repel(aes(label=name), size=3, colour="black")
  } else if (samplename == "text") {
    ggtemp <- ggtemp + geom_text(aes(label=name), size=3, vjust=2, colour="black") # previous size=4
  } else {
    ggtemp <- ggtemp
  }
  return(ggtemp)
}

# Input
# inputdata ... a vector with original id (has to be in dds); replace ID and filename
# ddstemp ... DESeq object
# conditions ... is a vector and entries should be column names from the colData(dds) or coldata
#                the order is import for the ggplot; it is: colour, shape, alpha
#                DON'T USE MORE THAN 3 CONDITIONS!!!
plotNormCounts <- function(ensemblid, geneid, ddstemp, conditions, filename) {
  require(DESeq2)
  require(ggplot2)
  require(ggrepel)
  theme_set(theme_bw(14))
  
  pp <- plotCounts(dds, ensemblid, intgroup = conditions, returnData = T)
  pp$Sample <- rownames(pp)
  ppplot <- ggplot(pp, aes_string(conditions[1], "count"))
  ppplot <- addStandards(ppplot, ggtitle(paste0("normalised counts for ", geneid, " (", ensemblid, ")")), F)
  
  if (length(conditions) > 1) {
    ppplot <- addConditions(ppplot, conditions, length(unique(pp[, conditions[2]])))  
  } else {
    ppplot <- addConditions(ppplot, conditions)  
  }

  ppplot <- ppplot + geom_text_repel(aes(label=Sample), size=3, colour="black") + ylab("normalised counts")
  if (is.character(filename)){
    ggsave(filename, ppplot)
  } else {
    ppplot
  }
}

# This PCA method plots a number of genes and shows their contribution of the sample's separation
# Input
# readtransform ... DESeqTransform object which comes from the transformation process (vsd, rld, ...)
# genelist ... data.frame which has two columns (Ensembl_ID and Gene_Symbol [from the anno function])
# nog ... number of genes used for pca calculation
# show ... number of genes used for showing in the pca bi plot
# condvec ... vector in the same order as the readtransform/correctedvsd to display the samples cond
# display_column ... string that explain the shown condition
# filename ... if provided plot is printed directly
# title ... title for pca
# correctedvsd ... matrix with batch corrected expression values
# see http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/#biplot
# above "Supplementary elements"
plotBiplot = function(readtransform, genelist = '', nog, show, condvec, display_column, filename = "", title = paste0("PCA - BiPlot (",nog, " genes, ", show, " genes labelled)"), correctedvsd = NULL) {
  require(factoextra)
  require(genefilter)
  require(ggplot2)
  
  if (is.matrix(correctedvsd)) {
    topVarGenes <- head(order(rowVars(correctedvsd),decreasing=TRUE), nog)
    mat <- correctedvsd[ topVarGenes, ]
  } else{
    topVarGenes <- head(order(rowVars(assay(readtransform)),decreasing=TRUE), nog)
    mat <- assay(readtransform)[ topVarGenes, ]
  }
  
  if (is.data.frame(genelist)) {
    genelist$Gene_Symbol <- make.unique(as.character(genelist$Gene_Symbol))
    mat <- merge(mat, genelist, by.x = "row.names", by.y = "Ensembl_ID")
    mat$Gene_Symbol <- ifelse(mat$Gene_Symbol == '', as.character(mat$Row.names), as.character(mat$Gene_Symbol))
    rownames(mat) <- mat$Gene_Symbol
    mat$Row.names <- NULL
    mat$Gene_Symbol <- NULL
  }
  
  pca <- prcomp(t(mat), center = TRUE, scale = TRUE)
  
  p <- fviz_pca_biplot(pca, geom.ind = c('point', 'text'), repel = T,
                       fill.ind = condvec, 
                       pointshape = 21, pointsize = 2,
                       invisible = "quali",
                       palette = "jco",
                       col.var = 'contrib',
                       select.var = list(contrib = show),
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  ) + labs(fill = display_column, color = "Contrib") + ggtitle(title) + theme(legend.position = 'bottom')
  
  if (filename == '') {
    p
  } else{
    ggsave(filename, p)
  }
}

# This PCA method plots a variable PCA. With comp1 and comp2 you can choose your principle component
# Input
# readtransform ... DESeqTransform object which comes from the transformation process (vsd, rld, ...)
# conditions ... is a vector and entries should be column names from the colData(dds) or coldata
#                the order is import for the ggplot; it is: colour, shape, alpha
#                DON'T USE MORE THAN 3 CONDITIONS!!!
# comp1 ... Principle Component e.g. 1 or 2 or 3
# comp2 ... Principle Component e.g. 1 or 2 or 3
# samplename ... print sample names in the plot can be "", "text" or "variable"
# filename ... if provided plot is printed directly
# ntop ... calculation is based on ntop most diverse genes
# title ... title for pca
plotPCAVariable = function(readtransform, conditions, comp1 = 1, comp2 = 2, samplename = "text", filename = "", ntop = 500, title = paste0("PCA (", comp1, "-", comp2, ") of top ", ntop, " most diverse genes"), correctedvsd = NULL) {
  require(DESeq2)
  require(ggplot2)
  require(genefilter)
  theme_set(theme_bw(14))
  
  if (!all(conditions %in% names(colData(readtransform)))) {
    stop("the argument 'conditions' should specify columns of colData(dds)")
  }
  
  if (is.matrix(correctedvsd)) {
    topVarGenes <- head(order(rowVars(correctedvsd),decreasing=TRUE), ntop)
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(correctedvsd[topVarGenes,]), center=TRUE, scale=TRUE)
  } else{
    topVarGenes <- head(order(rowVars(assay(readtransform)),decreasing=TRUE), ntop)
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(readtransform)[topVarGenes,]), center=TRUE, scale=TRUE)
  }
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  # add the intgroup factors together to create a new grouping factor
  intgroup.df <- as.data.frame(colData(readtransform)[, conditions, drop=FALSE])    
  
  group <- factor(apply( intgroup.df, 1, paste, collapse=" : "))
  
  # assembly the data for the plot
  pcadata <- data.frame(Var1=pca$x[,comp1], Var2=pca$x[,comp2], group=group, intgroup.df, name=colnames(readtransform))
  
  pca <- ggplot(pcadata, aes(Var1, Var2)) +
    xlab(paste0("PC", comp1, ": ",round(percentVar[comp1]*100),"% variance")) +
    ylab(paste0("PC", comp2, ": ",round(percentVar[comp2]*100),"% variance"))
  
  pca <- addStandards(pca, title)
  if (length(conditions) > 1) {
    pca <- addConditions(pca, conditions, length(unique(pcadata[, conditions[2]])))  
  } else {
    pca <- addConditions(pca, conditions)  
  }
  pca <- addSampleLabel(pca, samplename)
  
  if (filename == '') {
    pca
  } else {
    ggsave(filename, pca) 
  }
}

# Another plot, very similar to the PCA plot, can be made using the 
# multidimensional scaling (MDS) function in base R. This is useful
# when we don’t have a matrix of data, but only a matrix of distances.
# Here we compute the MDS plot for the distances calculated from the
# transformed counts (vsd, rlog)
# Input
# readtransform ... DESeqTransform which comes from the transformation process (vsd, rld, ...)
# conditions ... is a vector and entries should be column names from the colData(dds) or coldata
#                the order is import for the ggplot; it is: colour, shape, alpha
#                DON'T USE MORE THAN 3 CONDITIONS!!!
# samplename ... print sample names in the plot can be "", "text" or "variable"
# filename ... if provided plot is printed directly
# title ... title for plot
plotMDS <- function(readtransform, conditions, samplename = "text", filename = "", title = "MDS plot of Euclidean distance") {
  require(ggplot2)
  theme_set(theme_bw(14))
  
  if (!all(conditions %in% names(colData(readtransform)))) {
    stop("the argument 'conditions' should specify columns of colData(dds)")
  }
  
  sampleDists <- dist(t(assay(readtransform)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  mds <- data.frame(cmdscale(sampleDistMatrix))
  mds <- cbind(mds, as.data.frame(colData(readtransform)))
  mds$name <- rownames(mds)
  
  mdsplot <- ggplot(mds, aes(X1, X2))
  
  mdsplot <- addStandards(mdsplot, title)
  
  if (length(conditions) > 1) {
    mdsplot <- addConditions(mdsplot, conditions, length(unique(mds[, conditions[2]])))  
  } else {
    mdsplot <- addConditions(mdsplot, conditions)  
  }
  
  mdsplot <- addSampleLabel(mdsplot, samplename)
  
  ggsave(filename, mdsplot)
}

# Input
# deseqcounts ... DESeq2 countstable not normalised, counts(dds), has to come from ddstemp
# ddstemp ... DESeq2 object
# conditions ... is a vector and entries should be column names from the colData(dds) or coldata
#                the order is import for the ggplot; it is: colour, shape, alpha
#                DON'T USE MORE THAN 3 CONDITIONS!!!
# samplename ... print sample names in the plot can be "", "text" or "variable"
# filename ... if provided plot is printed directly
# title ... title for plot
plotPoissonDistance <- function(deseqcounts, ddstemp, conditions, samplename = "text", filename = "", title = "MDS plot of Poisson Distance"){
  require(PoiClaClu)
  require(ggplot2)
  theme_set(theme_bw(14))
  
  if (!all(conditions %in% names(colData(ddstemp)))) {
    stop("the argument 'conditions' should specify columns of colData(dds)")
  }
  
  poisd <- PoissonDistance(t(deseqcounts))
  samplePoisDistMatrix <- as.matrix(poisd$dd)
  mdspois <- data.frame(cmdscale(samplePoisDistMatrix))
  mdspois <- cbind(mdspois, as.data.frame(colData(ddstemp)))
  mdspois$name <- colnames(deseqcounts)
  
  mdspoisplot <- ggplot(mdspois, aes(X1, X2))
  
  mdspoisplot <- addStandards(mdspoisplot, title)
  
  if (length(conditions) > 1) {
    mdspoisplot <- addConditions(mdspoisplot, conditions, length(unique(mdspois[, conditions[2]])))  
  } else {
    mdspoisplot <- addConditions(mdspoisplot, conditions)  
  }
  
  mdspoisplot <- addSampleLabel(mdspoisplot, samplename)
  
  ggsave(filename, mdspoisplot)
}

# Function plots the heatmap for the samples using the pheatmap and RColorBrewer packages
# It follows the proposed heatmap from the rnaseqGene workflow
# http://www.bioconductor.org/help/workflows/rnaseqGene/#eda
# Input
# readtransform ... DESeqTransform which comes from the transformation process (vsd, rld, ...)
# conditions ... is a vector and entries should be column names from the colData(dds) or coldata
#               columns are the conditions which will be plotted in the heatmap
# filename ... if provided plot is printed directly
# title ... title for heatmap
plotHeatmap <- function(readtransform, conditions = "", cluster = FALSE, filename = "", title = "Heatmap of Euclidean distance", correctedvsd = NULL) {
  require(pheatmap)
  require(RColorBrewer)
  #get nice colours
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  if (is.matrix(correctedvsd)) {
    sampleDists <- dist(t(correctedvsd))
    sampleDistMatrix <- as.matrix(sampleDists)    
  } else {
    sampleDists <- dist(t(assay(readtransform)))
    sampleDistMatrix <- as.matrix(sampleDists) 
  }
  
  if (!is.data.frame(conditions)) {
    conditions <- NA
  } else {
    if (!all(names(conditions) %in% names(colData(readtransform)))) {
      stop("the argument 'conditions' should specify columns of colData(dds)")
    }
  }
  
  fontvec <- c(9,9,9)
  cellvec <- c(20,20)
  
  clust <- ifelse(cluster, TRUE, FALSE)
  cols <- ncol(readtransform)
  if (filename == '') {
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
             cluster_rows = clust, cluster_cols = clust,
             color=colors, annotation_col = conditions,
             main = title, border_color = NA,
             fontsize = fontvec[1],
             fontsize_row = fontvec[2],
             fontsize_col = fontvec[3],
             cellwidth = cellvec[1],
             cellheight = cellvec[2]
    )   
  } else {
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
             cluster_rows = clust, cluster_cols = clust,
             color=colors, annotation_col = conditions,
             main = title, border_color = NA,
             fontsize = fontvec[1],
             fontsize_row = fontvec[2],
             fontsize_col = fontvec[3],
             cellwidth = cellvec[1],
             cellheight = cellvec[2],
             filename = filename
    )
  }
}

# Function plots the heatmap for the samples using the pheatmap and RColorBrewer packages
# It follows the proposed heatmap from the rnaseqGene workflow
# http://www.bioconductor.org/help/workflows/rnaseqGene/#eda
#
# Another option for calculating sample distances is to use the Poisson Distance [@Witten2011Classification], 
# implemented in the PoiClaClu package. This measure of dissimilarity between counts also takes the inherent 
# variance structure of counts into consideration when calculating the distances between samples. 
# The PoissonDistance function takes the original count matrix (not normalized) with samples as rows instead
# of columns, so we need to transpose the counts in dds.
# Input
# deseqcounts ... DESeq2 countstable not normalised, counts(dds)
# conditions ... data.frame; rownames must be the same rownames as readtransform; 
#               columns are the conditions which will be plotted in the heatmap
# filename ... if provided plot is printed directly
# title ... title for heatmap
plotHeatmapPoiClaClu <- function(deseqcounts, conditions, cluster = FALSE, filename = "", title = "Heatmap of Poisson Distance") {
  require(PoiClaClu)
  require(RColorBrewer)
  require(DESeq2)
  require(pheatmap)
  #get nice colours
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  poisd <- PoissonDistance(t(deseqcounts))
  samplePoisDistMatrix <- as.matrix(poisd$dd)
  rownames(samplePoisDistMatrix) <- colnames(deseqcounts)
  colnames(samplePoisDistMatrix) <- colnames(deseqcounts)
  
  if (!is.data.frame(conditions)) {
    conditions <- NA
  } 
  
  fontvec <- c(9,9,9)
  cellvec <- c(20,20)
  
  clust <- ifelse(cluster, TRUE, FALSE)
  cols <- ncol(samplePoisDistMatrix)

  if (filename == '') {
    pheatmap(samplePoisDistMatrix,
             clustering_distance_rows=poisd$dd, clustering_distance_cols=poisd$dd,
             cluster_rows = clust, cluster_cols = clust,
             color=colors, annotation_col = conditions,
             main = title, border_color = NA,
             fontsize = fontvec[1],
             fontsize_row = fontvec[2],
             fontsize_col = fontvec[3],
             cellwidth = cellvec[1],
             cellheight = cellvec[2],
             display_numbers = F
    )
  } else {
    pheatmap(samplePoisDistMatrix,
             clustering_distance_rows=poisd$dd, clustering_distance_cols=poisd$dd,
             cluster_rows = clust, cluster_cols = clust,
             color=colors, annotation_col = conditions,
             main = title, border_color = NA,
             fontsize = fontvec[1],
             fontsize_row = fontvec[2],
             fontsize_col = fontvec[3],
             cellwidth = cellvec[1],
             cellheight = cellvec[2],
             display_numbers = F,
             filename = filename
    )
  }
}

# Function plots the spearman or pearson correlation using the pheatmap and RColorBrewer packages
# Input
# deseqmatrix ... counts table
# conditions ... data.frame; rownames must be the same rownames as deseqmatrix; 
#               columns are the conditions which will be plotted in the heatmap
# cormethod ... pearson or spearman
# fontsize ... size of the numbers in the cells
# cluster ... if provided rows and column are cluster according to values of correlation calculation
# filename ... if provided plot is printed directly
# titleappend ... add to Pearson Correlation an appendix like 'for fubar samples'
plotCorrelation <- function(deseqmatrix, conditions = "", cormethod = "pearson", cluster = FALSE, filename = "") {
  require(pheatmap)
  require(RColorBrewer)
  
  cormatrix <- cor(deseqmatrix, use = "pairwise.complete.obs", method = cormethod)
  cornumbers <- round(cormatrix, 3)
  # colnames(cormatrix) <- NULL
  if (cormethod == "pearson") {
    colors <- colorRampPalette(brewer.pal(8, "Blues"))(255)
    title <- "Pearson Correlation"  
  }
  if (cormethod == "spearman") {
    colors <- colorRampPalette(brewer.pal(8, "Reds"))(255)
    title <- "Spearman Correlation"
  }
  
  if (!is.data.frame(conditions)) {
    conditions <- NA
  }

  fontvec <- c(9,9,9)
  cellvec <- c(20,20)
  
  clust <- ifelse(cluster, TRUE, FALSE)
  cols <- ncol(cormatrix)
  if (cols < 19) {
    if (filename == '') {
      pheatmap(cormatrix,
               cluster_rows = clust, cluster_cols = clust,
               color=colors, annotation_col = conditions,
               main = title, border_color = NA,
               fontsize = fontvec[1],
               fontsize_number = fontvec[1] - 2,
               fontsize_row = fontvec[2],
               fontsize_col = fontvec[3],
               cellwidth = cellvec[1],
               cellheight = cellvec[2],
               display_numbers = cornumbers, number_color = "Black"
      )
    } else {
      pheatmap(cormatrix,
               cluster_rows = clust, cluster_cols = clust,
               color=colors, annotation_col = conditions,
               main = title, border_color = NA,
               fontsize = fontvec[1],
               fontsize_number = fontvec[1] - 2,
               fontsize_row = fontvec[2],
               fontsize_col = fontvec[3],
               cellwidth = cellvec[1],
               cellheight = cellvec[2],
               display_numbers = cornumbers, number_color = "Black",
               filename = filename
      )        
    }
  } else {
    if (filename == '') {
      pheatmap(cormatrix,
               cluster_rows = clust, cluster_cols = clust,
               color=colors, annotation_col = conditions,
               main = title, border_color = NA,
               fontsize = fontvec[1],
               fontsize_row = fontvec[2],
               fontsize_col = fontvec[3],
               cellwidth = cellvec[1],
               cellheight = cellvec[2],
               display_numbers = F
      )
    } else {
      pheatmap(cormatrix,
               cluster_rows = clust, cluster_cols = clust,
               color=colors, annotation_col = conditions,
               main = title, border_color = NA,
               fontsize = fontvec[1],
               fontsize_row = fontvec[2],
               fontsize_col = fontvec[3],
               cellwidth = cellvec[1],
               cellheight = cellvec[2],
               display_numbers = F,
               filename = filename
      )
    }
  }
}

# The heatmap becomes more interesting if we do not look at absolute expression strength
# but rather at the amount by which each gene deviates in a specific sample from the gene’s
# average across all samples. Hence, we center each genes’ values across samples, and plot
# a heatmap. We provide a data.frame that instructs the pheatmap function how to label
# the columns.
# Input:
# readtransform ... DESeqTransform which comes from the transformation process (vsd, rld, ...)
# genelist ... data.frame which has two columns (Ensembl_ID and Gene_Symbol [from the anno function])
# conditions ... is a data.frame and entries should be column names from the colData(dds)
#               columns are the conditions which will be plotted in the heatmap
# nog ... number of genes to be displayed
# clusterCol ... boolean if Columns should be clustered
# filename ... if provided plot is printed directly
# title ... title for heatmap
plotVarGenes <- function(readtransform, genelist = "", conditions = "", nog = 20, cluster = TRUE, filename = "", title = paste0("Deviation from the gene's average across all samples\nfor top ", nog, " of genes with the highest variance"), correctedvsd = NULL) {
  require(genefilter)
  require(pheatmap)
  require(DESeq2)
  
  if (is.matrix(correctedvsd)) {
    topVarGenes <- head(order(rowVars(correctedvsd),decreasing=TRUE), nog)
    mat <- correctedvsd[ topVarGenes, ]
    mat <- mat - rowMeans(mat)
  } else{
    topVarGenes <- head(order(rowVars(assay(readtransform)),decreasing=TRUE), nog)
    mat <- assay(readtransform)[ topVarGenes, ]
    mat <- mat - rowMeans(mat)
  }
  
  if (is.data.frame(genelist)) {
    genelist$Gene_Symbol <- make.unique(as.character(genelist$Gene_Symbol))
    mat <- merge(mat, genelist, by.x = "row.names", by.y = "Ensembl_ID")
    mat$Gene_Symbol <- ifelse(mat$Gene_Symbol == '', as.character(mat$Row.names), as.character(mat$Gene_Symbol))
    rownames(mat) <- mat$Gene_Symbol
    mat$Row.names <- NULL
    mat$Gene_Symbol <- NULL
  }
  
  if (!is.data.frame(conditions)){
    conditions <- NULL
  } else {
    if (!all(names(conditions) %in% names(colData(readtransform)))) {
      stop("the argument 'conditions' should specify columns of colData(dds)")
    }
  }

  clust <- ifelse(cluster, TRUE, FALSE)
  
  no <- nrow(mat)
  if (no <= 25) {
    rowsize <- 6
    h <- 5
    w <- 8
  } else if (no > 25 & no <100) {
    rowsize <- 6
    h <- 7
    w <- 8
  } else if (no >= 100 & no < 200) {
    rowsize <- 6
    h <- 11
    w <- 8
  } else if (nrow(mat)>=200) {
    rowsize <- 6
    h <- 19
    w <- 8
  }

  if (filename == ""){
    pheatmap(mat,
             annotation_col = conditions,
             cluster_cols = clust,
             main = title,
             fontsize = 8,
             fontsize_row = rowsize, # previous: 8
             cellheight = rowsize,
             fontsize_col = 8,
             treeheight_row = 30,
             treeheight_col = 30,
             height = h
    )
  } else {
    pheatmap(mat,
             annotation_col = conditions,
             cluster_cols = clust,
             main = title,
             fontsize = 8,
             fontsize_row = rowsize, # previous: 8
             cellheight = rowsize,
             fontsize_col = 8,
             filename = filename,
             treeheight_row = 30,
             treeheight_col = 30,
             height = h
    )
  }
}

# function plots an individual list and shows gene deviation, if it is set
# input:
# readtransform ... counttable or readtransformed table (assay(vsd), assay(rld)); only your genes of interest!!!
#                   in any case rownames should be gene ids
# genelist ... data.frame which has two columns (Ensembl_ID and Gene_Symbol [from the anno function])
#              if not provided, not taken
# conditions ... is a data.frame and entries should be column names from the colData(dds)
#               columns are the conditions which will be plotted in the heatmap
# rowMean ... boolean, either calculate rowMean or not
# clusterCol ... boolean if Columns should be clustered
# clusterRow ...boolean if Rows should be clustered
# colnames ... boolean if colnames are displayed
# rownames ... boolean if rownames are displayed
# title ... show title
# filename ... if provided plot is printed directly
plotVarGenesIndividual <- function(readtransform, genelist = "", conditions = "", rowMean = TRUE, clusterCol = TRUE, clusterRow = TRUE, colnames = TRUE, rownames = FALSE, title = TRUE, filename = "") {
  require(genefilter)
  require(pheatmap)
  # require(RColorBrewer)
  
  if (rowMean){
    mat <- readtransform - rowMeans(readtransform)
    if (title) {
      title = paste0("Gene deviation from the average across all samples for ", nrow(readtransform), " genes")      
    } else {
      title = ""
    }
  } else {
    mat <- readtransform
    if (title) {
      title = paste0("Gene expression of ", nrow(readtransform), " genes") 
    } else {
      title = ""
    }
  }
  
  # colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
  if (rownames) {
    if (is.data.frame(genelist)) {
      genelist$Gene_Symbol <- make.unique(as.character(genelist$Gene_Symbol))
      mat <- merge(mat, genelist, by.x = "row.names", by.y = "Ensembl_ID")
      mat$Gene_Symbol <- ifelse(mat$Gene_Symbol == '', as.character(mat$Row.names), as.character(mat$Gene_Symbol))
      rownames(mat) <- mat$Gene_Symbol
      mat$Row.names <- NULL
      mat$Gene_Symbol <- NULL
    }
  }
  
  if (!is.data.frame(conditions)){
    conditions <- NULL
  }
  
  no <- nrow(mat)
  
  if (no <= 25) {
    rowsize <- 6
    h <- 5
    w <- 8
  } else if (no > 25 & no <100) {
    rowsize <- 6
    h <- 7
    w <- 8
  } else if (no >= 100 & no < 200) {
    rowsize <- 6
    h <- 11
    w <- 8
  } else if (nrow(mat)>=200) {
    rowsize <- 6
    h <- 19
    w <- 8
  }
  
  if (filename == ""){
    pheatmap(mat,
             annotation_col = conditions,
             # color = colors,
             cluster_cols = clusterCol,
             cluster_rows = clusterRow,
             show_rownames = rownames,
             show_colnames = colnames,
             fontsize = 8,
             fontsize_row = rowsize,
             cellheight = rowsize,
             height = h,
             fontsize_col = 8,
             main = title
    )
  } else {
    pheatmap(mat,
             annotation_col = conditions,
             # color = colors,
             cluster_cols = clusterCol,
             cluster_rows = clusterRow,
             show_rownames = rownames,
             show_colnames = colnames,
             fontsize = 8,
             fontsize_row = rowsize,
             cellheight = rowsize,
             height = h,
             fontsize_col = 8,
             main = title,
             filename = filename
    )
  }
}

# The heatmap becomes more interesting if we do not look at absolute expression strength
# but rather at the amount by which each gene deviates in a specific sample from the gene’s
# average across all samples. Hence, we center each genes’ values across samples, and plot
# a heatmap. We provide a data.frame that instructs the pheatmap function how to label
# the columns.
# Input:
# readtransform ... DESeqTransform which comes from the transformation process (vsd, rld, ...)
# genelist ... data.frame which has two columns (Ensembl_ID and Gene_Symbol [from the anno function])
# conditions ... is a data.frame and entries should be column names from the colData(dds)
#               columns are the conditions which will be plotted in the heatmap
# nog ... number of genes to be displayed
# clusterCol ... boolean if Columns should be clustered
# filename ... if provided plot is printed directly
# title ... title for heatmap
plotHighExpressedGenes <- function(readtransform, genelist = "", conditions = "", nog = 20, clusterCol = TRUE, filename = "", title = paste0("Top ", nog, " of genes with highest mean expression across samples")) {
  require(genefilter)
  require(pheatmap)
  require(DESeq2)
  topHighGenes <- head(order(rowMeans(assay(readtransform)),decreasing=TRUE), nog)
  mat <- assay(readtransform)[ topHighGenes, ]
  
  if (is.data.frame(genelist)) {
    genelist$Gene_Symbol <- make.unique(as.character(genelist$Gene_Symbol))
    mat <- merge(mat, genelist, by.x = "row.names", by.y = "Ensembl_ID")
    mat$Gene_Symbol <- ifelse(mat$Gene_Symbol == '', as.character(mat$Row.names), as.character(mat$Gene_Symbol))
    rownames(mat) <- mat$Gene_Symbol
    mat$Row.names <- NULL
    mat$Gene_Symbol <- NULL
  }
  
  if (!is.data.frame(conditions)){
    conditions <- NULL
  } else {
    if (!all(names(conditions) %in% names(colData(readtransform)))) {
      stop("the argument 'conditions' should specify columns of colData(dds)")
    }
  }
  
  no <- nrow(mat)
  
  if (no <= 25) {
    rowsize <- 6
    h <- 5
    w <- 8
  } else if (no > 25 & no <100) {
    rowsize <- 6
    h <- 7
    w <- 8
  } else if (no >= 100 & no < 200) {
    rowsize <- 6
    h <- 11
    w <- 8
  } else if (nrow(mat)>=200) {
    rowsize <- 6
    h <- 19
    w <- 8
  }
  
  if (filename == ""){
    pheatmap(mat,
             annotation_col = conditions,
             cluster_cols = clusterCol,
             main = title,
             fontsize = 8,
             fontsize_row = rowsize, # previous: 8
             fontsize_col = 8,
             height = h,
             cellheight = rowsize
    )
  } else {  
    pheatmap(mat,
             annotation_col = conditions,
             cluster_cols = clusterCol,
             main = title,
             fontsize = 8,
             fontsize_row = rowsize, # previous: 8
             fontsize_col = 8,
             height = h,
             cellheight = rowsize,
             filename = filename
    )
  }
}

#plot the differentiell expressed genes
# input
# deseqresult ... DESeqResults object
# titleprefix ... start of the title like "KO vs WT (sample 2 removed)
# padj_value ... padj which is used for getting the DESeqResults object 
# ylow, ymax: min and max coordinates on y axis
# filename: save name
MAplot <- function(deseqresult, titleprefix, padj_value = 0.1, ylow = -4, ymax = 4, filename = "") {
  require(DESeq2)
  genes <- paste0(prettyNum(sum(deseqresult$padj < padj_value, na.rm=T), big.mark=",", decimal.mark="."), ' DEGs')
  fdr <- paste0("FDR ", padj_value*100, "%")
  title <- paste0(titleprefix, '\n with ', fdr, ' and ', genes)
  if(filename == "") {
    DESeq2::plotMA(deseqresult, ylim=c(ylow, ymax), main=title, alpha=padj_value, ylab=expression(paste(log["2"], " fold change")))
  } else {
    if (grepl('pdf$', filename)) {
      pdf(filename)
      DESeq2::plotMA(deseqresult, ylim=c(ylow, ymax), main=title, alpha=padj_value, ylab=expression(paste(log["2"], " fold change")))
      dev.off()
    } else if(grepl('png$', filename))  {
      png(filename)
      DESeq2::plotMA(deseqresult, ylim=c(ylow, ymax), main=title, alpha=padj_value, ylab=expression(paste(log["2"], " fold change")))
      dev.off()
    }
  }
}

doMAplotting <- function(detemp, runvectemp, bfxid, pngtemp, pdftemp, degene, padjstr, padjv) {
  namestring <- paste0(bfxid, '.maplot.', paste(runvectemp[6], collapse = '_'), '-', runvectemp[1], '-', runvectemp[4], '_vs_', runvectemp[5], '.', padjstr, '-', degene)
  filenamepng <- createfilename(pngtemp, namestring, runvectemp[7], 'png')
  filenamepdf <- createfilename(pdftemp, namestring, runvectemp[7], 'pdf')
  MAplot(detemp, paste0(runvectemp[1], ': ', runvectemp[4], ' vs ', runvectemp[5]), padjv, min(detemp$log2FoldChange, na.rm = T)-0.5, max(detemp$log2FoldChange, na.rm = T)+0.5, filenamepng)
  MAplot(detemp, paste0(runvectemp[1], ': ', runvectemp[4], ' vs ', runvectemp[5]), padjv, min(detemp$log2FoldChange, na.rm = T)-0.5, max(detemp$log2FoldChange, na.rm = T)+0.5, filenamepdf)
  return(filenamepng)
}

doMAGlobalplotting <- function(detemp, decomptemp, bfxid, pngtemp, pdftemp) {
  detemp$padj <- detemp$padjGlob; detemp$padjGlob <- NULL
  deformulavector <- getFormulavector(as.character(decomptemp[, 'Formula']))
  fact <- as.character(decomptemp[, 'Factor'])
  comp <- gsub(' ', '_', as.character(decomptemp[, 'Comparison']))
  gene <- decomptemp[, 'DEGenes']
  padjv <- decomptemp[, 'FDR']/100
  padjstr <- as.character(decomptemp[, 'FDR'])
  padjstr <- ifelse(nchar(padjstr) == 1, paste0('00', padjstr), paste0('0', padjstr))
  
  namestring <- paste0(bfxid, '.maplot.', paste(deformulavector, collapse = '_'), '-', fact, '-', comp, '.', padjstr, '-', gene)
  filenamepng <- createfilename(pngtemp, namestring, decomptemp[, 'App'], 'png')
  filenamepdf <- createfilename(pdftemp, namestring, decomptemp[, 'App'], 'pdf')
  MAplot(detemp, paste0(fact, ': ', as.character(decomptemp[, 'Comparison'])), padjv, min(detemp$log2FoldChange, na.rm = T)-0.5, max(detemp$log2FoldChange, na.rm = T)+0.5, filenamepng)
  MAplot(detemp, paste0(fact, ': ', as.character(decomptemp[, 'Comparison'])), padjv, min(detemp$log2FoldChange, na.rm = T)-0.5, max(detemp$log2FoldChange, na.rm = T)+0.5, filenamepdf)
  return(filenamepng)
}

doVolcanoPlotting <- function(detemp, runvectemp, bfxid, pngtemp, pdftemp, padjstr, padj_value, geneensembltab, fc_up, fc_down, showgenes){
  # filename preparation
  sign_up <- nrow(subset(detemp, padj < padj_value & log2FoldChange > fc_up))
  sign_down <- nrow(subset(detemp, padj < padj_value & log2FoldChange < fc_down))
  namestring <- paste0(bfxid, '.volcano_plot.', paste(runvectemp[6], collapse = '_'), '-', runvectemp[1], '-', runvectemp[4], '_vs_', runvectemp[5], '.', padjstr, '-', sign_up + sign_down)
  filenamepng <- createfilename(pngtemp, namestring, runvectemp[7], 'png')
  filenamepdf <- createfilename(pdftemp, namestring, runvectemp[7], 'pdf')
  fdr <- paste0("FDR ", padj_value*100, "%")
  title <-  paste0(runvectemp[1], ': ', runvectemp[4], ' vs ', runvectemp[5], '\nwith ', fdr, ' and ', prettyNum(sign_up + sign_down, big.mark=",", decimal.mark="."), ' DEGs')
  VolcanoPlot(detemp, geneensembltab, title, padj_value, fc_up, fc_down, showgenes, filenamepdf)
  VolcanoPlot(detemp, geneensembltab, title, padj_value, fc_up, fc_down, showgenes, filenamepng)
  return(filenamepng)
}

doVolcanoGlobalPlotting <- function(detemp, decomptemp, bfxid, pngtemp, pdftemp, geneensembltab, fc_up, fc_down, showgenes) {
  detemp$padj <- detemp$padjGlob; detemp$padjGlob <- NULL
  deformulavector <- getFormulavector(as.character(decomptemp[, 'Formula']))
  fact <- as.character(decomptemp[, 'Factor'])
  comp <- gsub(' ', '_', as.character(decomptemp[, 'Comparison']))
  gene <- decomptemp[, 'DEGenes']
  padjv <- decomptemp[, 'FDR']/100
  padjstr <- as.character(decomptemp[, 'FDR'])
  padjstr <- ifelse(nchar(padjstr) == 1, paste0('00', padjstr), paste0('0', padjstr))
  sign_up <- nrow(subset(detemp, padj < padjv & log2FoldChange > fc_up))
  sign_down <- nrow(subset(detemp, padj < padjv & log2FoldChange < fc_down))

  # filename preparation
  namestring <- paste0(bfxid, '.volcano_plot.', paste(deformulavector, collapse = '_'), '-', fact, '-', comp, '.', padjstr, '-', sign_up + sign_down)
  filenamepng <- createfilename(pngtemp, namestring, decomptemp[, 'App'], 'png')
  filenamepdf <- createfilename(pdftemp, namestring, decomptemp[, 'App'], 'pdf')
  fdr <- paste0("FDR ", padjv*100, "%")
  title <-  paste0(fact, ': ', as.character(decomptemp[, 'Comparison']), '\nwith ', fdr, ' and ', prettyNum(sign_up + sign_down, big.mark=",", decimal.mark="."), ' DEGs')
  VolcanoPlot(detemp, geneensembltab, title, padjv, fc_up, fc_down, showgenes, filenamepdf)
  VolcanoPlot(detemp, geneensembltab, title, padjv, fc_up, fc_down, showgenes, filenamepng)
  return(filenamepng)
}

# Input for Volcano Plot
# detemp ... DESeqResults S4 Class
# runvectemp ... custom vector for DeSeq analysis
# bfxid ... String
# pngtemp, pdftemp ... Directories
# padjstr ... String
# padj_value ... Float
# fc_up, fc_down ... Float
# showgenes ... integer; show gene names of x up and x down genes
# "
VolcanoPlot <- function(detemp, annodf, title, padj_value, fc_up, fc_down, showgenes, filename = '') {
  require(ggplot2)
  require(magrittr)
  require(ggrepel)
  theme_set(theme_bw(14))
  # modify the deseq result set
  detemp <- data.frame(detemp[!is.na(detemp$pvalue), ])
  detemp <- merge(detemp, annodf, by.x = 'row.names', by.y = 'Ensembl_ID')
  rownames(detemp) <- detemp$Row.names; detemp$Row.names <- NULL
  
  detemp$sign <- ifelse((detemp$padj < padj_value & detemp$log2FoldChange > fc_up), 'Up', 
                        ifelse((detemp$padj < padj_value & detemp$log2FoldChange < fc_down), 'Down', 'Not Sig'))
  sign_up <- sum(detemp$sign == 'Up'); sign_down <- sum(detemp$sign == 'Down'); not_sig <- sum(detemp$sign == 'Not Sig')
  
  factor_vector <- c(); colour_vector <- c(); gene_ids <- c()
  if (sign_up != 0) {
    factor_vector <- c(factor_vector, 'Up'); colour_vector <- c(colour_vector, 'darkgoldenrod2')# '#F8766D'
    deup <- subset(detemp, sign == 'Up')
    deup <- deup[order(deup$log2FoldChange, decreasing = T), ]
    if (nrow(deup) < showgenes) {
      gene_ids <- c(gene_ids, rownames(deup[1:nrow(deup), ]))
    } else{
      gene_ids <- c(gene_ids, rownames(deup[1:showgenes, ]))
    }
  }
  if (not_sig != 0) {
    factor_vector <- c(factor_vector, 'Not Sig'); colour_vector <- c(colour_vector, 'grey70')
    denotsign <- subset(detemp, sign == 'Not Sig')
  }
  if (sign_down != 0) {
    factor_vector <- c(factor_vector, 'Down'); colour_vector <- c(colour_vector,'steelblue1')# '#619CFF'
    dedown <- subset(detemp, sign == 'Down')
    dedown <- dedown[order(dedown$log2FoldChange), ]
    if (nrow(dedown) < showgenes) {
      gene_ids <- c(gene_ids, rownames(dedown[1:nrow(dedown), ]))
    } else{
      gene_ids <- c(gene_ids, rownames(dedown[1:showgenes, ]))
    }
  }
  detemp <- rbind(denotsign, deup, dedown)
  detemp$sign <- factor(detemp$sign, levels = factor_vector)
  
  ggobject <- ggplot(detemp , aes(log2FoldChange, -log10(pvalue), colour = sign)) + geom_point(size = 1) +
    theme(legend.title= element_blank(), legend.position = 'none') + 
    xlab('log2 Fold Change') + 
    ylab('-log10 p-value') +
    scale_colour_manual(values = colour_vector) + 
    ggtitle(title) +
    geom_text_repel(data=detemp[gene_ids,], aes(label=Gene_Symbol), color = 'grey30', force = 1) #'#619CFF')
  if (sign_down != 0){
    ggobject <- ggobject + annotate('label', x = -Inf, y = 0, label = paste0('log2 FC<', fc_down, ': ',  prettyNum(sign_down, big.mark=",", decimal.mark=".")), col= 'steelblue1', vjust = 0.5, hjust = -0.3, alpha = 0.5)
  }
  if (sign_up != 0){
    ggobject <- ggobject + annotate('label', x = Inf, y = 0, label = paste0('log2 FC>', fc_up, ': ',  prettyNum(sign_up, big.mark=",", decimal.mark=".")), col= 'darkgoldenrod2', vjust = 0.5, hjust = 1.4, alpha = 0.5)
  }

  if(filename == "") {
    ggobject
  } else {
    if (grepl('pdf$', filename)) {
      ggsave(filename, ggobject)
    } else if(grepl('png$', filename))  {
      ggsave(filename, ggobject)
    }
  }
}
