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


library(magrittr)
library(parallel)
library(BiocParallel)
library(data.table)
library(openxlsx)
library(DESeq2)

# Allows to use script on cluster without X11
options(bitmapType='cairo')

# Ask ggrepel to plot all labels
options(ggrepel.max.overlaps=Inf)

num_cores <- min(future::availableCores(), 12)
register(MulticoreParam(num_cores)) 

args <- commandArgs(TRUE)

configfilename <- args[1]
configfile <- read.table(configfilename, header = T, sep = '\t', quote = "", as.is = T)
workingdir <- subset(configfile, key == 'workingdir')$value
ensembldir <- subset(configfile, key == 'ensembldir')$value
logodir <- subset(configfile, key == 'logodir')$value
file.genes_remove <- subset(configfile, key == 'exclude_genes')$value

setwd(workingdir)

# source standard functions
modulefolder <- subset(configfile, key == 'modules')$value
source(file.path(modulefolder, 'modules-deseq-report.r'))
source(file.path(modulefolder, 'modules-ensembl.r'))
source(file.path(modulefolder, 'modules-plots.r'))
source(file.path(modulefolder, 'modules-html.r'))
vsdbatchwritten <- FALSE

appendix <- subset(configfile, key == 'appendix')$value #added to filenames
appendix <- ifelse(nchar(appendix)==0, 'de-expl', appendix)
bfxid <- subset(configfile, key == 'bfxid')$value #bfxid; for filenames
# create directories and set rights
pdfdir <- file.path(workingdir, 'pdf');mkdirectory(pdfdir)
pngdir <- file.path(workingdir, 'picture');mkdirectory(pngdir)
datadir <- file.path(workingdir, 'data');mkdirectory(datadir)
gseadir <- file.path(workingdir, 'data', 'gsea');mkdirectory(gseadir)

whereiswhat <- data.frame(what = character(), xls = character(), pdf = character(), png = character(), stringsAsFactors = FALSE)

#########
# Step 1: import the conditions
#########
print('READING CONDITION FILE')
conditionfilename <- subset(configfile, key == 'conditionfile')$value
coldata <- read.table(conditionfilename, header=T, sep='\t', check.names = F)
coldata <- changeContinuous(coldata)

#########
# STEP 2: read the gene count table. script will end here if a wrong filepath is given
#########
print('READING GENE COUNT TABLE')
countfilename <- subset(configfile, key == 'countfile')$value
counttab <- read.table(countfilename, header=T, sep='\t', skip=1, check.names = F, row.names = 1) #skip first row because it shows the command call of featureCounts
genelength <- counttab[!grepl('ERCC', rownames(counttab)), 'Length', drop = F] # saving the gene length for some tables
counttab <- counttab[, !(colnames(counttab) %in% c('Chr', 'Start', 'End', 'Strand', 'End', 'Length')), drop = F]
counttab <- counttab[!grepl('^ERCC', rownames(counttab)),] # remove ERCC
# this tryCatch checks if the same samples are in coldata$Sample and colnames(counttab)
tryCatch(counttab <- counttab[, as.character(coldata$Sample), drop = F],
         error = function(e){
           print(e)
           print('Error: Some sample names from the condition file are not sample names in the counts table')
           print('COLDATA sample names:'); print(coldata$Sample)
           print('COUNTTAB sample names:'); print(colnames(counttab))
           stop()
         }
)
allgenes <- rownames(counttab)
counttab <- removeLowCounts(counttab, 10)
rownames(coldata) <- as.vector(colnames(counttab))

print(rownames(coldata) == coldata$Sample)
if (sum(rownames(coldata) == coldata$Sample) != ncol(counttab)){
  warning('If you see FALSE in the values above it means samples in coldata and counttab are in the wrong order')
  stop('Check your code an rerun program again.')
}


#########
# STEP 3: retrieve the annotation, either from ensembl or local storage
#########
species <- subset(configfile, key == 'species')$value
ensembl <- subset(configfile, key == 'ensembl')$value
annotable <- buildAnnotation(ensembldir, species, ensembl, allgenes, genelength)
annotable <- subset(annotable, Ensembl_ID %in% rownames(counttab)) # only genes seen in experiment
annotable$Gene_Symbol <- ifelse(annotable$Gene_Symbol == '', as.character(annotable$Ensembl_ID), as.character(annotable$Gene_Symbol))
geneensembltab <- annotable[, c('Ensembl_ID', 'Gene_Symbol')] # table for the top var genes and highest expressed genes

#########
# Exploratory
#########
exploratory <- subset(configfile, key == 'exploratory')$value

if (nchar(exploratory) != 0) {
  exploratory <- gsub(' ', '', exploratory) %>% strsplit(',') %>% unlist()
  if (sum(exploratory %in% colnames(coldata)) != length(exploratory)) {
    print(paste0('Exploratory variable: ', paste(exploratory, collapse = ',')))
    print(paste0('Factors from coldata: ', paste(colnames(coldata), collapse = ',')))
    warning('Some entries in the exploratory variable are not listed in the condition file. Skipping Exploratory!')
    Sys.sleep(3)
  } else if (length(exploratory) > 3) {
    warning('Only less than 4 exploratory factors are accepted. Skipping Exploratory!')
  } else {
    dds <- DESeqDataSetFromMatrix(countData = counttab, colData = coldata, design = ~ 1)
    dds <- DESeq(dds, parallel=T, betaPrior = TRUE)
    
    coldatanames <- data.frame(colData(dds))[, exploratory, drop = F]
    trans.method <- ''
    if (ncol(counttab) < 30){
      vsd <- rlogTransformation(dds, blind = T)
      trans.method <- 'rlog'
    } else {
      if (exists('vst')){
        vsd <- vst(dds, blind = T) # DESeq2 1.12+ 
        trans.method <- 'vst'
      } else {
        vsd <- varianceStabilizingTransformation(dds, blind = T)
        trans.method <- 'vst'
      }  
    }
    
    #########
    # workbook for counts
    #########
    deseq.wb <- createWorkbook()
    namestring <- paste(bfxid, 'various_countTables', sep = '.')
    deseq.wb.filename <- createfilename(datadir, namestring, appendix, 'xlsx')
    
    #########
    # calculate the raw counts
    ######### 
    rawcounts <- merge(annotable, counts(dds, normalized = F), by.x = "Ensembl_ID", by.y = "row.names", all.y = T)
    rawcounts <- rawcounts[, c('Ensembl_ID', 'Gene_Symbol', sort(rownames(coldata)), 
                               c('Chromosome', 'Start', 'End', 'Length', 'Strand', 'Biotype', 'Description'))]
    # write data to workbook
    addWorksheet(deseq.wb, 'raw-counts')
    writeData(deseq.wb, 'raw-counts', rawcounts[, c(setdiff(colnames(rawcounts), colnames(counttab)), colnames(counttab))])
    rawcounts <- NULL
    
    #########
    # calculate the normalised counts
    #########
    normcountsrounded <- round(counts(dds, normalized = T), 1)
    normcountsrounded <- merge(annotable, normcountsrounded, by.x = "Ensembl_ID", by.y = "row.names", all.y = T)
    normcountsrounded <- normcountsrounded[, c('Ensembl_ID', 'Gene_Symbol', sort(rownames(coldata)), 
                                               c('Chromosome', 'Start', 'End', 'Length', 'Strand', 'Biotype', 'Description'))]
    
    # write data to workbook
    addWorksheet(deseq.wb, 'normalised-counts')
    writeData(deseq.wb, 'normalised-counts', normcountsrounded[, c(setdiff(colnames(normcountsrounded), colnames(counttab)), colnames(counttab))])
    normcountsrounded <- NULL
    
    #########
    # calculate the variance stabilizing transformation
    #########    
    vsd.df <- as.data.frame(assay(vsd))
    zscore.df <- vsd.df
    vsd.df$Ensembl_ID <- rownames(vsd.df)
    vsd.df <- merge(vsd.df, annotable, by = 'Ensembl_ID')
    # write data to workbook
    
    addWorksheet(deseq.wb, paste0(trans.method,'-values'))
    writeData(deseq.wb, paste0(trans.method,'-values'), vsd.df[, c(setdiff(colnames(vsd.df), colnames(counttab)), colnames(counttab))])
    
    #########
    # calculate the Z-score
    #########
    zscore.df <- build_zscore(zscore.df)
    zscore.df$Ensembl_ID <- rownames(zscore.df)
    zscore.df <- merge(zscore.df, annotable, by = 'Ensembl_ID')
    # write data to workbook
    addWorksheet(deseq.wb, 'zscore-values')
    writeData(deseq.wb, 'zscore-values', zscore.df[, c(setdiff(colnames(zscore.df), colnames(counttab)), colnames(counttab))])
    
    # dispersion
    filenamepdf <- createfilename(pdfdir, paste(bfxid, 'dispEst', sep = '.'), appendix, 'pdf')
    filenamepng <- createfilename(pngdir, paste(bfxid, 'dispEst', sep = '.'), appendix, 'png')
    pdf(filenamepdf);plotDispEsts(dds);dev.off()
    png(filenamepng);plotDispEsts(dds);dev.off()
    whereiswhat <- rbind(whereiswhat, data.frame(what = 'dispersion', xls = '', pdf = filenamepdf, png = filenamepng, stringsAsFactors = FALSE))

    # create vsd and countstable for exploratory analysis
    if ((file.exists(file.genes_remove)) & (file.info(file.genes_remove)$size != 0)) {
      genes_remove <- read.table(file.genes_remove)[, 1]
      counttab.explore <- counttab[!rownames(counttab) %in% genes_remove, ]
      anno_genes_remove <- subset(annotable, Gene_Symbol %in% genes_remove)$Ensembl_ID
      if (length(anno_genes_remove) != 0){
        counttab.explore <- counttab.explore[!rownames(counttab.explore) %in% anno_genes_remove, ]
      }
      dds <- DESeqDataSetFromMatrix(countData = counttab.explore, colData = coldata, design = ~ 1)
      dds <- DESeq(dds, parallel=T, betaPrior = TRUE)
      if (ncol(counttab) < 30){
        vsd <- rlogTransformation(dds, blind = T)
      } else {
        if (exists('vst')){
          vsd <- vst(dds, blind = T) # DESeq2 1.12+ 
        } else {
          vsd <- varianceStabilizingTransformation(dds, blind = T)
        }  
      }
    }
    
    counts.raw.df <- counts(dds)
    counts.norm.df <- counts(dds, normalized = T)
    
    ########
    # GSEA #
    ########
    # data formats: https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CLS:_Categorical_.28e.g_tumor_vs_normal.29_class_file_format_.28.2A.cls.29
    # data formats: https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GCT:_Gene_Cluster_Text_file_format_.28.2A.gct.29
    namestring <- paste(bfxid, 'normalisedCounts', 'GSEA', sep = '.')
    gsea.filename.gct <- createfilename(file.path(gseadir), namestring, appendix, 'gct')
    gsea.counts.df <- counts.norm.df
    gsea.header <- paste(c('#1.2', paste(nrow(gsea.counts.df), ncol(gsea.counts.df), sep = '\t')), sep = '\n')
    anno.temp <- annotable; anno.temp$Description <- paste(annotable$Gene_Symbol, annotable$Biotype, annotable$Description, sep = '; ')
    anno.temp <- anno.temp[, c('Ensembl_ID', 'Description')]
    gsea.counts.df <- merge(anno.temp, gsea.counts.df, by.x = "Ensembl_ID", by.y = "row.names")

    con <- file(gsea.filename.gct, open="wt")
    writeLines(gsea.header, con)
    write.table(gsea.counts.df, con, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    close(con)
    anno.temp <- NULL; gsea.counts.df <- NULL; gsea.filename.gct <- NULL

    # Phenotype data in CLS format
    for(c in colnames(coldata)[2:ncol(coldata)]) {
      v <-  coldata[, c, drop=TRUE]
      namestring <- paste(bfxid, 'phenotypeLabels', 'GSEA', c, sep = '.')
      gsea.filename.cls <- createfilename(gseadir, namestring, appendix, 'cls')
      con <- file(gsea.filename.cls, open="wt")
      writeLines(paste(length(v), length(unique(v)), 1), con)
      writeLines(paste(c("#", unique(v)), collapse=" "), con)
      writeLines(paste(v, collapse=" "), con)
      close(con)
    }
    if (ncol(coldata) > 2) {
      v <- gsub(' ', '_', do.call(paste, coldata[,colnames(coldata)[2:ncol(coldata)]]))
      c <- paste(colnames(coldata)[2:ncol(coldata)], collapse = '_')
      namestring <- paste(bfxid, 'phenotypeLabels', 'GSEA', c, sep = '.')
      gsea.filename.cls <- createfilename(gseadir, namestring, appendix, 'cls')
      con <- file(gsea.filename.cls, open="wt")
      writeLines(paste(length(v), length(unique(v)), 1), con)
      writeLines(paste(c("#", unique(v)), collapse=" "), con)
      writeLines(paste(v, collapse=" "), con)
      close(con)
    }
    
    #########
    # PCA
    #########
    createplot <- subset(configfile, key == 'pca')$value
    pcagenes <- subset(configfile, key == 'pcagenes')$value %>% as.integer()
    dopca <- FALSE
    
    pcadata <- getPCAdata(vsd, exploratory, pcagenes)
    addWorksheet(deseq.wb, 'pca-data')
    writeData(deseq.wb, 'pca-data', pcadata)
    
    # PCA 1-2
    if(createplot == '1'){
      dopca <- TRUE
      filenamepdf <- createfilename(pdfdir, paste(bfxid, 'pca_1-2', sep = '.'), appendix, 'pdf')
      filenamepng <- createfilename(pngdir, paste(bfxid, 'pca_1-2', sep = '.'), appendix, 'png')
      plotPCAVariable(vsd, exploratory, 1, 2, 'variable', filenamepdf, pcagenes)
      plotPCAVariable(vsd, exploratory, 1, 2, 'variable', filenamepng, pcagenes)
      whereiswhat <- rbind(whereiswhat, data.frame(what = 'pca12', xls = '', pdf = filenamepdf, png = filenamepng, stringsAsFactors = FALSE))
    }
    # PCA 2-3
    createplot <- subset(configfile, key == 'pca23')$value
    if(createplot=='1'){
      filenamepdf <- createfilename(pdfdir, paste(bfxid, 'pca_2-3', sep = '.'), appendix, 'pdf')
      filenamepng <- createfilename(pngdir, paste(bfxid, 'pca_2-3', sep = '.'), appendix, 'png')
      plotPCAVariable(vsd, exploratory, 2, 3, 'variable', filenamepdf, pcagenes)
      plotPCAVariable(vsd, exploratory, 2, 3, 'variable', filenamepng, pcagenes)
      whereiswhat <- rbind(whereiswhat, data.frame(what = 'pca23', xls = '', pdf = filenamepdf, png = filenamepng))
    }
    #########
    # BiPlot
    #########
    createplot <- subset(configfile, key == 'biplot')$value
    biplotgenes <- subset(configfile, key == 'biplotgenes')$value %>% as.integer()
    dobiplot <- FALSE
    
    if(createplot == '1' & biplotgenes > 0){
      dobiplot <- TRUE
      filenamepdf <- createfilename(pdfdir, paste(bfxid, 'biplot', sep = '.'), appendix, 'pdf')
      filenamepng <- createfilename(pngdir, paste(bfxid, 'biplot', sep = '.'), appendix, 'png')
      plotBiplot(vsd, geneensembltab, pcagenes, biplotgenes, coldata[, exploratory[1]], exploratory[1], filenamepdf)
      plotBiplot(vsd, geneensembltab, pcagenes, biplotgenes, coldata[, exploratory[1]], exploratory[1], filenamepng)
      whereiswhat <- rbind(whereiswhat, data.frame(what = 'biplot', xls = '', pdf = filenamepdf, png = filenamepng, stringsAsFactors = FALSE))
    }
    #########
    # MDS
    #########
    createplot <- subset(configfile, key == 'mds')$value
    if(createplot=='1'){
      namestring <- paste(bfxid, 'mds_plot', sep = '.')
      filenamepdf <- createfilename(pdfdir, namestring, appendix, 'pdf')
      filenamepng <- createfilename(pngdir, namestring, appendix, 'png')
      plotMDS(vsd, exploratory, samplename = "variable", filename = filenamepdf)
      plotMDS(vsd, exploratory, samplename = "variable", filename = filenamepng)
      whereiswhat <- rbind(whereiswhat, data.frame(what = 'mds', xls = '', pdf = filenamepdf, png = filenamepng, stringsAsFactors = FALSE))
    }
    #########
    # Poisson distance
    #########
    createplot <- subset(configfile, key == 'poissondist')$value
    if(createplot=='1'){
      namestring <- paste(bfxid, 'poisson_plot', sep = '.')
      filenamepdf <- createfilename(pdfdir, namestring, appendix, 'pdf')
      filenamepng <- createfilename(pngdir, namestring, appendix, 'png')
      plotPoissonDistance(counts.raw.df, dds, exploratory, samplename = "variable", filenamepdf)
      plotPoissonDistance(counts.raw.df, dds, exploratory, samplename = "variable", filenamepng)
      whereiswhat <- rbind(whereiswhat, data.frame(what = 'poisson', xls = '', pdf = filenamepdf, png = filenamepng, stringsAsFactors = FALSE))
    }
    #########
    # spearman and pearson correlation
    #########
    createplot <- subset(configfile, key == 'pearson')$value
    sorter <- ifelse(createplot %in% c('0', '1N'), FALSE, TRUE)
    if(createplot %in% c('1', '1N')){
      namestring <- paste(bfxid, 'pearson_plot', sep = '.')
      filenamepdf <- createfilename(pdfdir, namestring, appendix, 'pdf')
      filenamepng <- createfilename(pngdir, namestring, appendix, 'png')
      plotCorrelation(counts.norm.df, coldatanames, "pearson", sorter, filenamepdf)
      plotCorrelation(counts.norm.df, coldatanames, "pearson", sorter, filenamepng)
      whereiswhat <- rbind(whereiswhat, data.frame(what = 'pearson', xls = '', pdf = filenamepdf, png = filenamepng, stringsAsFactors = FALSE))
    }
    createplot <- subset(configfile, key == 'spearman')$value
    sorter <- ifelse(createplot %in% c('0', '1N'), FALSE, TRUE)
    if(createplot %in% c('1', '1N')){
      namestring <- paste(bfxid, 'spearman_plot', sep = '.')
      filenamepdf <- createfilename(pdfdir, namestring, appendix, 'pdf')
      filenamepng <- createfilename(pngdir, namestring, appendix, 'png')
      plotCorrelation(counts.norm.df, coldatanames, "spearman", sorter, filenamepdf)
      plotCorrelation(counts.norm.df, coldatanames, "spearman", sorter, filenamepng)
      whereiswhat <- rbind(whereiswhat, data.frame(what = 'spearman', xls = '', pdf = filenamepdf, png = filenamepng, stringsAsFactors = FALSE))
    }
    #########
    # heatmap
    #########
    createplot <- subset(configfile, key == 'heatmap')$value
    sorter <- ifelse(createplot %in% c('0', '1N'), FALSE, TRUE)
    doheatmap <- FALSE
    if(createplot %in% c('1', '1N')){
      doheatmap <- TRUE
      namestring <- paste(bfxid, 'heatmap_plot', sep = '.')
      filenamepdf <- createfilename(pdfdir, namestring, appendix, 'pdf')
      filenamepng <- createfilename(pngdir, namestring, appendix, 'png')
      plotHeatmap(vsd, coldatanames, sorter, filenamepdf)
      plotHeatmap(vsd, coldatanames, sorter, filenamepng)
      whereiswhat <- rbind(whereiswhat, data.frame(what = 'heatmap', xls = '', pdf = filenamepdf, png = filenamepng, stringsAsFactors = FALSE))
    }
    createplot <- subset(configfile, key == 'heatmappoisson')$value
    sorter <- ifelse(createplot %in% c('0', '1N'), FALSE, TRUE)
    if(createplot %in% c('1', '1N')){
      namestring <- paste(bfxid, 'heatmappoisson_plot', sep = '.')
      filenamepdf <- createfilename(pdfdir, namestring, appendix, 'pdf')
      filenamepng <- createfilename(pngdir, namestring, appendix, 'png')
      plotHeatmapPoiClaClu(counts.raw.df, coldatanames, sorter, filenamepdf)
      plotHeatmapPoiClaClu(counts.raw.df, coldatanames, sorter, filenamepng)
      whereiswhat <- rbind(whereiswhat, data.frame(what = 'heatmappoisson', xls = '', pdf = filenamepdf, png = filenamepng, stringsAsFactors = FALSE))
    }
    
    #########
    # vargenes
    #########
    createplot <- subset(configfile, key == 'varplot')$value
    vargenes <- subset(configfile, key == 'vargenes')$value
    vargenes <- suppressWarnings(subset(configfile, key == 'vargenes')$value %>% strsplit(',') %>% unlist() %>% as.integer())
    sorter <- ifelse(createplot %in% c('0', '1N'), FALSE, TRUE)
    dovargenes <- FALSE
    
    if(createplot %in% c('1', '1N')){
      if (!is.na(vargenes[1])) {
        filevecpdf <- c(); filevecpng <- c();dovargenes <- TRUE
        for (i in vargenes) {
          namestring <- paste0('varGenes_', i)
          filenamepdf <- createfilename(pdfdir, paste(bfxid, namestring, sep = '.'), appendix, 'pdf')
          filenamepng <- createfilename(pngdir, paste(bfxid, namestring, sep = '.'), appendix, 'png')
          plotVarGenes(vsd, geneensembltab, coldatanames, i, sorter, filenamepdf)
          plotVarGenes(vsd, geneensembltab, coldatanames, i, sorter, filenamepng)
          filevecpdf <- c(filevecpdf, filenamepdf)
          filevecpng <- c(filevecpng, filenamepng)
        }
        whereiswhat <- rbind(whereiswhat, data.frame(what = 'varGenes', xls = '', pdf = paste(filevecpdf, collapse = ','), png = paste(filevecpng, collapse = ',')))
      }
    }
    
    #########
    # highexpressedgenes
    #########
    createplot <- subset(configfile, key == 'highexpressedplot')$value
    highgenes <- suppressWarnings(subset(configfile, key == 'highgenes')$value %>% strsplit(',') %>% unlist() %>% as.integer())
    sorter <- ifelse(createplot %in% c('0', '1N'), FALSE, TRUE)
    
    if(createplot %in% c('1', '1N')){
      if (!is.na(highgenes[1])) {
        filevecpdf <- c(); filevecpng <- c()
        for (i in highgenes) {
          namestring <- paste0('highExpressedGenes_', i)
          filenamepdf <- createfilename(pdfdir, paste(bfxid, namestring, sep = '.'), appendix, 'pdf')
          filenamepng <- createfilename(pngdir, paste(bfxid, namestring, sep = '.'), appendix, 'png')
          plotHighExpressedGenes(vsd, geneensembltab, coldatanames, i, sorter, filenamepdf)
          plotHighExpressedGenes(vsd, geneensembltab, coldatanames, i, sorter, filenamepng)
          filevecpdf <- c(filevecpdf, filenamepdf)
          filevecpng <- c(filevecpng, filenamepng)
        }
        whereiswhat <- rbind(whereiswhat, data.frame(what = 'highExpressedGenes', xls = '', pdf = paste(filevecpdf, collapse = ','), png = paste(filevecpng, collapse = ','), stringsAsFactors = FALSE))
      }
    }
    
    #########
    # plot Genes
    #########
    genenamefile <- subset(configfile, key == 'plotcounts')$value
    if ((file.exists(genenamefile)) & (file.info(genenamefile)$size != 0)) {
      genenamesdf <- read.table(genenamefile, header = F, sep = '\t', col.names = c('Input'))
      genenamesdf.ensembl <- subset(geneensembltab, Ensembl_ID %in% genenamesdf$Input)
      genenamesdf.gene <- subset(geneensembltab, Gene_Symbol %in% genenamesdf$Input)
      genenamesdf <- rbind(genenamesdf.ensembl, genenamesdf.gene)
      if (nrow(genenamesdf)!=0) {
        filevecpdf <- c(); filevecpng <- c()
        for (i in 1:nrow(genenamesdf)) {
          if (genenamesdf[i, 1] %in% rownames(counts.norm.df)) {
            namestring <- paste0('normalisedCounts_', genenamesdf[i, 2])
            filenamepdf <- createfilename(pdfdir, paste(bfxid, namestring, sep = '.'), appendix, 'pdf')
            filenamepng <- createfilename(pngdir, paste(bfxid, namestring, sep = '.'), appendix, 'png')
            plotNormCounts(as.character(genenamesdf[i, 1]), as.character(genenamesdf[i, 2]), dds, exploratory, filenamepdf)
            plotNormCounts(as.character(genenamesdf[i, 1]), as.character(genenamesdf[i, 2]), dds, exploratory, filenamepng)
            filevecpdf <- c(filevecpdf, filenamepdf)
            filevecpng <- c(filevecpng, filenamepng)
          }
        }
        whereiswhat <- rbind(whereiswhat, data.frame(what = 'normalisedGenes', xls = '', pdf = paste(filevecpdf, collapse = ','), png = paste(filevecpng, collapse = ','), stringsAsFactors = FALSE))
      }
    }
    ###########
    # Batch
    ###########
    # that is for doing a batch correction
    correctionvec <- subset(configfile, key == 'correction')$value
    correctionvec <- gsub(' ', '', correctionvec) %>% strsplit(',') %>% unlist()
    if ((length(correctionvec) >= 1) & (sum(correctionvec %in% exploratory) == length(correctionvec))) {
      correctdf <- data.frame(colData(dds))[, correctionvec, drop = F]
      correctedvsdmatrix <- removeBatchEffectWrapper(vsd, correctdf)
      
      vsd.corrected <- as.data.frame(correctedvsdmatrix)
      zscore.corrected <- vsd.corrected
      vsd.corrected$Ensembl_ID <- rownames(vsd.corrected)
      vsd.corrected <- merge(vsd.corrected, annotable, by = 'Ensembl_ID')
      addWorksheet(deseq.wb, paste0(trans.method, '-values-batch-corr'))
      writeData(deseq.wb, paste0(trans.method, '-values-batch-corr'), vsd.corrected[, c(setdiff(colnames(vsd.corrected), colnames(counttab)), colnames(counttab))])
      vsdbatchwritten <- TRUE
      
      zscore.corrected <- build_zscore(zscore.corrected)
      zscore.corrected$Ensembl_ID <- rownames(zscore.corrected)
      zscore.corrected <- merge(zscore.corrected, annotable, by = 'Ensembl_ID')
      addWorksheet(deseq.wb, 'z-score-batch-corr')
      writeData(deseq.wb, 'z-score-batch-corr', zscore.corrected[, c(setdiff(colnames(zscore.corrected), colnames(counttab)), colnames(counttab))])
      
      #########
      # close workbook for counts
      #########
      saveWorkbook(deseq.wb, deseq.wb.filename, overwrite = T)
      
      if (dopca) {
        filenamepdf <- createfilename(pdfdir, paste(bfxid, 'batch_pca_1-2', paste0('batchcorrection-', paste(correctionvec, collapse = '-')), sep = '.'), appendix, 'pdf')
        filenamepng <- createfilename(pngdir, paste(bfxid, 'batch_pca_1-2', paste0('batchcorrection-', paste(correctionvec, collapse = '-')), sep = '.'), appendix, 'png')
        title <- paste0("PCA (1-2) of top ", pcagenes, " most diverse genes (correction for ", paste(correctionvec, collapse = ', '), ")")
        plotPCAVariable(vsd, exploratory, 1, 2, 'variable', filenamepdf, pcagenes, title, correctedvsdmatrix)
        plotPCAVariable(vsd, exploratory, 1, 2, 'variable', filenamepng, pcagenes, title, correctedvsdmatrix)
        whereiswhat <- rbind(whereiswhat, data.frame(what = 'batchpca12', xls = '', pdf = filenamepdf, png = filenamepng, stringsAsFactors = FALSE))
      }
      if (dobiplot) {
        filenamepdf <- createfilename(pdfdir, paste(bfxid, 'batch_biplot', paste0('batchcorrection-', paste(correctionvec, collapse = '-')), sep = '.'), appendix, 'pdf')
        filenamepng <- createfilename(pngdir, paste(bfxid, 'batch_biplot', paste0('batchcorrection-', paste(correctionvec, collapse = '-')), sep = '.'), appendix, 'png')
        title <- paste0("PCA - BiPlot (", biplotgenes, " genes)", paste0("\n(correction for ", paste(correctionvec, collapse = ', '),")"))
        plotBiplot(vsd, geneensembltab, pcagenes, biplotgenes, coldata[, exploratory[1]], exploratory[1], filenamepdf, title, correctedvsdmatrix)
        plotBiplot(vsd, geneensembltab, pcagenes, biplotgenes, coldata[, exploratory[1]], exploratory[1], filenamepng, title, correctedvsdmatrix)
        whereiswhat <- rbind(whereiswhat, data.frame(what = 'batchbiplot', xls = '', pdf = filenamepdf, png = filenamepng))
      }
      if (dovargenes){
        createplot <- subset(configfile, key == 'varplot')$value
        sorter <- ifelse(createplot %in% c('0', '1N'), FALSE, TRUE)
        filevecpdf <- c(); filevecpng <- c()
        for (i in vargenes) {
          namestring <- paste(bfxid, paste0('batch_varGenes_', i), paste0('batchcorrection-', paste(correctionvec, collapse = '-')), sep = '.')
          filenamepdf <- createfilename(pdfdir, namestring, appendix, 'pdf')
          filenamepng <- createfilename(pngdir, namestring, appendix, 'png')
          title = paste0("Deviation from the gene's average across all samples (top ", i, " of genes)\n", paste0("(correction for ", paste(correctionvec, collapse = ', '),")"))
          plotVarGenes(vsd, geneensembltab, coldatanames, i, sorter, filenamepdf, title, correctedvsdmatrix)
          plotVarGenes(vsd, geneensembltab, coldatanames, i, sorter, filenamepng, title, correctedvsdmatrix)
          filevecpdf <- c(filevecpdf, filenamepdf)
          filevecpng <- c(filevecpng, filenamepng)
        }
        whereiswhat <- rbind(whereiswhat, data.frame(what = 'batchvarGenes', xls = '', pdf = paste(filevecpdf, collapse = ','), png = paste(filevecpng, collapse = ','), stringsAsFactors = FALSE))
      }
      if(doheatmap){
        createplot <- subset(configfile, key == 'heatmap')$value
        sorter <- ifelse(createplot %in% c('0', '1N'), FALSE, TRUE)
        namestring <- paste(bfxid, 'batch_heatmap_plot', paste0('batchcorrection-', paste(correctionvec, collapse = '-')), sep = '.')
        filenamepdf <- createfilename(pdfdir, namestring, appendix, 'pdf')
        filenamepng <- createfilename(pngdir, namestring, appendix, 'png')
        title <- paste0("Heatmap of Euclidean distance ", paste0("(correction for ", paste(correctionvec, collapse = ', '),")"))
        plotHeatmap(vsd, coldatanames, sorter, filenamepdf, title, correctedvsdmatrix)
        plotHeatmap(vsd, coldatanames, sorter, filenamepng, title, correctedvsdmatrix)
        whereiswhat <- rbind(whereiswhat, data.frame(what = 'batchheatmap', xls = '', pdf = filenamepdf, png = filenamepng, stringsAsFactors = FALSE))
      }
    } else if(length(correctionvec) == 0) {
      print('Skipping Correction!')
      correctedvsdmatrix <- NULL
    } else {
      print(paste0('Correction variable: ', paste0(correctionvec, collapse = ',')))
      print(paste0('Factors from exploratory: ', paste(exploratory, collapse = ',')))
      print('Some entries in the correction variable are not listed in the exploratory vector. Skipping Correction!')
      correctedvsdmatrix <- NULL
    }
    
    #########
    # close workbook for counts
    #########
    print('PREPARING OUTPUT FOR COUNTS TABLE')
    saveWorkbook(deseq.wb, deseq.wb.filename, overwrite = T)
    if (vsdbatchwritten == FALSE) {
      whereiswhat <- rbind(whereiswhat, data.frame(what = 'countstable', xls = deseq.wb.filename, pdf = '', png = '0'))
    } else{
      whereiswhat <- rbind(whereiswhat, data.frame(what = 'countstable', xls = deseq.wb.filename, pdf = '', png = '1'))
    }
    
    #####################
    # Data for morpheus #
    #####################
    namestring <- paste(bfxid, 'column_annotation', 'Morpheus', sep = '.')
    morpheus.filename <- createfilename(file.path(datadir), namestring, appendix, 'xlsx')
    openxlsx::write.xlsx(coldata, morpheus.filename, rowNames=FALSE)
    whereiswhat <- rbind(whereiswhat, data.frame(what = 'morpheus', xls = morpheus.filename, pdf = '', png = ''))
  }
}

#########
# expression analysis
#########
contrastfilename <- subset(configfile, key == 'contrastfile')$value
contrastdf <- read.table(contrastfilename, header = T, sep = '\t', as.is = T, fill = T)

contrastdf$appendix <- ifelse(is.na(contrastdf$appendix), 'none', contrastdf$appendix)
contrastdf$comment <- ifelse(is.na(contrastdf$comment), '', contrastdf$comment)

doglobalpadj <- subset(configfile, key == 'doglobalpadj')$value
useIHW <- subset(configfile, key == 'useIHW')$value
dohtml <- subset(configfile, key == 'htmlreport')$value

if(nrow(contrastdf)!=0) {contrastdf$formula <- gsub(' ', '', contrastdf$formula)}
padjvalues <- subset(configfile, key == 'padjvalues')$value

# building the combine factor, if given
combinevalues <- subset(configfile, key == 'combine')$value %>% gsub(' ', '', .) %>% strsplit(',') %>% unlist()
if (length(combinevalues) != 0) {
  coldata$Comb <- interaction(coldata[, combinevalues], drop = T, sep = '')
}

decomp <- data.frame(Comparison = character(), Factor = character(), Formula = character(), FDR = character(), DEGenes = integer(), 
                     UpRegulated = integer(), DownRegulated = integer(), MA = character(), Volcano = character(), Count = integer(),
                     Comment = character(), App = character(), Num = character(), Denom = character())
deresults <- vector('list'); countresults <- vector('list'); countresultsrounded <- vector('list')

if (padjvalues != 'none') {
  padjvalues <- suppressWarnings(strsplit(padjvalues, ',') %>% unlist() %>% as.numeric())
  counter <- 1
  for (deformula in unique(contrastdf$formula)) {
    if (checkFormula(deformula, colnames(coldata))) {next} # check if the formula has valid factors
    deformulavector <- getFormulavector(deformula)
    coldatatemp <- removeNoneSampleFromColdata(deformulavector, coldata, combinevalues) # remove samples which are marked as none in the colums of coldata
    counttabtemp <- counttab[, as.character(coldatatemp$Sample)] # reduce the counttab to the selected samples
    print(paste0('FORMULA: ', deformula, ' - COMPARISON: ', nrow(subset(contrastdf, formula == deformula))))
    deseqformula <- as.formula(deformula) # create formula from string object
    ddsrun <- DESeqDataSetFromMatrix(countData = counttabtemp, colData = coldatatemp, design = deseqformula)
    ddsrun <- DESeq(ddsrun, parallel=T, betaPrior = TRUE)
    factors <- unique(subset(contrastdf, formula == deformula)$factor) # get the unique factors
    for (fact in factors) {
      comparisontab <- subset(contrastdf, factor == fact & formula == deformula)
      for (compindex in 1:nrow(comparisontab)) {
        if (checkFactorCorrectness(comparisontab, compindex, coldata)) {next} # checking if factors, conditions are correct
        runvec <- getFactorNumDenomComment(comparisontab, compindex)
        print(paste0("FORMULA: '", deformula, "' - COMPARISON - FACTOR: '", fact))
        for(padjvalue in padjvalues) {
          padjstring <- getpadjstring(padjvalue) # format the padj to a string
          # here it is time to add ihw
          if (useIHW == '1') {
            require(IHW)
            tempderesult <- results(ddsrun, contrast = c(runvec[1], runvec[2], runvec[3]), alpha = padjvalue, format = "DataFrame", parallel = T, filterFun = ihw)
          } else {
            tempderesult <- results(ddsrun, contrast = c(runvec[1], runvec[2], runvec[3]), alpha = padjvalue, format = "DataFrame", parallel = T)
          }
          deresults[[counter]] <- tempderesult # store the de-analysis, can be retrieved by checking decomp data frame for Index
          tempderesult_anno <- prepResultCounttable(tempderesult, ddsrun, runvec, annotable, doglobalpadj) # add counts, sort by padj, remove genes which are na
          tempderesult_anno = tempderesult_anno[, c(colnames(annotable),setdiff(colnames(tempderesult_anno), colnames(annotable)))]
          countresultsrounded[[counter]] <- roundResults(tempderesult_anno, 2) # round all numbers
          countresults[[counter]] <- tempderesult_anno
          if (doglobalpadj == '0') {
            de <- nrow(subset(tempderesult, padj < padjvalue)); deup <- nrow(subset(tempderesult, padj < padjvalue & log2FoldChange > 0)); dedown <- nrow(subset(tempderesult, padj < padjvalue & log2FoldChange < 0))
            filenamema <- doMAplotting(tempderesult, runvec, bfxid, pngdir, pdfdir, de, padjstring, padjvalue)
            filenamevolcano <- doVolcanoPlotting(tempderesult, runvec, bfxid, pngdir, pdfdir, padjstring, padjvalue, geneensembltab, 0, 0, 15)
            decomp <- rbind(decomp, data.frame(Comparison = c(paste0(runvec[4], ' vs ', runvec[5])), Factor = c(fact), Formula = c(deformula),
                                               FDR = c(padjvalue*100), DEGenes = c(de), UpRegulated = c(deup), DownRegulated = c(dedown),
                                               MA=filenamema, Volcano = filenamevolcano, Count = c(counter), Comment = c(runvec[8]), App = c(runvec[7]),
                                               Num = runvec[2], Denom = runvec[3]))
          } else {
            decomp <- rbind(decomp, data.frame(Comparison = c(paste0(runvec[4], ' vs ', runvec[5])), Factor = c(fact), Formula = c(deformula),
                                               FDR = c(padjvalue*100), DEGenes = c(0), UpRegulated = c(0), DownRegulated = c(0),
                                               MA=c('no'), Volcano = c('no'), Count = c(counter), Comment = c(runvec[8]), App = c(runvec[7]),
                                               Num = runvec[2], Denom = runvec[3]))
          }
          
          counter <- counter + 1
          # rm(tempderesult, tempderesult_anno)
        }
      }
    }
  }
  
  # global padj of the results is applied here
  if (doglobalpadj == '1' & nrow(contrastdf) > 0) {
    decomp$MA <- as.character(decomp$MA)
    decomp$Volcano <- as.character(decomp$Volcano)
    for (padjvalue in padjvalues) {
      loi <- subset(decomp, FDR == 100 * padjvalue)$Count # get the resultobj with the same fdr
      if (length(loi) == 0) {next}
      pvaluevec <- lapply(loi, function(x){countresults[[x]]$pvalue}) %>% unlist() %>% p.adjust(method = 'BH') # build a vector and apply BH
      pvaluevec <- split(pvaluevec, rep(1:length(loi), each = nrow(counttab))) # split it in the number of initial resultobj
      pvaluevecraw <- lapply(loi, function(x){deresults[[x]]$pvalue}) %>% unlist() %>% p.adjust(method = 'BH') # build a vector and apply BH
      pvaluevecraw <- split(pvaluevecraw, rep(1:length(loi), each = nrow(counttab))) # split it in the number of initial resultobj
      pvaluecount <- 1
      for (i in loi) {
        deresults[[i]] <- addGlobalPadj(deresults[[i]], pvaluevecraw[[pvaluecount]]) # calculate new p adjusted
        countresults[[i]] <- addGlobalPadj(countresults[[i]], pvaluevec[[pvaluecount]]) # calculate new p adjusted
        countresultsrounded[[i]] <- addGlobalPadj(countresultsrounded[[i]], pvaluevec[[pvaluecount]]) # calculate new p adjusted
        decomp[i, 'DEGenes'] <- nrow(subset(deresults[[i]], padjGlob < padjvalue))
		decomp[i, 'UpRegulated'] <- nrow(subset(deresults[[i]], padjGlob < padjvalue & log2FoldChange > 0))
		decomp[i, 'DownRegulated'] <- nrow(subset(deresults[[i]], padjGlob < padjvalue & log2FoldChange < 0))  
		
        filenamema <- doMAGlobalplotting(deresults[[i]], decomp[i, ], bfxid, pngdir, pdfdir) # new ma
        decomp[i, 'MA'] <- filenamema; # decomp[i, 'CSV'] <- as.character(filenamecsv)
		
		filenamevolcano <- doVolcanoGlobalPlotting(deresults[[i]],  decomp[i, ], bfxid, pngdir, pdfdir, geneensembltab, 0, 0, 15)
		decomp[i, 'Volcano'] <- filenamevolcano

        pvaluecount <- pvaluecount + 1
      }
    }
  }
  # create and write the combined results files
  if (nrow(contrastdf) > 1) {
    tabnames <- c('Summary', paste0('FDR ', padjvalues*100))
    combinedresultslist <- lapply(padjvalues, buildCombinedComparisionDataFrame, countresultsrounded, decomp, annotable)
    combinedresultslist <- lapply(combinedresultslist, function(df) {
      return(df[, c(colnames(annotable),setdiff(colnames(df), colnames(annotable)))])
    })    

    # function builds the overview of the list which is created by the previous called function
    # lapply returns a list and ldply makes it into a data.frame
    combinedresultsinfo <- as.data.frame(rbindlist(lapply(seq_along(padjvalues), function(x) {buildSubsetFromSummaryTable(padjvalues[x], decomp, x)})))
    combinedresultlist <- append(list(combinedresultsinfo), combinedresultslist)
    combinefilenamexlsx <- createfilename(datadir, paste(bfxid, 'deseq-results', 'combined', appendix, sep = '.'), '', 'xlsx')
    writeXLSX(combinedresultlist, tabnames, combinefilenamexlsx)
    whereiswhat <- rbind(whereiswhat, data.frame(what = 'combineresultsDE', xls = combinefilenamexlsx, pdf = '', png = '', stringsAsFactors = FALSE))
  }
  # all results in separate sheets
  allresultsfilenamexlsx <- createfilename(datadir, paste(bfxid, 'deseq-results', 'separate', appendix, sep = '.'), '', 'xlsx')
  allresultstabnames <- as.vector(paste(decomp$Count, decomp$Factor, decomp$Num,'v', decomp$Denom, paste0('F', decomp$FDR), sep = '_'))
  decomptemp <- decomp[, c('Count', 'Comparison', 'Factor', 'Formula', 'FDR', 'DEGenes', 'UpRegulated', 'DownRegulated', 'Comment')]
  colnames(decomptemp)[1] <- 'TabIndex'
  allresultslist <- append(list(decomptemp), countresultsrounded)
  tabnames <- c('Summary', allresultstabnames)
  if (sum(nchar(tabnames)>31) == 0) {
    writeXLSX(allresultslist, tabnames, allresultsfilenamexlsx)
    whereiswhat <- rbind(whereiswhat, data.frame(what = 'allresultsDE', xls = allresultsfilenamexlsx, pdf = '', png = '', stringsAsFactors = FALSE))
  } else {
    stop(paste('NOT POSSIBLE TO WRITE ALL DE-RESULTS IN ONE EXCEL FILE BECAUSE CONDITION NAMES ARE TOO LONG! Only 31 characters can be used, not more:', paste(tabnames, collapse=','), '. Shorten your condition names and values.'))
  }
}

desummaryfilename <- createfilename(workingdir, paste(bfxid, 'deseq-results', 'summary', sep = '.'), '', 'csv')
printResultSummary(decomp[, c('Comparison', 'Factor', 'Formula', 'FDR', 'DEGenes', 'UpRegulated', 'DownRegulated', 'App', 'Comment')], desummaryfilename)
writeLines(capture.output(sessionInfo()), file.path(workingdir, "sessionInfo.txt"))
save.image()
invisible(if (file.exists(file.path(workingdir, 'Rplots.pdf'))) {file.remove(file.path(workingdir, 'Rplots.pdf'))})
