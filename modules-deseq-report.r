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

mkdirectory <- function(dirname) {
  dir.create(dirname, showWarnings = FALSE)
}

changeContinuous <- function(df){
  for (i in 1:ncol(df)){
    if (is.numeric(df[, i])) {
      df[, i] <- factor(df[, i], levels = sort(unique(df[, i])))
    }
  }
  return(df)
}

# function creates filename out of 4 parts
# input
# dir ... storage directory
# firstpart ... string
# appendix ... string
# extension ... file suffix
createfilename <- function(dir, firstpart, appendix, extension) {
  if (appendix != '') {
    filename <- file.path(dir, paste(firstpart, appendix, extension, sep = '.'))
  } else {
    filename <- file.path(dir, paste(firstpart, extension, sep = '.'))
  }
  return(filename)
}

# creates the name for the expression tables
# input
# worktemp ... data storage directory
# runvectemp ... helper vector having useful strings for the name
# bfxid ... bfxid
# app ... appendix for the file
# degene ... number of degs
# padjstr ... string representing the padj value
createCSVname <- function(worktemp, runvectemp, bfxid, degene, padjstr) {
  namestring <- paste0(bfxid, '.deseq-results.', paste(runvectemp[6], collapse = '_'), '-', runvectemp[1], '-', runvectemp[4], '_vs_', runvectemp[5], '.', padjstr, '-', degene)
  resultstablecsvname <- createfilename(worktemp, namestring, runvectemp[7], 'csv')
  return(resultstablecsvname)
}

# similar function as above but for expression tables with global padj
# input
# worktemp ... data storage directory
# decomptemp ... data frame having usefule strings
# bfxid ... bfxid
# app ... appendix for the file
createCSVglobalname <- function(worktemp, decomptemp, bfxid) {
  deformulavector <- getFormulavector(as.character(decomptemp[, 'Formula']))
  fact <- as.character(decomptemp[, 'Factor'])
  comp <- gsub(' ', '_', as.character(decomptemp[, 'Comparison']))
  gene <- decomptemp[, 'DEGenes']
  padjstr <- as.character(decomptemp[, 'FDR'])
  padjstr <- ifelse(nchar(padjstr) == 1, paste0('00', padjstr), paste0('0', padjstr))
  
  namestring <- paste0(bfxid, '.deseq-results.', paste(deformulavector, collapse = '_'), '-', fact, '-', comp, '.', padjstr, '-', gene)
  resultstablecsvname <- createfilename(worktemp, namestring, decomptemp[, 'App'], 'csv')
  return(resultstablecsvname)
}

# creates a workbook and writes the inputlist back to separate sheets
# input
# inputlist ... list of data frames
# tabnames ... names for singe spread sheets
# filename ... name of the excel file
writeXLSX <- function(inputlist, tabnames, filename) {
  require(openxlsx)
  wb <- createWorkbook()
  invisible(lapply(seq_along(inputlist), function(i){
    addWorksheet(wb, sheetName = tabnames[i])
    writeData(wb, sheet = tabnames[i], keepNA = T, x = inputlist[[i]], rowNames = F, colNames = T)
    setColWidths(wb, sheet = tabnames[i], widths = 'auto', cols = 1:ncol(inputlist[[i]]))
  }))
  saveWorkbook(wb, filename, overwrite = T)
}

writeRNK <- function(inputlist, tabnames, bfxid, appendix, directory) {
  invisible(lapply(seq_along(inputlist), function(i){
    filename <-file.path(directory, paste(bfxid, tabnames[i], appendix, 'rnk', sep = '.'))
    temp.df <- inputlist[[i]][, c('Ensembl_ID', 'stat')]
    temp.df <- temp.df[order(temp.df$stat, decreasing = T), ]
    write.table(temp.df, file = filename, sep = '\t', quote = F, row.names = F, col.names = F)
  }))
}

# write the combined results as separate csv files 
# input
# combinelist ... list of the combined result tables
# combineinfor ... info overview of the combined result tables
# workingdir ... data storage directory
# bfxid ... bfx id
# return ... vector of strings
writeCombineCSV <- function(combinelist, combineinfo, workingdir, bfxid) {
  names <- character()
  for (fdr in unique(combineinfo$FDR)) {
    table <- unique(subset(combineinfo, FDR == fdr)[, 'TabIndex'])
    padjstr <- ifelse(nchar(fdr) == 1, paste0('00', fdr), paste0('0', fdr))
    namestring <- paste(bfxid, 'deseq-results', paste('combined', padjstr, sep = '_'), sep = '.')
    filename <- createfilename(workingdir, namestring, '', 'csv')
    write.table(combinelist[[table]], filename, sep = '\t', row.names = F, col.names = T, quote = F)
    names <- append(names, filename)
  }
  return(names)
}

# Function which removes rows that have no reads
# Input
# datatable ... data.frame or matrix with counts
# fragmentcount ... filter for the number of fragments
# columns ... columns which are used for the counting
removeLowCounts <- function(datatable, fragmentcount = 0, columns=c(1:ncol(datatable))) {
  print(paste0("Number of genes: ", nrow(datatable)))
  datatable <- datatable[rowSums(datatable) > fragmentcount, ]
  print(paste0("Number of genes (non-zero): ", nrow(datatable)))
  return(datatable)
}

removeBatchEffectWrapper <- function(readtransform, inputdf) {
  if (ncol(inputdf) > 2) {
    warning('Can only correct for up to two conditions')
    return(NULL)
  }
  require(edgeR)
  inputmx <- assay(readtransform)

  if (ncol(inputdf) == 2) {
    inputmx <- removeBatchEffect(inputmx, inputdf[, 1], inputdf[, 2])
  } else if (ncol(inputdf)  == 1) {
    inputmx <- removeBatchEffect(inputmx, inputdf[, 1])  
  }
  return(inputmx)
}

# Function for getting a matrix with the pca data
# Input
# readtransform ... DESeqTransform which comes from the transformation process (vsd, rld, ...)
# conditions ... is a vector and entries should be column names from the colData(dds) or coldata
#                one is enough because we are only interested in the table
# ntop ... number of most diverse genes
getPCAdata <- function(readtransform, conditions, ntop=500){
  require(DESeq2)
  require(genefilter)
  
  if (!all(conditions %in% names(colData(readtransform)))) {
    stop("the argument 'conditions' should specify columns of colData(dds)")
  }
  
  topVarGenes <- head(order(rowVars(assay(readtransform)),decreasing=TRUE), ntop)
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(readtransform)[topVarGenes,]))
  # the contribution to the total variance for each component
  
  # add the intgroup factors together to create a new grouping factor
  intgroup.df <- data.frame(colData(readtransform)[, conditions, drop=FALSE])
  intgroup.df[] <- lapply(intgroup.df, as.character)
  intgroup.df <- rbind(intgroup.df, rep('NA', ncol(intgroup.df)))
  percentVar <- round(pca$sdev^2 / sum( pca$sdev^2 ) * 100)
  
  pcadata <- as.data.frame(rbind(pca$x, percentVar))
  pcadata$Sample <- c(colnames(readtransform), 'VARIANCE')
  pcadata <- cbind(pcadata, intgroup.df)
  return(pcadata)
}


# Function for getting a matrix with the pca data
# Input
# readtransform ... DESeqTransform which comes from the transformation process (vsd, rld, ...)
# conditions ... is a vector and entries should be column names from the colData(dds) or coldata
#                one is enough because we are only interested in the table
# ntop ... number of most diverse genes
getPCAdataVariable <- function(readtransform, conditions, comp1 = 1, comp2 = 2, ntop=500){
  require(DESeq2)
  require(genefilter)
  
  if (!all(conditions %in% names(colData(readtransform)))) {
    stop("the argument 'conditions' should specify columns of colData(dds)")
  }
  
  topVarGenes <- head(order(rowVars(assay(readtransform)),decreasing=TRUE), ntop)
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(readtransform)[topVarGenes,]))
  # the contribution to the total variance for each component
  
  # add the intgroup factors together to create a new grouping factor
  intgroup.df <- as.data.frame(colData(readtransform)[, conditions, drop=FALSE])
  group <- factor(apply( intgroup.df, 1, paste, collapse=" : "))
  
  # assembly the data for the plot
  pcadata <- data.frame(var1 = pca$x[, comp1], var2 = pca$x[,comp2], group=group, intgroup.df, name=colnames(readtransform))
  names(pcadata)[names(pcadata) == "var1"] <- paste0('PC', comp1)
  names(pcadata)[names(pcadata) == "var2"] <- paste0('PC', comp2)
  return(pcadata)
}

#function for getting the percent variation
# readtransform ... DESeqTransform which comes from the transformation process (vsd, rld, ...)
# conditions ... is a vector and entries should be column names from the colData(dds) or coldata
#               the order is import for the ggplot; it is: colour, shape, alpha
# ntop ... number of most diverse genes
getPCAdataPercentVarVariable <- function(readtransform, conditions, comp1 = 1, comp2 = 2, ntop=500){
  require(DESeq2)
  require(genefilter)
  
  if (!all(conditions %in% names(colData(readtransform)))) {
    stop("the argument 'conditions' should specify columns of colData(dds)")
  }
  
  topVarGenes <- head(order(rowVars(assay(readtransform)),decreasing=TRUE), ntop)
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(readtransform)[topVarGenes,]))
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  return(c(round(percentVar[comp1]*100), round(percentVar[comp2]*100)))
}

# Function adds counts to your resulttable
# Input
# resulttable ... resulttable either data.frame or matrix
# dds_copy ... dds object from DESeq2
# cond ... column index of your condition which you want to select for samples; from colData(dds_copy)
# cond1 ... first condition which is appended to the sample
#           as vector or string (usually the denominator in the expression comparison)
# cond2 ... second condition which is appended to the sample
#           as vector or string (usually the denominator in the expression comparison)
# isnormalized ... read count normalised or not
# samplescounts ... print counts for each sample of the conditions
# meansample... print the mean read counts for the condition
addCounts <- function(resulttable, dds_copy, cond, cond1, cond2, cond1name=c(), cond2name=c(), isnormalized=TRUE, samplescounts=TRUE, meansample=TRUE) {
  require(DESeq2)
  singlecolumn_s1 <- FALSE
  singlecolumn_s2 <- FALSE
  resulttable <- as.data.frame(resulttable)
  # that'a test to check wheter it's only a single sample condition, then only this count will be shown and mean is not calculated
  column_s1 <- as.data.frame(counts(dds_copy, normalized=TRUE)[, colData(dds_copy)[, cond] %in% cond1])
  column_s2 <- as.data.frame(counts(dds_copy, normalized=TRUE)[, colData(dds_copy)[, cond] %in% cond2])
  
  if (ncol(column_s1)==1) {
    singlecolumn_s1 <- TRUE
    meansample <- FALSE
  }
  if (ncol(column_s2)==1) {
    singlecolumn_s2 <- TRUE
    meansample <- FALSE
  }
  
  if(samplescounts == TRUE) {
    if(isnormalized == TRUE) {
      cond1_counts <- as.data.frame(counts(dds_copy, normalized=TRUE)[, colData(dds_copy)[, cond] %in% cond1])
      cond2_counts <- as.data.frame(counts(dds_copy, normalized=TRUE)[, colData(dds_copy)[, cond] %in% cond2])
      if (singlecolumn_s1) {colnames(cond1_counts)[1] <- ifelse(length(cond1name)==0, c(cond1), cond1name)}
      if (singlecolumn_s2) {colnames(cond2_counts)[1] <- ifelse(length(cond2name)==0, c(cond2), cond2name)}
      samplescounttable <- merge(cond1_counts, cond2_counts, by.x="row.names", by.y="row.names")
      rownames(samplescounttable) <- c(samplescounttable$Row.names)
      samplescounttable$Row.names <- NULL      
      samplescounttable <- round(samplescounttable,1)
    } else {
      cond1_counts <- as.data.frame(counts(dds_copy, normalized=FALSE)[, colData(dds_copy)[, cond] %in% cond1])
      cond2_counts <- as.data.frame(counts(dds_copy, normalized=FALSE)[, colData(dds_copy)[, cond] %in% cond2])
      if (singlecolumn_s1) {colnames(cond1_counts)[1] <- ifelse(length(cond1name)==0, c(cond1), cond1name)}
      if (singlecolumn_s2) {colnames(cond2_counts)[1] <- ifelse(length(cond2name)==0, c(cond2), cond2name)}
      samplescounttable <- merge(cond1_counts, cond2_counts, by.x="row.names", by.y="row.names")      
      rownames(samplescounttable) <- c(samplescounttable$Row.names)
      samplescounttable$Row.names <- NULL
    }
    resulttable <- merge(resulttable, samplescounttable, by.x="row.names", by.y="row.names")
    rownames(resulttable) <- resulttable$Row.names
    resulttable$Row.names <- NULL
  }
  if(meansample == TRUE) {
    if(isnormalized==TRUE) {
      cond1_counts <- as.data.frame(rowMeans(counts(dds_copy, normalized=TRUE)[, colData(dds_copy)[, cond] %in% cond1]))
      colnames(cond1_counts) <- ifelse(length(cond1name)==0, c(cond1), cond1name)
      cond2_counts <- as.data.frame(rowMeans(counts(dds_copy, normalized=TRUE)[, colData(dds_copy)[, cond] %in% cond2]))
      colnames(cond2_counts) <- ifelse(length(cond2name)==0, c(cond2), cond2name)
    } else {
      cond1_counts <- as.data.frame(rowMeans(counts(dds_copy, normalized=TRUE)[, colData(dds_copy)[, cond] %in% cond1]))
      colnames(cond1_counts) <- ifelse(length(cond1name)==0, c(cond1), cond1name)
      cond2_counts <- as.data.frame(rowMeans(counts(dds_copy, normalized=TRUE)[, colData(dds_copy)[, cond] %in% cond2]))
      colnames(cond2_counts) <- ifelse(length(cond2name)==0, c(cond2), cond2name)
    }
    cond1_cond2 <- merge(cond1_counts, cond2_counts, by.x="row.names", by.y="row.names")
    rownames(cond1_cond2) <- c(cond1_cond2$Row.names)
    cond1_cond2$Row.names <- NULL
    cond1_cond2 <- round(cond1_cond2, 1)
    resulttable <- merge(resulttable, cond1_cond2, by.x="row.names", by.y="row.names")
    rownames(resulttable) <- resulttable$Row.names
    resulttable$Row.names <- NULL
  }
  return(resulttable)
}

getFactorNumDenomComment <- function(contrasttemp, index) {
  fact <- as.character(contrasttemp[index, 'factor']) # factor for the comparison
  num <- as.character(contrasttemp[index, 'numerator']) # condition
  denom <- as.character(contrasttemp[index, 'denominator']) # condition to compare against
  numhead <- as.character(contrasttemp[index, 'numerator_name']) # string which is used for plots and column names
  denomhead <- as.character(contrasttemp[index, 'denominator_name']) # string which is used for plots and column names
  comment <- as.character(contrasttemp[index, 'comment'])
  appendix <- as.character(contrasttemp[index, 'appendix'])
  appendix <- ifelse(nchar(appendix) == 0, 'none', appendix)
  formula <- paste(getFormulavector(as.character(contrasttemp[index, 'formula'])), collapse = '_')
  numheadreplace <- gsub(' ', '_', numhead); denomheadreplace <- gsub(' ', '_', denomhead)
  return(c(fact, num, denom, numheadreplace, denomheadreplace, formula, appendix, comment))
}

# function creates a string representation of a fdr value
# Input
# padj ... padj value
getpadjstring <- function(padj){
  padj <- gsub('\\.', '', as.character(padj))
  padj <- ifelse(nchar(padj)==2, paste0(padj, '0'), padj)
  return(padj)
}

addGlobalPadj <- function(dedf, globvec) {
  colno <- ncol(dedf)
  dedf$padjGlob <- globvec
  padjindex <- which(colnames(dedf) == 'padj')
  if (padjindex != colno) {
    dedf <- dedf[, c(1:padjindex, colno+1, (padjindex+1):colno)]
  }
  if ('padjGlob' %in% colnames(dedf)) {dedf[, 'padjGlob'] <- suppressWarnings(as.numeric(format(dedf[, 'padjGlob'], scientific = T, digits = 7)))}
  return(dedf)
}

# function removes samples which are marked as 'none' in selected columns of the data frame
# input
# levvec ... vector of levels
# coldatadf ... dataframe showing sample factor relation
removeNoneSampleFromColdata <- function(levvec, coldatadf, combinevalues) {
  if ('Combine' %in% levvec) {
    levvec <- union(intersect(c('Combine'), levvec), combinevalues)
  }
  for (lev in levvec) {
    levindex <- which(lev == colnames(coldatadf))
    coldatadf <- coldatadf[which('none' != coldatadf[, levindex]),]
  }
  coldatadf <- droplevels(coldatadf)
  return(coldatadf)
}

# small wrapper to produce a df for temp table; add counts and annotation and remove genes which have 0 counts
# also sort the table
# temp ... initial deseq results ob or a data.frame
# dds ... deseq object
# runvectemp ... vector having the factor levels and replacement names
# dosubset ... boolean filters rows which have NA in the log2FoldChange
# sortme ... boolean sorts table by padj, ascending
# return ... data.frame
prepResultCounttable <- function(temp, dds, runvectemp, annotable, globalpadj = '0', sortme = TRUE, nafold = TRUE, napadj = FALSE, napvalue = FALSE) {
  require(dplyr)
  require(DESeq2)
  temptable <- addCounts(temp, dds, cond=runvectemp[1], cond1=runvectemp[3], cond2=runvectemp[2], cond1name=runvectemp[5], cond2name=runvectemp[4])
  temptable$Ensembl_ID <- rownames(temptable)
  temptable <- merge(temptable, annotable, by.x = "Ensembl_ID", by.y = "Ensembl_ID")
  temptable <- temptable[, c('Ensembl_ID', 'Gene_Symbol', setdiff(colnames(temptable), c('Ensembl_ID', 'Gene_Symbol')))]
  if (globalpadj == '0') {
    if(nafold) {temptable <- subset(temptable, !(is.na(log2FoldChange)))}
    if(napadj) {temptable <- subset(temptable, !(is.na(padj)))}
    if(napvalue) {temptable <- subset(temptable, !(is.na(pvalue)))}
  }
  if(sortme) {temptable <- arrange(temptable, padj)}
  return(temptable)
}

# round all the given columns in a df by number
# padj and pvalue are not rounded anymore
# by default 9(8) which is padj is rounded to 6 positions after the . (via scientific printing)
# by default 8(7) which is p-value is rounded to 3 positions after the . (via scientific printing)
# input
# temp ... data frame
# number ... number of signs after the .
# return ... data frame
roundResults <- function(temp, number) {
  if ('baseMean' %in% colnames(temp)) {temp[, 'baseMean'] <- round(temp[, 'baseMean'], 1)}
  if ('log2FoldChange' %in% colnames(temp)) {temp[, 'log2FoldChange'] <- round(temp[, 'log2FoldChange'], number)}
  if ('lfcMLE' %in% colnames(temp)) {temp[, 'lfcMLE'] <- round(temp[, 'lfcMLE'], number)}
  if ('lfcSE' %in% colnames(temp)) {temp[, 'lfcSE'] <- round(temp[, 'lfcSE'], number)}
  # if ('stat' %in% colnames(temp)) {temp[, 'stat'] <- round(temp[, 'stat'], number)}
  if ('pvalue' %in% colnames(temp)) {temp[, 'pvalue'] <- suppressWarnings(as.numeric(format(temp[, 'pvalue'], scientific = T, digits = 7)))}
  if ('padj' %in% colnames(temp)) {temp[, 'padj'] <- suppressWarnings(as.numeric(format(temp[, 'padj'], scientific = T, digits = 7)))}
  if ('padjGlob' %in% colnames(temp)) {temp[, 'padjGlob'] <- suppressWarnings(as.numeric(format(temp[, 'padjGlob'], scientific = T, digits = 7)))}
  if ('weight' %in% colnames(temp)) {temp[, 'weight'] <- round(temp[, 'weight'], 3)}
  return(temp)
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

# check if the factors given in the formula are in your condition file
# input
# formula ... string which is an formula to calculate expression
# conditionnames ... names of the conditions
# return ... boolean
checkFormula <- function(formula, conditionnames) {
  require(magrittr)
  formulachanged <- getFormulavector(formula)
  if (sum(formulachanged %in% conditionnames) != length(formulachanged)) {
    warning(paste0("Check formula '", formula, "', has factors which are not column names in the condition table! skipping ...")); return(TRUE)
    return(TRUE)
  }
  return(FALSE)
}

# methods checks each row of the contrast file. it detects factors, conditions which are not condition file
# and prints a warning
# input
# contrasttemp ... subset of rows of the contrastfile
# index ... exact row in the subset
# coldatatemp ... the information from the condition file
# return ... 0 or vector of strings
checkFactorCorrectness <- function(contrasttemp, index, coldatatemp) {
  returnmsg <- TRUE
  fact <- as.character(contrasttemp[index, 'factor']) # factor for the comparison
  num <- as.character(contrasttemp[index, 'numerator']) # condition
  denom <- as.character(contrasttemp[index, 'denominator']) # condition to compare against
  # just warnings if someone didn't use the correct factors or conditions
  if(!(fact %in% colnames(coldatatemp))) {
    warning(paste0("'", fact, "' not a factor in the condition table! skipping ...'")); return(returnmsg)
  }
  if(!(num %in% unique(coldatatemp[, fact]))) {
    warning(paste0("'", num, "' not in '", fact), "' in the condition table! skipping ..."); return(returnmsg)
  }
  if(!(denom %in% unique(coldatatemp[, fact]))) {
    warning(paste0("'", denom, "' not in '", fact), "' in the condition table! skipping ..."); return(returnmsg)
  }
  formula <- getFormulavector(as.character(contrasttemp[index, 'formula'])) # factor for the comparison
  if (sum(fact %in% formula) == 0) {
    warning(paste0("'", fact, "' not part of formula '", contrasttemp[index, 'formula'], "'! skipping ...")); return(returnmsg)
  }
  return(FALSE)
}

# function reduces a result table to a columns of interest
# new df is only the mean counts of the two Conditions plus foldchange, unshrunken fold change (if available), pvalue and padj
# input
# df ... data frame
# resultcounter ... added to the columnames
# return ... data frame
getReducedCombineTable <- function(df, decomphelp, resultcounter) {
  require(magrittr)
  levelnames <- as.character(decomphelp[, 1]) %>% strsplit(' vs ') %>% unlist()
  if ('padjGlob' %in% colnames(df)) {
    newdf <- df[, c(levelnames[2], levelnames[1], 'log2FoldChange', 'pvalue', 'padj', 'padjGlob')]
  } else {
    newdf <- df[, c(levelnames[2], levelnames[1], 'log2FoldChange', 'pvalue', 'padj')]
  }
  colnames(newdf) <- paste(colnames(newdf), resultcounter, sep = "_")
  newdf$Ensembl_ID <- df$Ensembl_ID
  return(newdf)
}

# function builds a data frame out of many data frame by the given padj value
# annotation is added as well
# input
# padj ... padj value
# deresults ... list of data frames with results
# summary ... data frame which represents the deresults list
# return ... data frame
buildCombinedComparisionDataFrame <- function(padj, deresults, summary, annotable) {
  padj <- ifelse(padj<1, padj*100, padj)
  # select the data frames which have the required padj values
  indexentries <- subset(summary, FDR==padj)$Count
  # go through the selected data frames by using indexentries and reduce the frames
  allresults <- lapply(seq_along(indexentries), function(x) {getReducedCombineTable(deresults[[indexentries[x]]], summary[summary$Count == indexentries[x], ], x)})
  # merge them into a single data frame
  # not sure why by Ensembl_ID is magically the first column
  allresults <- Reduce(function(x, y) merge(x, y, by.x = 'Ensembl_ID', by.y = 'Ensembl_ID', all = T), allresults)
  # add annotation and order it
  allresults <- merge(allresults, annotable, by.x = "Ensembl_ID", by.y = "Ensembl_ID")
  allresults <- allresults[, c('Ensembl_ID', 'Gene_Symbol', setdiff(colnames(allresults), c('Ensembl_ID', 'Gene_Symbol')))]
  # allresults <- allresults[, c(1,ncol(allresults)-7, seq(2, ncol(allresults)-8), seq(ncol(allresults)-6, ncol(allresults)))]
  return(allresults)
}

# subset a table by FDR and add index counter
# input
# summary ... data frame which represents the deresults list
# padj ... padj value for subset
# tabindex ... is added for each row
# return ... data.frame
buildSubsetFromSummaryTable <- function(padj, summary, tabindex) {
  padj <- ifelse(padj<1, padj*100, padj)
  temp <- subset(summary, FDR == padj)
  temp$TabIndex <- tabindex
  temp$TableIndex <- seq(1:nrow(temp))
  temp <- temp[, c('TabIndex', 'TableIndex', 'Comparison', 'Factor', 'Formula', 'FDR', 'DEGenes', 'UpRegulated', 'DownRegulated', 'Comment')]
  return(temp)
}

# small function to print the summary of the results
# input
# df ... data frame
# filename ... filename to write to
printResultSummary <- function(df, filename){
  df[,'DEGenes'] <- prettyNum(df[,'DEGenes'], format="d", big.mark=',', decimal.mark = ".")
  df[,'UpRegulated'] <- prettyNum(df[,'UpRegulated'], format="d", big.mark=',', decimal.mark = ".")
  df[,'DownRegulated'] <- prettyNum(df[,'DownRegulated'], format="d", big.mark=',', decimal.mark = ".")
  write.table(df, filename, quote=F, sep="\t", row.names=F, col.names=T)
}


funcname2string <- function(v1) {
  return(deparse(substitute(v1)))
}



#clusterprofiler version
getEntrez <- function(x, from, to, db){
  require(clusterProfiler)
  return(bitr(x, fromType = from, toType = to, annoDb = db))
}

build_zscore <- function(inputdf){
  m = apply(inputdf, 1, mean, na.rm = T)
  s = apply(inputdf, 1, sd, na.rm = T)
  return((inputdf - m) / s)
}
