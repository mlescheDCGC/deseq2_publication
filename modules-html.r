# function creates the basic structure for the html report
# input
# dedir ... working directory for the de analysis
# bfxid ... string of the bfx id
# decompdf ... data frame with the results of the single comparisons
# whereiswhatdf ... data frame with locations of pdfs, pngs, csv, xlsx
setuppublishHTML <- function(dedir, logodir, bfxid, decompdf = '', whereiswhatdf = '') {
  require(R.utils)
  reportdir <- file.path(dedir, 'de_report')
  folders <- c(reportdir, file.path(reportdir, 'picture'), file.path(reportdir, 'data'),
               file.path(reportdir, 'pdf'), file.path(reportdir, 'report'), file.path(reportdir, 'data', 'gsea'))
  if(dir.exists(folders[1])) {
    unlink(folders[1], recursive = T)
  }
  invisible(sapply(folders, mkdirectory))
  # copy logos
  file.copy(file.path(logodir, 'CRTD/crtd-farbig.png'), file.path(folders[2], 'crtd.png'))
  file.copy(file.path(logodir, 'BIOTEC/biotec-farbig.png'), file.path(folders[2], 'biotec.png'))
  file.copy(file.path(logodir, 'dcgc/Logo_DDC_GenomeCenter.png'), file.path(folders[2], 'group.png'))
  file.copy(file.path(logodir, 'joint_technology_platform/jtp-farbig.png'), file.path(folders[2], 'jtp.png'))
  
  # copy ma and volcano plots
  if (nrow(decompdf)>0) {
    # apply(decompdf, 1, function(x) {gzip(as.character(x[9]), file.path(folders[3], basename(paste0(as.character(x[9]), '.gz'))), remove = F)})
    invisible(file.copy(as.character(decompdf$MA), folders[2]))
    invisible(file.copy(gsub("png", "pdf", gsub("picture", "pdf", as.character(decompdf$MA))), folders[4]))
    invisible(file.copy(as.character(decompdf$Volcano), folders[2]))
    invisible(file.copy(gsub("png", "pdf", gsub("picture", "pdf", as.character(decompdf$Volcano))), folders[4]))
  }
  if (is.data.frame(whereiswhatdf)) {
    if (nrow(subset(whereiswhatdf, what == 'countstable')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'countstable')
      file.copy(as.character(tmpdf$xls), folders[3])
      # cls files
      cls.files <- list.files(file.path(dirname(tmpdf$xls), 'gsea'), "*.cls", full.names=T)
      file.copy(cls.files, folders[6])
      # gct files
      gct.files <- list.files(file.path(dirname(tmpdf$xls), 'gsea'), "*.gct", full.names=T)
      file.copy(gct.files, folders[6])
    }
    if (nrow(subset(whereiswhatdf, what == 'batchvsd')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'batchvsd')
      file.copy(as.character(tmpdf$xls), folders[3])
    }
    if (nrow(subset(whereiswhatdf, what == 'combineresultsDE')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'combineresultsDE')
      file.copy(as.character(tmpdf$xls), folders[3])
    }
    if (nrow(subset(whereiswhatdf, what == 'morpheus')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'morpheus')
      file.copy(as.character(tmpdf$xls), folders[3])
    }
    if (nrow(subset(whereiswhatdf, what == 'allresultsDE')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'allresultsDE')
      file.copy(as.character(tmpdf$xls), folders[3])
      # rnk files
      # rnk.files <- list.files(file.path(dirname(tmpdf$xls), 'gsea'), "*.rnk", full.names=T)
      # file.copy(rnk.files, folders[6])
    }
    if (nrow(subset(whereiswhatdf, what == 'pca12')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'pca12')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'batchpca12')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'batchpca12')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'biplot')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'biplot')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'batchbiplot')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'batchbiplot')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'pca23')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'pca23')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'mds')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'mds')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'poisson')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'poisson')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'heatmap')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'heatmap')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'batchheatmap')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'batchheatmap')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'heatmappoisson')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'heatmappoisson')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'pearson')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'pearson')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'spearman')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'spearman')
      file.copy(as.character(tmpdf$png), folders[2])
      file.copy(as.character(tmpdf$pdf), folders[4])
    }
    if (nrow(subset(whereiswhatdf, what == 'varGenes')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'varGenes')
      tmpfiles <- unlist(strsplit(as.character(tmpdf$png) , ','))
      invisible(sapply(tmpfiles, function(x){file.copy(x, folders[2])}))
      tmpfiles <- unlist(strsplit(as.character(tmpdf$pdf) , ','))
      invisible(sapply(tmpfiles, function(x){file.copy(x, folders[4])}))
    }
    if (nrow(subset(whereiswhatdf, what == 'batchvarGenes')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'batchvarGenes')
      tmpfiles <- unlist(strsplit(as.character(tmpdf$png) , ','))
      invisible(sapply(tmpfiles, function(x){file.copy(x, folders[2])}))
      tmpfiles <- unlist(strsplit(as.character(tmpdf$pdf) , ','))
      invisible(sapply(tmpfiles, function(x){file.copy(x, folders[4])}))
    }
    if (nrow(subset(whereiswhatdf, what == 'highExpressedGenes')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'highExpressedGenes')
      tmpfiles <- unlist(strsplit(as.character(tmpdf$png) , ','))
      invisible(sapply(tmpfiles, function(x){file.copy(x, folders[2])}))
      tmpfiles <- unlist(strsplit(as.character(tmpdf$pdf) , ','))
      invisible(sapply(tmpfiles, function(x){file.copy(x, folders[4])}))
    }
    if (nrow(subset(whereiswhatdf, what == 'normalisedGenes')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'normalisedGenes')
      tmpfiles <- unlist(strsplit(as.character(tmpdf$png) , ','))
      invisible(sapply(tmpfiles, function(x){file.copy(x, folders[2])}))
      tmpfiles <- unlist(strsplit(as.character(tmpdf$pdf) , ','))
      invisible(sapply(tmpfiles, function(x){file.copy(x, folders[4])}))
    }
  }
}

# basic function to produce a single html report page
# it plots a data frame, adds title, and if available other plots and ma plot
# df ... data frame
# wd ... working directory
# shortname ... name of the html file
# reporttitle ... header of the html file
# ma ... filename of a ma plot
# addedf ... list of data frames
# return ... site in html reports class
# buildReportPage <- function(df, wd, shortname, reporttitle, maplot = '', addeddf = NA, addedLines = c()) {
#   require(ReportingTools)
#   require(hwriter)
#   reportname <- 'de_report'
#   reportfolder <- 'report'
#   picdir <- 'picture'
#   tempreport <- HTMLReport(shortName = shortname, title = '', reportDirectory = file.path(reportname, reportfolder), file.path(wd))
#   crtdlogo <- hwriteImage(file.path('..', picdir, 'crtd.png'), width=150)
#   bioteclogo <- hwriteImage(file.path('..', picdir, 'biotec.png'), width=150)
#   grouplogo <- hwriteImage(file.path('..', picdir, 'group.jpg'), width=120)
#   publish(hwrite(c(crtdlogo, bioteclogo,grouplogo), br=T, center=T, border=0), tempreport)
#   publish(hwrite(reporttitle, heading=2, center = T), tempreport)
#   if(maplot!='') {
#     ma <- hwriteImage(file.path('..', 'picture', basename(maplot)), width = 500)
#     publish(hwrite(ma, br = TRUE, center=T, border=0), tempreport)
#   }
#   if(length(addedLines)!=0) {
#     sapply(addedLines, function(x) {publish(hwrite(x, br = T, border=0), tempreport)})
#   }
#   if(is.list(addeddf)) {
#     sapply(addeddf, publish, tempreport)
#   }
#   publish(df, tempreport)
#   finish(tempreport)
#   return(tempreport)
# }

# function formats the appendix to add it to the title for the report page
# input
# appendix ... string
# return ... string
getAdjAppendixString <- function(appendix) {
  if (grepl('^remove', appendix)) {
    tempvec <- unlist(strsplit(appendix, '_'))
    stringback <- paste0('(', paste(tempvec[2:length(tempvec)], collapse = ', '), ' removed)')
  } else if(grepl('^add', appendix)) {
    tempvec <- unlist(strsplit(appendix, '_'))
    stringback <- paste0('(', paste(tempvec[2:length(tempvec)], collapse = ', '), ' added)')
  } else if(grepl('^only', appendix)){
    tempvec <- unlist(strsplit(appendix, '_'))
    stringback <- paste0('(only ', paste(tempvec[2:length(tempvec)], collapse = ', '), ')') 
  } else if(grepl('^comb', appendix)){
    tempvec <- unlist(strsplit(appendix, '_'))
    stringback <- paste0('(', paste(tempvec[2:length(tempvec)], collapse = ', '), ' combined)') 
  } else {
    stringback <- paste0('(', gsub('_', ' ', appendix), ')')
  }
  return(stringback)
}


# complete method to build the html report
# input
# dedir ... working directory
# bfxid .. bfx id
# appendix ... appendix for the documents
# decompdf ... summary table for all single comparisons
# whereiswhatdf ... data frame with locations of pdfs, pngs, csv, xlsx
publishHTML <- function(dedir, bfxid, decompdf, whereiswhatdf, trans.method) {
  require(ReportingTools)
  require(hwriter)
  require(magrittr)
    
  picdir <- 'picture'; datadir <- 'data'; reportname <- 'de_report'
  grouplogo <- hwriteImage(file.path(picdir, 'group.png'), width=210)
  jtplogo <- hwriteImage(file.path(picdir, 'jtp.png'), width=200)
  crtdlogo <- hwriteImage(file.path(picdir, 'crtd.png'), width=185)
  bioteclogo <- hwriteImage(file.path(picdir, 'biotec.png'), width=185)

  
  indexpage <- HTMLReport(shortName = 'report', title = '', reportDirectory = reportname, basePath = dedir)
  publish(hwrite(c(grouplogo, jtplogo, crtdlogo, bioteclogo), br=T, center = T, border = 0), indexpage)
  publish(hwrite(paste0('RNA-seq analysis report for Project ', bfxid), heading=2, center = T), indexpage)
  # as overview the pca is show and if that's not avaiable fall back to heatmap
  if (nrow(subset(whereiswhatdf, what == 'pca12'))==1) {
    tmpdf <- subset(whereiswhatdf, what == 'pca12')
    preimg  <- hwriteImage(file.path(picdir, basename(as.character(tmpdf$png))), width = 700)
    publish(hwrite(c(preimg), br=T, center = T, border = 0), indexpage)
  } else if (nrow(subset(whereiswhatdf, what == 'heatmap')) == 1) {
    tmpdf <- subset(whereiswhatdf, what == 'heatmap')
    preimg  <- hwriteImage(file.path(picdir, basename(as.character(tmpdf$png))), width = 700)
    publish(hwrite(c(preimg), br=T, center = T, border = 0), indexpage)
  }
  # publish the single reports and attach download links
  if (nrow(decompdf)>0) {
    publish(hwrite('Overview Reports', heading=3, center = T), indexpage)
    decompdf$Download <- apply(decompdf, 1, function(x) {hwrite('Download', link = file.path(datadir, paste0(basename(as.character(x[9])), '.gz')))})
    # publish(decompdf[, c('Comparison', 'Factor', 'Formula', 'FDR', 'DEGenes', 'UpRegulated', 'DownRegulated', 'Comment', 'Download')], indexpage)
    decompdf$MA <- hwrite('Image', link = file.path(picdir, basename(as.character(decompdf$MA))))
    decompdf$Volcano <- hwrite('Image', link = file.path(picdir, basename(as.character(decompdf$Volcano))))
	
	# AP: fixes a strange BUG which causes incorrect numbers
	decompdf$FDR <- as.character(decompdf$FDR)
	decompdf$DEGenes <- as.character(decompdf$DEGenes)
	decompdf$UpRegulated <- as.character(decompdf$UpRegulated)
	decompdf$DownRegulated <- as.character(decompdf$DownRegulated)
    publish(decompdf[, c('Comparison', 'Factor', 'Formula', 'FDR', 'DEGenes', 'UpRegulated', 'DownRegulated', 'MA', 'Volcano', 'Comment')], indexpage)
  }
  publish(hwrite('Overview Downloads', heading=3, center = T), indexpage)
  publish(hwrite("Download the report (as pdf)", link = file.path("report", paste0(bfxid, ".pdf"))),indexpage)
  if (is.data.frame(whereiswhatdf)) {
    if (nrow(subset(whereiswhatdf, what == 'allresultsDE')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'allresultsDE')
      publish(hwrite("Download all analyses in separate sheets (as xlsx)", link = file.path(datadir, basename(as.character(tmpdf$xls)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'combineresultsDE')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'combineresultsDE')
      publish(hwrite("Download combined analyses (as xlsx)", link = file.path(datadir, basename(as.character(tmpdf$xls)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'countstable')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'countstable')
      # publish(hwrite("Download normalised counts (as csv)", link = file.path(datadir, paste0(basename(as.character(tmpdf$csv)), ".gz"))), indexpage)
      if (tmpdf$png == '0') {
        input <- paste0("Download count tables (raw counts, normalised counts, pca-data, ", trans.method, "-values [similar to log2 scale], z-score (from vst-values)) (as xlsx)")
        publish(hwrite(input, link = file.path(datadir, basename(as.character(tmpdf$xls)))), indexpage)
      } else if (tmpdf$png == '1') {
        input <- paste0("Download count tables (raw counts, normalised counts, pca-data, ", trans.method, "-values [similar to log2 scale], z-score (from vst-values), ", trans.method, "-values [after batch correction], z-score [after batch correction]) (as xlsx)")
        publish(hwrite(input, link = file.path(datadir, basename(as.character(tmpdf$xls)))), indexpage)
      }
      
    }
    if (nrow(subset(whereiswhatdf, what == 'pca12')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'pca12')
      publish(hwrite("Show PCA (1-2)", link = file.path(picdir, basename(as.character(tmpdf$png)))),indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'batchpca12')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'batchpca12')
      publish(hwrite("Show PCA (1-2) after batch correction", link = file.path(picdir, basename(as.character(tmpdf$png)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'pca23')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'pca23')
      publish(hwrite("Show PCA (2-3)", link = file.path(picdir, basename(as.character(tmpdf$png)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'biplot')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'biplot')
      publish(hwrite("Show PCA - BiPlot (The further a vector is away from PC origin, the more contribution it has)", link = file.path(picdir, basename(as.character(tmpdf$png)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'batchbiplot')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'batchbiplot')
      publish(hwrite("Show PCA - BiPlot after batch correction", link = file.path(picdir, basename(as.character(tmpdf$png)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'mds')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'mds')
      publish(hwrite("Show MDS plot of Euclidean distance", link = file.path(picdir, basename(as.character(tmpdf$png)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'poisson')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'poisson')
      publish(hwrite("Show MDS plot of Poisson Distance", link = file.path(picdir, basename(as.character(tmpdf$png)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'heatmap')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'heatmap')
      publish(hwrite("Show Heatmap", link = file.path(picdir, basename(as.character(tmpdf$png)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'batchheatmap')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'batchheatmap')
      publish(hwrite("Show Heatmap after batch correction", link = file.path(picdir, basename(as.character(tmpdf$png)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'heatmappoisson')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'heatmappoisson')
      publish(hwrite("Show Heatmap of Poisson Distance", link = file.path(picdir, basename(as.character(tmpdf$png)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'pearson')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'pearson')
      publish(hwrite("Show Pearson correlation", link = file.path(picdir, basename(as.character(tmpdf$png)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'spearman')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'spearman')
      publish(hwrite("Show Spearman correlation", link = file.path(picdir, basename(as.character(tmpdf$png)))), indexpage)
    }
    if (nrow(subset(whereiswhatdf, what == 'normalisedGenes')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'normalisedGenes')
      tmpfiles <- as.character(tmpdf$png) %>% strsplit(',') %>% unlist() %>% basename()
      sapply(tmpfiles, function(x){
        genes <- unlist(strsplit(unlist(strsplit(x, '\\.'))[2], '_'))[2]
        publish(hwrite(paste0("Normalised counts plot for ", genes), link = file.path(picdir, basename(as.character(x)))), indexpage)
      })
    }
    if (nrow(subset(whereiswhatdf, what == 'varGenes')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'varGenes')
      tmpfiles <- as.character(tmpdf$png) %>% strsplit(',') %>% unlist() %>% basename()
      sapply(tmpfiles, function(x){
        genes <- unlist(strsplit(unlist(strsplit(x, '\\.'))[2], '_'))[2]
        publish(hwrite(paste0("Show Deviation from gene's average accross all samples for top ", genes, " genes (having highest variance)"), link = file.path(picdir, basename(as.character(x)))), indexpage)
      })
    }
    if (nrow(subset(whereiswhatdf, what == 'batchvarGenes')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'batchvarGenes')
      tmpfiles <- as.character(tmpdf$png) %>% strsplit(',') %>% unlist() %>% basename()
      sapply(tmpfiles, function(x){
        genes <- unlist(strsplit(unlist(strsplit(x, '\\.'))[2], '_'))[3]
        publish(hwrite(paste0("Show Deviation from gene's average accross all samples for top ", genes, " genes (having highest variance) (after batch correction)"), link = file.path(picdir, basename(as.character(x)))), indexpage)
      })
    }
    if (nrow(subset(whereiswhatdf, what == 'highExpressedGenes')) == 1) {
      tmpdf <- subset(whereiswhatdf, what == 'highExpressedGenes')
      tmpfiles <- as.character(tmpdf$png) %>% strsplit(',') %>% unlist() %>% basename()
      sapply(tmpfiles, function(x){
        genes <- unlist(strsplit(unlist(strsplit(x, '\\.'))[2], '_'))[2]
        publish(hwrite(paste0("Show top ", genes, " genes (having highest expression)"), link = file.path(picdir, basename(as.character(x)))), indexpage)
      })
    }
  }
  publish(hwrite(c('If you have questions regarding the analysis please write an ', hwrite('email',link = paste0('mailto: genomecenter@tu-dresden.de?subject=analysis ',bfxid))), center = T, border = 0, table = T, cellpadding = 2), indexpage)
  # publish(hwrite(c('Download', hwrite('7-Zip', link = 'http://www.7-zip.org/', target='_blank'), 'for unpacking the csv files'), center = T, border = 0, table = T, cellpadding = 2, br = T), indexpage)
  invisible(finish(indexpage))
}
