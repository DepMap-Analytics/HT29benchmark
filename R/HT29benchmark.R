

# exported

### documented
HT29R.download_ref_dataset<-function(whatToDownload='FCs',
                                     destFolder='./',
                                     dataRepoURL = "https://cog.sanger.ac.uk/cmp/downloads/crispr_cas9_benchmark/",
                                     expNames=c("HT29_c903","HT29_c904","HT29_c905","HT29_c906","HT29_c907","HT29_c908")){
  nexp<-length(expNames)

  for (i in 1:nexp){
    if (whatToDownload != 'rawCounts' &
        whatToDownload != 'FCs'){
      stop('Available data types (to be specified in whatToDownload): normalised sgRNA depletion fold-changes ("FCs") and rawCounts ("rawCounts")')
    }

    print(paste('Downloading data for ',expNames[i],' experiment (',i,' of ',nexp,')',sep=''))

    if (whatToDownload == 'FCs'){
      download.file(url = paste(dataRepoURL,'01_normalised_and_FCs/',expNames[i],'.tsv_foldChanges.Rdata',sep=''),
                    destfile = paste(destFolder,'/',expNames[i],'_foldChanges.Rdata',sep=''))
    }

    if (whatToDownload == 'rawCounts'){
      download.file(url = paste(dataRepoURL,'00_rawCounts%20assembled/',expNames[i],'.tsv',sep=''),
                    destfile = paste(destFolder,'/',expNames[i],'.tsv',sep=''))
    }
  }
}


### non documented


# non exported
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)

  if(r<0.68){COL='red'}
  else{COL='black'}

  text(0.5, 0.5, txt, cex = cex.cor * r,col=COL)
}


HT29R.replicateCorr<-function(refDataDir='./'){




}
