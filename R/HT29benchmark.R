

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

HT29R.replicateCorr_Pscore<-function(refDataDir='./',resDir='./',userFCs=NULL){


  fn<-dir(refDataDir)
  fn<-grep('_foldChanges.Rdata',fn,value=TRUE)

  if (length(fn)==0){
    stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')
  }

  data(HT29R.prSCORE_rCorr_Reprod)
  data(HT29R.reproducible_GeneGuides)

  toPlot<-list(Pr.Score_bg=density(HT29R.prSCORE_rCorr_Reprod$BGscores),
               Pr.Score_repCor=density(HT29R.prSCORE_rCorr_Reprod$REPscores))


  pdf(paste(resDir,'/RepCor_Vs_PrScore.pdf',sep=''),5,15)
  layout(matrix(c(1,1,2,3,4,5,6,7)))

  XLIM<-c(0.3,0.95)
  ccr.multDensPlot(toPlot,XLIMS = XLIM,TITLE='',COLS=c('gray','darkgreen'),
                   LEGentries = c('Project Score background',
                                  'Project Score replicates'),XLAB = 'R')

  if (length(userFCs)>0){
    fc<-userFCs
    nr<-ncol(fc)-2

    fc<-fc[,3:ncol(fc)]

    colnames(fc)<-paste('UserData_R',1:ncol(fc),sep='')

    rownames(fc)<-userFCs$sgRNA

    cc<-c(as.dist(cor(fc[HT29R.reproducible_GeneGuides,])))

    points(cc,rep(0,length(cc)),cex=3,pch=21,
           bg=rgb(0,0,255,maxColorValue = 255,alpha = 120))

    abline(h=0,col='gray')

    legend('topleft',legend = 'User data',inset = c(0,0.15),pt.cex = 3,pch=21,pt.bg=rgb(0,0,255,maxColorValue = 255,alpha = 120),bty = 'n')

    print(as.dist(cor(fc[HT29R.reproducible_GeneGuides,])))

    vv<-c(as.dist(cor(fc[HT29R.reproducible_GeneGuides,])))

    print(paste(length(which(vv>=0.68)),' pair-wise replicate comparisons (out of ',length(vv),') yield correlation scores greater or equal than the QC threshold',sep=''))

  }

  abline(v=0.68)

  lapply(fn,function(x){
    load(paste(refDataDir,'/',x,sep=''))

    nr<-ncol(foldchanges)-2
    fc<-foldchanges[,3:ncol(foldchanges)]
    rownames(fc)<-foldchanges$sgRNA

    cc<-c(as.dist(cor(fc[HT29R.reproducible_GeneGuides,])))

    nn<-str_split(x,'_foldChanges.Rdata')[[1]][1]
    plot(cc,rep(0,length(cc)),xlim=c(0.3,0.95),cex=3,pch=21,
         bg=rgb(0,0,255,maxColorValue = 255,alpha = 120),main=nn,frame.plot = FALSE,
         xaxt='n',yaxt='n',xlab='',ylab='',ylim=c(-1,1),col='white')
    abline(h=0,col='gray')
  }
  )

  dev.off()
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



