

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

HT29R.evaluate_reps<-function(refDataDir='./',resDir='./',userFCs=NULL, geneLevel=TRUE){

  fn<-dir(refDataDir)
  fn<-grep('_foldChanges.Rdata',fn,value=TRUE)

  if (length(fn)==0){
    stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')
  }

  if(!geneLevel){
    data(HT29R.prSCORE_rCorr)
    ScoreCorrs<-HT29R.prSCORE_rCorr
    sigTH<-0.55
  }else{
    data(HT29R.GL_prSCORE_rCorr)
    ScoreCorrs<-HT29R.GL_prSCORE_rCorr
    sigTH<-0.68
    data(KY_Library_v1.0)
  }

  data(HT29R.reproducible_GeneGuides)

  toPlot<-list(Pr.Score_bg=density(ScoreCorrs$BGscores),
               Pr.Score_repCor=density(ScoreCorrs$REPscores))

  pdf(paste(resDir,'/RepCor_Vs_PrScore.pdf',sep=''),6,13)
  layout(matrix(c(1,1,2,3,4,5,6,7)))

  XLIM<-c(0.3,0.95)
  ccr.multDensPlot(toPlot,XLIMS = XLIM,TITLE='',COLS=c('gray','darkgreen','black'),
                   LEGentries = c('Project Score background',
                                  'Project Score replicates','Reproducibility Threshold'),XLAB = 'R')

  if (length(userFCs)>0){
    fc<-userFCs
    nr<-ncol(fc)-2

    fc<-fc[,3:ncol(fc)]

    colnames(fc)<-paste('UserData_R',1:ncol(fc),sep='')

    rownames(fc)<-userFCs$sgRNA

    dataSET<-fc[HT29R.reproducible_GeneGuides,]

    if(geneLevel){
      dataSET<-apply(dataSET,MARGIN = 2,function(x){
          cdata<-x
          names(cdata)<-rownames(dataSET)
          ccr.geneMeanFCs(cdata,KY_Library_v1.0)
        })
      }


    cc<-c(as.dist(cor(dataSET)))

    points(cc,rep(0,length(cc)),cex=3,pch=21,
           bg=rgb(200,0,255,maxColorValue = 255,alpha = 120))

    abline(h=0,col='gray')

    if(geneLevel){
      abline(v=0.68,lwd=2)
    }else{
      abline(v=0.55,lwd=2)
    }

    legend('topleft',legend = 'User data',inset = c(0,0.25),pt.cex = 3,pch=21,pt.bg=rgb(200,0,255,maxColorValue = 255,alpha = 120),bty = 'n')


    print(as.dist(cor(dataSET)))

    vv<-c(as.dist(cor(dataSET)))

    cat(blue(paste(length(which(vv>=0.68)),' pair-wise replicate comparisons (out of ',length(vv),
              ') yield correlation scores greater or equal than Project Score QC threshold\n',sep='')))

  }


  lapply(fn,function(x){
    load(paste(refDataDir,'/',x,sep=''))

    nr<-ncol(foldchanges)-2
    fc<-foldchanges[,3:ncol(foldchanges)]
    rownames(fc)<-foldchanges$sgRNA

    dataSET<-fc[HT29R.reproducible_GeneGuides,]

    if(geneLevel){
      dataSET<-apply(dataSET,MARGIN = 2,function(x){
        cdata<-x
        names(cdata)<-rownames(dataSET)
        ccr.geneMeanFCs(cdata,KY_Library_v1.0)
      })
    }

    cc<-c(as.dist(cor(dataSET)))

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
HT29R.FC_dist_properties<-function(refDataDir='./',resDir='./',userFCs=NULL){

  fn<-dir(refDataDir)
  fn<-grep('_foldChanges.Rdata',fn,value=TRUE)

  if (length(fn)==0){
    stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')
  }

  nf<-length(fn)

  n <- nf
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  COL<-col_vector[1:n]
  names(COL)<-fn

  data(KY_Library_v1.0)
  cguides<-rownames(KY_Library_v1.0)

  all_sgRNAs<-lapply(fn,function(x){
    load(paste(refDataDir,'/',x,sep=''))
    return(foldchanges$sgRNA)
  })

  cguides<-Reduce(intersect, all_sgRNAs)

  if(length(userFCs)>1){
    cguides<-intersect(cguides,userFCs$sgRNA)
  }

  firstEl<-fn[1]

  par(mfrow=c(1,2))

  par(mar=c(5,4,2,0))
  allAvgFc<-
    lapply(fn,function(x){
      cat(paste('loading',x,'\n'))
      load(paste(refDataDir,'/',x,sep=''))
      nr<-ncol(foldchanges)-2

      fc<-foldchanges[match(cguides,foldchanges$sgRNA),3:ncol(foldchanges)]
      rownames(fc)<-cguides
      nreps<-ncol(fc)-2

      currentCol<-makeTransparent(COL[x])

      for (j in 1:nreps){
        if(x==firstEl & j==1){
          plot(density(fc[,j]),frame.plot=FALSE,xlim=c(-9,3),lwd=7,col=currentCol,main='',xlab='sgRNA fc')

        }else{
          par(new=TRUE)
          plot(density(fc[,j]),frame.plot=FALSE,xlim=c(-9,3),lwd=7,xaxt='n',yaxt='n',xlab='',ylab='',col=currentCol,
               main='')
        }
      }
      fc<-rowMeans(fc)
      return(fc)
    }
    )

  if(length(userFCs)>0){
    nr<-ncol(userFCs)-2

    for (i in 1:nr){
      par(new=TRUE)

      plot(density(userFCs[,2+i]),frame.plot=FALSE,xlim=c(-9,3),lwd=1,xaxt='n',yaxt='n',xlab='',ylab='',
           main='',col='black')
    }

    legend('topleft','User data',lwd = 1)
  }

  GlobalFC<-do.call(cbind,allAvgFc)
  colnames(GlobalFC)<-unlist(lapply(str_split(fn,'_foldChanges.Rdata'),function(x){x[[1]][1]}))

  if(length(userFCs)>0){
    nr<-ncol(userFCs)-2
    GlobalFC<-cbind(GlobalFC,rowMeans(userFCs[match(cguides,userFCs$sgRNA),3:(2+nr)]))
    colnames(GlobalFC)[ncol(GlobalFC)]<-'User data'
    COL<-c(COL,'black')
    nf<-nf+1
    NC<-ncol(GlobalFC)-1
  }else{
    NC<-ncol(GlobalFC)
  }

  par(mar=c(5,4,0,0))

  tmp<-boxplot(GlobalFC,col=makeTransparent(COL),frame.plot=FALSE,ylab='sgRNA Avg FC',pch=16,
               pars=list(outcol=makeTransparent(COL)),las=2,cex.axis=0.9,xaxt='n')

  axis(1, at=1:nf, labels=FALSE)
  text(x=1:nf, y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
       labels=colnames(GlobalFC), srt=45, adj=1, xpd=TRUE)

  abline(h=median(tmp$stats[1,1:NC]),col='darkgray',lty=2)
  abline(h=median(tmp$stats[5,1:NC]),col='darkgray',lty=2)

  abline(h=median(tmp$stats[2,1:NC]),col='darkgray')
  abline(h=median(tmp$stats[4,1:NC]),col='darkgray')

  abline(h=median(tmp$stats[3,1:NC]),lwd=3,col='darkgray')

  abline(h=median(apply(GlobalFC[,1:NC],2,'min')),col='darkgray')
  abline(h=median(apply(GlobalFC[,1:NC],2,'max')),col='darkgray')

  cat('\n\nReference average sgRNA FC range:\n')

  ranges<-apply(GlobalFC[,1:NC],2,range)
  avgRanges<-apply(ranges[,1:NC],1,mean)


  SE<-round(apply(ranges,1,sd),digits = 2)

  cat(blue(c(paste(c('min = ','max = '),round(apply(ranges,1,'mean'),digits=2),' ± ',SE,c(',',''),sep=''))))

  if(length(userFCs)>0){
    usFC<-GlobalFC[,ncol(GlobalFC)]

    cat('\n\nUser sgRNA FC range:\n')


  }




  # print('range se')
  # sd(apply(GlobalFC,2,range)[1,], na.rm=TRUE)/sqrt(6)
  # sd(apply(GlobalFC,2,range)[2,], na.rm=TRUE)/sqrt(6)
  #
  # print('average median')
  # mean(tmp$stats[3,])
  # print('se')
  # sd(tmp$stats[3,])/sqrt(6)
  #
  # print('average inerquartile range')
  # apply(tmp$stats[c(2,4),],1,mean)
  #
  # print('range se')
  # apply(tmp$stats[c(2,4),],1,sd)/sqrt(6)
  #
  # print('average 10-90 perc. range')
  # apply(apply(GlobalFC,2,quantile,c(0.10,0.90)),1,'mean')
  #
  # print('range se')
  # apply(apply(GlobalFC,2,quantile,c(0.10,0.90)),1,'sd')/sqrt(6)
  #
  #
  # print('range se')
  # apply(tmp$stats[c(2,4),],1,sd)/sqrt(6)
  #
  # print('average skeweness')
  # mean(apply(GlobalFC,2,skewness))
  #
  # print('range se')
  # sd(apply(GlobalFC,2,skewness))/sqrt(6)
  #
  # print('average kurtosis')
  # mean(apply(GlobalFC,2,kurtosis))
  #
  # print('range se')
  # sd(apply(GlobalFC,2,kurtosis))/sqrt(6)



}

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
makeTransparent <- function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}







