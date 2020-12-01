

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

  if(!geneLevel){
    pdf(paste(resDir,'/RepCor_Vs_PrScore.pdf',sep=''),6,13)
  }else{
    pdf(paste(resDir,'/GL_RepCor_Vs_PrScore.pdf',sep=''),6,13)
  }

  layout(matrix(c(1,1,2,3,4,5,6,7)))

  XLIM<-c(0.3,0.95)
  ccr.multDensPlot(toPlot,XLIMS = XLIM,TITLE='',COLS=c('gray','darkgreen','black'),
                   LEGentries = c('Project Score background',
                                  'Project Score replicates','Reproducibility Threshold'),XLAB = 'R')

  abline(v=sigTH)

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
HT29R.exp_similarity<-function(refDataDir='./',resDir='./',userFCs=NULL){

  data(KY_Library_v1.0)
  data(HT29R.prSCORE_bkgr_screen_similarity)
  bgCorr<-HT29R.prSCORE_bkgr_screen_similarity

  if (length(userFCs)>0){
    gg<-userFCs$sgRNA
    userFCs<-userFCs[,3:ncol(userFCs)]
    rownames(userFCs)<-gg
    userFCs<-apply(userFCs,MARGIN = 1,'mean')
    userFCs<-ccr.geneMeanFCs(userFCs,KY_Library_v1.0)
    cgenes<-names(userFCs)
  }else{
    cgenes<-unique(KY_Library_v1.0$GENES)
  }

  fn<-dir(refDataDir)
  fn<-grep('_foldChanges.Rdata',fn,value=TRUE)

  if (length(fn)==0){
    stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')
  }

  ref_fcs<-lapply(fn,function(x){
      load(paste(refDataDir,'/',x,sep=''))
      nr<-ncol(foldchanges)-2
      fc<-foldchanges[,3:ncol(foldchanges)]
      rownames(fc)<-foldchanges$sgRNA
      fc<-apply(fc,MARGIN = 1,FUN = 'mean')
      fc<-ccr.geneMeanFCs(fc,libraryAnnotation = KY_Library_v1.0)
      })

  cgenes<-intersect(cgenes,Reduce('intersect',lapply(ref_fcs,'names')))

  ref_fcs<-do.call(cbind,lapply(ref_fcs,function(x){x[cgenes]}))

  fn<-unlist(lapply(str_split(fn,'_foldChanges.Rdata'),function(x){x[1]}))

  colnames(ref_fcs)<-fn

  obsCorr<-c(as.dist(cor(ref_fcs)))

  toPlot<-list(Score_bg=density(bgCorr),
               RefScreen_sim=density(obsCorr))

  pdf(paste(resDir,'/Screen_sim.pdf',sep=''),5.23,7.03)

  layout(mat = c(1,1,2))

  ccr.multDensPlot(toPlot,XLIMS = c(0.7,1),TITLE='Screen similarity',COLS=c('gray','darkgreen'),
                                      LEGentries = c('Project Score background',
                                  'Ref.pair-wise HT29 screens'),XLAB = 'R')

  userCorrs<-NULL

  if(length(userFCs)>0){

    ref_fcs<-cbind(userFCs,ref_fcs)
    userCorrs<-cor(ref_fcs)
    userCorrs<-userCorrs[1,2:ncol(userCorrs)]
    points(userCorrs,rep(0,length(userCorrs)),cex=1.5,pch=21,
           bg=rgb(200,0,255,maxColorValue = 255,alpha = 120))
    legend('topleft',legend = 'User data',inset = c(0,0.10),
           pt.cex = 1.5,pch=21,pt.bg=rgb(200,0,255,maxColorValue = 255,alpha = 120),bty = 'n')


    colnames(ref_fcs)[1]<-'User data'

  }

  tres<-t.test(bgCorr,obsCorr)


  mtitle<-paste('bg Vs ref = ',format(tres$p.value,scientific = TRUE,digits = 2))

  if (length(userFCs)){
    tres<-t.test(obsCorr,userCorrs)
    mtitle<-paste(mtitle,', user Vs ref = ',format(tres$p.value,scientific = TRUE, digits= 2,sep=''))
  }

  boxplot(bgCorr,obsCorr,userCorrs,horizontal = TRUE,ylim=c(0.7,1),col=c('gray','darkgreen',rgb(200,0,255,maxColorValue = 255,alpha = 120)),
          main=mtitle)

  dev.off()

  pdf(paste(resDir,'/Sreen_pair_cor.pdf',sep=''),15,15)
  pairs(ref_fcs,lower.panel = panel.cor,upper.panel = my.panelSmooth)
  dev.off()
}
HT29R.PhenoIntensity<-function(refDataDir='./',resDir='./',userFCs=NULL, geneLevel=TRUE){

  fn<-dir(refDataDir)
  fn<-grep('_foldChanges.Rdata',fn,value=TRUE)

  if (length(fn)==0){
    stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')
  }

  pdf(paste(resDir,'/allScreens_PhenoIntensity.pdf',sep=''),11.88,9.71)
  par(mfrow=c(3,3))

  RES<-do.call(rbind,lapply(fn,function(x){
    load(paste(refDataDir,x,sep=''))
    currentFc<-rowMeans(foldchanges[,3:ncol(foldchanges)])
    names(currentFc)<-foldchanges$sgRNA
    HT29R.singleScreen_PhenoIntensity(currentFc,geneLevel = geneLevel,
                                      expName = unlist(lapply(str_split(x,'_foldChanges.Rdata'),
                                                              function(y){y[1]})))
    }))


  if(length(userFCs)>0){
    currentFc<-rowMeans(userFCs[,3:ncol(userFCs)])
    names(currentFc)<-userFCs$sgRNA
    usrRes<-HT29R.singleScreen_PhenoIntensity(currentFc,geneLevel = geneLevel,
                                      expName = 'User Data')
  }

  dev.off()

  pdf(paste(resDir,'/PhenoIntensity_reference.pdf',sep=''),3,7)

  par(mar=c(12,4,1,1))
  boxplot(list(RES$GD_ribProt,RES$GD_essential),
          names = c('Ribosomal protein genes','Other essential genes'),
          ylab='Glass Delta',col=c('darkblue','blue'),border = 'darkgrey',lwd=2,las=2)

  if (length(userFCs)>0){
    points(c(1,2),c(usrRes$GD_ribProt,usrRes$GD_essential),
           cex=2,col=rgb(200,0,255,maxColorValue = 255,alpha = 255),pch=16)

  }
  dev.off()
}
HT29R.singleScreen_PhenoIntensity<-function(FCprofile, geneLevel=TRUE,expName=NULL){
  data(EssGenes.ribosomalProteins)
  data(BAGEL_essential)
  data(BAGEL_nonEssential)

  data(KY_Library_v1.0)

  if(geneLevel){
    FCprofile<-ccr.geneMeanFCs(sgRNA_FCprofile = FCprofile,libraryAnnotation = KY_Library_v1.0)

  }else{
    EssGenes.ribosomalProteins<-rownames(KY_Library_v1.0)[match(EssGenes.ribosomalProteins,KY_Library_v1.0$GENES)]
    EssGenes.ribosomalProteins<-EssGenes.ribosomalProteins[!is.na(EssGenes.ribosomalProteins)]

    BAGEL_essential<-rownames(KY_Library_v1.0)[match(BAGEL_essential,KY_Library_v1.0$GENES)]
    BAGEL_essential<-BAGEL_essential[!is.na(BAGEL_essential)]

    BAGEL_nonEssential<-rownames(KY_Library_v1.0)[match(BAGEL_nonEssential,KY_Library_v1.0$GENES)]
    BAGEL_nonEssential<-BAGEL_nonEssential[!is.na(BAGEL_nonEssential)]

    }

  ribProtFC<-FCprofile[intersect(names(FCprofile),EssGenes.ribosomalProteins)]
  BagEssFC<-FCprofile[intersect(names(FCprofile),BAGEL_essential)]
  BagNonEssFC<-FCprofile[intersect(names(FCprofile),BAGEL_nonEssential)]

  toPlot<-list(NonEss=density(BagNonEssFC),
               Ess=density(BagEssFC),
               ribProt=density(ribProtFC))

  GD_essential<-abs(mean(BagEssFC)-mean(BagNonEssFC))/sd(BagEssFC)
  GD_ribProt<-abs(mean(ribProtFC)-mean(BagNonEssFC))/sd(ribProtFC)


  ccr.multDensPlot(toPlot,COLS=c('gray','blue','darkblue'),XLIMS = range(FCprofile),XLAB = 'depletion fold-change',
                   TITLE=paste(expName,'\nGD ess. = ',format(GD_essential,digits = 3),', GD rib.prot. = ',format(GD_ribProt,digits=3)),
                   LEGentries = c('Non Essential','Essential','Ribosomal Proteins'))
  abline(v=mean(BagNonEssFC),col='gray')
  abline(v=mean(BagEssFC),col='blue')
  abline(v=mean(ribProtFC),col='darkblue')

  res<-data.frame(GD_essential=GD_essential,GD_ribProt=GD_ribProt,stringsAsFactors = FALSE)
  rownames(res)<-expName

  return(res)

  }

HT29R.ROCanalysis<-function(refDataDir='./',resDir='./',positives,negatives,userFCs=NULL, geneLevel=TRUE){

  fn<-dir(refDataDir)
  fn<-grep('_foldChanges.Rdata',fn,value=TRUE)

  if (length(fn)==0){
    stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')
  }

  pdf(paste(resDir,'/allScreens_ROC_c.pdf',sep=''),18,6)
  par(mfrow=c(length(fn)/3,6))

  RES<-lapply(fn,function(x){
    load(paste(refDataDir,x,sep=''))
    currentFc<-rowMeans(foldchanges[,3:ncol(foldchanges)])
    names(currentFc)<-foldchanges$sgRNA
    ename<-str_split(x,'_foldChanges.Rdata')[[1]][1]
    HT29R.individualROC(currentFc,positives,negatives,geneLevel = FALSE,expName=ename)
  })

  dev.off()

  if(length(userFCs)>0){
     currentFc<-rowMeans(userFCs[,3:ncol(userFCs)])
     names(currentFc)<-userFCs$sgRNA
     if (geneLevel){
       data(KY_Library_v1.0)
       currentFc<-ccr.geneMeanFCs(currentFc,KY_Library_v1.0)
     }
     pdf(paste(resDir,'/UserData_ROC_c.pdf',sep=''),10,5)
     par(mfrow=c(1,2))
     usrRes<-HT29R.individualROC(currentFc,positives,negatives,geneLevel = FALSE,expName='User Data')
     dev.off()
  }


  plot(0,0,xlim=c(1,0),ylim=c(0,1),col=NA,xlab='TNR',ylab='Recall')
  abline(1,-1,col='darkgray')
  lapply(1:length(RES),function(x){
    lines(RES[[x]]$ROCres$curve[,1],RES[[x]]$ROCres$curve[,2],lwd=5,col=rgb(0,0,255,alpha = 100,maxColorValue = 255))
    abline(h=RES[[x]]$ROCres$Recall,lty=2,col=rgb(0,0,255,alpha = 100,maxColorValue = 255),lwd=2)
  })

  if(length(userFCs)>0){
    lines(usrRes$ROCres$curve[,1],usrRes$ROCres$curve[,2],lwd=2,col=rgb(255,0,255,alpha = 255,maxColorValue = 255))
    abline(h=usrRes$ROCres$Recall,lty=2,col=rgb(255,0,255,alpha = 255,maxColorValue = 255),lwd=2)
  }



  pdf(paste(resDir,'/PhenoIntensity_reference.pdf',sep=''),4,3)

  par(mar=c(12,4,1,1))
  boxplot(list(RES$GD_ribProt,RES$GD_essential),
          names = c('Ribosomal protein genes','Other essential genes'),
          ylab='Glass Delta',col=c('darkblue','blue'),border = 'darkgrey',lwd=2,las=2)

  if (length(userFCs)>0){
    points(c(1,2),c(usrRes$GD_ribProt,usrRes$GD_essential),
           cex=2,col=rgb(200,0,255,maxColorValue = 255,alpha = 255),pch=16)

  }
  dev.off()
}

HT29R.individualROC<-function(FCprofile, positives,negatives,geneLevel=TRUE,expName=NULL){

  data(KY_Library_v1.0)

  if(geneLevel){
    FCprofile<-ccr.geneMeanFCs(sgRNA_FCprofile = FCprofile,libraryAnnotation = KY_Library_v1.0)

  }else{

    positives<-rownames(KY_Library_v1.0)[match(positives,rownames(KY_Library_v1.0))]
    positives<-positives[!is.na(positives)]

    negatives<-rownames(KY_Library_v1.0)[match(negatives,rownames(KY_Library_v1.0))]
    negatives<-negatives[!is.na(negatives)]

  }

  ROCres<-ccr.ROC_Curve(FCprofile,positives = positives,negatives = negatives,FDRth = 0.05,expName = paste(expName,'ROC'))
  PRRCres<-ccr.PrRc_Curve(FCprofile,positives = positives,negatives = negatives,FDRth = 0.05,expName = paste(expName,'PrRc'))

  return(list(ROCres=ROCres,PRRCres=PRRCres))

}




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
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

my.panelSmooth<-function (x, y, col = 'blue', bg = NA, pch = 16,
          cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)
{
  par(new = TRUE)

  smoothScatter(x, y, pch = pch, col = col, colramp = colorRampPalette(c('white',col)), cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
          col = col.smooth, lwd=3,...)
}

makeTransparent <- function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}







