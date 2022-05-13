# function definition
plot.DensityPrjSCORE <- function(refDataDir='./', userFCs=NULL, geneLevel=TRUE){
  
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
  
  XLIM<-c(0.3,0.95)
  
  ccr.multDensPlot(toPlot,
                   XLIMS = XLIM,
                   TITLE='',
                   COLS=c('gray','darkgreen','black'),
                   LEGentries = c('Project Score background',
                                  'Project Score replicates',
                                  'Reproducibility Threshold'),
                   XLAB = 'R')
  
  abline(v=sigTH)
  
  # plot pair-wise user replicate R scores onto multDensPlot
  
  if (length(userFCs)>0){
    
    fc<-userFCs
    nr<-ncol(fc)-2
    fc<-fc[,3:ncol(fc)]
    colnames(fc) <- paste('UserData_R', 1:ncol(fc), sep='')
    rownames(fc) <- userFCs$sgRNA
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
    
    legend('topleft',legend = 'User data',inset = c(0,0.25), pt.cex = 3,
           pch=21, pt.bg=rgb(200, 0, 255, maxColorValue = 255, alpha = 120),bty = 'n')
    
  }
}

plot.HT29RScores <- function(refDataDir='./', geneLevel=TRUE){
  
  layout(matrix(1:6, nrow=2, ncol=3, byrow=TRUE))
  
  for (f in ref_fn) {
    
    load(paste(refDataDir,f,sep=""))
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
    nn<-str_split(f,'_foldChanges.Rdata')[[1]][1]
    
    plot(cc,
         rep(0,length(cc)),
         xlim=c(0.3,0.95),
         cex=3,
         pch=21,
         bg=rgb(0,0,255,maxColorValue = 255,alpha = 120),
         main=nn, 
         frame.plot = FALSE,
         yaxt='n', 
         xlab='R',
         ylab='',
         ylim=c(-1,1),
         col='white')
    
    abline(h=0,col='gray')
    
  }
  
}

