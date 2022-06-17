
HT29R.downloadRefData <- function(whatToDownload='FCs',
                                  destFolder='./',
                                  dataRepoURL = "https://cog.sanger.ac.uk/cmp/downloads/crispr_cas9_benchmark/",
                                  expNames=c("HT29_c903","HT29_c904","HT29_c905","HT29_c906","HT29_c907","HT29_c908")){
  nexp <- length(expNames)

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

HT29R.FCdistributions <- function(refDataDir='./',
                                  resDir='./',
                                  userFCs=NULL, 
                                  saveToFig=TRUE, 
                                  display=FALSE){

    fn <- dir(refDataDir)
    fn <- grep('_foldChanges.Rdata', fn, value=TRUE)
    
    if (length(fn)==0) {
        stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')}
    
    nf <- length(fn)
    n <- nf
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    COL <- col_vector[1:n]
    names(COL) <- fn

    data(KY_Library_v1.0)
    cguides <- rownames(KY_Library_v1.0)

    all_sgRNAs <- lapply(fn,function(x){
            load(paste(refDataDir,'/',x,sep=''))
            return(foldchanges$sgRNA)
    })

    cguides <- Reduce(intersect, all_sgRNAs)

    if(!is.null(userFCs)){
        cguides <- intersect(cguides,userFCs$sgRNA)
    }

    firstEl <- fn[1]

    if(saveToFig){
        display <- TRUE
        pdf(paste(resDir,'/QC_FCdistproperties.pdf',sep=''),10,6)
        }
    
    if(display) {
        
        par(mfrow=c(1,2))

        allAvgFc <- lapply(fn,function(x){
            load(paste(refDataDir,'/',x,sep=''))
            nr <- ncol(foldchanges)-2
            fc <- foldchanges[match(cguides,foldchanges$sgRNA),3:ncol(foldchanges)]
            rownames(fc)<-cguides
            nreps<-ncol(fc)-2
            currentCol<-makeTransparent(COL[x])
            for (j in 1:nreps){
                if (x==firstEl & j==1) {
                    plot(density(fc[,j]),
                        frame.plot=FALSE,
                        xlim=c(-9,3),
                        lwd=7,
                        col=currentCol,
                        main='',
                        xlab='sgRNA log fold-change')
                }else{
                    par(new=TRUE)
                    plot(density(fc[,j]),
                        frame.plot=FALSE,
                        xlim=c(-9,3),
                        lwd=7,
                        xaxt='n',
                        yaxt='n',
                        xlab='',
                        ylab='',
                        col=currentCol,
                        main='')
                    }
                } 
    
            fc<-rowMeans(fc)
            return(fc)
        }
        )
    
        if(!is.null(userFCs)){
            nr <- ncol(userFCs)-2
            for (i in 1:nr){
                par(new=TRUE)
                plot(density(userFCs[,2+i]),
                    frame.plot=FALSE,
                    xlim=c(-9,3),
                    lwd=1,
                    xaxt='n',
                    yaxt='n',
                    xlab='',
                    ylab='',
                    main='',
                    col='black')
                }
            legend('topleft','User data',lwd = 3, bty="n")
        }


        GlobalFC <- do.call(cbind,allAvgFc)
        colnames(GlobalFC) <- unlist(lapply(str_split(fn,'_foldChanges.Rdata'),
                                    function(x){x[[1]][1]}))

        if(!is.null(userFCs)){
            nr <- ncol(userFCs)-2
            GlobalFC <- cbind(GlobalFC,rowMeans(userFCs[match(cguides,userFCs$sgRNA),3:(2+nr)]))
            colnames(GlobalFC)[ncol(GlobalFC)] <- 'User data'
            COL <- c(COL,'black')
            nf <- nf+1
            NC <- ncol(GlobalFC)-1
        }else{
            NC <- ncol(GlobalFC)
        }

        tmp <- vioplot(GlobalFC,
                    col=makeTransparent(COL),
                    frame.plot=FALSE,
                    ylab='Avg sgRNA log fold-change',
                    pch=16,
                    pars=list(outcol=makeTransparent(COL)),
                    las=2,
                    cex.axis=0.9,
                    xaxt='n')

        axis(1, at=1:nf, labels=FALSE)
        text(x=1:nf,y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
        labels=colnames(GlobalFC),srt=45, adj=1, xpd=TRUE)

        abline(h=median(tmp$median), col='darkgray')
        abline(h=median(tmp$q1),col='darkgray',lty=2)
        abline(h=median(tmp$q3),col='darkgray', lty=2)
        abline(h=median(tmp$lower),col='darkgray', lty=2)
        abline(h=median(tmp$upper),col='darkgray', lty=2)
        abline(h=median(apply(GlobalFC,2,'min')),col='darkgray')
        abline(h=median(apply(GlobalFC,2,'max')),col='darkgray')

        
    }
    
    if(saveToFig) {
        dev.off()
    }
    
    return(GlobalFC)
}

HT29R.evaluateReps <- function(refDataDir='./',
                               resDir='./',
                               userFCs=NULL, 
                               geneLevel=TRUE,
                               saveToFig=TRUE, 
                               display=FALSE) {

    data(HT29R.reproducible_GeneGuides)
    
    if(!geneLevel){
        data(HT29R.prSCORE_rCorr)
        ScoreCorrs <- HT29R.prSCORE_rCorr
        sigTH <- 0.55
    }else{
        data(HT29R.GL_prSCORE_rCorr)
        ScoreCorrs <- HT29R.GL_prSCORE_rCorr
        sigTH <- 0.68
        data(KY_Library_v1.0)
    }
    
    fn <- dir(refDataDir)
    fn <- grep('_foldChanges.Rdata',fn,value=TRUE)
    
    if (length(fn)==0){
        stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')
    }
    
    toPlot<-list(Pr.Score_bg=density(ScoreCorrs$BGscores), # all the possible pairwise R corr between replicates across different experiment (i.e. PJSC BACKGROUND)
           Pr.Score_repCor=density(ScoreCorrs$REPscores)) # pairwise R corr between replicates of the same experiment

    if(saveToFig){     
        display <- TRUE     
        if(!geneLevel){
            pdf(paste(resDir,'/QC_REPRODUCIBILITY_sgRNAlevel.pdf',sep=''),6,8)
        }else{
            pdf(paste(resDir,'/QC_REPRODUCIBILITY_GENElevel.pdf',sep=''),6,8)
        }
    }

    if(display){

        layout(matrix(c(1,1,1,1,
                        2,3,4,5,
                        6,7)), widths = lcm(12), heights = lcm(6))
                

        XLIM <- c(0,1.0)

        ccr.multDensPlot(toPlot,
                        XLIMS = XLIM,
                        TITLE='',
                        COLS=c('gray','darkgreen','darkblue'),
                        LEGentries = c('Project Score background',
                              'Project Score replicates',
                              'Reproducibility Threshold'),
                        XLAB = 'R')
                                  
        abline(v=sigTH, col = "darkblue", lwd = 3, lty = 2)

        if (!is.null(userFCs)){
            fc <- userFCs
            nr <- ncol(fc)-2
            fc <- fc[,3:ncol(fc)]
            colnames(fc) <- paste('UserData_R', 1:ncol(fc), sep='')
            rownames(fc) <- userFCs$sgRNA
            dataSET <- fc[HT29R.reproducible_GeneGuides,]
    
            if(geneLevel){
                dataSET <- apply(dataSET,MARGIN = 2,function(x){
                    cdata <- x
                    names(cdata)<-rownames(dataSET)
                    ccr.geneMeanFCs(cdata,KY_Library_v1.0)
                })}

            cc<-c(as.dist(cor(dataSET)))

            points(cc,
                   rep(0,length(cc)),
                   cex=1.5,
                   pch=21,
                   bg=rgb(200,0,255,maxColorValue = 255)) 
        
            legend(x=0,
                   y=max(toPlot$Pr.Score_repCor$y) - max(toPlot$Pr.Score_repCor$y)/8,
                   legend = 'User data',
                   inset = c(0,0.25),
                   pt.cex = 2,
                   pt.bg=rgb(200,0,255,maxColorValue = 255,alpha = 255),
                   pch=21,            
                   bty = 'n')

            cat(blue(paste(length(which(cc>=sigTH)),' pair-wise replicate comparisons (out of ',length(cc),
            ') yield correlation scores greater or equal than Project Score QC threshold\n',sep='')))
    
        }


        lapply(fn,function(x){ 
            load(paste(refDataDir,'/',x,sep=''))
            nr <- ncol(foldchanges)-2     
            fc <- foldchanges[,3:ncol(foldchanges)]
            rownames(fc) <- foldchanges$sgRNA
            dataSET <- fc[HT29R.reproducible_GeneGuides,]

            if(geneLevel){
                dataSET<-apply(dataSET,MARGIN = 2,function(x){
                        cdata<-x
                        names(cdata)<-rownames(dataSET)
                        ccr.geneMeanFCs(cdata,KY_Library_v1.0)})
            }

            cc <- c(as.dist(cor(dataSET)))

            nn <- str_split(x,'_foldChanges.Rdata')[[1]][1]
    
            par(mar = c(1, 4, 1, 2))  
            plot(cc,
                rep(0,length(cc)),
                xlim=c(0,1),
                cex=3,
                pch=21,
                bg=rgb(0,0,255,maxColorValue = 255,alpha = 120),
                frame.plot = FALSE,
                xaxt='n',
                yaxt='n',
                xlab='',
                ylab='',
                ylim=c(-1,1),
                col='white')

            abline(h=0,col='gray')
    
            legend("topleft",legend=nn,bty = 'n' )
        })
    }
    
    if(saveToFig){
        dev.off()
        }
    
}

HT29R.expSimilarity <- function(refDataDir='./',
                                resDir='./',
                                userFCs=NULL, 
                                geneGuides=c("All","HI"), 
                                geneLevel=TRUE,
                                Rscores=TRUE,
                                saveToFig=TRUE,
                                display=FALSE) {

    data(KY_Library_v1.0)
    data(HT29R.prSCORE_bkgr_screen_similarity)
    data(HT29R.prSCORE_bkgr_screen_similarity_HI)
    data(HT29R.prSCORE_bkgr_screen_similarity_sgRNA)
    data(HT29R.prSCORE_bkgr_screen_similarity_sgRNA_HI)
    data(HT29R.reproducible_GeneGuides)
    data(HT29R.consensus_GeneGuides)

    geneGuides <- match.arg(geneGuides)

    fn <- dir(refDataDir)
    fn <- grep('_foldChanges.Rdata', fn, value=TRUE)

    if (length(fn)==0) {
        stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')
    }

    if (!is.null(userFCs)) {
    
        cguides <- userFCs$sgRNA
    
        if (geneGuides == "HI") {
            cguides <- intersect(cguides, HT29R.reproducible_GeneGuides)
        } else {
            cguides <- intersect(cguides, HT29R.consensus_GeneGuides) 
        }
    
        userFCs <- userFCs[match(cguides, userFCs$sgRNA),3:ncol(userFCs)]
        rownames(userFCs) <- cguides
        userFCs <- apply(userFCs, MARGIN = 1, 'mean') 
    
        if (geneLevel) {
            userFCs <- ccr.geneMeanFCs(userFCs,KY_Library_v1.0)
            cgenes <- names(userFCs)
        }

    } else {

        if (geneGuides == "HI") {
            cguides <- HT29R.reproducible_GeneGuides
            if (geneLevel) {
                cgenes <- unique(KY_Library_v1.0[which(rownames(KY_Library_v1.0) %in% HT29R.reproducible_GeneGuides), "GENES"])
            }
        } else {
            cgenes <- unique(KY_Library_v1.0$GENES) 
        }
    }

    ref_fcs <- lapply(fn, function(x) {
        load(paste(refDataDir,'/',x,sep=''))
        nr <- ncol(foldchanges)-2
        fc <- foldchanges[, 3:ncol(foldchanges)]
        rownames(fc) <- foldchanges$sgRNA
        fc <- apply(fc, MARGIN = 1, FUN = 'mean')
        if (geneLevel) {
            fc <- ccr.geneMeanFCs(fc, KY_Library_v1.0)
        }
        return(fc)
    })

    if (geneLevel) {
        cgenes <- intersect(cgenes,Reduce('intersect',lapply(ref_fcs,'names')))
        ref_fcs <- do.call(cbind,lapply(ref_fcs,function(x){x[cgenes]}))
    } else {
        cguides <- intersect(cguides,Reduce('intersect',lapply(ref_fcs,'names')))
        ref_fcs <- do.call(cbind,lapply(ref_fcs,function(x){x[cguides]}))
    }

    fn <- unlist(lapply(str_split(fn,'_foldChanges.Rdata'),function(x){x[1]}))

    colnames(ref_fcs) <- fn
    obsCorr <- c(as.dist(cor(ref_fcs)))

    if (geneGuides == "All") {
        if (geneLevel) {
            bgCorr <- HT29R.prSCORE_bkgr_screen_similarity 
        } else {
            bgCorr <- HT29R.prSCORE_bkgr_screen_similarity_sgRNA
        }
    } else {
        if (geneLevel) {
            bgCorr <- HT29R.prSCORE_bkgr_screen_similarity_HI 
        } else {
            bgCorr <- HT29R.prSCORE_bkgr_screen_similarity_sgRNA_HI
        }
    }


    toPlot <- list(Score_bg=density(bgCorr),
                RefScreen_sim=density(obsCorr))

    userCorrs <- NULL

    if(!is.null(userFCs)) {
        if (geneLevel) {
            ref_fcs <- cbind(userFCs[cgenes], ref_fcs[cgenes,])
        } else { 
            ref_fcs <- cbind(userFCs[cguides],ref_fcs[cguides,])
        }
        colnames(ref_fcs)[1] <- 'User data'
        userCorrs <- cor(ref_fcs)
        userCorrs <- userCorrs[1, 2:ncol(userCorrs)]  
    }

    if(saveToFig){
        display <- TRUE
        if(geneLevel){
            pdf(paste(resDir,geneGuides,'_QC_SIMILARITY_GENElevel.pdf',sep=''),5.23,7.03)
        }else{
            pdf(paste(resDir,geneGuides,'_QC_SIMILARITY_sgRNAlevel.pdf',sep=''),5.23,7.03)
        }
    }

    if(display){

        layout(mat = c(1,1,2))
        ccr.multDensPlot(toPlot,
                        XLIMS = c(0.3,1), 
                        TITLE='Screen similarity',
                        COLS=c('gray','darkgreen'),
                        LEGentries = c('Project Score background','Ref.pair-wise HT29 screens'),
                        XLAB = 'R')

        if(!is.null(userFCs)) {
            
            points(userCorrs,
                    rep(0,length(userCorrs)),
                    cex=1.5,
                    pch=21,
                    bg=rgb(200,0,255,maxColorValue = 255,alpha = 120))

            legend('topleft',
                    legend = 'User data',
                    inset = c(0,0.10),
                    pt.cex = 1.5,
                    pch=21,
                    pt.bg=rgb(200,0,255,maxColorValue = 255,alpha = 120),
                    bty = 'n')

        }

        tres <- t.test(bgCorr, obsCorr)

        mtitle <- paste('PRJ SCORE BACKGROUND vs REFERENCE = ', format(tres$p.value, scientific = TRUE, digits = 2))

        if (!is.null(userFCs)) {
            tres <- t.test(obsCorr,userCorrs)
            mtitle <- paste(mtitle,'\n USER-DATA vs REFERENCE = ', format(tres$p.value, scientific = TRUE, digits= 2, sep=''))
        }

        boxplot(bgCorr, 
                obsCorr, 
                userCorrs, 
                horizontal = TRUE,
                ylim=c(0.3,1),
                col=c('gray','darkgreen',rgb(200,0,255,maxColorValue = 255,alpha = 120)),
                main=mtitle)

        
        if(Rscores) {
            
            if (saveToFig) {
                dev.off()
                if(geneLevel){
                    pdf(paste(resDir, geneGuides,'_SCATTERPLOT_R_MATRIX_GENElevel.pdf',sep=''),15,15)
                } else {
                    pdf(paste(resDir, geneGuides,'_SCATTERPLOT_R_MATRIX_sgRNAlevel.pdf',sep=''),15,15)
                }
            }
            
            pairs(ref_fcs,
                    lower.panel = panel.cor, 
                    upper.panel = my.panelSmooth)
            }
    
    }

    if(saveToFig){
        dev.off()
    }

    return(ref_fcs)
}


HT29R.PhenoIntensity <- function(refDataDir='./',
                                 resDir='./',
                                 userFCs=NULL, 
                                 geneLevel=TRUE,
                                 saveToFig=TRUE,
                                 display=FALSE){

    fn<-dir(refDataDir)
    fn<-grep('_foldChanges.Rdata',fn,value=TRUE)

    if (length(fn)==0){
        stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')    
    }

    if(saveToFig){
        display <- TRUE
        pdf(paste(resDir,'/ALLSCREENS_PHENOINTENSITY.pdf',sep=''),11.88,9.71)
    }

    if(display) {

        par(mfrow=c(3,3))

        RES <- do.call(rbind,lapply(fn, function(x){
            load(paste(refDataDir,x,sep=''))
            currentFc<-rowMeans(foldchanges[,3:ncol(foldchanges)])
            names(currentFc)<-foldchanges$sgRNA
            HT29R.singleScreen_PhenoIntensity(currentFc,
                                          geneLevel = geneLevel,
                                          expName = unlist(lapply(str_split(x,'_foldChanges.Rdata'),
                                                              function(y){y[1]})))
            }))


        if(!is.null(userFCs)){
            currentFc <- rowMeans(userFCs[,3:ncol(userFCs)])
            names(currentFc) <- userFCs$sgRNA
            usrRes <- HT29R.singleScreen_PhenoIntensity(currentFc,
                                              geneLevel = geneLevel,
                                              expName = 'User Data')
        }

        if(saveToFig){
            dev.off()
            pdf(paste(resDir,'/REFERENCE_PHENOINTENSITY.pdf',sep=''),3,7)
        }

        par(mar=c(12,4,1,1))
        boxplot(list(RES$GD_ribProt,RES$GD_essential),
                names = c('Ribosomal protein genes','Other essential genes'),
                ylab='Glass Delta',
                col=c('darkblue','blue'),
                border = 'darkgrey',
                lwd=2,
                las=2)

        if (!is.null(userFCs)){
            points(c(1,2),
            c(usrRes$GD_ribProt,usrRes$GD_essential),
            cex=2,
            col=rgb(200,0,255,maxColorValue = 255,alpha = 255),
            pch=16)
        }

        if(saveToFig){
            dev.off()
            }
    }

}

HT29R.singleScreen_PhenoIntensity <- function(FCprofile, geneLevel=TRUE, expName=NULL){
  
  data(EssGenes.ribosomalProteins)
  data(BAGEL_essential)
  data(BAGEL_nonEssential)

  data(KY_Library_v1.0)

  if(geneLevel){
    FCprofile <- ccr.geneMeanFCs(sgRNA_FCprofile = FCprofile, libraryAnnotation = KY_Library_v1.0)

  }else{
    EssGenes.ribosomalProteins <- rownames(KY_Library_v1.0)[match(EssGenes.ribosomalProteins,KY_Library_v1.0$GENES)]
    EssGenes.ribosomalProteins <- EssGenes.ribosomalProteins[!is.na(EssGenes.ribosomalProteins)]

    BAGEL_essential <- rownames(KY_Library_v1.0)[match(BAGEL_essential,KY_Library_v1.0$GENES)]
    BAGEL_essential <- BAGEL_essential[!is.na(BAGEL_essential)]

    BAGEL_nonEssential <- rownames(KY_Library_v1.0)[match(BAGEL_nonEssential,KY_Library_v1.0$GENES)]
    BAGEL_nonEssential <- BAGEL_nonEssential[!is.na(BAGEL_nonEssential)]

    }

  ribProtFC <- FCprofile[intersect(names(FCprofile),EssGenes.ribosomalProteins)]
  BagEssFC <- FCprofile[intersect(names(FCprofile),BAGEL_essential)]
  BagNonEssFC <- FCprofile[intersect(names(FCprofile),BAGEL_nonEssential)]

  toPlot<-list(NonEss=density(BagNonEssFC),
               Ess=density(BagEssFC),
               ribProt=density(ribProtFC))

  GD_essential <- abs(mean(BagEssFC)-mean(BagNonEssFC))/sd(BagEssFC)
  GD_ribProt <- abs(mean(ribProtFC)-mean(BagNonEssFC))/sd(ribProtFC)


  ccr.multDensPlot(toPlot,COLS=c('gray','blue','darkblue'),XLIMS = range(FCprofile),XLAB = 'depletion fold-change',
                   TITLE=paste(expName,'\nGD ess. = ',format(GD_essential,digits = 3),', GD rib.prot. = ',format(GD_ribProt,digits=3)),
                   LEGentries = c('Non Essential','Essential','Ribosomal Proteins'))
  
  abline(v=mean(BagNonEssFC),col='gray')
  abline(v=mean(BagEssFC),col='blue')
  abline(v=mean(ribProtFC),col='darkblue')
  

  res <- data.frame(GD_essential=GD_essential,GD_ribProt=GD_ribProt,stringsAsFactors = FALSE)
  rownames(res) <- expName
  
  return(res)
  
}

HT29R.ROCanalysis <- function(refDataDir='./',
                              resDir='./',
                              positives,
                              negatives,
                              userFCs=NULL, 
                              geneLevel=TRUE,
                              saveToFig=TRUE,
                              display=FALSE){
    
    fn <- dir(refDataDir)
    fn <- grep('_foldChanges.Rdata',fn,value=TRUE)

    if (length(fn)==0){
        stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')
    }

    if (geneLevel) {
        warning(strwrap("Please make sure that your Positive and Negative controls 
                        are character vector of gene names when geneLevel parameter is set
                        to TRUE",prefix = "\n", initial = "\n"))
    } else {
        warning(strwrap("Please make sure that your Positive and Negative controls
                        are character vector of sgRNAs when geneLevel parameter is set
                        to FALSE",prefix = "\n", initial = "\n"))
    }   

    if(saveToFig){
        display <- TRUE
        pdf(paste(resDir,'/ALLSCREENS_ROCs.pdf',sep=''),18,6)
    }

    if(display) {
        
        par(mfrow=c(length(fn)/3,6))
        
        RES <- lapply(fn,function(x){
            load(paste(refDataDir,x,sep=''))
            currentFc <- rowMeans(foldchanges[,3:ncol(foldchanges)])
            names(currentFc) <- foldchanges$sgRNA
            ename <- str_split(x,'_foldChanges.Rdata')[[1]][1]
            HT29R.individualROC(currentFc,
                                positives,
                                negatives,
                                geneLevel = geneLevel,
                                expName=ename)
        })

        if(saveToFig){
            dev.off()
        }

        if(!is.null(userFCs)){
            currentFc<-rowMeans(userFCs[,3:ncol(userFCs)])
            names(currentFc)<-userFCs$sgRNA

            if(saveToFig) {
                display <- TRUE
                pdf(paste(resDir,'/USER_ROCs.pdf',sep=''),10,5)
            }
        
            par(mfrow=c(1,2))
            usrRes <- HT29R.individualROC(currentFc,
                                      positives,
                                      negatives,
                                      geneLevel = geneLevel,
                                      expName='User Data')
            
            if(saveToFig) {
                dev.off()
            }
        }

        plot(0,0,xlim=c(1,0),ylim=c(0,1),col=NA,xlab='TNR',ylab='Recall')
        abline(1,-1,col='darkgray')
        lapply(1:length(RES),function(x){
            lines(RES[[x]]$ROCres$curve[,1],RES[[x]]$ROCres$curve[,2],lwd=5,col=rgb(0,0,255,alpha = 100,maxColorValue = 255))
            abline(h=RES[[x]]$ROCres$Recall,lty=2,col=rgb(0,0,255,alpha = 100,maxColorValue = 255),lwd=2)
        })

        if(!is.null(userFCs)){
            lines(usrRes$ROCres$curve[,1],usrRes$ROCres$curve[,2],lwd=2,col=rgb(255,0,255,alpha = 255,maxColorValue = 255))
            abline(h=usrRes$ROCres$Recall,lty=2,col=rgb(255,0,255,alpha = 255,maxColorValue = 255),lwd=2)
        }
    }
}

HT29R.individualROC <- function(FCprofile, 
                                positives, 
                                negatives,
                                geneLevel=TRUE,
                                expName=NULL){

  data(KY_Library_v1.0)

  if(geneLevel){
    FCprofile <- ccr.geneMeanFCs(sgRNA_FCprofile = FCprofile,libraryAnnotation = KY_Library_v1.0)
    
  }else{

    positives <- rownames(KY_Library_v1.0)[match(positives,rownames(KY_Library_v1.0))]
    positives <- positives[!is.na(positives)]

    negatives <- rownames(KY_Library_v1.0)[match(negatives,rownames(KY_Library_v1.0))]
    negatives <- negatives[!is.na(negatives)]

  }

  ROCres <- ccr.ROC_Curve(FCprofile,positives,negatives, FDRth = 0.05,expName = paste(expName,'ROC'))
  PRRCres <- ccr.PrRc_Curve(FCprofile,positives,negatives, FDRth = 0.05,expName = paste(expName,'PrRc'))

  return(list(ROCres=ROCres,PRRCres=PRRCres))
}

HT29R.FDRconsensus <- function(refDataDir="./", 
                                resDir="./", 
                                userFCs=NULL, 
                                distance=c("GlDelta","Cohen's"), 
                                FDRth=0.05,
                                saveToFig=TRUE,
                                display=FALSE) {

    distance <- match.arg(distance)

    data(KY_Library_v1.0)
    data(BAGEL_essential)
    data(BAGEL_nonEssential)

    TruePositives <- BAGEL_essential
    TrueNegatives <- BAGEL_nonEssential

    fn <- dir(refDataDir)
    fn <- grep("_foldChanges.Rdata", fn, value=TRUE)

    if (!is.null(userFCs)) {
        gg <- userFCs$sgRNA
        userFCs <- userFCs[, 3:ncol(userFCs)]
        rownames(userFCs) <- gg
        userFCs <- apply(userFCs, MARGIN = 1, 'mean')
        userFCs <- ccr.geneMeanFCs(userFCs, KY_Library_v1.0)
        } 

    cgenes <- unique(KY_Library_v1.0$GENES)

    ref_fcs <- lapply(fn, function(x) {
        load(paste(refDataDir, x, sep=""))
        nr <- ncol(foldchanges)-2
        fc <- foldchanges[, 3:ncol(foldchanges)]
        rownames(fc) <- foldchanges$sgRNA
        fc <- apply(fc, MARGIN = 1, FUN = "mean")
        fc <- ccr.geneMeanFCs(fc,libraryAnnotation = KY_Library_v1.0)})

    cgenes <- intersect(cgenes, Reduce("intersect", lapply(ref_fcs, "names")))

    ref_fcs <- do.call(cbind, lapply(ref_fcs, function(x){x[cgenes]}))
    fn <- unlist(lapply(str_split(fn, "_foldChanges.Rdata"), function(x){x[1]}))
    colnames(ref_fcs) <- fn

    predictions <- ref_fcs[intersect(c(TruePositives, TrueNegatives),row.names(ref_fcs)),]
    observations <- is.element(row.names(predictions), TruePositives)+0

    RES <- lapply(1:ncol(ref_fcs), function(x) { 
        currentFc <- predictions[, x]
        names(observations) <- names(currentFc)
        res <- roc(observations, currentFc, direction = ">", quiet = TRUE) 
        COORS <- coords(res, "all", ret = c("threshold", "ppv", "sensitivity"), transpose=TRUE) 
        FDRpercTh <- max(COORS['threshold', which(COORS['ppv',] >= (1 - FDRth))])
        threshold <- COORS["threshold", min(which(COORS["threshold",] <= FDRpercTh))]
        list(TH=threshold)  
        })
    
    FDR5_Positives <- Reduce(intersect, lapply(1:length(RES),function(x){names(ref_fcs[ref_fcs[,x] <= RES[[x]]$TH, x])}))
    FDR5_Negatives <- Reduce(intersect, lapply(1:length(RES),function(x){names(ref_fcs[ref_fcs[,x] > RES[[x]]$TH, x])}))

    consensus <- ref_fcs[intersect(c(FDR5_Positives, FDR5_Negatives), row.names(ref_fcs)),]
    colnames(consensus) <- colnames(ref_fcs)

    if(!is.null(userFCs)) {
        consensus <- cbind(as.numeric(userFCs[rownames(consensus)]), consensus)
        colnames(consensus)[1] <- "User data" 
        }

    C_dist <- function(x,y) {
        mu_x <- mean(x)
        mu_y <- mean(y)
        n1 <- length(x)
        n2 <- length(y)
        sigma_pooled <- sqrt(((n1-1)*sd(x)^2 + (n2-1)*sd(y)^2) / (n1 + n2 - 2))
        Size <- abs(mu_x - mu_y) / sigma_pooled
        return(Size)
        }

    dist <- list()

    for (i in 1:ncol(consensus)) {
        x1 <- as.numeric(consensus[FDR5_Positives,i])
        x2 <- as.numeric(consensus[FDR5_Negatives,i])
        if (distance == "Cohen's") {
            dist[i] <- C_dist(x1,x2)
        } else if (distance == "GlDelta") {
            dist[i] <- abs(mean(x1)-mean(x2))/sd(x1)
        } 
    }

    dist <- unlist(dist)
    names(dist) <- colnames(consensus)

    # graphics 
    group <- rep(c("Positives","Negatives"), c(length(FDR5_Positives), length(FDR5_Negatives)))
    consensus <- cbind.data.frame(consensus, group)
    consensus$group <- factor(group)
    resh_consensus <- melt(consensus)

    if(saveToFig){
        display <- TRUE
        pdf(paste(resDir,'/FDR_CONSENSUS_DIST.pdf',sep=''),6,10)
    }

    if(display) {

        layout(mat = matrix(c(1, 2, 0, 0), nrow = 2, ncol = 1),
            heights = c(2, 1),    
            widths = c(10, 5))     

        par(mar=c(5,7,2,2))

        tmp <- boxplot(value~group*variable,
                data = resh_consensus,
                frame.plot = FALSE,
                yaxt = 'n',
                ylab = '',
                xlab = "log fold-change",
                ylim = c(-7,3),
                las=1,
                col=c("#117733","#CC6677"),
                pch=14,
                pars=list(outcol=c("#117733","#CC6677")),
                horizontal = TRUE,
                cex.axis=0.9)

        if(!is.null(userFCs)){
            NL <- 7
        } else {
            NL <- 6
        }

        axis(2, at=(1:(NL*2))[-2*(1:NL)], labels=FALSE)
        text(y=(1:(NL*2))[-2*(1:NL)], 
             x=-7.55,
             labels=colnames(consensus)[1:(ncol(consensus)-1)], 
             srt=45, 
             adj=1, 
             xpd=TRUE)
        
        #legend(1 ,15, legend=c("FDR5% pos", "FDR5% neg"), col=c("#CC6677","#117733"), box.lty = 0, pch=22, cex=1.0)

        par(mar=c(4,2,2,2))

        ylim <- c(min(dist), max(dist))
        
        if(!is.null(userFCs)){
            userDist <- dist[["User data"]]
            dist <- dist[2:length(dist)]
            if (as.numeric(userDist) < min(as.numeric(dist))) {
                ylim <- c(as.numeric(userDist)-0.2, max(as.numeric(dist)))
            } else {
                ylim <- c(min(as.numeric(dist)),max(as.numeric(dist)))
            }
        }

        boxplot(dist, 
                type="o",
                lwd=1,
                xlab="Distance",
                ylim=ylim,
                main="",
                outline=FALSE,
                range=0,
                border="black",
                boxlwd=4,
                whiskcol = "black",
                whisklty = 2,
                horizontal=TRUE,
                col="white")

        if(!is.null(userFCs)) {
            points(userDist,1,cex=1.5,pch=21,bg=rgb(200,0,255,maxColorValue = 255))  
            legend("topleft", legend="User data", pt.cex = 1.5,pch=21,pt.bg=rgb(200,0,255,maxColorValue = 255),bty = 'n')
        }
        
        if(saveToFig){
            dev.off()
        }
    }
    return(list(POS=FDR5_Positives, NEG=FDR5_Negatives, Universe=cgenes, DF=ref_fcs))
}

HT29R.sgRNAFCStats <- function(x, userFCs=NULL) {
  
  stats <- describe(x, quant=c(.1,.25,.75,.90))
  
  if(!is.null(userFCs)) {
    userStats <- stats['User data', 5:ncol(stats)]
    NC <- ncol(x)-1
    } else {
    NC <- ncol(x)
    }
  
  AvgStats <- apply(stats[1:NC,5:ncol(stats)],2,'mean')
  
  SE <- lapply(stats[1:NC,5:ncol(stats)],function(x){sd(x)/sqrt(NC)})
  
  cat("HT29 sgRNAs logFCs statistics:\n\n")
  cat(paste(c('Avg. Range: ','; '), round(c(AvgStats[4],AvgStats[5]),digits=2),
            '±',round(c(SE$min,SE$max),digits=3),sep=''),'\n')
  cat(paste('Avg. Median: ', 
            round(AvgStats[1],digits=3),'±',
            round(SE$median,digits=3),sep=''),'\n')
  cat(paste(c('Avg. IQR range: ','; '), 
            round(c(AvgStats[11],AvgStats[12]),digits=2),
            '±',round(c(SE$Q0.25,SE$Q.075),digits=2),sep=''),'\n')
  cat(paste(c('Avg. 10-90th perc range: ','; '), 
            round(c(AvgStats[10],AvgStats[13]),digits=2),
            '±',round(c(SE$Q0.1,SE$Q0.9),digits=2),sep=''),'\n')
  cat(paste('Avg. Skewness: ', 
            round(AvgStats[7],digits=2),
            '±',round(SE$skew,digits=2),sep=''),'\n')
  cat(paste('Avg. Kurtosis: ', 
            round(AvgStats[8],digits=2),
            '±',round(SE$kurtosis,digits=2),sep=''),'\n')
  
  if(!is.null(userFCs)) {
    cat('\nUser screen sgRNA logFCs statistics:\n\n')
    cat(paste(c('Range min: ','; Range max: '),
              round(c(userStats[,4], userStats[,5]), digits=3),sep=""),'\n')
    cat(paste('Median: ', round(userStats[,1], digits=3), sep=""),'\n')
    cat(paste(c('IQR min: ','; IQR max: '),
              round(c(userStats[,11], userStats[,12]), digits=3), sep=""),'\n')
    cat(paste(c('10th perc: ','; 90th perc: '),round(c(userStats[,10],userStats[,13]),digits=3),sep=""),'\n')
    cat(paste('Skewness: ', round(userStats[,7], digits=3), sep=""),'\n')
    cat(paste('Kurtosis: ', round(userStats[,8], digits=3), sep=""),'\n')
    }

}

HT29R.runCrisprQC_Analysis <- function(refData = "FCs", userFCs=NULL, outdir="./", positives, negatives, FDRth=0.05) {

    cat(paste("Downloading HT-29", refData, "data in", outdir, "...\n"), sep=" ")
    HT29R.downloadRefData(whatToDownload = refData, destFolder = outdir)
    Sys.sleep(1)
    cat("...Done!\n\n")

    dir.create(paste(outdir, "PLOTS/",sep=""))
    pathToDir <- paste(outdir,"PLOTS/", sep="")

    cat(paste(" === QC PLOTS WILL BE STORED IN:", pathToDir, " ===\n"), sep="HI")
    cat("Saving sgRNA average fold-changes distribution across replicates for each HT-29 experiment...\n")
    
    if(!is.null(userFCs)){
        cat("Saving sgRNA average fold-changes distribution across replicates for User-provided data...\n")
        }

    STATS <- HT29R.FCdistributions(refDataDir = outdir, resDir = pathToDir, userFCs = NULL)
    
    if(file.exists(paste(pathToDir, "/QC_FCdistproperties.pdf", sep=""))) {
        Sys.sleep(1)
        cat("...Done!\n\n")
        }
    
    cat("Average parameters and confidence intervals:\n\n")
    HT29R.sgRNAFCStats(STATS, userFCs = NULL)

    cat("\n1) LOW-LEVEL QC USING THE DOWNLOADED REFERENCE CELL LINES:\n")
    cat("Select the sgRNA set to use - i.e, type \"HI\" for the highly-informatives ones and \"All\" for the entire library:\n")
    var1 <- readline()
    var1 <- as.character(var1)

    cat("Select at what level you wish to assess reproducibility - i.e, type TRUE if at gene-level (by default) OR FALSE if at sgRNA-level - :\n")
    var2 <- readline()
    var2 <- as.logical(var2)

    cat("Evaluating reproducibility across replicates...\n")
    RES1 <- HT29R.evaluateReps(refDataDir = outdir, resDir = pathToDir, geneLevel=var2, userFCs = NULL)
    Sys.sleep(1)
    cat("...Done!\n\n")

    cat("Evaluating similarity across averaged replicates...\n")
    HT29R.expSimilarity(refDataDir = outdir, resDir = pathToDir, geneGuides = var1, geneLevel = var2, userFCs = NULL)
    Sys.sleep(1)
    cat("...Done!\n\n")
    
    cat("2) HIGH-LEVEL QC USING THE DOWNLOADED REFERENCE CELL LINES:\n")
    cat("Evaluating phenotype intensity...\n")
    HT29R.PhenoIntensity(refDataDir = outdir, resDir = pathToDir, userFCs=NULL, geneLevel = var2)
    Sys.sleep(1)
    cat("...Done!\n\n")
    
    cat("Performing ROC-PrRc analysis...\n")
    HT29R.ROCanalysis(refDataDir = outdir, resDir = pathToDir, positives, negatives, userFCs = NULL, geneLevel = var2)
    Sys.sleep(1)
    cat("...Done!\n\n")

    cat("3) COMPUTING HT-29-SPECIFICS GENES AT 5% FDR....\n")
    cat("Select the distance to be computed - i.e type \"Cohen\'s\" or \"GlDelta\"")
    var3 <- readline()
    var3 <- as.character(var3)

    RES2 <- HT29R.FDRconsensus(refDataDir = outdir, resDir = pathToDir, userFCs = NULL, distance=var3)
    
    write.csv(RES2$FDR5_Positives, paste(pathToDir, "HT-29-specific_genes.csv", sep=""))
    Sys.sleep(2)
    cat("...Done!\n\n")

}


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



