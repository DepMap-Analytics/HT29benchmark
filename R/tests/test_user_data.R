# CREATE USER_DATA.TSV

dir <- "/Users/raffaele.iannuzzi/Documents/Project-HT29/OneDrive_1_5-19-2022"
fn <- dir(dir)

userData <- lapply(fn,function(x){paste(dir,x,sep='/')})

DF <- lapply(userData, function(x) {
    if (x == paste(dir, fn[[1]], sep="/")) {
        read.table(x, header = TRUE, stringsAsFactors = FALSE) }
    else {
        read.table(x, header = TRUE, stringsAsFactors = FALSE, colClasses = c(gene="NULL",ERS717283.plasmid="NULL"))}
})

cguides <- rownames(KY_Library_v1.0)
cguides <- intersect(cguides, Reduce("intersect", lapply(1:length(DF), function(x){DF[[x]]$sgRNA})))

DF <- do.call(cbind, lapply(1:length(DF), function(x){DF[[x]][which(cguides %in% DF[[x]]$sgRNA),]}))
DF[,c(5,7,9,11,13)] <- NULL
COLS <- colnames(DF)[5:ncol(DF)]
DF <- DF[,c("sgRNA", "gene", "ERS717283.plasmid", "X3980STDY6393091.sample", COLS)]

colnames(DF)[4:ncol(DF)] <- c("HT29_c920R1", "HT29_c920R2", "HT29_c920R3",
                                 "HT29_c920R4","HT29_c920R5", "HT29_c920R6")

write.table(DF, file="User_data.tsv", row.names=TRUE, sep="\t")


NormDF <- ccr.NormfoldChanges(Dframe = DF,  
                        min_reads = 30, 
                        EXPname = "", 
                        libraryAnnotation = KY_Library_v1.0, 
                        display = TRUE,
                        outdir = resultsDir)

UserFCs <- NormDF$logFCs

# CREATE BACKGROUND SIMILARITY AND DATA

## project score background similarity 
# sgRNA level (all) 
# sgRNA level (most-informative) 
# gene level (all) 
# gene level (most-informative) 

# load data
fn <- dir(filesDir)
fn <- grep('HT29R.',fn,value=TRUE)
for (f in fn) {
  cat(paste('loading',f,'\n'))
  load(paste(filesDir,f,sep=""))
  }

PathToDownload <- "/Users/raffaele.iannuzzi/Downloads/"

fn <- "EssentialityMatrices/00_logFCs.tsv" # genes x cells
fn <- 'Sanger_corrected_sgRNA_logFCs.tsv' # sgRNAs x cells


PrjSc <- read.table(paste(PathToDownload, fn,sep=""), header=TRUE, row.names = 1)
PrjSc <- na.omit(PrjSc)
dim(PrjSc)

# ALL , GENELEVEL TRUE
df <- sapply(PrjSc, as.numeric)
df <- df[,2:ncol(df)]
HT29R.prSCORE_bkgr_screen_similarity <- c(as.dist(cor(df)))
save(HT29R.prSCORE_bkgr_screen_similarity, file="HT29R.prSCORE_bkgr_screen_similarity.RData")

#HI, GENELEVEL TRUE


# ALL, GENELEVEL FALSE
PrjSc <- PrjSc[HT29R.consensus_GeneGuides, ]
df <- sapply(PrjSc, as.numeric)
df <- df[,2:ncol(df)]
HT29R.prSCORE_bkgr_screen_similarity_sgRNAs <- c(as.dist(cor(df)))
save(HT29R.prSCORE_bkgr_screen_similarity_sgRNAs, file="HT29R.prSCORE_bkgr_screen_similarity_sgRNA.RData")
toPlot<-list(Similarity=density(HT29R.prSCORE_bkgr_screen_similarity_sgRNAs),
                Similarity_previous=density(HT29R.prSCORE_bkgr_screen_similarity))

ccr.multDensPlot(toPlot,
                    XLIMS = c(0.2,1), 
                    TITLE='Screen similarity',
                    COLS=c("grey","black"),
                    LEGentries = c("sgRNA","genes"),
                    XLAB = 'R')


# HI, GENELEVEL FALSE
PrjSc <- PrjSc[HT29R.reproducible_GeneGuides, ]
df <- sapply(PrjSc, as.numeric)
df <- df[,2:ncol(df)]
HT29R.prSCORE_bkgr_screen_similarity_sgRNAs_HI <- c(as.dist(cor(df)))
save(HT29R.prSCORE_bkgr_screen_similarity_sgRNAs_HI, file="HT29R.prSCORE_bkgr_screen_similarity_sgRNA_HI.RData")



HT29R.consensus_similarity <- function(refDataDir="./", resDir="./", userFCs=NULL, FDRth=0.05) {

    data(KY_Library_v1.0)
    data(BAGEL_essential)
    data(BAGEL_nonEssential)

    fn <- dir(refDataDir)
    fn <- grep("_foldChanges.Rdata", fn, value=TRUE)

    if (!is.null(userFCs)) {
        gg <- userFCs$sgRNA
        userFCs <- userFCs[, 3:ncol(userFCs)]
        rownames(userFCs) <- gg
        userFCs <- apply(userFCs, MARGIN = 1, 'mean')
        userFCs <- ccr.geneMeanFCs(userFCs, KY_Library_v1.0)} 

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

    TruePositives <- BAGEL_essential
    TrueNegatives <- BAGEL_nonEssential

    predictions <- ref_fcs[intersect(c(TruePositives, TrueNegatives),row.names(ref_fcs)),]
    observations <- is.element(row.names(predictions), TruePositives)+0

    RES <- lapply(1:ncol(ref_fcs), function(x) { 
        currentFc <- predictions[, x]
        names(observations) <- names(currentFc)
        res <- roc(observations, currentFc, direction = ">", quiet = TRUE) 
        COORS <- coords(res, "all", ret = c("threshold", "ppv", "sensitivity"), transpose=TRUE) 
        FDRpercTh <- max(COORS['threshold', min(which(COORS['ppv',] >= 1 - FDRth ))])
        threshold <- COORS["threshold", min(which(COORS["threshold",] <= FDRpercTh))]
        list(TH=threshold)  
        })

    FDR5_CONSENSUS_P <- Reduce(intersect, lapply(1:length(RES),function(x){
        names(ref_fcs[ref_fcs[,x] <= RES[[x]]$TH, x])}))

    consensus <- ref_fcs[intersect(FDR5_CONSENSUS_P, row.names(ref_fcs)),]
    colnames(consensus) <- colnames(ref_fcs)

    C_dist <- function(x,y) {
        mu_x <- mean(x)
        mu_y <- mean(y)
        n1 <- length(x)
        n2 <- length(y)
        sigma_pooled <- sqrt(((n1-1)*sd(x)^2 + (n2-1)*sd(y)^2) / (n1 + n2 - 2))
        Size <- abs(mu_x - mu_y) / sigma_pooled
        return(Size)
        }

    effect_size <- function(x) {
        M <- matrix(0,ncol(x),ncol(x))
        for (i in 1:ncol(x)) {
            for (j in 1:ncol(x)) {
                x1 <- as.numeric(x[,i])
                x2 <- as.numeric(x[,j])
                M[i,j] <- C_dist(x1,x2)
                }
            }
        return(M)
        }

    if(!is.null(userFCs)) {
        consensus <- cbind(as.numeric(userFCs[FDR5_CONSENSUS_P]), consensus)
        colnames(consensus)[1] <- "User data" 
        }

    Cohens <- effect_size(consensus)
    colnames(Cohens) <- colnames(consensus)

    pdf(paste(resDir,'/PositiveConsensus_distributions.pdf',sep=''),5,5)

    res <- lapply(1:ncol(consensus), function(x){
        currentCol<-makeTransparent(COL[x])
        if (x == 1) {
            plot(density(consensus[,x]), xlim=c(-7,1), lwd=5, xlab='gene-level depletion log FC', ylab='density', main='5% FDR Positive Controls', col=currentCol)
        } else {
            par(new=TRUE)
            plot(density(consensus[,x]), xlim=c(-7,1),lwd=5, xaxt='n',yaxt='n',xlab='',ylab='', main='', col=currentCol)}
        })

    if(!is.null(userFCs)){
        par(new=TRUE)
        plot(density(consensus[,1]),frame.plot=FALSE,xlim=c(-7,1),lwd=3,xaxt='n',yaxt='n',xlab='',ylab='',
        main='',col='black')}
    
    dev.off() 

    # Dist panel
    pdf(paste(resDir,'/Cohens_matrix.pdf',sep=''),10,10)

    if (!is.null(userFCs)) {
        layout_matrix <- matrix(c(1:49), nrow=7, ncol=7, byrow=TRUE)
        layout_matrix[lower.tri(layout_matrix)] <- c(1:21)
        layout_matrix <- t(layout_matrix)
        layout_matrix[lower.tri(layout_matrix)] <- c(22:42)
        diag(layout_matrix) <- c(43:49)
        layout(layout_matrix) 
        } else { 
            layout_matrix <- matrix(c(1:36), nrow=6, ncol=6, byrow=TRUE)
            layout_matrix[lower.tri(layout_matrix)] <- c(1:15)
            layout_matrix <- t(layout_matrix)
            layout_matrix[lower.tri(layout_matrix)] <- c(16:30)
            diag(layout_matrix) <- c(31:36)
            layout(layout_matrix)
            }
    
    par(mar=c(1, 1, 1, 1))
    for (i in 1:(ncol(consensus)-1)) {
        for (j in (i+1):ncol(consensus)) {
            plot(density(consensus[,i]), xlim=c(-7,1), lwd=2, xlab='', ylab='density', main='', col="darkgreen")
            par(new=TRUE)
            plot(density(consensus[,j]), xlim=c(-7,1),lwd=2, xaxt='n',yaxt='n',xlab='',ylab='', main='', col="blue")
            }
        }

    par(mar=c(0, 0, 0, 0))
    
    C <- as.dist(Cohens)

    for (val in C){
        plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, labels = round(val, digits=2), 
            cex = 2.2, col="magenta")
        }

    par(mar=c(0, 0, 0, 0))

    for (expName in colnames(consensus)){
        plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, labels = expName, 
            cex = 2, col = "black")
    }   
    dev.off()
}

HT29R.consensus_similarity(HT29FCsDir, resultsDir, userFCs)



res <- HT29R.FDRconsensus(HT29FCsDir,resultsDir, UserFCs, distance="Cohen's")


install.packages("VennDiagram")
library(VennDiagram)
myCol <- brewer.pal(3, "Pastel1")
Colors <- c("#CC6677", myCol[2])
Colors <- c("#117733", myCol[3])

grid.newpage()

v <- venn.diagram(
  x = list(res$POS, BAGEL_essential),
  category.names = c("Positive Consensus" , "Bagel Essential"),
  filename = paste(resultsDir,'_Venn.png',sep=''),
  output=FALSE,
  fill=Colors,
  lwd = 2,
  lty = 'blank',        
  cex = 1.8,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.8,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.pos = c(-21, 10)
  )
grid::grid.draw(v)

# Fisher test (DO IT ALSO ON negative consensus)
# take the intersection between 

TP <- length(intersect(res$Pos, BAGEL_essential)) # TP
FP <- length(intersect(res$Pos, setdiff(res$Universe, BAGEL_essential))) # FP setdiff(UNIVERSE, BAGEL_essential) instead of BAGEL_nonEssential
FN <- length(intersect(setdiff(res$Universe, res$Pos), BAGEL_essential)) # FN
TN <- length(intersect(setdiff(res$Universe, res$Pos),setdiff(res$Universe, BAGEL_essential))) # TN setdiff(UNIVERSE, BAGEL_essential) instead of BAGEL_nonEssential


dat <- data.frame(
    "col_1" = c(TP, FP),
    "col_2" = c(FN, TN),
    row.names = c("Is", "Is not"),
    stringsAsFactors = FALSE
  )
colnames(dat) <- c("consensus", "Other")
dat
fisher.test(dat)$p.value



# gsea on FDR5_POSITIVE AND FDR5_NEGATIVE
library(topGO)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)

genes <- res$Pos
genes <- res$Neg

# USING cluster Profiler
install.packages("clusterProfiler")
library(clusterProfiler)
library(enrichplot)

BPmapping <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol") 

genesUniverse <- unique(unlist(BPmapping))

GO <- enrichGO(gene = genes, keyType = "SYMBOL", universe = genesUniverse, ont="BP", OrgDb = "org.Hs.eg.db")

png(paste(resultsDir,"GO.png"), width=8*300, height=6*300, res=300, pointsize=8)
dotplot(GO, showCategory=10)
dev.off()
##########

# MAPPING
BPmapping <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol") 

genesUniverse <- unique(unlist(BPmapping))
geneList <- factor(as.integer(genesUniverse%in%genes))
names(geneList) <- genesUniverse
selection <- function(allScore) { 
  return(allScore < 0.05)
}

# TESTING
GO_BP <- new("topGOdata", 
             ontology = "BP", 
             allGenes = geneList,
             geneSel = selection,
             nodeSize = 5, 
             annot = annFUN.org, 
             mapping = "org.Hs.eg.db", 
             ID = "symbol") 


resultFisher.classic_BP <- runTest(GO_BP, algorithm = "classic", statistic = "fisher")
resultFisher.elim_BP <- runTest(GO_BP, algorithm = "elim", statistic = "fisher")

Res_BP <- GenTable(GO_BP, classicFisher = resultFisher.classic_BP, elimFisher = resultFisher.elim_BP, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 500, numChar = 100)
Res_BP <- Res_BP[as.numeric(Res_BP$classicFisher)<0.05,]
fisher_BP <- Res_BP

# PLOTTING
fisher_BP$FractionGenes <- as.numeric((fisher_BP$Significant/fisher_BP$Annotated)*100)
fisher_BP$classicFisher <- as.double(fisher_BP$classicFisher)
fisher_BP_top <- fisher_BP %>% filter(classicFisher < 0.01) %>% arrange(desc(classicFisher))

# save significant genes per GO terms
AnnotatedGenes <- lapply(fisher_BP_top$GO.ID, function(x) as.character(unlist(genesInTerm(object = GO_BP, whichGO = x))))  
SignificantGenes <- lapply(AnnotatedGenes, function(x) intersect(x, genes))
lst <- c()
myterms <- fisher_BP_top$Term

for (i in 1:length(SignificantGenes)) { 
  myterm <- myterms[i]
  mygenesforterm <- SignificantGenes[[i]]
  mygenesforterm <- paste(mygenesforterm, collapse=',')
  lst[i] <- paste(myterm,":",mygenesforterm)
}


png(paste(resultsDir,"topGO.png"), width=8*300, height=6*300, res=300, pointsize=8)
ggplot(fisher_BP_top, aes(x=FractionGenes, y=reorder(Term, -(classicFisher)), limit(0,10))) +
  geom_point(data = head(fisher_BP_top,30), aes(color=classicFisher, size=as.factor(Significant))) +
  scale_color_gradientn(colours = c("red", "blue", "lightblue"), name="p-Value\n") +
  scale_size_discrete(name="Number of\nenriched\ngenes") +
  scale_x_continuous(name="Fraction of genes (%)") +
  scale_y_discrete(name="") +
  #ggtitle("")+
  theme_light()

dev.off()


userFCs <- UserFCs

HT29R.FC_distributions <- function(refDataDir='./',resDir='./',userFCs=NULL){

    fn<-dir(refDataDir)
    fn<-grep('_foldChanges.Rdata',fn,value=TRUE)
    
    if (length(fn)==0){
        stop('No normalised sgRNA depletion fold-changes in a suitable format found in the indicated directory')    }
    
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

    pdf(paste(resDir,'/QC_FCdistproperties.pdf',sep=''),10,6)
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

    dev.off()
}

HT29R.FC_distributions(HT29FCsDir,resultsDir, UserFCs)

stats <- describe(GlobalFC, quant=c(.1,.25,.75,.90))

if(!is.null(userFCs)) {
    userStats <- stats['User data', 5:ncol(stats)]
    NC <- ncol(GlobalFC)-1
} else {
    NC <- ncol(GlobalFC)
}

AvgStats <- apply(stats[1:NC,5:ncol(stats)],2,'mean')
SE <- lapply(stats[1:NC,5:ncol(stats)],function(x){sd(x)/sqrt(NC)})

cat('\nUser screen sgRNA logFCs statistics:\n\n')
cat(paste(c('Range min: ','; Range max: '),round(c(userStats[,4], userStats[,5]), digits=3),sep=""),'\n')
cat(paste('Median: ', round(userStats[,1], digits=3), sep=""),'\n')
cat(paste(c('IQR min: ','; IQR max: '),round(c(userStats[,11], userStats[,12]), digits=3), sep=""),'\n')
cat(paste(c('10th perc: ','; 90th perc: '),round(c(userStats[,10],userStats[,13]),digits=3),sep=""),'\n')
cat(paste('Skewness: ', round(userStats[,7], digits=3), sep=""),'\n')
cat(paste('Kurtosis: ', round(userStats[,8], digits=3), sep=""),'\n')


#par(mar = c(1,1,1,3))

#plot(1,ann = F,bty = "n",type = "n")
#text(x = 1,y = 1.4, "HT29 sgRNAs logFCs statistics:",font=4)
#text(x = c(0.91,1.10), y = 1.3,paste(c('Avg. Range: ','; '),round(c(AvgStats[4],AvgStats[5]),digits=2),'±',round(c(SE$min,SE$max),digits=3),sep=''),font=3)
#text(x = 0.91, y = 1.2, paste('Avg. Median: ', round(AvgStats[1],digits=3),'±', round(SE$median,digits=3),sep=''),font=3)
#text(x = c(0.92,1.12),y = 1.1, paste(c('Avg. IQR range: ','; '), round(c(AvgStats[11],AvgStats[12]),digits=2),'±',round(c(SE$Q0.25,SE$Q.075),digits=2),sep=''),font=3)
#text(x = c(0.96,1.20),y = 1 ,paste(c('Avg. 10-90th perc range: ','; '), round(c(AvgStats[10],AvgStats[13]),digits=2),'±',round(c(SE$Q0.1,SE$Q0.9),digits=2),sep=''),font=3)
#text(x = 0.91, y=0.9, paste('Avg. Skewness: ', round(AvgStats[7],digits=2),'±',round(SE$skew,digits=2),sep=''),font=3)
#text(x =0.90, y=0.8, paste('Avg. Kurtosis: ', round(AvgStats[8],digits=2),'±',round(SE$kurtosis,digits=2),sep=''),font=3)


#if(!is.null(userFCs)) {
 #   par(mar = c(1,1,1,3))
 #   plot(1,ann = F,bty = "n",type = "n")
 #   text(x=1, y=1.4,'User screen sgRNA logFCs statistics:',font=4)
 #   text(x=c(0.94,1.15), y=1.3, paste(c('Range min: ','; Range max: '),round(c(userStats[,4], userStats[,5]), digits=3),sep=""),font=3)
 #   text(x=0.91, y=1.2, paste('Median: ', round(userStats[,1], digits=3), sep=""),font=3)
 #   text(x=c(0.91,1.08), y=1.1, paste(c('IQR min: ','; IQR max: '),round(c(userStats[,11], userStats[,12]), digits=3), sep=""),font=3)
 #   text(x=c(0.91,1.10), y=1.0, paste(c('10th perc: ','; 90th perc: '),round(c(userStats[,10],userStats[,13]),digits=3),sep=""),font=3)
 #   text(x=0.91, y=0.9, paste('Skewness: ', round(userStats[,7], digits=3), sep=""),font=3)
    #text(x=0.91, y=0.8, paste('Kurtosis: ', round(userStats[,8], digits=3), sep=""),font=3)
 #   }

#dev.off()






