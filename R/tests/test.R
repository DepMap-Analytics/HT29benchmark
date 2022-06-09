# install packages
install.packages('dplyr')
install.packages('moments')
install.packages('psych')
install.packages("devtools")
install_github("francescojm/CRISPRcleanR")
install.packages("vioplot")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MAGeCKFlute")

# libraries
library(crayon)
library(dplyr)
library(moments)
library(psych)
library(devtools)
library(CRISPRcleanR)
library(RColorBrewer)
library(MAGeCKFlute)
library(vioplot)
library(data.table)
library(reshape2)

# set path
filesDir <- "/Users/raffaele.iannuzzi/Documents/Project-HT29/HT29benchmark/data/"
HT29FCsDir <- "/Users/raffaele.iannuzzi/Documents/Project-HT29/user_data/HT29FCsDATA/"
resultsDir <- "/Users/raffaele.iannuzzi/Documents/Project-HT29/resDir/"

# load data
fn <- dir(filesDir)
fn <- grep('HT29R.',fn,value=TRUE)
for (f in fn) {
  cat(paste('loading',f,'\n'))
  load(paste(filesDir,f,sep=""))
  }

# download data
HT29R.download_ref_dataset(HT29FCsDir)

# load KY lib
pathToKYLib <- paste(system.file('data', package = 'CRISPRcleanR'),'/KY_Library_v1.0.Rdata', sep='')
load(pathToKYLib)

# test HT29.FC_dist_properties
test <- HT29R.FCdistributions(refDataDir = HT29FCsDir,
                          resDir = resultsDir,
                          userFCs = UserFCs,
                          saveToFig = TRUE)

# test HT29R.evaluate_reps function
HT29R.evaluateReps(refDataDir = HT29FCsDir,
                    resDir = resultsDir,
                    userFCs = UserFCs,
                    geneLevel = FALSE, 
                    saveToFig = TRUE) 

# test HT29R.exp_similarity function

# all and FALSE
HT29R.expSimilarity(refDataDir = HT29FCsDir, 
                    resDir = resultsDir,
                    geneGuides = "All",
                    geneLevel = FALSE,
                    userFCs = userFCs,
                    saveToFig = TRUE)

# all and TRUE
RES <- HT29R.expSimilarity(refDataDir = HT29FCsDir, 
                    resDir = resultsDir,
                    geneGuides = "All",
                    geneLevel = TRUE,
                    Rscore=TRUE,
                    userFCs = UserFCs,
                    saveToFig = TRUE)

# HI and FALSE
HT29R.expSimilarity(refDataDir = HT29FCsDir, 
                    resDir = resultsDir,
                    geneGuides = "HI",
                    geneLevel = FALSE,
                    userFCs = UserFCs,
                    saveToFig = TRUE)

# HI and TRUE
RES <- HT29R.expSimilarity(refDataDir = HT29FCsDir, 
                    resDir = resultsDir,
                    geneGuides = "HI",
                    geneLevel = TRUE,
                    userFCs = UserFCs,
                    saveToFig = TRUE)


# test HT29R.PhenoIntensity
# geneLevel FALSE
HT29R.PhenoIntensity(refDataDir = HT29FCsDir,
                    resDir = resultsDir,
                    userFCs = UserFCs, 
                    geneLevel=FALSE)
# geneLevel TRUE
HT29R.PhenoIntensity(refDataDir = HT29FCsDir,
                    resDir = resultsDir,
                    userFCs = UserFCs, 
                    geneLevel=TRUE)

# test HT29R.ROCanalysis function
#sgRNA
positives <-ccr.genes2sgRNAs(KY_Library_v1.0, BAGEL_essential)
negatives <-ccr.genes2sgRNAs(KY_Library_v1.0, BAGEL_nonEssential)

HT29R.ROCanalysis(refDataDir = HT29FCsDir,
                  resDir = resultsDir,
                  positives = positives,
                  negatives = negatives,
                  userFCs = UserFCs,
                  geneLevel = FALSE,
                  saveToFig = TRUE)


# genes
positives <- BAGEL_essential
negatives <- BAGEL_nonEssential

HT29R.ROCanalysis(refDataDir = HT29FCsDir,
                  resDir = resultsDir,
                  positives = positives,
                  negatives = negatives,
                  userFCs = UserFCs,
                  geneLevel = TRUE,
                  saveToFig = TRUE)


res <- HT29R.FDRconsensus(refDataDir = HT29FCsDir,
                    resDir = resultsDir,
                    userFCs = UserFCs,
                    distance = "Cohen's",
                    FDRth = 0.05, 
                    saveToFig = FALSE)


data(EssGenes.ribosomalProteins)
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.HISTONES)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)
data(BAGEL_essential)
data(BAGEL_nonEssential)


EssGenes <- list(EssGenes.ribosomalProteins,
               BAGEL_essential,
               EssGenes.DNA_REPLICATION_cons,
               EssGenes.HISTONES,
               EssGenes.KEGG_rna_polymerase,
               EssGenes.PROTEASOME_cons,
               EssGenes.SPLICEOSOME_cons)

names(EssGenes) <- c("Ribosomal Proteins genes",
                      "BAGEL essential",
                      "DNA replication genes",
                      "Histones genes",
                      "KEGG_RNA_polymerase",
                      "Proteasome genes",
                      "Spliceosome genes")

POS <- RES[which(rownames(RES) %in% res$POS),]
group <- rep("HT-29-specific genes", nrow(POS))
POS <- cbind.data.frame(POS, group)

RIBO <- RES[which(rownames(RES) %in% EssGenes.ribosomalProteins),]
group <- rep("Ribosomal Proteins genes", nrow(RIBO))
RIBO <- cbind.data.frame(RIBO, group)

ESS <- RES[which(rownames(RES) %in% BAGEL_essential),]
group <- rep("BAGEL essential", nrow(ESS))
ESS <- cbind.data.frame(ESS, group)

NNESS <- RES[which(rownames(RES) %in% BAGEL_nonEssential),]
group <- rep("BAGEL non-essential", nrow(NNESS))
NNESS <- cbind.data.frame(NNESS, group)

DF <- rbind(POS,RIBO,ESS,NNESS)
group <- DF$group
DF$group <- factor(group)
DF <- melt(DF)


pdf(paste(resultsDir,'/Genesets_logFCs.pdf',sep=''),20,8)
par(mar=c(5,5,0,1))
tmp <- vioplot(value~group*variable,
                data = DF,
                frame.plot = FALSE,
                xaxt = 'n',
                xlab='',
                ylab = "log fold-change",
                ylim = c(-7,3),
                cex.axis=1.2,
                las=1,
                col=c("#B3CDE3","#888888","#CC6677","darkblue"),
                pch=14,
                pars=list(outcol=c("#B3CDE3","#888888","#CC6677","darkblue")),
                horizontal = FALSE,
                at=seq(1,40)[-c(5,6,11,12,17,18,23,24,29,30,35,36)]
                )

axis(1, at=(1:47)[c(2,8,14,20,26,32,38)], labels=FALSE)
text(x=(1:47)[c(2,8,14,20,26,32,38)], 
             y=-7.7,
             labels=unique(DF$variable), 
             srt=45, 
             adj=1, 
             xpd=TRUE,cex=1.2)

abline(h=median(DF[which(DF$group == "HT-29-specific genes"),]$value), col="#CC6677",lty=2,lwd=2)
abline(h=median(DF[which(DF$group == "Ribosomal Proteins genes"),]$value), col="darkblue",lty=2,lwd=2)
abline(h=median(DF[which(DF$group == "BAGEL essential"),]$value), col="#B3CDE3",lty=2,lwd=2)
abline(h=median(DF[which(DF$group == "BAGEL non-essential"),]$value), col="#888888",lty=2,lwd=2)
dev.off()


FDRth <- 0.25

TrueNegatives <- BAGEL_nonEssential
TruePositives <- EssGenes.SPLICEOSOME_cons

predictions <- RES[intersect(c(TruePositives, TrueNegatives),row.names(RES)),]
observations <- is.element(row.names(predictions), TruePositives)+0

res <- lapply(1:ncol(RES), function(x) {
  currentFc <- predictions[, x]
  names(observations) <- names(currentFc)
  res <- roc(observations, currentFc, direction = ">", quiet = TRUE) 
  COORS <- coords(res, "all", ret = c("threshold", "ppv", "sensitivity"), transpose=TRUE) 
  FDRpercTh <- max(COORS['threshold', which(COORS['ppv',] >= (1 - FDRth))])
  threshold <- COORS["threshold", min(which(COORS["threshold",] <= FDRpercTh))]
  recall <- COORS["sensitivity", min(which(COORS["threshold",] <= threshold))]
  list(RC=recall)
  })

res

median(unlist(res))
quantile(unlist(res))

signGenes <- lapply(1:ncol(RES),function(x){names(RES[RES[,x] <= res[[x]]$TH, x])})
NsignGenes <- lapply(1:length(signGenes), function(x){length(signGenes[[x]])})
NsignGenes
median(unlist(NsignGenes))
quantile(unlist(NsignGenes))


df1 <- data.frame(
  FDR="1%",
  TH=as.numeric(unlist(res))
)
df2 <- data.frame(
  FDR="5%",
  TH=as.numeric(unlist(res))
)
df3 <- data.frame(
  FDR="10%",
  TH=as.numeric(unlist(res))
)
df4 <- data.frame(
  FDR="25%",
  TH=as.numeric(unlist(res))
)


df <- rbind(df1,df2,df3,df4)
group <- factor(df$FDR,levels =c("1%","5%","10%","25%"))
df$FDR <- NULL
df <- cbind.data.frame(group, df)
df

pdf(paste(resultsDir,'/recall_Spliceosome.pdf',sep=''),6,7)
par(mar=c(4,8,2,8))
boxplot(TH~FDR,
        data = df,
        ylab="% Recall",
        col=c("#332288", "#AA4499","#44AA99", "#999933"),
        pch=14,
        pars=list(outcol=c("#332288", "#AA4499","#44AA99", "#999933")),
        las=1)
dev.off()

m <- apply(RES[res$POS,], MARGIN = 1, FUN = mean, na.rm = TRUE)
o <- m[order(m,decreasing = FALSE)]

pdf(paste(resultsDir,'/top50_logFCs.pdf',sep=''),6,7)
par(mar=c(4,4,1,4))
boxplot.matrix(RES[names(o),][1:50,], col="white", border="black", whiskcol="black", use.cols = FALSE,ylab="log fold-change",xaxt="n",xlab="Rank",las=1,cex.axis=1.2)
dev.off()





pdf(paste(resultsDir,'/test_density_user_Replicates.pdf',sep=''),10,10)
par(mfrow=c(1,1))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
COL = tail(col_vector,ncol(userFCs)-2)
names(COL) <- tail(colnames(userFCs), ncol(userFCs)-2) 
currentCol <- makeTransparent(COL)
idx <- 1
legendspace <- 0.02
x_lim <- c(-9,3)
y_lim <- c(0, 0.8)
for (j in 3:ncol(userFCs)){
  if (j==3) {
    plot(density(userFCs[,j]),frame.plot=FALSE, xlim=x_lim, ylim=y_lim, col=currentCol[idx],lwd=7,main='',xlab='sgRNA fc')
    legend(x=x_lim[1],y=y_lim[2],names(COL)[idx],lwd = 3, bty='n', col = currentCol[idx])
    } else { 
      idx <- idx + 1
      legendspace <- legendspace * 2
      par(new=TRUE)
      plot(density(userFCs[,j]),frame.plot=FALSE,xlim=x_lim,col=currentCol[idx],lwd=7,xaxt='n',yaxt='n',xlab='',ylab='',
      main='')
      legend(x=x_lim[1],y=y_lim[2]-legendspace,names(COL)[idx],lwd = 3, bty='n', col = currentCol[idx])
    }
}      

dev.off()


### compute effect size 

C_dist <- function(x,y) {
  mu_x <- mean(x)
  mu_y <- mean(y)
  n1 <- length(x)
  n2 <- length(y)
  sigma_pooled <- sqrt(((n1-1)*sd(x)^2 + (n2-1)*sd(y)^2) / (n1 + n2 - 2))
  Size <- abs(mu_x - mu_y) / sigma_pooled
  return(Size)
 }

effectSize <- function(FCsprofile, positives) {
  dat <- FCsprofile[positives,]
  m <- matrix(0,ncol(dat),ncol(dat))
  for (i in 1:ncol(dat)) {
    for (j in 1:ncol(dat)) {
      x <- as.numeric(dat[,i])
      y <- as.numeric(dat[,j])
      m[i,j] <- C_dist(x,y)
     }
   }
   return(m)
 }

if(!is.null(userFCs)) {
  ref_data <- cbind(as.numeric(userFCs[FDR5_CONSENSUS_P]), ref_data)
  colnames(ref_data)[1] <- "User data"
}

Cohens <- effectSize(ref_data, positives = FDR5_CONSENSUS_P)

pdf(paste(resultsDir,'/Cohens_FDR5PositiveConsensus.pdf',sep=''),10,10)

res <- lapply(1:ncol(ref_data), function(x){
  currentCol<-makeTransparent(COL[x])
  if (x == 1) {
    plot(density(ref_data[,x]), xlim=c(-7,1), lwd=5, xlab='5% FDR Positive Controls (log FC)', ylab='density', main='title', col=currentCol)
    } else {
        par(new=TRUE)
        plot(density(ref_data[,x]), xlim=c(-7,1),lwd=5, xaxt='n',yaxt='n',xlab='',ylab='', main='', col=currentCol)}
  }
)

if(!is.null(userFCs)){
  par(new=TRUE)
  #uData <- as.numeric(userFCs[FDR5_CONSENSUS_P])
  plot(density(ref_data[,1]),frame.plot=FALSE,xlim=c(-7,1),lwd=3,xaxt='n',yaxt='n',xlab='',ylab='',
           main='',col='black') 
  } 

# Dist panel

N1 <- 6
N2 <- 15
N3 <- 30

layout_matrix <- matrix(c(1:36), nrow=6, ncol=6, byrow=TRUE)
layout_matrix[lower.tri(layout_matrix)] <- c(1:15)
layout_matrix <- t(layout_matrix)
layout_matrix[lower.tri(layout_matrix)] <- c(16:30)
diag(layout_matrix) <- c(31:36)
layout(layout_matrix)

if (!is.null(userFCs)) {
  layout_matrix <- matrix(c(1:49), nrow=7, ncol=7, byrow=TRUE)
  layout_matrix[lower.tri(layout_matrix)] <- c(1:21)
  layout_matrix <- t(layout_matrix)
  layout_matrix[lower.tri(layout_matrix)] <- c(22:42)
  diag(layout_matrix) <- c(43:4)
  layout(layout_matrix)
}

par(mar=c(1, 1, 1, 1))
for (i in 1:(ncol(ref_data)-1)) {
  for (j in (i+1):ncol(ref_data)) {
      plot(density(ref_data[,i]), xlim=c(-7,1), lwd=2, xlab='', ylab='density', main='', col="darkgreen")
      par(new=TRUE)
      plot(density(ref_data[,j]), xlim=c(-7,1),lwd=2, xaxt='n',yaxt='n',xlab='',ylab='', main='', col="blue")
      
    }
  }

par(mar=c(0, 0, 0, 0))
m <- as.dist(Cohens)
for (val in m){
  plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, labels = round(val,digits=2), 
     cex = 2, col="magenta")
  }

par(mar=c(0, 0, 0, 0))
for (expName in colnames(ref_data)){
  plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, labels = expName, 
     cex = 2, col = "black")
}

dev.off()


# MAGECK

fn <- dir(HT29FCsDir)
fn <- grep('.tsv',fn,value=TRUE)

RES <- lapply(fn, function(x) { 
  normANDfcs <- ccr.NormfoldChanges(
    paste(HT29FCsDir,x,sep=''),
    method = 'MedRatios',
    display = FALSE,
    min_reads=30,
    EXPname=unlist(lapply(str_split(x, '.tsv'), function(y){y[1]})),
    libraryAnnotation = KY_Library_v1.0)
  
  fname <- ccr.PlainTsvFile(
    sgRNA_count_object = normANDfcs$norm_counts,
    fprefix = unlist(lapply(str_split(x,'.tsv'),
    function(y){y[1]})),
    path = resultsDir)
  
  ccr.ExecuteMageck(fname,
                    expName = unlist(lapply(str_split(x, '.tsv'), 
                    function(y){y[1]})),
                    normMethod = 'none', 
                    outputPath = resultsDir)
    })



if (method == "fixedFDR") {
  RES <- lapply(seq(ncol(ref_fcs)), function(x) { 
      currentFc <- predictions[, x]
      names(observations) <- names(currentFc)
      roc(observations, currentFc, direction = ">", quiet = TRUE)
      })
    }

ROCvalues <- sapply(1:length(RES), function(x){unlist(lapply(RES[[x]]$auc, unique))})
names(ROCvalues) <- colnames(ref_fcs)
x <- sort(ROCvalues,decreasing=TRUE)
plot(x, frame.plot=FALSE, type="o", col="blue",lwd=2,ylab= "AUROC", xlab="",xaxt="n", main="")
axis(1, at=1:ncol(ref_fcs), labels=FALSE)
text(x=1:ncol(ref_fcs), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
      labels=names(x), srt=45, adj=1, xpd=TRUE)


NC <- ncol(ref_fcs)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
COL <- col_vector[1:NC]
currentCol<-makeTransparent(COL[x])

par(mfrow=c(1,ncol(ref_fcs)))
firstEl <- ref_fcs[names(RES[[1]]$AUC),1]

for (i in 1:ncol(ref_fcs)) {
  x <- ref_fcs[names(RES[[i]]$AUC),i]
  if (x == firstEl) {
    boxplot(x,frame.plot = FALSE,ylim=c(-9,2),lwd=2,col=makeTransparent('black'),ylab= "log fold-changes",main="")
  } else {
    boxplot(x,frame.plot = FALSE,lwd=2,col=makeTransparent(COL[i]),xaxt="n",yaxt="n",xlab="",ylab="")
    }

}


