\name{HT29R.FDRconsensus}
\alias{HT29R.FDRconsensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Defining a positive consensus of essential genes across the HT-29 reference screens.
}
\description{
For each reference HT-29 screen a set of genes significantly depleted at a specified FDR (False Discovery Rate) threshold and its complement, i.e. a set of genes not significantly depleted, are defined using reference sets of essential and non-essential genes [1]. 

HT-29-specific essential genes are obtained as the intersection between all lists of depleted genes across each HT-29 reference screen. Conversely, HT-29-specific non-essential genes are defined as the intersection between all lists of genes not significantly depleted across each HT-29 reference screen.
}
\usage{
HT29R.FDRconsensus(refDataDir="./", 
                   resDir="./", 
                   userFCs=NULL, 
                   distance=c("GlDelta","Cohen's"), 
                   FDRth=0.05,
                   saveToFig=TRUE,
                   display=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refDataDir}{
      Reference HT-29 dataset directory: a string specifying the path of the processed HT-29 reference dataset.
}
  \item{resDir}{
      Output directory: a string specifying the directory where the output of this function will be saved.
}
  \item{userFCs}{
      Data from a user-provided screen: a data frame with the first two columns containing sgRNAs' identifiers and HGNC symbols of the targeted genes (headers = sgRNA and gene, respectivley), followed by one column per screen replicate containing sgRNAs' log fold-changes.
      }
  \item{distance}{
  Distance to measure the difference between the two consensus lists: "Cohen's" for Cohen's d or "GlDelta" for Glass's Delta.
  }
  \item{FDRth}{
   False discovery rate threshold at which fixed recall will be computed. Default value = 0.05.
  }
  \item{saveToFig}{
      Boolean specifying whether to save all figures. Default value = TRUE.
  }
  \item{display}{
      Boolean specifying whether to display graphics on the R console. Default value = FALSE.
  }
}
\details{
This function generates a PDF file (named \code{'FDR_CONSENSUS_DIST.pdf'} with multiple plots showing the depletion fold-changes (FCs) of the positive and negative consensus of genes for each HT-29 screen plus for the user screen, if provided. A second plot displays the distribution distances between known essential and non-essential genes in each screen (see the Vignette for details).
}
\value{
}
\references{
[1] Hart Traver and Moffat Jason. BAGEL: a computational framework for identifying essential genes from pooled library screens. BMC Bioinformatics. 2016;17:164.
}
\author{
Raffaele Iannuzzi & Francesco Iorio (\email{fi1@sanger.ac.uk})
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{HT29R.PhenoIntensity},
\code{HT29R.ROCanalysis},
\code{BAGEL_essential},
\code{BAGEL_nonEssential}
}
\examples{
  ## Creating a temporary folder to store the HT29 reference dataset
  ## and the pdf created by this function
  dir.create('tmpDir')

  ## Downloading the HT29 reference dataset in the temporary folder
  HT29R.downloadRefData(destFolder = 'tmpDir')

  ## Loading CRISPRcleanR library to use example screen data
  library(CRISPRcleanR)

  ## Deriving the path of the file with the example dataset,
  ## from the mutagenesis of the HT-29 colorectal cancer cell line
  fn<-paste(system.file('extdata', package = 'CRISPRcleanR'),'/HT-29_counts.tsv',sep='')

  ## Loading library Annotation
  data('KY_Library_v1.0')

  ## Loading, median-normalizing and computing fold-changes for the example dataset
  normANDfcs<-ccr.NormfoldChanges(fn,min_reads=30,EXPname='ExampleScreen',
                                libraryAnnotation = KY_Library_v1.0,
                                display = FALSE)
  
  ExampleScreen<-normANDfcs$logFCs
 
  ## Compute 5 percent FDR consensus 
  RES <- HT29R.FDRconsensus(refDataDir = 'tmpDir',
                       resDir = 'tmpDir',
                       userFCs = ExampleScreen,
                       distance = "Cohen's",
                       FDRth = 0.05,
                       saveToFig = TRUE,
                       display=FALSE)
  
  ## Checking results
  system2('open', args = 'tmpDir/FDR_CONSENSUS_DIST.pdf', wait = FALSE)

  ## Removing Example dataset processed files
  file.remove('ExampleScreen_foldChanges.Rdata')
  file.remove('ExampleScreen_normCounts.Rdata')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
\keyword{benchmarking}% __ONLY ONE__ keyword per line
