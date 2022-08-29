\name{HT29R.FDRconsensus}
\alias{HT29R.FDRconsensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Defining a positive consensus of essential gene across the HT29 screens using the FDR ranking method.
}
\description{
The FDR method is applied for each HT29 screen, then a strict intersection of genes classified as "essential" at a defined FDR threshold (0.05 by default) is taken. The positive consensus of genes identified is used as a probe to measure the distance (Cohen's or Glass's Delta) between the depletion fold-changes (FCs) of the positive and the negative consensus of genes for each screen, which is also computed by the same function. The function uses a set of prior-known essential and non essential genes (respectively the \code{BAGEL_essential} and \code{BAGEL_nonEssential}) as reference. 
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
      Reference HT29 dataset directory: a string specifying the location of the processed HT29 reference dataset.
}
  \item{resDir}{
      Output directory: a string specifying the directory where the output of this function (a pdf file with multiple plots) should be saved.
}
  \item{userFCs}{
      Data from a user performed screen: A data frame with the same format of the R ojbects composing the reference dataset, i.e. first two columns containing sgRNAs' identifiers and HGNC symbols of targeted gene (headers = sgRNA and gene, respectivley), followed by one column per screen replicate containing sgRNAs' fold-changes.
      }
  \item{distance}{
  Type of distance the user want to measure between the two consensus: A character string indicating "Cohen's" for the Cohen's d or "GlDelta" for Glass's Delta.
  }
  \item{FDRth}{
   	If different from NULL, will be a numerical value >=0 and <=1 specifying the false discovery rate threshold at which fixed recall will be computed. Default value = 0.05.
  }
  \item{saveToFig}{
      Boolean parameter specifying whether to save all figures in PDF format. Default value = TRUE.
  }
  \item{display}{
      Boolean parameter specifying whether to display graphics on R console. Default value = FALSE.
  }
}
\details{
This function produces a pdf file (named \code{'FDR_CONSENSUS_DIST.pdf'} in the specified directory with multiple plots showing the depletion fold-changes (FCs) of the positive and negative consensus of genes for each HT29 screen plus for the user defined screen if provided through the \code{userFCs} parameter. A second plot shows the distribution of the distances computed for each screen. If data are provided by the user, it allows to measure how far the HT29-specific genes log fold-changes are from the reference.  
}
\value{
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
Raffaele Iannuzzi & Francesco Iorio (fi1@sanger.ac.uk)
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
