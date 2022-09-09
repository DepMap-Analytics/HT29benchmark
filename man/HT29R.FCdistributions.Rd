\name{HT29R.FCdistributions}
\alias{HT29R.FCdistributions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Screens depletion fold-changes (FCs) average statistics and distributions.
}
\description{
  The function computes sgRNA logFCs statistics and distributions for the reference HT-29 screening obtained with reagent and experimental settings described in [1], and for an optional user-provided screen.
} 
\usage{
HT29R.FCdistributions(refDataDir='./',
                      resDir='./',
                      userFCs=NULL,
                      stats=TRUE,
                      saveToFig=TRUE,
                      display=FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refDataDir}{
      Reference HT-29 dataset directory: a string specifying the path of the processed HT-29 reference dataset.
  }
  \item{resDir}{
      Output directory: a string specifying the directory where the output plots should be saved.
  }
  \item{userFCs}{
      Data from a user-provided screen: a data frame with the first two columns containing sgRNAs' identifiers and HGNC symbols of targeted genes (headers = sgRNA and gene, respectively), followed by one column per screen replicate containing sgRNAs' log fold-changes.
  }
  \item{stats}{
      Boolean specifying whether to show average statistics. Default is TRUE.
  }
  \item{saveToFig}{
      Boolean specifying whether to save all figures. Default is TRUE.
  }
  \item{display}{
      Boolean specifying whether to display graphics on the R console. Default is FALSE.
  }
}
\details{
  The function computes sgRNAs logFCs’ statistics such as average range, median, inter-quartile range, 10th-90th percentile range, skewness and kurtosis across the six HT-29 screens plus the user-defined screen, if provided.
}

\value{
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
Raffaele Iannuzzi, Ichcha Manipur, Francesco Iorio (\email{fi1@sanger.ac.uk})
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{H29R.downloadRefData}
}
\examples{
  ## Creating a temporary folder to store the HT29 reference dataset
  ## and the pdf created by this function
  dir.create('tmpDir')

  ## Downloading the HT29 reference dataset in the temporary folder
  HT29R.downloadRefData(destFolder = 'tmpdir')

  ## Loading CRISPRcleanR library to use example screen data
  library(CRISPRcleanR)

  ## Deriving the path of the file with the example dataset,
  ## from the mutagenesis of the HT-29 colorectal cancer cell line
  fn <- dir('tmpDir')
  fn <- grep('_foldChanges.Rdata',fn,value=TRUE)

  ## Loading library Annotation
  data('KY_Library_v1.0')

  ## Loading, median-normalizing and computing fold-changes for the example dataset
  normANDfcs<-ccr.NormfoldChanges(fn,min_reads=30,EXPname='ExampleScreen',
                                libraryAnnotation = KY_Library_v1.0,
                                display = FALSE)
  
  ExampleScreen<-normANDfcs$logFCs

  ## Evaluating screen reproducibility of HT29 reference and user defined data
  ## both, using Project Score criteria
  HT29R.FCdistributions(refDataDir='tmpDir',
                        resDir='tmpDir',
                        userFCs=ExampleScreen,
                        stats=TRUE,
                        saveToFig=TRUE, 
                        display=FALSE)

  ## Checking results
  system2('open', args = 'tmpDir/QC_FCdistproperties.pdf', wait = FALSE)
  
  
  ## Removing Example dataset processed files
  file.remove('ExampleScreen_foldChanges.Rdata')
  file.remove('ExampleScreen_normCounts.Rdata')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
\keyword{benchmarking}% __ONLY ONE__ keyword per line
