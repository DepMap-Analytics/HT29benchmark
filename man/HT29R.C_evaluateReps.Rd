\name{HT29R.evaluateReps}
\alias{HT29R.evaluateReps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Low-level quality control (QC) analysis: reproducibility across technical replicates.
}
\description{
Assessing CRISPR-Cas9 genome-wide pooled screens' reproducibility.
}
\usage{
HT29R.evaluateReps(refDataDir='./',
                   resDir='./',
                   userFCs=NULL, 
                   geneLevel=TRUE,
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
      Data from a user-provided screen: a data frame with the first two columns containing sgRNAs' identifiers and HGNC symbols of targeted genes (headers = sgRNA and gene, respectively), followed by one column per screen replicate containing sgRNAs' log fold-changes.
      }
  \item{geneLevel}{
      Boolean specifying whether the sgRNA profiles of depletion log fold-changes should be collapsed by averaging them on a targeted gene basis. Default value = TRUE.
  }
  \item{saveToFig}{
      Boolean specifying whether to save all figures. Default value = TRUE.
  }
  \item{display}{
      Boolean specifying whether to display graphics on the R console. Default value = FALSE.
  }
}
\details{
This function computes intra-screen replicate Pearson's correlation scores for the HT-29 reference screen, as well as for the user-provided screen. This is performed as for the Project Score data [1]: considering only the set of the 838 most informative sgRNAs (in the \code{HT29R.reproducible_GeneGuides object}).

Particularly, it generates a plot of Pearson’s correlation scores computed between sgRNA depletion log fold-change profiles of replicates for each of the six HT-29 reference screens, and between replicates of the user-provided screen (see the Vignette for details).
}
\value{
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
Ichcha Manipur, Raffaele Iannuzzi, Francesco Iorio (\email{fi1@sanger.ac.uk})
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{H29R.downloadRefData},
\code{HT29R.reproducible_GeneGuides},
\code{HT29R.prSCORE_rCorr},
\code{HT29R.GL_prSCORE_rCorr}
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

  ## Evaluating screen reproducibility of HT29 reference and user defined data
  ## both, using Project Score criteria
  HT29R.evaluateReps(refDataDir='tmpDir',
                     resDir='tmpDir',
                     userFCs=ExampleScreen, 
                     geneLevel=TRUE,
                     saveToFig=TRUE, 
                     display=FALSE)

  ## Checking results
  system2('open', args = 'tmpDir/QC_REPRODUCIBILITY_GENElevel.pdf',sep=""), wait = FALSE)

  ## Removing Example dataset processed files
  file.remove('ExampleScreen_foldChanges.Rdata')
  file.remove('ExampleScreen_normCounts.Rdata')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
\keyword{benchmarking}% __ONLY ONE__ keyword per line
