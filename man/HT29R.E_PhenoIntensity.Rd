\name{HT29R.PhenoIntensity}
\alias{HT29R.PhenoIntensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  High-level quality controls (QC) using Glass’s Delta scores. 
}
\description{
Assessing screens’ performance in recovering gene essentiality profiles using sets of known essential (True Positives, TP) and known non-essential (True Negatives, TN) genes.
}
\usage{
HT29R.PhenoIntensity(refDataDir='./',
                     resDir='./',
                     userFCs=NULL, 
                     geneLevel=TRUE,
                     saveToFig=TRUE,
                     display=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refDataDir}{
      Reference HT-29 dataset directory: a string specifying the location of the processed HT-29 reference dataset.
}
  \item{resDir}{
      Output directory: a string specifying the directory where the output of this function will be saved.
}
  \item{userFCs}{
      Data from a user-provided screen: data frame with the first two columns containing sgRNAs' identifiers and HGNC symbols of targeted genes (headers = sgRNA and gene, respectivley), followed by one column per screen replicate containing sgRNAs' log fold-changes.
      }
  \item{geneLevel}{
      Boolean specifying whether the sgRNA profiles of depletion fold-changes should be collapsed by averaging them on a targeted gene basis. Default value = TRUE.
  }
  \item{saveToFig}{
      Boolean specifying whether to save all figures in PDF format. Default value = TRUE.
  }
  \item{display}{
      Boolean specifying whether to display graphics on R console. Default value = FALSE.
  }
}
\details{
This function computes Glass’s Delta scores as a metric to quantify screen quality by leveraging the intensity of the measured essentiality of prior-known essential (E), ribosomal (R) and non-essential (N) genes (see the Vignette for details).

This function generates a PDF file (named \code{'ALLSCREENS_PHENOINTENSITY.pdf'} in the specified directory with multiple plots for each HT-29 screen. If provided through the \code{userFCs} parameter, the function produces a second PDF file (named \code{'REFERENCE_PHENOINTENSITY.pdf'}) for the user-provided data.
}
\value{
}
\author{
Ichcha Manipur & Francesco Iorio (\email{fi1@sanger.ac.uk})
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{EssGenes.ribosomalProteins},
\code{BAGEL_essential},
\code{BAGEL_nonEssential}
}
\examples{
  ## Creating a temporary folder to store the HT29 reference dataset
  ## and the pdf created by this function
  dir.create('tmpdir')

  ## Downloading the HT29 reference dataset in the temporary folder
  HT29R.downloadRefData(destFolder = 'tmpdir')

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

  ## Evaluating phenotype intensity of HT29 reference and user defined data
  HT29R.PhenoIntensity(refDataDir = 'tmpDir',
                       resDir = 'tmpDir',
                       userFCs = ExampleScreen, 
                       geneLevel=TRUE,
                       saveToFig = TRUE,
                       display=FALSE)
  
  ## Checking results
  system2('open', args = 'tmpDir/ALLSCREENS_PHENOINTENSITY.pdf', wait = FALSE)
  system2('open', args = 'tmpDir/REFERENCE_PHENOINTENSITY.pdf', wait = FALSE)

  ## Removing Example dataset processed files
  file.remove('ExampleScreen_foldChanges.Rdata')
  file.remove('ExampleScreen_normCounts.Rdata')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
\keyword{benchmarking}% __ONLY ONE__ keyword per line
