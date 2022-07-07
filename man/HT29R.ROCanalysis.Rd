\name{HT29R.ROCanalysis}
\alias{HT29R.ROCanalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  High-level quality control (QC) analysis using AUROC and AUPRRC metrics.
}
\description{
Assessing screens’ performance in recovering gene essentiality profiles using sets of known essential (True Positives, TP) and known non-essential (True Negatives, TN) genes based on their sgRNA logFC rank positions while computing and visualising ROC and PrRC curves .
}
\usage{
HT29R.ROCanalysis(refDataDir='./',
                  resDir='./',
                  positives,
                  negatives,
                  userFCs=NULL, 
                  geneLevel=TRUE,
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
  \item{positives}{
      A vector of strings containing a reference set of positive cases: HGNC symbols of essential genes or identifiers of their targeting sgRNAs. This must be a subset disjointed from negatives. 
  }
  \item{negatives}{
      A vector of strings containing a reference set of negative cases: HGNC symbols of essential genes or identifiers of their targeting sgRNAs. This must be a subset disjointed from positives. 
  }
  \item{userFCs}{
      Data from a user performed screen: A data frame with the same format of the R ojbects composing the reference dataset, i.e. first two columns containing sgRNAs' identifiers and HGNC symbols of targete gene (headers = sgRNA and gene, respectivley), followed by one column per screen replicate containing sgRNAs' fold-changes.
      }
  \item{geneLevel}{
      Boolean parameter specifying whether the sgRNA profiles of depletion fold-changes should be collapsed by averaging them on a targeted gene basis. Default value = TRUE.
  }
  \item{saveToFig}{
      Boolean parameter specifying whether to save all figures in PDF format. Default value = TRUE.
  }
  \item{display}{
      Boolean parameter specifying whether to display graphics on R console. Default value = FALSE.
  }
}
\details{
This function computes AUROC and AUPRRC in order to assess the ability of each HT29 screen plus for the user defined screen (if provided through the \code{userFCs} parameter) to detect prior-known essential (positives) and non-essential (negatives) set of genes using ranked depletion fold-changes (FCs) values as predictors. 
This function produces a pdf file (named \code{'/ALLSCREENS_ROCs.pdf'} in the specified directory with multiple plots showing ROC and PrRc curves for each HT29 screen. If provided through the \code{userFCs} parameter, the function produces a second pdf file (named \code{'/USER_ROCs.pdf'}) with ROC and PrRc curves of the user defined data. If \code{display} is set to TRUE, a plot showing the overlayed curves' profiles will be displayed on R console.

}
\value{
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
Ichcha Manipur & Francesco Iorio (fi1@sanger.ac.uk)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{HT29R.PhenoIntensity},
\code{HT29R.FDRconsensus},
\code{ccr.ROC_Curve},
\code{ccr.PrRc_Curve},
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
  
  ## deriving sgRNAs targeting essential and non-essential genes (respectively)
  BAGEL_essential_sgRNAs<-ccr.genes2sgRNAs(KY_Library_v1.0,BAGEL_essential)
  BAGEL_nonEssential_sgRNAs<-ccr.genes2sgRNAs(KY_Library_v1.0,BAGEL_nonEssential)
  
  ## Computing classification performances at the sgRNA level
  HT29R.ROCanalysis(refDataDir = 'tmpDir',
                       resDir = 'tmpDir',
                       positives = BAGEL_essential_sgRNAs
                       negatives = BAGEL_nonEssential_sgRNAs
                       userFCs = ExampleScreen, 
                       geneLevel=FALSE,
                       saveToFig = TRUE,
                       display=FALSE)
  
  ## Checking results
  system2('open', args = 'tmpDir/ALLSCREENS_ROCs.pdf', wait = FALSE)
  system2('open', args = 'tmpDir/USER_ROCs.pdf', wait = FALSE)

  ## Removing Example dataset processed files
  file.remove('ExampleScreen_foldChanges.Rdata')
  file.remove('ExampleScreen_normCounts.Rdata')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
\keyword{benchmarking}% __ONLY ONE__ keyword per line
