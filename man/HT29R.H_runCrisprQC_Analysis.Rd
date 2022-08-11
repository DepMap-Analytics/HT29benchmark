\name{HT29R.runCrisprQC_Analysis}
\alias{HT29R.runCrisprQC_Analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Wrapper function: download reference data and run all QC analysis from low-level to high-level.
}
\description{
}
\usage{
HT29R.runCrisprQC_Analysis(data = c("FCs","rawCounts") , 
                           userFCs=NULL,
                           positives,
                           negatives,
                           outdir="./")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
      String parameter specifying what type of data to dowload.
      Possible values are "rawCounts" for plain .tsv files containing raw sgRNA countsare or "FCs" (default) for R objects containing sgRNA normalised depletion fold-changes: data frames in which the first two columns contain sgRNAs' identifiers and HGNC symbols of targeted gene, followed by one column per screen replicate containing sgRNAs' fold-changes.
}
  \item{userFCs}{
    Data from a user performed screen: A data frame with the same format of the R ojbects composing the reference dataset, i.e. first two columns containing sgRNAs' identifiers and HGNC symbols of targete gene (headers = sgRNA and gene, respectivley), followed by one column per screen replicate containing sgRNAs' fold-changes.
}
  \item{positives}{
     A vector of strings containing a reference set of positive cases: HGNC symbols of essential genes or identifiers of their targeting sgRNAs. This must be a subset disjointed from negatives. 
      }
  \item{negatives}{
    A vector of strings containing a reference set of negative cases: HGNC symbols of essential genes or identifiers of their targeting sgRNAs. This must be a subset disjointed from positives.
  }
  \item{outdir}{
       Output directory: a string specifying the directory where the output of this function (a pdf file with multiple plots) should be saved
       }
}
\details{
The function create a folder named \code{"PLOTS"} in which the figures will be stored. The code will ask the user to provide two parameters whose are necessary for the low-level QC analysis. The first define the subset of guides the user intends to use (i.e., "All" for the sgRNAs targeting all genes in the genome-wide library and "HI" for the sgRNAs targeting only the most-informative and highly-reproducible ones as described in [1]), the second define the level of inspection (i.e., TRUE if inspecting at the gene-level or FALSE if at the guide-level). 

}
\value{
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
Raffaele M. Iannuzzi & Francesco Iorio (fi1@sanger.ac.uk)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{HT29R.downloadRefData},
\code{HT29R.evaluateReps},
\code{HT29R.expSimilarity},
\code{HT29R.PhenoIntensity},
\code{HT29R.ROCanalysis},
\code{HT29R.FDRconsensus},
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
 
  ## Run all with interactive wrapper 
  
  HT29R.runCrisprQC_Analysis(data = "FCs",
                             userFCs = ExampleScreen,
                             positives = BAGEL_essential,
                             negatives = BAGEL_nonEssential,
                             outdir = 'tmpDir')
  
  ## Checking results
  system2('open', args = 'tmpDir/PLOTS/', wait = FALSE)

  ## Removing Example dataset processed files
  file.remove('ExampleScreen_foldChanges.Rdata')
  file.remove('ExampleScreen_normCounts.Rdata')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
\keyword{benchmarking}% __ONLY ONE__ keyword per line
