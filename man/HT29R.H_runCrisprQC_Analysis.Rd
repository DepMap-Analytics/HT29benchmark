\name{HT29R.runCrisprQC_Analysis}
\alias{HT29R.runCrisprQC_Analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Wrapper function: download the HT-29 reference data and run all quality control (QC) analyses.
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
      String specifying which data to dowload.
      Possible values are:
      \itemize{
      \item "rawCounts" to download .tsv files containing raw sgRNA counts; 
      \item "FCs" (default) to download .RData objects containing data.frames with normalised sgRNA depletion log fold-changes. The first two columns contain sgRNAs' identifiers and HGNC symbols of targeted gene, followed by one column per screen replicate with sgRNAs' log fold-changes;}
}
}
  \item{userFCs}{
    User-provided data: data frame with the first two columns containing sgRNAs' identifiers and HGNC symbols of targeted genes (headers = sgRNA and gene, respectivley), followed by one column per screen replicate containing sgRNAs' log fold-changes.
}
  \item{positives}{
     A vector with known essential genes (HGNC symbols or identifiers of their targeting sgRNAs) There must be no overlapping genes/sgRNAs between positives and negatives. 
      }
  \item{negatives}{
A vector with known non-essential genes (HGNC symbols or identifiers of their targeting sgRNAs) There must be no overlapping genes/sgRNAs between positives and negatives.}
  \item{outdir}{
       Output directory: directory where the output of this function will be saved
       }
}
\details{
The function creates a folder named \code{"PLOTS"} where the figures will be stored. The user is asked to provide two parameters necessary for low-level QC analysis. The first defines the subset of guides the user intends to use (i.e., "All" for all sgRNAs and "HI" for the most informative sgRNAs, selected as described in [1]); the second parameter defines the inspection level (i.e., TRUE if inspecting at the gene level or FALSE if at the sgRNA level). 
}
\value{
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
Raffaele M. Iannuzzi & Francesco Iorio (\email{fi1@sanger.ac.uk})
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
