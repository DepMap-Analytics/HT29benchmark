\name{HT29R.evaluateReps}
\alias{HT29R.evaluateReps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Low-level quality control (QC) analysis: reproducibility across technical replicates.
}
\description{
Assessing CRISPR-Cas9 genome-wide pooled screens' reproducibility using Project Score metrics and criteria.
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
      Reference HT29 dataset directory: a string specifying the location of the processed HT29 reference dataset.
}
  \item{resDir}{
      Output directory: a string specifying the directory where the output of this function (a pdf file with multiple plots) should be saved.
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
This function computes correlation scores between each pair of screen replicates for the HT29 reference dataset as well as for used defined data.
This is performed as for the Project Score data [1]: considering only a set of 838 most informative sgRNAs (in the \code{HT29R.reproducible_GeneGuides object}), defined as those targeting the same genes and with an average pairwise Pearson’s correlation > 0.6 between corresponding patterns of depletion fold-changes (FCs) across hundreds of screened cell lines. Per construction, the depletion patterns of these sgRNAs are both reproducible and informative (as they involve genes carrying an actual fitness signal). Computing correlation scores between replicates of the same screen on the domain of these sgRNAs only allows estimating a null distribution of replicate correlations and computing a reproducibility threshold defined as the minimal correlation score that should be observed between replicates of the same screen (R = 0.68, using the
\code{HT29R.GL_prSCORE_rCorr_Reprod} object).
This is performed as genome-wide correlation scores computed between replicates of the same CRISPR-Cas9 pooled genome-wide viability screen are generally always very high and indistiguishable from expectation due to only a small percentage of genes exerting an effect on cellular fitness upon knock-out.

This function produces a pdf file showing the plot of the expected/observed distributions of replicate pair-wise correlation scores and reproducibility threshold from the Project Score dataset. Below this plots, the pair-wise correlation scores computed between replicates of the HT29 reference dataset are also plotted. If the user provides its own data (trough the \code{userFCs} parameter) then pair-wise correlation scores computed between replicates of this screens are overlaied on the kernel distributions. In this case, a specification of the number of replicates reaching threshold is printed in the console.
}
\value{
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
Ichcha Manipur, Raffaele Iannuzzi, Francesco Iorio (fi1@sanger.ac.uk)
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
