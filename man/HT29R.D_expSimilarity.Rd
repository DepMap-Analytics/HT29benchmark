\name{HT29R.expSimilarity}
\alias{HT29R.expSimilarity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Low-level quality control (QC) analysis: quantifying screens' similarity. 
}
\description{
Assessing CRISPR-Cas9 genome-wide pooled screens' comparability using Project Score data.
}
\usage{
HT29R.expSimilarity(refDataDir='./',
                    resDir='./',
                    userFCs=NULL, 
                    geneGuides=c("All","HI"), 
                    geneLevel=TRUE,
                    Rscores=TRUE,
                    saveToFig=TRUE,
                    display=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refDataDir}{
      Reference HT-29 dataset directory: a string specifying the path of the processed HT-29 reference dataset.
}
  \item{resDir}{
      Output directory: a string specifying the path where the output of this function will be saved.
}
  \item{userFCs}{
      User-provided data: data frame with the first two columns containing sgRNAs' identifiers and HGNC symbols of targeted genes (headers = sgRNA and gene, respectivley), followed by one column per screen replicate containing sgRNAs' log fold-changes.
}
  \item{geneGuides}{
      String specifying whether to use all sgRNAs (genome-wide) or just the most reproducible and informative sgRNAs (n=838)
}
  \item{geneLevel}{
      Boolean specifying whether the sgRNA profiles of depletion log fold-changes should be averaged gene-wise. Default value = TRUE.
  }
  \item{Rscores}{
      Boolean specifying whether to compute correlation scores matrix between averaged screen replicates for the HT-29 reference dataset as well as for user-provided data. Default value = TRUE.
  }
  \item{saveToFig}{
      Boolean specifying whether to save all figures. Default value = TRUE.
  }
  \item{display}{
      Boolean specifying whether to display graphics on the R console. Default value = FALSE.
  }
}
\details{
This function compares inter-screen replicate Pearson's correlation scores for the HT-29 reference screen, as well as for the user-provided screen. This is performed as for the Project Score data [1]: considering only the set of the 838 most informative sgRNAs (in the \code{HT29R.reproducible_GeneGuides object}).

Particularly, it generates a plot of Pearson’s correlation scores computed between sgRNA depletion log fold-change profiles of replicates for each of the six HT-29 reference screens, and between replicates of the user-provided screen (see the Vignette for details).

This function computes correlation scores between each averaged screens' replicate across the HT29 reference dataset and between each averaged HT29 screens' replicate and user defined data (through the \code{userFCs} parameter) in order to measure depletion fold-changes (FCs) similarities across the HT29 screens and the user defined depletion fold-changes patterns. The analysis can be performed considering the set of 838 most informative sgRNAs only (named \code{HT29R.reproducible_GeneGuides} as seen for the \code{HT29R.evaluateReps} function, or considering the sgRNAs consensus (named \code{HT29R.consensus_GeneGuides}) across the reference and for the user defined screen. If the \code{geneLevel} parameter is set to TRUE, HGNC genes' symbols will be used instead of sgRNAs identifiers. The background distributions will be chosen based on the combination of parameters selected (e.g., if \code{geneGuides} parameter is set to "All" and \code{geneLevel} is set to TRUE, \code{HT29R.prSCORE_bkgr_screen_similarity} will be selected instead of \code{HT29R.prSCORE_bkgr_screen_similarity_sgRNA}). This function produces two pdf files (named '_QC_SIMILARITY_GENElevel.pdf' or '_QC_SIMILARITY_sgRNAlevel.pdf' and '_SCATTERPLOT_R_MATRIX_GENElevel.pdf' or '_SCATTERPLOT_R_MATRIX_sgRNAlevel.pdf') in the specified directory. The first pdf contains a plot of the kernel distribution of Pearson's correlation scores computed using the selected sgRNAs (through the \code{geneGuides} parameter) between pairs of HT29 screens and between reference screens versus the user defined data (if provided through the \code{userFCs} parameter). The Project Score background distribution is represented in grey. The same distributions are represented as boxplots coupled with two-sided t-test results for each comparison considered (Background vs Reference and, if provided, User-data vs Reference). If \code{Rscores} parameter is set to TRUE, the function saves scatterplot correlation matrix of the pair-wise Pearson’s correlation scores computed within HT-29 references and, if provided, between user versus HT-29 references.
}

\value{
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
Raffaele Iannuzzi, Ichcha Manipur, Francesco Iorio (fi1@sanger.ac.uk)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{HT29R.downloadRefData},
\code{HT29R.reproducible_GeneGuides},
\code{HT29R.consensus_GeneGuides},
\code{HT29R.prSCORE_bkgr_screen_similarity},
\code{HT29R.prSCORE_bkgr_screen_similarity_HI},
\code{HT29R.prSCORE_bkgr_screen_similarity_sgRNA},
\code{HT29R.prSCORE_bkgr_screen_similarity_sgRNA_HI}
}
\examples{

  ## Creating a temporary folder to store the HT29 reference dataset
  ## and the pdf created by this function
  dir.create('tmpDir')

  ## Downloading the HT29 reference dataset in the temporary folder
  HT29R.downloadRefData(destFolder = 'tmpDir')

  ## Loading CRISPRcleanR library to use example screen data
  library(CRISPRcleanR)

  ## Loading library Annotation
  data('KY_Library_v1.0')
  
  ## Deriving the path of the file with the example dataset,
  ## from the mutagenesis of the HT-29 colorectal cancer cell line
  fn<-paste(system.file('extdata', package = 'CRISPRcleanR'),'/HT-29_counts.tsv',sep='')
  
  ## Loading, median-normalizing and computing fold-changes for the example dataset
  normANDfcs<-ccr.NormfoldChanges(fn,min_reads=30,EXPname='ExampleScreen',
                                libraryAnnotation = KY_Library_v1.0,
                                display = FALSE)
  
  ExampleScreen<-normANDfcs$logFCs

  ## Evaluating screen similarity of HT29 reference and user defined data
  ## both, using Project Score criteria
  
  HT29R.expSimilarity(refDataDir='tmpDir', 
                      resDir='tmpDir',
                      userFCs=ExampleScreen,
                      geneGuides='All',
                      geneLevel=TRUE,
                      Rscore=TRUE,
                      saveToFig=TRUE,
                      display=FALSE)

  ## Checking results
  system2('open', args = 'tmpDir/All_QC_SIMILARITY_GENElevel.pdf', wait = FALSE)
  system2('open', args = 'tmpDir/All_SCATTERPLOT_R_MATRIX_GENElevel.pdf', wait = FALSE)
  
  ## Removing Example dataset processed files
  file.remove('ExampleScreen_foldChanges.Rdata')
  file.remove('ExampleScreen_normCounts.Rdata')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
\keyword{benchmarking}% __ONLY ONE__ keyword per line
