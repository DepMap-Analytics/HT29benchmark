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
      String specifying whether all sgRNAs (genome-wide, geneGuides="All") or just the most reproducible and informative sgRNAs (n=838, geneGuides="HI") will be used
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
This function computes inter-screen Pearson's correlation scores for the HT-29 reference screen (averaging replicates), as well as for the user-provided screen. This is performed as for the Project Score data [1].

Particularly, it generates a plot of Pearson’s correlation scores computed between sgRNA depletion log fold-change profiles between pairs of the six HT-29 reference screens, and between reference screens and user-provided data (see the Vignette for details).
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
