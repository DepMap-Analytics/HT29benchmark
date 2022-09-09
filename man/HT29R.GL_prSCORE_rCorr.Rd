\name{HT29R.GL_prSCORE_rCorr}
\alias{HT29R.GL_prSCORE_rCorr}
\docType{data}
\title{
  Gene level Pearson's correlations between screen replicates and background correlations from Project Score.
}
\description{
   Similarity between individual replicates of CRISPR-Cas9 screens from the Project Score [1], performed with the Kosuke Yusa library (AddGene: 67989) [2], using depletion log fold-changes only for a set of informative/reproducible sgRNAs (averaged on a targeted gene basis) [1].
}
\usage{data("HT29R.GL_prSCORE_rCorr")}
\format{
  A list of two numerical vectors:
  \describe{
    \item{\code{BGscores}}{a numeric vector with 441,387 entries, containing the pairwise Pearson's correlation scores between all possible pairs of individual replicates across experiments (i.e. background correlation), employing depletion log fold-changes for a set of informative/reproducible sgRNAs (averaged on a targeted gene basis).}
    \item{\code{REPscores}}{a numeric vector with 883 entries, containing the correlation scores obtained by comparing profiles of KY-Library[1] specific informative/reproducible sgRNAs depletion fold-changes (avereged on a targeted gene basis) between replicates of the same experiment.}
  }
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, et al. A CRISPR Dropout Screen Identifies Genetic Vulnerabilities and Therapeutic Targets in Acute Myeloid Leukemia. Cell Rep. 2016;17:1193–205.
}
\author{
Ichcha Manipur, Raffaele Iannuzzi, Francesco Iorio (\email{fi1@sanger.ac.uk})
}
\examples{

library(CRISPRcleanR)
data(HT29R.GL_prSCORE_rCorr)

ccr.multDensPlot(
  list(density(HT29R.GL_prSCORE_rCorr$BGscores),
       density(HT29R.GL_prSCORE_rCorr$REPscores)),
  XLIMS = c(0,1),
  TITLE = 'Project Score Observed vs. Expected replicate correlations\n (Gene level depletion fc, informative sgRNAs only)',
  COLS = c('gray','darkgreen'),LEGentries = c('expected','observed'),XLAB='R')
}
\keyword{datasets}
