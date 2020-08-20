\name{HT29R.prSCORE_rCorr}
\alias{HT29R.prSCORE_rCorr}
\docType{data}
\title{
  Pair-wise screen replicate correlations and background correlations from Project Score.
}
\description{
  Correlation scores obtained by comparing sgRNA essentiality profiles, considering KY-Library[1] specific informative/reproducible sgRNAs depletion fold-changes between replicates of the same experiment, and between all possible pairs of individual replicates across experiments, from Project Score[2]].
}
\usage{data("HT29R.prSCORE_rCorr")}
\format{
  A list of two numerical vectors:
  \describe{
    \item{\code{BGscores}}{a numeric vector with 441387 entries, containing the correlation scores obtained by comparing profiles of KY-Library[1] specific informative/reproducible sgRNAs depletion fold-changes between all possible pairs of individual replicates across experiments, i.e. background correlation.}
    \item{\code{REPscores}}{a numeric vector with 883 entries, containing the correlation scores obtained by comparing profiles of KY-Library[1] specific informative/reproducible sgRNAs depletion fold-changes between replicates of the same experiment.}
  }
}
\references{
[1] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, et al. A CRISPR Dropout Screen Identifies Genetic Vulnerabilities and Therapeutic Targets in Acute Myeloid Leukemia. Cell Rep. 2016;17:1193–205.

[2] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\examples{

library(CRISPRcleanR)
data(HT29R.prSCORE_rCorr)

ccr.multDensPlot(
  list(density(HT29R.prSCORE_rCorr$BGscores),
       density(HT29R.prSCORE_rCorr$REPscores)),
  XLIMS = c(0,1),
  TITLE = 'Project Score Observed vs. Expected replicate correlations\n (informative sgRNAs only)',
  COLS = c('gray','darkgreen'),LEGentries = c('expected','observed'),XLAB='R')
}
\keyword{datasets}
