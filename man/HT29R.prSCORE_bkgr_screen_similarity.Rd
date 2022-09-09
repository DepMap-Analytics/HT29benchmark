\name{HT29R.prSCORE_bkgr_screen_similarity}
\alias{HT29R.prSCORE_bkgr_screen_similarity}
\docType{data}
\title{
  Pairwise screen replicate correlations (i.e., background correlations) from Project Score.
}
\description{
  Pearson's correlation scores obtained by comparing genes' essentiality profiles (i.e., genes' depletion log fold-changes) between all possible pairs of sample (325 x 325), considering the Project Score [1] data matrix for 17,995 genes x 325 samples (= 324 cell lines + HT29v1.1).
}
\usage{data("HT29R.prSCORE_bkgr_screen_similarity")}
\format{
  A numerical vector of Pearson's correlation coefficents.
}
\details{
  A numerical vector with 52,650 entries, containing the correlation scores obtained by comparing profiles of 17,995 genes' depletion log fold-changes between all possible pairs of individual experiments (i.e., the "All" gene-level background correlation). See function \code{HT29R.expSimilarity} description in the reference manual for further details. 
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
Ichcha Manipur & Francesco Iorio (\email{fi1@sanger.ac.uk})
}
\examples{

library(HT29benchmark)
data(HT29R.prSCORE_bkgr_screen_similarity)

plot(density(HT29R.prSCORE_bkgr_screen_similarity), 
  xlim=c(0,1), 
  main="Project Score background distribution\n (17995 targeted genes)",
  col="grey",
  lwd=2,
  xlab="R")
}


