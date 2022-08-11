\name{HT29R.prSCORE_bkgr_screen_similarity_sgRNA}
\alias{HT29R.prSCORE_bkgr_screen_similarity_sgRNA}
\docType{data}
\title{
  Pair-wise screen replicate correlations (i.e., background correlations) from Project Score.
}
\description{
  Correlation scores obtained by comparing genes' essentiality profiles (i.e., sgRNAs' depletion log fold-changes) between all possible pairs of sample (325 x 325), considering the Project Score[1] data matrix for 86875 KY-library[1] specific sgRNAs x 325 samples (= 324 cell lines + HT29v1.1).
}
\usage{data("HT29R.prSCORE_bkgr_screen_similarity_sgRNA")}
\format{
  A numerical vector of Pearson's correlation coefficents.
}
\details{
  A numerical vector with 52650 entries, containing the correlation scores obtained by comparing profiles of 86875 sgRNAs' depletion fold-changes between all possible pairs of individual experiment (i.e., the "All" guide-level background correlation). See function \code{HT29R.expSimilarity} description in the reference manual for further details.
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\examples{

library(HT29benchmark)
data(HT29R.prSCORE_bkgr_screen_similarity_sgRNA)

plot(density(HT29R.prSCORE_bkgr_screen_similarity_sgRNA), 
  xlim=c(0,1), 
  main="Project Score background distribution\n (86875 sgRNAs)",
  col="grey",
  lwd=2,
  xlab="R")

}

