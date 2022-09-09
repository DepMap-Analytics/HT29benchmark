\name{HT29R.prSCORE_bkgr_screen_similarity_sgRNA}
\alias{HT29R.prSCORE_bkgr_screen_similarity_sgRNA}
\docType{data}
\title{
  Pairwise screen replicate correlations (i.e., background correlations) from Project Score.
}
\description{
  Pearson's correlation scores obtained by comparing genes' essentiality profiles (i.e., sgRNAs' depletion log fold-changes) between all possible pairs of samples (325 x 325) from the Project Score [1] with Kosuke Yusa library [2] sgRNAs x 325 samples (= 324 cell lines + HT29v1.1).
}
\usage{data("HT29R.prSCORE_bkgr_screen_similarity_sgRNA")}
\format{
  A numerical vector of Pearson's correlation coefficents.
}
\details{
  A numerical vector with 52,650 entries, containing the correlation scores obtained by comparing profiles of 86,875 sgRNAs' depletion log fold-changes between all possible pairs of individual experiments (i.e., the "All" guide-level background correlation). See function \code{HT29R.expSimilarity} description in the reference manual for further details.
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
[2] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, et al. A CRISPR Dropout Screen Identifies Genetic Vulnerabilities and Therapeutic Targets in Acute Myeloid Leukemia. Cell Rep. 2016;17:1193–205.
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

