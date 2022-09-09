\name{HT29R.prSCORE_bkgr_screen_similarity_sgRNA_HI}
\alias{HT29R.prSCORE_bkgr_screen_similarity_sgRNA_HI}
\docType{data}
\title{
  Pairwise correlations between screens replicates (i.e., background correlations) from Project Score, using the most informative sgRNAs.
}
\description{
  Pearson's correlation scores obtained by comparing genes' essentiality profiles (i.e., genes' depletion log fold-changes) between all possible pairs of sample (325 x 325) from the Project Score [1] with the Kosuke Yusa library [2] informative/reproducible sgRNAs x 325 samples (= 324 cell lines + HT29v1.1).
}
\usage{data("HT29R.prSCORE_bkgr_screen_similarity_sgRNA_HI")}
\format{
  A numerical vector of Pearson's correlation coefficents.
}
\details{
  A numerical vector with 52,650 entries, containing the correlation scores obtained by comparing profiles of 838 sgRNAs' depletion log fold-changes between all possible pairs of individual experiments (i.e., the "HI" sgRNA-level background correlation). See function \code{HT29R.expSimilarity} and \code{HT29R.reproducible_GeneGuides} descriptions in the reference manual for further details.
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
[2] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, et al. A CRISPR Dropout Screen Identifies Genetic Vulnerabilities and Therapeutic Targets in Acute Myeloid Leukemia. Cell Rep. 2016;17:1193–205.
}
\examples{

library(HT29benchmark)
data(HT29R.prSCORE_bkgr_screen_similarity_sgRNA_HI)

plot(density(HT29R.prSCORE_bkgr_screen_similarity_sgRNA_HI), 
  xlim=c(-0.2,1), 
  main="Project Score background distribution\n (838 most-informative sgRNAs)",
  col="grey",
  lwd=2,
  xlab="R")
}
