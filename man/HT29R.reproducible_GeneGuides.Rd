\name{HT29R.reproducible_GeneGuides}
\alias{HT29R.reproducible_GeneGuides}
\docType{data}
\title{
Library-specific informative/reproducible sgRNAs.
}
\description{
838 informative/reproducible sgRNAs (targeting 308 genes) from the Kosuke Yusa library [1] for evaluating CRISPR-Cas9 pooled genome-wide viability screen replicates.
}
\usage{data(HT29R.reproducible_GeneGuides)}
\format{
  A vector of strings with entries corresponding to sgRNAs' identifiers.
}

\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
[1] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, et al. A CRISPR Dropout Screen Identifies Genetic Vulnerabilities and Therapeutic Targets in Acute Myeloid Leukemia. Cell Rep. 2016;17:1193–205.
}
\examples{
data(HT29R.reproducible_GeneGuides)
head(HT29R.reproducible_GeneGuides)
}
\keyword{datasets}
