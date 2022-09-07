\name{HT29R.consensus_GeneGuides}
\alias{HT29R.consensus_GeneGuides}
\docType{data}
\title{
Library specific sgRNAs' consensus across HT-29 screens.
}
\description{
  86875 KY-library[1] specific sgRNAs for evaluating CRISPR-Cas9 pooled genome-wide viability screen replicates.
}
\usage{data(HT29R.consensus_GeneGuides)}
\format{
  A vector of strings with entries corresponding to sgRNAs' identifiers.
}
\details{
  This file contain a list of sgRNAs strictly shared by each HT29 genome-wide CRISPR-Cas9 depletion screen and it's used by the \code{HT29R.expSimilarity} function to evaluate the User-provided screen (i.e., the intersection between consensus and User's sgRNAs is considered for comparing screens).
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
[1] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, et al. A CRISPR Dropout Screen Identifies Genetic Vulnerabilities and Therapeutic Targets in Acute Myeloid Leukemia. Cell Rep. 2016;17:1193–205.

[2] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\examples{
data(HT29R.consensus_GeneGuides)
head(HT29R.consensus_GeneGuides)
}
\keyword{datasets}
