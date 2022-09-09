\name{HT29R.consensus_GeneGuides}
\alias{HT29R.consensus_GeneGuides}
\docType{data}
\title{
Library specific sgRNAs' consensus across HT-29 screens.
}
\description{
  86,875 sgRNAs in the Kosuke Yusa library (AddGene: 67989) [1] employedfor evaluating CRISPR-Cas9 pooled genome-wide viability-reduction screen replicates.
}
\usage{data(HT29R.consensus_GeneGuides)}
\format{
  A vector of sgRNAs' identifiers.
}
\details{
  This file contains a list of sgRNAs strictly shared by all HT-29 genome-wide CRISPR-Cas9 screens. 
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\author{
Raffaele Iannuzzi, Ichcha Manipur, Francesco Iorio (\email{fi1@sanger.ac.uk})
}
\references{
[1] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, et al. A CRISPR Dropout Screen Identifies Genetic Vulnerabilities and Therapeutic Targets in Acute Myeloid Leukemia. Cell Rep. 2016;17:1193–205.
}
\examples{
data(HT29R.consensus_GeneGuides)
head(HT29R.consensus_GeneGuides)
}
\keyword{datasets}
