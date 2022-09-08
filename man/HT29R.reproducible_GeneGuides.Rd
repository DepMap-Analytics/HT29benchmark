\name{HT29R.reproducible_GeneGuides}
\alias{HT29R.reproducible_GeneGuides}
\docType{data}
\title{
Library specific informative/reproducible sgRNAs.
}
\description{
838 KY-library[1] specific informative/reproducible sgRNAs (targeting 308 genes) for evaluating CRISPR-Cas9 pooled genome-wide viability screen replicates.
}
\usage{data(HT29R.reproducible_GeneGuides)}
\format{
  A vector of strings with entries corresponding to sgRNAs' identifiers.
}
\details{
Genome-wide correlation scores computed between replicates of the same CRISPR-Cas9 pooled genome-wide viability screen are generally always very high and indistiguishable from expectation due to only a small percentage of genes exerting an effect on cellular fitness upon knock-out. In [2] we have selected a set of 838 most informative sgRNAs, defined as those targeting the same genes and with an average pairwise Pearson’s correlation > 0.6 between corresponding patterns of depletion fold-changes (FCs) across hundreds of screened cell lines.
Per construction, the depletion patterns of these sgRNAs are both reproducible and informative (as they involve genes carrying an actual fitness signal).
Computing correlation scores between replicates of the same screen on the domain of these sgRNAs only allowed the estimation of a null distribution of replicate correlations and computing a reproducibility threshold defined as the minimal correlation score that should be observed between replicates of the same screen (R = 0.68).
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
[1] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, et al. A CRISPR Dropout Screen Identifies Genetic Vulnerabilities and Therapeutic Targets in Acute Myeloid Leukemia. Cell Rep. 2016;17:1193–205.

[2] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\examples{
data(HT29R.reproducible_GeneGuides)
head(HT29R.reproducible_GeneGuides)
}
\keyword{datasets}
