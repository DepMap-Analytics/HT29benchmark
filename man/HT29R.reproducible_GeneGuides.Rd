\name{HT29R.reproducible_GeneGuides}
\alias{HT29R.reproducible_GeneGuides}
\docType{data}
\title{
Library specific informative/reproducible sgRNAs.
}
\description{
838 KY-library specific [1] informative/reproducible sgRNAs (targeting 308 genes) for evaluating CRISPR-Cas9 pooled genome-wide viability screen replicates.
Genome-wide correlation scores computed between replicates of the same CRISPR-Cas9 pooled genome-wide viability screen are generally always very high and indistiguishable from expectation due to only a small percentage of genes exerting an effect on cellular fitness upon knock-out. In [2] the authors selected this set of 838 most informative sgRNAs, defined as those targeting the same genes and with an average pairwise Pearson’s correlation > 0.6 between corresponding patterns of depletion fold-changes (FCs) across hundreds of screened cell lines.
Per construction, these sgRNAs are both reproducible and informative (as they involve genes carrying an actual fitness signal).
Computing correlation scores between replicates of the same screen on the domain of these sgRNAs only for 308 genes targeted by these informative sgRNAs for each individual technical replicate, and then computed all possible pairwise Pearson’s correlation scores between the resulting profiles. This allowed the estimation of a null distribution of replicate correlations (plotted in gray in Extended Data Fig. 1d).
}
\usage{data(HT29R.reproducible_GeneGuides)}
\format{
  A vector of strings with entries corresponding to sgRNA identifiers.
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(HT29R.reproducible_GeneGuides)
## maybe str(HT29R.reproducible_GeneGuides) ; plot(HT29R.reproducible_GeneGuides) ...
}
\keyword{datasets}
