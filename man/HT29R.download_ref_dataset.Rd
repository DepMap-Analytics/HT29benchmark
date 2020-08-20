\name{HT29R.download_ref_dataset}
\alias{HT29R.download_ref_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Download reference HT-29 screens data
}
\description{
  This function allows downloading reference datasets from high-quality CRISPR-Cas9 pooled screens of the HT-29 cell line performed with the KY sgRNA library[1]. This data has been generated through the experimental pipeline described in [2] and it is also public aviable on the Project Score web-site (https://score.depmap.sanger.ac.uk/downloads), part of the Cancer Dependency Map portfolio of tools and resources at the Wellcome Sanger Institute (https://depmap.sanger.ac.uk/).
}
\usage{
H29R.download_ref_dataset(
      whatToDownload = "FCs",
      destFolder = "./",
      dataRepoURL = "https://cog.sanger.ac.uk/cmp/downloads/crispr_cas9_benchmark/",
      expNames=c("HT29_c903","HT29_c904","HT29_c905",
                 "HT29_c906","HT29_c907","HT29_c908"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{whatToDownload}{
      String parameter specifying what type of data to dowload.
      Possible values are "rawCounts" for plain .tsv files containing raw sgRNA countsare or "FCs" (default) for R objects containing sgRNA normalised depletion fold-changes: data frames in which the first two columns contain sgRNAs' identifiers and HGNC symbols of targete gene, followed by one column per screen replicate containing sgRNAs' fold-changes;
}
  \item{destFolder}{
      String specifying where the dataset should be saved;
}
  \item{dataRepoURL}{
      The URL of the data repository;
}
  \item{expNames}{
      A vector of strings specifying the experiment names for the dataset to download.
}
}
\references{
[1] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, et al. A CRISPR Dropout Screen Identifies Genetic Vulnerabilities and Therapeutic Targets in Acute Myeloid Leukemia. Cell Rep. 2016;17:1193–205.

[2] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
Ichcha Manipur & Francesco Iorio (fi1@sanger.ac.uk)
}
\examples{
  #### creating a temporary directory
  dir.create('tempDir')

  #### downloading reference sgRNA depletion fold-changes from high-quality
  #### HT-29 screens into the temporary directory
  HT29R.download_ref_dataset(destFolder = 'tempDir')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
\keyword{data management}% __ONLY ONE__ keyword per line
