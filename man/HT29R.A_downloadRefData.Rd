\name{HT29R.downloadRefData}
\alias{HT29R.downloadRefData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Download reference HT-29 screens data.
}
\description{
  This function allows downloading reference datasets from high-quality CRISPR-Cas9 pooled screens of the HT-29 cell line performed with the employing a commercially available genome-wide library of single guide RNAs (the Kosuke Yusa library) (AddGene: 67989) [1] and setting described in [2]. These data are also publicly aviable on the Project Score website (\url{https://score.depmap.sanger.ac.uk/downloads}) as part of the Cancer Dependency Map portfolio of tools and resources at the Wellcome Sanger Institute (\url{https://depmap.sanger.ac.uk/}).
}
\usage{
H29R.downloadRefData(
      whatToDownload = "FCs",
      destFolder = "./",
      dataRepoURL = "https://cog.sanger.ac.uk/cmp/downloads/crispr_cas9_benchmark/",
      expNames=c("HT29_c903","HT29_c904","HT29_c905",
                 "HT29_c906","HT29_c907","HT29_c908"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{whatToDownload}{
      String specifying which data to dowload.
      Possible values are:
      \itemize{
      \item "rawCounts" to download .tsv files containing raw sgRNA counts; 
      \item "FCs" (default) to download .RData objects containing data.frames with normalised sgRNA depletion log fold-changes. The first two columns contain sgRNAs' identifiers and HGNC symbols of targeted gene, followed by one column per screen replicate with sgRNAs' log fold-changes;}
}
  \item{destFolder}{
      String specifying the path where the dataset will be saved;
}
  \item{dataRepoURL}{
      The URL of the data repository;
}
  \item{expNames}{
      A vector specifying experiment names for the downloaded data.
}
}
\references{
[1] Tzelepis K, Koike-Yusa H, De Braekeleer E, Li Y, Metzakopian E, Dovey OM, et al. A CRISPR Dropout Screen Identifies Genetic Vulnerabilities and Therapeutic Targets in Acute Myeloid Leukemia. Cell Rep. 2016;17:1193–205.

[2] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
Ichcha Manipur & Francesco Iorio (\email{fi1@sanger.ac.uk})
}
\examples{
  #### creating a temporary directory
  dir.create('tempDir')

  #### downloading reference sgRNA depletion fold-changes from high-quality
  #### HT-29 screens into the temporary directory
  HT29R.downloadRefDataset(destFolder = 'tempDir')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
\keyword{data management}% __ONLY ONE__ keyword per line
