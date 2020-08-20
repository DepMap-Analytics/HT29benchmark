#### creating a temporary directory
dir.create('tempDir')

#### downloading reference sgRNA depletion fold-changes from high-quality
#### HT-29 screens into the temporary directory
HT29R.download_ref_dataset(destFolder = 'tempDir')

## Loading CRISPRcleanR library to use example screen data
library(CRISPRcleanR)

## Deriving the path of the file with the example dataset,
## from the mutagenesis of the HT-29 colorectal cancer cell line
fn<-paste(system.file('extdata', package = 'CRISPRcleanR'),'/HT-29_counts.tsv',sep='')

## Loading library Annotation
data('KY_Library_v1.0')

## Loading, median-normalizing and computing fold-changes for the example dataset
normANDfcs<-ccr.NormfoldChanges(fn,min_reads=30,EXPname='ExampleScreen',
                                libraryAnnotation = KY_Library_v1.0,
                                display = FALSE)
ExampleScreen<-normANDfcs$logFCs

## Evaluating screen reproducibility of HT29 reference and user defined data
## both, using Project Score criteria
HT29R.replicateCorr_Pscore(refDataDir = 'tempDir',resDir = 'tempDir',userFCs = ExampleScreen)

## Checking results
system2('open', args = 'tempdir/RepCor_Vs_PrScore.pdf', wait = FALSE)

## Removing Example dataset processed files
file.remove('ExampleScreen_foldChanges.Rdata')
file.remove('ExampleScreen_normCounts.Rdata')





### Other figures

load('tempDir/HT29_c903_foldChanges.Rdata')

HT29R.download_ref_dataset(destFolder = 'tempDir',whatToDownload = 'rawCounts')

counts1<-read.table('tempDir/HT29_c903.tsv',stringsAsFactors = FALSE,row.names = 1,header=TRUE)

plot(log10(counts1$HT29_c903R1+0.05),log10(counts1$HT29_c903R2+0.05))


