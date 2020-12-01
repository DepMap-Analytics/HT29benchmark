### creating a temporary directory
dir.create('../../tmpFolder/')

### downloading reference sgRNA depletion fold-changes from high-quality
### HT-29 screens into the temporary directory
HT29R.download_ref_dataset(destFolder = '../../tmpFolder/')

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

## adding some noise to the data
length <- dim(ExampleScreen)[1] * 3
noise <- matrix(runif(length, -0.5, 0.5), dim(ExampleScreen)[1])
ExampleScreen[,3:ncol(ExampleScreen)] <- ExampleScreen[,3:ncol(ExampleScreen)] + noise

## Evaluating screen reproducibility of HT29 reference and user defined data
## both, using Project Score criteria, at the sgRNA level
HT29R.evaluate_reps(refDataDir = '../../tmpFolder/',resDir = '../../resFolder/',userFCs = ExampleScreen,geneLevel = FALSE)
## Checking results
system2('open', args = '../../resFolder/RepCor_Vs_PrScore.pdf', wait = FALSE)

## Evaluating screen reproducibility of HT29 reference and user defined data
## both, using Project Score criteria, at the Gene level
HT29R.evaluate_reps(refDataDir = '../../tmpFolder/',resDir = '../../resFolder/',userFCs = ExampleScreen)

## Checking results
system2('open', args = '../../resFolder/GL_RepCor_Vs_PrScore.pdf', wait = FALSE)


## Evaluating similarity of the HT29 reference screens and user defined data
HT29R.exp_similarity(refDataDir = '../../tmpFolder/',resDir = '../../resFolder/',userFCs = ExampleScreen)

## Checking results
system2('open', args = '../../resFolder/Screen_sim.pdf', wait = FALSE)
system2('open', args = '../../resFolder/Sreen_pair_cor.pdf', wait = FALSE)

## Evaluating phenotype intensity
HT29R.PhenoIntensity(refDataDir = '../../tmpFolder/',resDir = '../../resFolder/',userFCs = ExampleScreen)

## Checking results
system2('open', args = '../../resFolder/allScreens_PhenoIntensity.pdf', wait = FALSE)
system2('open', args = '../../resFolder/PhenoIntensity_reference.pdf', wait = FALSE)

## Detection of depleted sgRNAs targeting essential genes

## Loading reference sets of essential/non-essential genes from the CRISPRcleanR package
data(BAGEL_essential)
data(BAGEL_nonEssential)

## Loading KY library annotation from the CRISPRcleanR package
data(KY_Library_v1.0)

## deriving sgRNAs targeting essential/non-essential genes
BAGEL_essential_sgRNAs<-ccr.genes2sgRNAs(libraryAnnotation = KY_Library_v1.0,genes = BAGEL_essential)
BAGEL_nonEssential_sgRNAs<-ccr.genes2sgRNAs(libraryAnnotation = KY_Library_v1.0,genes = BAGEL_nonEssential)



## Removing Example dataset processed files
file.remove('ExampleScreen_foldChanges.Rdata')
file.remove('ExampleScreen_normCounts.Rdata')



