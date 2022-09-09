# HT29benchmark
R package for benchmarking genome-wide CRISPR-Cas9 knock-out viability screening pipelines making use of a reference dataset from six independent high-quality screens of the HT-29 cell line (30 replicates, overall).

## Description
The R package accompanying the HT-29 CRISPR screens allows assessing the quality of reference and user-provided CRISPR data. This R package implements several CRISPRcleanR’s utilities wrapped in novel ad-hoc functions, providing a powerful and easy-to-use tool to:  

* Download the HT-29 reference datasets. 

* Visualize sgRNAs depletion logFCs distributions of each screen.  

* Evaluate screens’ reproducibility both at the sgRNA and gene-level. 

* Evaluate screens’ similarity both at the sgRNA and gene-level. 

* Evaluate screens’ performance in recovering gene essentiality profiles using sets of known essential (True Positives, TP) and known non-essential (True Negatives, TN) genes. 

* Visualize depletion logFCs distributions for these gene sets observed both at the sgRNA and gene-level and provide a Glass’s Δ score for the difference between TP and TN classes. 

* Compute ROC and PrRc curves at the 5% False Discovery Rate (FDR). 

* Derive genome-wide consensus of positive - i.e., significantly depleted genes at 5% FDR - and negative controls in order to derive novel HT-29-specific gene sets and evaluate Cohen’s d or Glass’s Δ distances between the resulting depletion logFCs of the positive and negative consensus of genes for each screen.  
