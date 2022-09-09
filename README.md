# HT29benchmark
R package for benchmarking genome-wide CRISPR-Cas9 knock-out viability screening pipelines making use of a reference dataset from six independent high-quality screens of the HT-29 cell line (30 replicates, overall).

## Description
The R package accompanying the HT-29 CRISPR screens allows assessing the quality of reference and user-provided CRISPR data. This R package implements several CRISPRcleanR’s utilities wrapped in novel ad-hoc functions, providing a powerful and easy-to-use tool to compute evaluation metrics and generate corresponding plots for visual inspection.  

Specifically, it allows to:

* Download the HT-29 reference datasets. 

* Visualise sgRNAs depletion logFCs distributions of each screen.  

* Evaluate intra-screen reproducibility comparing depletion logFC (log Fold Change) profiles both at the sgRNA and gene-level. 

* Evaluate inter-screen similarity both at the sgRNA and gene-level. 

* Evaluate screens’ performance in recovering gene essentiality profiles of known essential (True Positives, TP) and known non-essential (True Negatives, TN) genes. 

* Visualise depletion logFCs distributions for TP and TN gene sets (as well as their targeting sgRNAs) and provide a Glass’s Δ score for the difference between TP and TN classes. 

* Derive HT-29-specific essential/non-essential genes, by analysing all screens in the reference dataset jointly, then use these sets as positive/negative controls to estimate to what extent a user-provided screen meets expectations, using the metric listed above.

Detailed examples of package usage can be found in the Vignette, which reproduces the full evaluation pipeline for both reference HT-29 screens and user-provided data.