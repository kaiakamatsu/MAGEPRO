# Benchmarking MAGEPRO against PRS-CSx

## PRS-CSx
https://github.com/getian107/PRScsx.git
- PRS-CSx is a method to integrate GWAS summary statistics from multiple populations to create PRS models for understudied populations
- they take summary statistics from multiple cohorts, apply a shrinkage prior
- integrate with a linear regression, evaluate with 5-fold cross validation
- phi = 1e-4 (for less polygenic traits)
- download LD reference panels and SNP information files as outlined in the PRS-CSx github

## Prune and Threshold Summary Stats
- use 80% of individuals to train summary stats (linear regression per snp)
- use plink --clump using AFR GTEX ldref files
- use these P+T weights to test on 20%


