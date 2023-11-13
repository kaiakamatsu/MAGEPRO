# MAGEPRO (Multi-Ancestry Gene Expression Prediction Regularized Optimization) <img src = IMAGES/veigar.png width="50" height="50"> <img src = IMAGES/wand.png width="50" height="50">
MAGEPRO is a predictive genetics method to create powerful gene expression prediction models. We leverage eQTL summary statistics from diverse ancestries to optimize our gene models and enable the identification of novel gene-trait associations through transcriptome-wide association studies. 

# USING MAGEPRO 

## Dependencies 

R version 4.2.3
Packages:
- data.table_1.14.8
- glmnet_4.1-7
- optparse_1.7.3
- plink2R_1.1
  > install.packages("devtools")
  
  > library(devtools)
  
  > devtools::install_github("carbocation/plink2R/plink2R", ref="carbocation-permit-r361")

GCTA 
Download gcta from either ...
- https://github.com/gusevlab/fusion_twas.git
- https://yanglab.westlake.edu.cn/software/gcta/#Overview

Plink 1.9 
Download from ...
- https://www.cog-genomics.org/plink/1.9/

## Preparing datasets

MAGEPRO utilizes external eQTL summary statistics to improve gene models trained on a specific population. See PROCESS_DATASET directory for more information on how to prepare summary statistics used in our analysis.

## Command-line options for the full pipeline (RUN_MAGEPRO_PIPELINE.R)

**required flags are bolded**

| Flags | Description |
| -------- | -------- |
| **--bfile** | Path to PLINK binary input file prefix (this pipeline will only use sample IDs with both GE and GENOTYPE data) |
| **--out** | Path to save output files |
| **--scratch** | Path to scratch directory (temporary files will be stored here) |
| **--ge** | Path to individual-level, normalized gene expression data in matrix format (see input formats below for more details) |
| --covar | Optional path to quantitative covariates (PLINK format) |
| --num_covar | Number of covariate to use (number of rows to extract from --covar file). "ALL" to use all covariates available. Default assumes gtex covariate file (first 5 genotype PC, first 15 gene expression inferredcov, pcr, platform, sex: ideal for N < 150, see GTEx pipeline for more information). |
| --num_batches | Number of batch jobs to split the genes into. Default 20. |
| --rerun | Boolean to indicate if the pipeline is being reran (TRUE = skip creation of plink files per gene) |
| --intermed_dir | Directory to store intermediate files |
| --subset_genes | Path to file with genes of interest in one column |
| --sumstats_dir | Path to external datasets (required if using MAGEPRO or META models) |
| --sumstats | Comma-separated list of external datasets to include (required if using MAGEPRO or META models) |
| --model | Comma-separated list of models to use. Options: "SINGLE" "META" and "MAGEPRO" (default SINGLE,META,MAGEPRO) |
| --ss | Comma-separated list of sample sizes of sumstats in the same order as --sumstats (required if using "META" model or --cell_type_meta) |
| --cell_meta  | Comma-separated list of prefixes of eqtl datasets to cell type meta-analyze (--ss required) |
| --PATH_plink | Path to plink executable (default "plink") |
| --PATH-gcta | Path to gcta executable (default "gcta_nr_robust") |
| --resid | Also regress the covariates out of the genotypes (default FALSE) |
| --hsq_p | Minimum heritability p-value for which to compute weights (default 0.01, significantly heritable) |
| --lassohsq | Backup heritability value to use in lasso regression if heritability is not calculated with gcta or gcta fails to produce a reasonable estimate (h2 > 0) (default 0.05) |
| --hsq_set | Skip heritability estimation and set hsq estimate to this value (optional, heritability computed with gcta otherwise) |
| --crossval | Number of cross-validations (0 or 1 to skip) |
| --verbose | How much chatter to print: 0=nothing; 1=minimal; 2=all (default 1) |
| --noclean | Do not delete any temporary files (default FLASE) |
| --save_hsq | Save heritability results even if weights are not computed (default FALSE) |

## Output format
MAGEPRO saves the output as a RData file which can be loaded into another R script with...
> load("<file path>")

The following variables are saved in this output file: 

| name | description | 
| ---- | ----- | 
| wgt.matrix | final gene model weights, one column for each model type used (SINGLE, META, MAGEPRO) | 
| snps | information about snps used in the gene model, in plink .bim format | 
| cv.performance | matrix with a column for each model type used and a row for each of r-squared and p-value. these are the results of cross-validation | 
| hsq | gcta-estimated heritability. set to NA if gcta fails to converge | 
| hsq.pv | gcta-estimated heritability p-value. set to NA if gcta fails to converge | 
| N.tot | sample size | 
| wgtmagepro | datasets used in the MAGEPRO model | 
| cf_total | mixing weights (alphas from regression) used to combine the datasets from wgtmagepro | 
| avg_training_r2_single | average r-squared of SINGLE model on training cohort | 
| avg_training_r2_meta | average r-squared of META model on training cohort | 
| avg_training_r2_magepro | average r-squared of MAGEPRO model on training cohort | 

> see "PROCESS_RESULTS" directory on how to format the results from cv.performance across all genes into a tab-delimited dataframe. 

## Input data format
> the MAGEPRO pipeline is designed to handle files that are formatted like the gene expression and covariate files made available by GTEx.

- GTEx portal: https://www.gtexportal.org/home/downloads/adult-gtex#qtl
- gene expression data from GTEx: https://storage.cloud.google.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_expression_matrices.tar
- covariates file from GTEx: https://storage.cloud.google.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_covariates.tar.gz 

### Gene Expression data format 
- Normalized gene expression measurements, inverse-normal-transformed across samples
- Tab-delimited matrix (gzipped) with the following columns:

| #chr | start | end | gene_id | SAMPLE1 | (the rest of the columns are sample IDs) | 
| ---- | ----- | ---- | ----- | ---- | ---- |
| 1 | 3465 | 4650 | ENSGXXXX | 0.20120726288521493 | (gene expression level of this gene for the other samples) |

Example: 

<img src = IMAGES/GTEx_GE.png>

### Covariates data format 
- Tab-delimited matrix (.txt) with the first row containing sample IDs and all other rows for covariates

| ID | SAMPLE1 | SAMPLE2 | SAMPLE3 | (the rest of the columns are sample IDs) | 
| ---- | ----- | ---- | ----- | ---- |
| PC1 | XXXX | XXXX | XXXX | (PC1 for rest of samples) | 
| PC2 | XXXX | XXXX | XXXX | (PC2 for rest of samples) | 

Example: 

<img src = IMAGES/GTEx_cov.png>

### Genotype data format
- plink bed/bim/fam files
- described here: https://www.cog-genomics.org/plink/1.9/formats

## Typical application of MAGEPRO on GTEx Tissues

## Command-line options for the computing gene-models one gene at a time

If you would like to run MAGEPRO on one gene, it is possible to run MAGEPRO.R separately with the following flags

**required flags are bolded**

| Flags | Description |
| -------- | -------- |
| **--gene** | ENSG ID of gene |
| **--bfile** | Path to PLINK binary input file prefix |
| **--out** | Path to save output files |
| **--tmp** | Path to store temporary files |
| --sumstats_dir | Path to external datasets (required if using MAGEPRO or META models) |
| --sumstats | Comma-separated list of external datasets to include (required if using MAGEPRO or META models) |
| --model | Comma-separated list of models to use. Options: "SINGLE" "META" and "MAGEPRO" (default SINGLE,META,MAGEPRO) |
| --ss | Comma-separated list of sample sizes of sumstats in the same order as --sumstats (required if using "META" model or --cell_type_meta) |
| --cell_meta  | Comma-separated list of prefixes of eqtl datasets to cell type meta-analyze (--ss required) |
| --pheno | Path to molecular phenotype file in PLINK format (taken from bfile otherwise) |
| --PATH_plink | Path to plink executable (default "plink") |
| --PATH-gcta | Path to gcta executable (default "gcta_nr_robust") |
| --covar | Path to quantitative covariates in PLINK format (optional) |
| --resid | Also regress the covariates out of the genotypes (default FALSE) |
| --hsq_p | Minimum heritability p-value for which to compute weights (default 0.01, significantly heritable) |
| --lassohsq | Backup heritability value to use in lasso regression if heritability is not calculated with gcta or gcta fails to produce a reasonable estimate (h2 > 0) |
| --hsq_set | Skip heritability estimation and set hsq estimate to this value (optional, heritability computed with gcta otherwise) |
| --crossval | Number of cross-validations (0 or 1 to skip) |
| --verbose | How much chatter to print: 0=nothing; 1=minimal; 2=all (default 1) |
| --noclean | Do not delete any temporary files (default FLASE) |
| --save_hsq | Save heritability results even if weights are not computed (default FALSE) |

## Directories 

| Directory | Description |
| -------- | -------- |
| BENCHMARK | Scripts to benchmark MAGEPRO against other models |
| GALA_SAGE | Scripts used to apply MAGEPRO to GALA II and SAGE data |
| GTEx_v8_pipeline | Example of applying MAGEPRO to GTEx v8 data |
| IMAGES | Images used in github |
| MAGEPRO_PIPELINE | Intermediate scripts that are a part of the main MAGEPRO pipeline |
| PLOTS | R scripts used to visualize data |
| PROCESS_DATASET | Scripts used to extract and process eQTL summary statistics from various studies |
| PROCESS_RESULTS | Scripts to process and summarize gene model performance |
| SIMULATIONS | Scripts for testing MAGEPRO in simulation
| SUPPLEMENTAL_ANALYSIS | All other supplementary analysis |
| TWAS | Transcriptome-wide association studies using MAGEPRO gene models |
| VALIDATE_MODELS | Validating MAGEPRO gene models | 
| test_1gene | Scripts for testing |

# Support 
> contact kakamatsu@ucsd.edu
