# MAGEPRO (Multi-Ancestry Gene Expression Prediction Regularized Optimization) <img src = IMAGES/veigar.png width="50" height="50"> <img src = IMAGES/wand.png width="50" height="50">
MAGEPRO is a predictive genetics method to create powerful gene expression prediction models. We leverage eQTL summary statistics from diverse ancestries to optimize our gene models and enable the identification of novel gene-trait associations through transcriptome-wide association studies. 

# USING MAGEPRO 

## Dependencies 

R version 4.2.3
Packages:
- data.table_1.14.8
- glmnet_4.1-7
- plink2R_1.1
- optparse_1.7.3

GCTA 
Download gcta from either ...
- https://github.com/gusevlab/fusion_twas.git
- https://yanglab.westlake.edu.cn/software/gcta/#Overview

Plink 1.9 
Download from ...
- https://www.cog-genomics.org/plink/1.9/

## Command-line options for the full pipeline (RUN_MAGEPRO_PIPELINE.R)

**required flags are bolded**

| Flags | Description |
| -------- | -------- |
| **--bfile** | Path to PLINK binary input file prefix (this pipeline will only use sample IDs with both GE and GENOTYPE data) |
| **--out** | Path to save output files |
| **--scratch** | Path to scratch directory (temporary files will be stored here) |
| **--ge** | Path to individual-level, normalized gene expression data in matrix format (see input formats below for more details) |
| --covar | Optional path to quantitative covariates (PLINK format) |
| --num_covar | Number of covariate to use (number of rows to extract from --covar file). Default ALL rows below sample ids. |
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
| --lassohsq | Backup heritability value to use in lasso regression if heritability is not calculated with gcta or gcta fails to produce a reasonable estimate (h2 > 0) |
| --hsq_set | Skip heritability estimation and set hsq estimate to this value (optional, heritability computed with gcta otherwise) |
| --crossval | Number of cross-validations (0 or 1 to skip) |
| --verbose | How much chatter to print: 0=nothing; 1=minimal; 2=all (default 1) |
| --noclean | Do not delete any temporary files (default FLASE) |
| --save_hsq | Save heritability results even if weights are not computed (default FALSE) |

## Output format


## Input data format


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
