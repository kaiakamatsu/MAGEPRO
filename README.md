# MAGEPRO (Multi-Ancestry Gene Expression Prediction Regularized Optimization) <img src = IMAGES/veigar.png width="50" height="50"> <img src = IMAGES/wand.png width="50" height="50">
MAGEPRO is a predictive genetics method to create powerful gene expression prediction models. We leverage eQTL summary statistics from diverse ancestries to optimize our gene models and enable the identification of novel gene-trait associations through transcriptome-wide association studies. 

# FIXES 
- in 1_CollectSamples.R, also subset to individuals in covar file, and create a set of indivdiuasl present in ge, genotype and covars
- --num_covar flag is hard to use for gtex formatted covars with sex and age in the last row
- --ldref_pt by default uses bfile

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

## Extra Dependencies 

Our tool enables users to build predictive models of gene expression using PRS and fine-mapping methods that have not been applied to gene expression prediction before.  
Although this may be a useful option, it require some extra dependencies:

### PRS-CSx 
> https://github.com/getian107/PRScsx  

**Python**
- scipy (https://www.scipy.org/)
- h5py (https://www.h5py.org/)

### Bridge-PRS
> https://www.bridgeprs.net/  

**R**
- BEDMatrix, boot, data.table, doMC, glmnet, MASS, optparse, parallel, and R.utils
> install.packages(c("BEDMatrix","boot","data.table","doMC","glmnet","MASS","optparse","parallel","R.utils"))  

**Python3**
- matplotlib

### SuSiE 
> https://stephenslab.github.io/susieR/index.html

**R**  
> install.packages("susieR")

## Preparing datasets

MAGEPRO utilizes external eQTL summary statistics to improve gene models trained on a specific population. See PROCESS_DATASET directory for more information on how to prepare summary statistics used in our analysis.

## Command-line options for the full pipeline (RUN_MAGEPRO_PIPELINE.R)

**required flags are bolded**

| Flags | Description |
| -------- | -------- |
| **--bfile** | Path to PLINK binary input file prefix (this pipeline will only use sample IDs with both GE and GENOTYPE data) |
| **--out** | Path to save output files |
| **--scratch** | Path to scratch directory (gene-specific files will be stored here) |
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
> load("file path")

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
| var_cov | variance explained in gene expression by covariates | 

> see "PROCESS_RESULTS" directory on how to format the results from cv.performance across all genes into a tab-delimited dataframe. 

## Input data format
> the MAGEPRO pipeline is designed to handle files that are formatted like the gene expression and covariate files made available by GTEx.

- GTEx portal: https://www.gtexportal.org/home/downloads/adult-gtex#qtl
- gene expression data from GTEx: https://storage.cloud.google.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_expression_matrices.tar
- covariates file from GTEx: https://storage.cloud.google.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_covariates.tar.gz 

### Gene Expression data format 
- Normalized gene expression measurements, inverse-normal-transformed across samples
- Tab-delimited matrix (gzipped) with the following columns:
- Sorted. First by increasing chromosome number. Then by increasing start position. 

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

## QUICKSTART: Typical application of MAGEPRO on GTEx Tissues

1. clone this repository
> git clone https://github.com/kaiakamatsu/MAGEPRO.git

2. create directories where intermediate files and eQTL summary statistics will be stored
> mkdir GE_PREDICTION
> 
> cd GE_PREDICTION
> 
> mkdir DATASETS
> 
> mkdir GTEX
> 
> mkdir MODELS
>
> mkdir SCRATCH

3. download and process eQTL summary statistics (**we will only use MESA eQTL summary statistics for this tutorial.** see PROCESS_DATASET directory)
> cd DATASETS
> 
> wget "https://www.dropbox.com/sh/f6un5evevyvvyl9/AAA3sfa1DgqY67tx4q36P341a?dl=1"
>
> mkdir mesahis
>
> mkdir mesaafr
>
> unzip MESA_cis-eQTLs.zip
> 
> pwd #copy this path
>
> cd ../../MAGEPRO/PROCESS_DATASET/mesa
>
> vim run_splitFilter.sh
>
> - paste path copied above to replace "PASTE PATH HERE"
> - replace "PATH TO SNPS HERE" with the path to a bim file containing all SNPs of interest (we used HM3 SNPs in GTEx)
>
> vim submitSplit.sh #edit this file to prepare for job submission on your HPC cluster
>
> sbatch submitSplit.sh
> - 'bash run_splitFilter.sh' #if not on cluster 
>
> cd ../../../GE_PREDICTION/DATASETS/
> 
> ls
> - after this job runs, 'mesaafr' and 'mesahis' directories will contain gene-specific eQTL summary statistics data files

4. download GTEx Gene Expression and Covariates data
> cd ../GTEX
>
> wget "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_expression_matrices.tar" # gene expression
>
> tar -xvf GTEx_Analysis_v8_eQTL_expression_matrices.tar
> 
> wget "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_covariates.tar.gz" # covariates
>
> tar -xzf GTEx_Analysis_v8_eQTL_covariates.tar.gz
>
> - individual-level GTEx genotypes are protected access. this should be in this same directory, plink format (bed/bim/fam) and split by chromosome
> - genotype files should contain individuals from the ancestry of interest (ex. only AFR individuals)

5. run MAGEPRO
> cd ../../MAGEPRO
>
> vim MAGEPRO_PIPELINE/5_RunJobs.sh # edit to prepare for job submission on your HPC cluster
> 
> #example with whole blood below, computing models for only significantly heritable
> 
> Rscript RUN_MAGEPRO_PIPELINE.R --bfile ../GE_PREDICTION/GTEX/"prefix of genotype files, preceeding chromosome number" --ge ../GE_PREDICTION/GTEX/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz --covar ../GE_PREDICTION/GTEX/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt --out ../GE_PREDICTION/MODELS --scratch ../GE_PREDICTION/SCRATCH --intermed_dir ../GE_PREDICTION --PATH_plink "path to plink" --PATH_gcta "path to gcta" --sumstats_dir ../GE_PREDICTION/DATASETS --sumstats mesahis,mesaafr --models SINGLE,META,MAGEPRO --ss 352,233 --verbose 2 --num_covar ALL
>
> cd ../GE_PREDICTION/MODELS #this directory will contain all output files after the job finishes running

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
