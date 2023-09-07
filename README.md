# MAGEPRO (Multi-Ancestry Gene Expression Prediction Regularized Optimization) <img src = IMAGES/veigar.png width="50" height="50"> <img src = IMAGES/wand.png width="50" height="50">
MAGEPRO is a predictive genetics method to create powerful gene expression prediction models for populations that are underrepresented in genetic studies. We leverage eQTL summary statistics from diverse ancestries to optimize our gene models and enable the identification of novel gene-trait associations through transcriptome-wide association studies. 

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


## Command-line options for the full pipeline

## Output format

## Input data format

## Typical application of MAGEPRO on GTEx Tissues

## Command-line options for the computing gene-models one gene at a time

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
