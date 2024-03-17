#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 03:00:00
#SBATCH -J MAGEPROsubmit
#SBATCH -A ddp412
#SBATCH -o ../../working_err/MAGEPROsubmit.%j.%N.out
#SBATCH -e ../../working_err/MAGEPROsubmit.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

source ~/.bashrc
conda activate r_env

Rscript RUN_GENE_MODEL_TEST.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_EUR/maf_hwe_rate_relatedness_inGEU/GTEx_EUR_QC_chr --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/LCL/Cells_EBV-transformed_lymphocytes.v8.normalized_expression.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/LCL/Cells_EBV-transformed_lymphocytes.v8.covariates.txt --out /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_APPLY_GTEX_EUR_LCL_VALIDATION/weights_IMPACT_10%_enet_important_insample_ciswindow --scratch /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_APPLY_GTEX_EUR_LCL_VALIDATION/scratch_insample_window --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_APPLY_GTEX_EUR_LCL_VALIDATION/intermediate_insample_window --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --hsq_p 1 --verbose 2 --genemodel /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_LCL/weights_IMPACT_10%_enet_important_b38_QCsnps_endTSS --num_covar ALL --insample_bim /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_LCL/scratch_b38_QC_endTSS/plink_gene
