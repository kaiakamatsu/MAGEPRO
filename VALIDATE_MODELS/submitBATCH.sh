#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 03:00:00
#SBATCH -J MAGEPROsubmit
#SBATCH -A csd832
#SBATCH -o ../../working_err/MAGEPROsubmit.%j.%N.out
#SBATCH -e ../../working_err/MAGEPROsubmit.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

source ~/.bashrc
conda activate r_py

Rscript RUN_GENE_MODEL_TEST.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_AFR/maf_hwe_rate_relatedness/GTExAFR_chr --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/Whole_Blood/Whole_Blood.v8.normalized_expression.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/Whole_Blood/Whole_Blood.v8.covariates.txt --out /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExEUR_WHOLEBLOOD_crossA_AFR/weights_5.12_benchmark --scratch /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExEUR_WHOLEBLOOD_crossA_AFR/scratch_redo --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExEUR_WHOLEBLOOD_crossA_AFR/intermediate_redo --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --hsq_p 1 --verbose 2 --genemodel /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExEUR_WHOLEBLOOD/weights_5.12_benchmark --insample_bim /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExEUR_WHOLEBLOOD/scratch_filtered/plink_gene
