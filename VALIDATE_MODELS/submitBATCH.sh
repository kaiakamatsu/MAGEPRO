#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 07:00:00
#SBATCH -J MAGEPROsubmit
#SBATCH -A ddp412
#SBATCH -o ../../working_err/MAGEPROsubmit.%j.%N.out
#SBATCH -e ../../working_err/MAGEPROsubmit.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

module purge
module load cpu/0.15.4
module load gcc/9.2.0
module load mvapich2/2.3.6
module load slurm
source ~/.bashrc
conda activate r_env

Rscript RUN_GENE_MODEL_TEST.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/1KG_geno/plink_HM3_variants_EUR/GEUVADIS_EUR_chr --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/GE/GEUVADIS_normalized_processed_ge_EUR.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/covar/GEUVADIS_EUR_covariates.txt --out /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTEX_EUR_LCL_APPLY_GEUVADIS_EUR_VALIDATION/weights --scratch /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTEX_EUR_LCL_APPLY_GEUVADIS_EUR_VALIDATION/scratch --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTEX_EUR_LCL_APPLY_GEUVADIS_EUR_VALIDATION/intermediate --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --hsq_p 1 --verbose 2 --genemodel /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExEUR_LCL/weights_renamed --num_covar ALL
