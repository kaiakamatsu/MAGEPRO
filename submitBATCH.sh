#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 04:00:00
#SBATCH -J MAGEPROsubmit
#SBATCH -A csd832
#SBATCH -o ../working_err/MAGEPROsubmit.%j.%N.out
#SBATCH -e ../working_err/MAGEPROsubmit.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

source ~/.bashrc
conda activate r_env

Rscript RUN_MAGEPRO_PIPELINE_IMPACT.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/1KG_geno/plink_HM3_variants_EUR/maf_hwe_rate_relatedness_HM3_ingenoa_b38_redo/GEUVADIS_EUR --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/GE/GEUVADIS_normalized_processed_ge_v26_b38_EUR_start_end.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/covar/GEUVADIS_EUR_covariates.txt --out /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_LCL/weights_IMPACT_10%_enet_important_ingenoa_redo --scratch /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_LCL/scratch_ingenoa_redo --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_LCL/intermediate_ingenoa_redo --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --sumstats_dir /expanse/lustre/projects/ddp412/kakamatsu/MAGEPRO_datasets_allinfo --sumstats eqtlgen,mesahis,mesaafr --models SINGLE,META,MAGEPRO --ss 31684,352,233 --hsq_p 1 --verbose 2 --num_covar 36 --IMPACT /expanse/lustre/projects/ddp412/kakamatsu/IMPACT_TRACKS/EUR_processed/means_LCL/IMPACT_707_LCL_HM3_EUR_chr --crossval 5
