#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 07:00:00
#SBATCH -J MAGEPROsubmit
#SBATCH -A csd832
#SBATCH -o ../working_err/MAGEPROsubmit.%j.%N.out
#SBATCH -e ../working_err/MAGEPROsubmit.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

source ~/.bashrc
conda activate r_py

Rscript RUN_MAGEPRO_PIPELINE.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_EUR/maf_hwe_rate_relatedness/GTEx_EUR_chr --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/GTEx_Analysis_v8_eQTL_ge/GTEx_Analysis_v8_eQTL_expression_matrices/Lung.v8.normalized_expression.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/GTEx_Analysis_v8_eQTL_covariates/Lung.v8.covariates.txt --out /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExEUR_LUNG/weights --scratch /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExEUR_LUNG/scratch --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExEUR_LUNG/intermediate --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --sumstats_dir /expanse/lustre/projects/ddp412/kakamatsu/SuSiE --sumstats eqtlgen,mesahis,mesaafr,genoa --models SINGLE,META,PT,SuSiE,PRSCSx,MAGEPRO_fullsumstats,MAGEPRO --ss 31684,352,233,1032 --hsq_p 1 --verbose 2 --ldref_pt /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_EUR/maf_hwe_rate_relatedness/GTEx_EUR_chr --ldref_PRSCSx /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/benchmark/LD_ref --dir_PRSCSx /expanse/lustre/projects/ddp412/kakamatsu/MAGEPRO/BENCHMARK/PRScsx --phi_shrinkage_PRSCSx 1e-7 --pops EUR,AMR,AFR,AFR --susie_pip 8 --susie_beta 9 --susie_cs 10
