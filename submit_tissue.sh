#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 02:00:00
#SBATCH -J MAGEPROsubmit
#SBATCH -A csd832
#SBATCH -o ../working_err/MAGEPROsubmit.%j.%N.out
#SBATCH -e ../working_err/MAGEPROsubmit.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

source ~/.bashrc
conda activate r_py

tissue=$1

cd /expanse/lustre/projects/ddp412/kakamatsu/MAGEPRO/

#--- MAKE DIRECTORIES HERE FOR WEIGHTS, SCRATCH FILES, AND INTERMEDIATE FILES 
mkdir /expanse/lustre/projects/ddp412/kakamatsu/BIMM182/${tissue}_weights
mkdir /expanse/lustre/projects/ddp412/kakamatsu/BIMM182/${tissue}_scratch
mkdir /expanse/lustre/projects/ddp412/kakamatsu/BIMM182/${tissue}_intermediate

#--- CHANGE THE --scratch, --intermed_dir, --out TO THE DIRECTORIES CREATED ABOVE
#--- OTHER FLAGS TO CHANGES 
	# --bfile (to downsampled genotype files)
	# NOTE THAT MY PIPELINE AUTOMATICALLY TAKES THE INTERSECTION OF --bfile and --ge so you only have to downsample in the plink genotype files! 
Rscript /expanse/lustre/projects/ddp412/kakamatsu/MAGEPRO/RUN_MAGEPRO_PIPELINE.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_EUR/maf_hwe_rate_relatedness/GTEx_EUR_chr --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/GTEx_Analysis_v8_eQTL_ge/GTEx_Analysis_v8_eQTL_expression_matrices/${tissue}.v8.normalized_expression.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/GTEx_Analysis_v8_eQTL_covariates/${tissue}.v8.covariates.txt --out /expanse/lustre/projects/ddp412/kakamatsu/BIMM182/${tissue}_weights --scratch /expanse/lustre/projects/ddp412/kakamatsu/BIMM182/${tissue}_scratch --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/BIMM182/${tissue}_intermediate --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --sumstats_dir /expanse/lustre/projects/ddp412/kakamatsu/SuSiE --sumstats eqtlgen,mesahis,mesaafr,genoa --models SINGLE,META,SuSiE,MAGEPRO --ss 31684,352,233,1032 --hsq_p 1 --verbose 2 --ldref_pt /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_EUR/maf_hwe_rate_relatedness/GTEx_EUR_chr --ldref_PRSCSx /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/benchmark/LD_ref --dir_PRSCSx /expanse/lustre/projects/ddp412/kakamatsu/MAGEPRO/BENCHMARK/PRScsx --phi_shrinkage_PRSCSx 1e-7 --pops EUR,AMR,AFR,AFR --susie_pip 8 --susie_beta 9 --susie_cs 10 --num_covar ALL
