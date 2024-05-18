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
conda activate r_py

Rscript RUN_MAGEPRO_PIPELINE.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/MESA/genotypes/consents_combined/HIS/maf_hwe_rate_filtered_HM3_filtered/MESA_HISchr --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/MESA/gene_expression/processed/MESA_HIS_ge_normalized.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/MESA/covars/ready/MESA_HIS_covariates.txt --out /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_HIS_Monocyte/weights_5.12_benchmark --scratch /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_HIS_Monocyte/scratch --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_HIS_Monocyte/intermediate --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --sumstats_dir /expanse/lustre/projects/ddp412/kakamatsu/SuSiE --sumstats eqtlgen,eurgtex,mesaafr,genoa --models SINGLE,META,PT,SuSiE,PRSCSx,MAGEPRO_fullsumstats,MAGEPRO --ss 31684,574,233,1032 --hsq_p 1 --verbose 2 --num_covar 42 --ldref_pt /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/MESA/genotypes/consents_combined/HIS/maf_hwe_rate_filtered_HM3_filtered/MESA_HISchr --ldref_PRSCSx /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/benchmark/LD_ref --dir_PRSCSx /expanse/lustre/projects/ddp412/kakamatsu/MAGEPRO/BENCHMARK/PRScsx --phi_shrinkage_PRSCSx 1e-7 --pops EUR,EUR,AFR,AFR --susie_pip 8 --susie_beta 9 --susie_cs 10
