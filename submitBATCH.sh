#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 04:00:00
#SBATCH -J MAGEPROsubmit
#SBATCH -A ddp412
#SBATCH -o ../working_err/MAGEPROsubmit.%j.%N.out
#SBATCH -e ../working_err/MAGEPROsubmit.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

module purge
module load cpu/0.15.4
module load gcc/9.2.0
module load mvapich2/2.3.6
module load slurm
source ~/.bashrc
conda activate r_env

Rscript RUN_MAGEPRO_PIPELINE.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_AFR/GTEx_v8_genotype_AFR_HM3_exclude_dups. --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/Whole_Blood/Whole_Blood.v8.normalized_expression.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/Whole_Blood/Whole_Blood.v8.covariates.txt --out /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExAFR_WHOLEBLOOD_MAGEPRO/weights --scratch /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExAFR_WHOLEBLOOD_MAGEPRO/scratch --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExAFR_WHOLEBLOOD_MAGEPRO/intermediate --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --sumstats_dir /expanse/lustre/projects/ddp412/kakamatsu/MAGEPRO_datasets --sumstats eqtlgen,genoa,mesahis,eurgtex,mesaafr,ota_CD16p_Mono,ota_CL_Mono,ota_LDG,ota_mDC,ota_Mem_CD4,ota_Mem_CD8,ota_Naive_B,ota_Naive_CD4,ota_Naive_CD8,ota_Neu,ota_NK,ota_pDC,ota_Plasmablast --models SINGLE,META,MAGEPRO --cell_meta ota --ss 31684,1031,352,574,233,416,416,416,416,416,416,416,416,416,416,416,416,416 --hsq_p 1 --verbose 2 --subset_genes /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExAFR_WHOLEBLOOD_MAGEPRO/intermediate/gene_subset.txt
