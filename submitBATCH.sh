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

source ~/.bashrc
conda activate r_env

Rscript RUN_MAGEPRO_PIPELINE_IMPACT.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GENOA/GENOTYPES/phg000927.v1.NHLBI_GENOA_GWAS_EA.genotype-calls-matrixfmt.c1/EA_genotypes_combined_affymetrix_Illumina660K/maf_hwe_rate_relatedness_HM3_filtered_people_with_GE_inGEU_EUR/GENOA_EA_chr --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GENOA/GE/GENOA_EA_ge_normalized.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GENOA/COVARS/EA_GENOA_covariates_ready.txt --out /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GENOA_EA_LCL/weights_IMPACT_enet --scratch /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GENOA_EA_LCL/scratch --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GENOA_EA_LCL/intermediate --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --sumstats_dir /expanse/lustre/projects/ddp412/kakamatsu/MAGEPRO_datasets_allinfo --sumstats eqtlgen,mesahis,mesaafr --models SINGLE,META,MAGEPRO --ss 31684,352,233 --hsq_p 1 --verbose 2 --num_covar 38 --IMPACT /expanse/lustre/projects/ddp412/kakamatsu/IMPACT_TRACKS/EUR_processed/means_LCL/IMPACT_707_LCL_HM3_EUR_chr --crossval 5
