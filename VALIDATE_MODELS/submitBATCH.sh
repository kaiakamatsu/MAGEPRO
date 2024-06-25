#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 06:00:00
#SBATCH -J MAGEPROsubmit
#SBATCH -A csd832
#SBATCH -o ../../working_err/MAGEPROsubmit.%j.%N.out
#SBATCH -e ../../working_err/MAGEPROsubmit.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

source ~/.bashrc
conda activate r_py

Rscript RUN_GENE_MODEL_TEST.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/MESA/genotypes/consents_combined/AA/maf_hwe_rate_filtered_HM3_filtered/MESA_AAchr --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/MESA/gene_expression/processed/MESA_AA_ge_normalized.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/MESA/covars/ready/MESA_AA_covariates.txt --out /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_HIS_Monocyte_APPLY_AA/weights_final_benchmark --scratch /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_HIS_Monocyte_APPLY_AA/scratch --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_HIS_Monocyte_APPLY_AA/intermediate --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --hsq_p 1 --verbose 2 --genemodel /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_HIS_Monocyte/weights_final_benchmark --insample_bim /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_HIS_Monocyte/scratch/plink_gene

Rscript RUN_GENE_MODEL_TEST.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/MESA/genotypes/consents_combined/AA/maf_hwe_rate_filtered_HM3_filtered/MESA_AAchr --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/MESA/gene_expression/processed/MESA_AA_ge_normalized.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/MESA/covars/ready/MESA_AA_covariates.txt --out /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_EUR_Monocyte_APPLY_AA/weights_final_benchmark --scratch /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_EUR_Monocyte_APPLY_AA/scratch --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_EUR_Monocyte_APPLY_AA/intermediate --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --hsq_p 1 --verbose 2 --genemodel /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_EUR_Monocyte/weights_final_benchmark --insample_bim /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/MESA_EUR_Monocyte/scratch/plink_gene
