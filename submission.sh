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
conda activate r_env

tissue=$1

cd /expanse/lustre/projects/ddp412/sgolzari/MAGEPRO/

#--- MAKE DIRECTORIES HERE FOR WEIGHTS, SCRATCH FILES, AND INTERMEDIATE FILES 
mkdir /expanse/lustre/projects/ddp412/sgolzari/testMAGEPRO/${tissue}_weights
mkdir /expanse/lustre/projects/ddp412/sgolzari/testMAGEPRO/${tissue}_scratch
mkdir /expanse/lustre/projects/ddp412/sgolzari/testMAGEPRO/${tissue}_intermediate

#--- CHANGE THE --scratch, --intermed_dir, --out TO THE DIRECTORIES CREATED ABOVE
#--- OTHER FLAGS TO CHANGES 
        # --bfile (to downsampled genotype files)
        # NOTE THAT MY PIPELINE AUTOMATICALLY TAKES THE INTERSECTION OF --bfile and --ge so you only have to downsample in the plink genotype files! 


Rscript RUN_MAGEPRO_PIPELINE.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/1KG_geno/plink_HM3_variants_EUR/maf_hwe_rate_relatedness_HM3_ingenoa_b38_redo/GEUVADIS_EUR --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/GE/GEUVADIS_normalized_processed_ge_v26_b38_EUR_start_end.bed.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/covar/GEUVADIS_EUR_covariates.txt --out /expanse/lustre/projects/ddp412/sgolzari/testMAGEPRO/_weights --scratch /expanse/lustre/projects/ddp412/sgolzari/testMAGEPRO/_scratch --intermed_dir /expanse/lustre/projects/ddp412/sgolzari/testMAGEPRO/_intermediate --PATH_plink /expanse/lustre/projects/ddp412/sgolzari/plink/plink --PATH_gcta /expanse/lustre/projects/ddp412/sgolzari/MAGEPRO/fusion_twas-master/gcta_nr_robust --sumstats_dir /expanse/lustre/projects/ddp412/sgolzari/fine_mapping/Kai --sumstats genoa --models MAGEPRO --ss 1032 --hsq_p 1 --verbose 2 --num_covar 36 --ldref_pt /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/1KG_geno/plink_HM3_variants_EUR/maf_hwe_rate_relatedness_HM3_ingenoa_b38_redo/GEUVADIS_EUR --pops EUR,EUR,AMR,AFR --susie_pip 8 --susie_beta 9 --susie_cs 10 --ldref_dir  /expanse/lustre/projects/ddp412/sgolzari/ldref --ldrefs GENOA_AA --num_batches 1 --subset_genes subset.txt