#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 08:00:00
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

Rscript RUN_GENE_MODEL_TEST.R --bfile /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/genotypes_1kg_hapmap/1000G_YRI_HM3_chr --ge /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/GEUVADIS_normalized_processed_ge.txt.gz --covar /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/covar/GEUVADIS_covar.txt --out /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_VALIDATION/weights --scratch /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_VALIDATION/scratch --intermed_dir /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_VALIDATION/intermediate --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --rerun FALSE --hsq_p 1 --verbose 2 --genemodel /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExAFR_WHOLEBLOOD_MAGEPRO/weights --num_covar 4
