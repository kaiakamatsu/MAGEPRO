#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH -t 02:00:00
#SBATCH -J eqtlgensplit
#SBATCH -A ddp412
#SBATCH -o eqtlgen.%j.%N.out
#SBATCH -e eqtlgen.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

#Compute Weights
module purge
module load cpu/0.15.4
module load gcc/9.2.0
module load mvapich2/2.3.6
module load slurm
source ~/.bashrc
conda activate r_env

bash eqtlgen_main.sh /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/eQTLGEN/smr-1.3.1-linux-x86_64 /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/eQTLGEN/ /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/eQTLGEN/chr /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/GTEx_plink/GTEx_v8_genotype_AFR_HM3_exclude_dups.allchr.bim /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/eQTLGEN/genes
