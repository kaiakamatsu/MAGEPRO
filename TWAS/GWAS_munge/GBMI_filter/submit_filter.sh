#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH -t 6:00:00
#SBATCH -J filter_gwas
#SBATCH -A csd832
#SBATCH -o ../../working_err/filter_gwas.%j.%N.out
#SBATCH -e ../../working_err/filter_gwas.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

source ~/.bashrc
conda activate r_env

pop=$1

Rscript filter_hm3.R $pop
