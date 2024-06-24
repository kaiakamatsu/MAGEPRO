#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH -t 22:00:00
#SBATCH -J runsimsLD
#SBATCH -A ddp412
#SBATCH -o runsimsLD.%j.%N.out
#SBATCH -e runsimsLD.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

source ~/.bashrc
conda activate r_py

numafr=$1 #comma separated list of 5 different sample sizes to test
heritability=$2
eur_geno_prefix=$3
afr_geno_prefix=$4
amr_geno_prefix=$5
out=$6

echo "running simulations with $numafr AFR individuals and a gene with heritability $heritability"

bash run_highLD_simulation_1000.sh $numafr $heritability $eur_geno_prefix $afr_geno_prefix $amr_geno_prefix $out
