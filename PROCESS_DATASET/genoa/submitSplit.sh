#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH -J genoasplit
#SBATCH -A ddp412
#SBATCH -o genoasplit.%j.%N.out
#SBATCH -e genoasplit.%j.%N.err
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


bash run_splitFilter.sh
