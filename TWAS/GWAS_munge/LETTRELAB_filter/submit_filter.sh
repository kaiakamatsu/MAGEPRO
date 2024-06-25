#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH -t 01:30:00
#SBATCH -J filter
#SBATCH -A ddp412
#SBATCH -o filter.%j.%N.out
#SBATCH -e filter.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

#Compute Weights
source ~/.bashrc
conda activate r_env

pop=(Trans)

for p in "${pop[@]}"; do
	bash filter_loop.sh $p	
done
