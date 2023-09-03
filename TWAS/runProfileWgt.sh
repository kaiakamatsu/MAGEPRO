#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 01:00:00
#SBATCH -J ComputeGeneModelProfiles
#SBATCH -A ddp412
#SBATCH -o ComputeGeneModelProfiles.%j.%N.out
#SBATCH -e ComputeGeneModelProfiles.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

module purge
module load cpu/0.15.4
module load gcc/9.2.0
module load mvapich2/2.3.6
module load slurm
source ~/.bashrc
conda activate r_env

Rscript FUSION.profile_wgt_multipop.R magepro_weights_paths.txt
