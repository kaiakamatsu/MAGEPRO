#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 01:00:00
#SBATCH -J MAGEPRO_process
#SBATCH -A ddp412
#SBATCH -o ../../workingerr/MAGEPRO_process.%j.%N.out
#SBATCH -e ../../workingerr/MAGEPRO_process.%j.%N.err
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

output=$1 #/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO_gtexEUR/weights
genes=$2 #/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO_gtexEUR/intermediate/genes_assign_Whole_Blood.txt

Rscript processWeights.R $output $genes
