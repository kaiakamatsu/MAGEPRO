#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH -t 20:00:00
#SBATCH -J runsims
#SBATCH -A ddp412
#SBATCH -o runsims.%j.%N.out
#SBATCH -e runsims.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

module purge
module load cpu/0.15.4
module load gcc/9.2.0
module load mvapich2/2.3.6
module load slurm
source ~/.bashrc
conda activate r_py

numafr=$1 #comma separated list of 5 different sample sizes to test
heritability=$2 #preset heritability

echo "running simulations with $numafr AFR individuals and a gene with heritability $heritability"

bash run_simulation_1000genes.sh $numafr $heritability
