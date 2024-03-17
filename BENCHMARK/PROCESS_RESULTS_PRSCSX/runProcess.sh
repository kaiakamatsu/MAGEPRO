#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 01:00:00
#SBATCH -J PRSCSX_process
#SBATCH -A ddp412
#SBATCH -o ../../working_err/PRSCSX_process.%j.%N.out
#SBATCH -e ../../working_err/PRSCSX_process.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

source ~/.bashrc
conda activate r_env

output=$1 #output dir
weightsdir=$2 #weights dir
genes=$3 #file with genes
models=$4 #models used SINGLE,META,MAGEPRO

Rscript processWeights.R $output $weightsdir $genes $models
