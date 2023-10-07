#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH -J merge_split
#SBATCH -A ddp412
#SBATCH -o merge_split.%j.%N.out
#SBATCH -e merge_split.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

module purge
module load cpu/0.15.4
module load gcc/9.2.0
module load mvapich2/2.3.6
module load slurm
source ~/.bashrc
conda activate r_env

basedir=$1 #directory where the parquet files and snp reference files are stored
outputgene=$2

# Rscript split_merge.R <eqtl full sumstats parquet file> <snps reference file> <output>  

for c in {1..22}
do
	Rscript split_merge.R ${basedir}/Whole_Blood.v8.EUR.allpairs.chr${c}.parquet ${basedir}/${c}snpMAP.txt $outputgene 
done
