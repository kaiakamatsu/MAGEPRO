#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 02:00:00
#SBATCH -J TWASprocess
#SBATCH -A csd832
#SBATCH -o ../working_err/TWASprocess.%j.%N.out
#SBATCH -e ../working_err/TWASprocess.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

#script to iteratively run combineOutput.R on all phenotypes - compile all Z scores from all chromosomes into one file

#phenos="Asthma Diabetes Hypertension MonocyteCount RedBloodCellCount SittingHeight BMI EosinophilCount LymphocyteCount PlateletCount RedBloodCellDistributionWidth WhiteBloodCellCount FEV1 FVC BasophilCount"

source ~/.bashrc
conda activate r_env

dataset=$1
anc=$2
anc_gbmi=$3

phenos="BAS HGB MCHC MPV RBC EOS LYM MCV NEU RDW HCT MCH MON PLT WBC"
#phenos="LYM"

for pheno in $phenos; do
	echo $pheno
	Rscript combineOutput.R $pheno LETTER $anc $dataset
done

phenos_gbmi="Asthma COPD Gout HF IPF Stroke VTE"

for pheno in $phenos_gbmi; do
	echo $pheno
	Rscript combineOutput.R $pheno GBMI $anc_gbmi $dataset
done
