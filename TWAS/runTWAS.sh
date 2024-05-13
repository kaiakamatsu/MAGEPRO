#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 10:00:00
#SBATCH -J RUNTWASmagepro
#SBATCH -A ddp412
#SBATCH -o ../../working_err/RUNTWASmagepro.%j.%N.out
#SBATCH -e ../../working_err/RUNTWASmagepro.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

module purge
module load cpu/0.15.4
module load gcc/9.2.0
module load mvapich2/2.3.6
module load slurm
source ~/.bashrc
conda activate r_env

home=/expanse/lustre/projects/ddp412/kakamatsu

pheno=$1
gwas=${home}/gwas/AA_BCX_LETTRELAB/munged/BCX2_${pheno}_AA_GWAMA_munged.sumstats.gz # path to GWAS summary statistics, munged
echo $gwas

wgt=${home}/fusion_multipop/magepro_weight.pos # path to .pos file 
wgtdir=/expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/<path to gene models>/ # path to gene models

out=${home}/fusion_multipop/outputMAGEPRO_BCX_LETTER/${pheno} # output directory 

ld=${home}/GENOTYPES_GE_data/<path to ld reference file> # ld reference file

#--force_model SINGLE

for chr in {1..22}; do \
	Rscript FUSION.assoc_test_multipop.R --chr $chr --ref_ld_chr $ld --sumstats $gwas --weights $wgt --weights_dir $wgtdir --out ${out}.${chr}.dat; \
done
