#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 10:00:00
#SBATCH -J RUNTWASafr
#SBATCH -A ddp412
#SBATCH -o RUNTWASafr.%j.%N.out
#SBATCH -e RUNTWASafr.%j.%N.err
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
gwas=${home}/gwas/AA_BCX_LETTRELAB/munged/BCX2_${pheno}_AA_GWAMA_munged.sumstats.gz
echo $gwas

wgt=${home}/fusion_multipop/magepro_weight.pos
wgtdir=/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/weights_MAGEPRO_fullsumstats/Whole_Blood/

out=${home}/fusion_multipop/outputMAGEPRO_BCX_LETTER/${pheno}

ld=${home}/eQTLsummary/1000g/1000G_AFR_gtexSNPs/1000G_AFR_

#--force_model AFRONLY

for chr in {1..22}; do \
	Rscript ${home}/fusion_multipop/FUSION.assoc_test_multipop.R --chr $chr --ref_ld_chr $ld --sumstats $gwas --weights $wgt --weights_dir $wgtdir --out ${out}.${chr}.dat; \
done
