#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 05:00:00
#SBATCH -J benchmark
#SBATCH -A ddp412
#SBATCH -o ../../working_err/benchmark.%j.%N.out
#SBATCH -e ../../working_err/benchmark.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

#Compute Weights
module purge
module load cpu/0.15.4
module load gcc/9.2.0
module load mvapich2/2.3.6
module load slurm
source ~/.bashrc
conda activate r_py

dataset=$3
tissue=$1
batch=$2

tmpdir=/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/tmp_${tissue}
mkdir $tmpdir
weights=/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/weights_benchmark_confirm/$tissue
mkdir -p $weights
wd=/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/AFR_$tissue
mkdir $wd
plinkdir=/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/plink_cissnps_AFR #stays the same for all tissues, just the cis SNPs 
home_dir=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles
batchfile=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/genes_assign_Whole_Blood2.txt #change this before redoing. 
genes=$(awk '$2 == '${batch}' {print $1}' $batchfile) #all genes with the batch number passed in to the script 
geneids=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/Genes_Expressed_in_${tissue}.txt #chrX gone, not gone in gefile. but it's at the end so ok. 
gefile=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/${tissue}.v8.normalized_expression.bed.gz
#grabbing just the individuals we want (pop) from all individuals - in this case AFR for this analysis
alldonors=$(zcat $gefile | head -n 1) #all individuals in Whole Blood GTEx data 
ind=$home_dir/All_Individuals_${tissue}.txt #all individuals of interest (AFR) - with both genotype and ge data 
colind=$(echo $alldonors | sed 's| |\n|g' | nl | grep -f $ind | awk 'BEGIN {ORS=","} {print $1}') #DONORS MAY BE TAB DELIMITED OR SPACE - ODDLY CHANGES
colind2=${colind%,} #remove last comma - at this point we have a comma delimited list of column numbers of individuals we want to extract from the ge data 

for gene in $genes
do
filecheck=$plinkdir/$gene.bed
if test -f "$filecheck"
then
    plink2 --bfile $plinkdir/$gene --make-bed --keep $ind --out $wd/$gene 
    rm $wd/$gene.log 
    rowid=$(cat $geneids | nl | grep $gene | awk '{print $1 + 1}') #+1 for col header in ge file - row number of the gene in the ge file 
    ge_donors=$(zcat $gefile | head -n $rowid | tail -n 1 | cut -f $colind2)  #gene expression from the individuals of interest
    paste --delimiters='\t' <(cut -f1-5 $wd/$gene.fam) <(echo $ge_donors | sed 's/ /\n/g') > $wd/$gene.mod.fam #modifying fam file with ge data #GE DONORS MAY BE TAB DELIMITED OR SPACE - ODDLY CHANGES 
    mv $wd/$gene.mod.fam $wd/$gene.fam #confirm that gene expression matches correct individual. Confirmed. Even though list of individuals are out of order. Need to make sure covariates are properly handled below in fusion. 
    TMP=$tmpdir/${gene}_${tissue} 
    OUT=$weights/${tissue}.$gene #CHANGE OUT TO MY DIRECTORY AFTER TESTING
    Rscript benchmark.R --gene $gene --tissue $tissue --bfile $wd/$gene --covar $home_dir/Covar_All_${tissue}.txt --hsq_p 1 --tmp $TMP --out $OUT --models lasso --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --gemmaout ${gene}_${tissue} --PATH_gemma /expanse/lustre/projects/ddp412/kakamatsu/gemma-0.98.5-linux-static-AMD64 --verbose 2 --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --datasets $dataset --crossval 5 --ldref /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/GTEx_plink
    rm $wd/$gene.* 

fi

done
