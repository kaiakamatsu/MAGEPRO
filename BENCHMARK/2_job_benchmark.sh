#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 10:00:00
#SBATCH -J benchmark
#SBATCH -A ddp412
#SBATCH -o ../../working_err/benchmark.%j.%N.out
#SBATCH -e ../../working_err/benchmark.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

source ~/.bashrc
conda activate r_py

batch=$1
dataset=$2
home_dir=$3
ldref=$4
gefile=$5

tmpdir=/scratch/$USER/job_$SLURM_JOBID/tmp
mkdir $tmpdir
wd=/scratch/$USER/job_$SLURM_JOBID/wd
mkdir $wd
weights=${home_dir}/weights_PRSCSx_alphasonly_eqtlgen_mesahis_mesaafr_10cv
mkdir $weights

#--- define file paths 
plinkdir=$home_dir/scratch/plink_gene
batchfile=$home_dir/intermediate/Genes_Assigned.txt 
geneids=$home_dir/intermediate/All_Genes_Expressed.txt 
ind=$home_dir/intermediate/All_Individuals.txt #all individuals with both genotype and ge data 
samples=$home_dir/intermediate/Sample_IDs.txt 

#--- get column numbers (in ge data) of individuals we are interested in
genes=$(awk '$2 == '${batch}' {print $1}' $batchfile) #all genes with the batch number passed in to the script 
#grabbing just the individuals we want
alldonors=$(zcat $gefile | head -n 1)
colind=$(echo $alldonors | sed 's| |\n|g' | nl | grep -f $samples | awk 'BEGIN {ORS=","} {print $1}') 
colind2=${colind%,} #remove last comma - at this point we have a comma delimited list of column numbers of individuals we want to extract from the ge data 

for gene in $genes
do
filecheck=$plinkdir/$gene.bed
if test -f "$filecheck"
then
    /expanse/lustre/projects/ddp412/kakamatsu/plink --bfile $plinkdir/$gene --allow-no-sex --make-bed --keep $ind --out $wd/$gene
    rm $wd/$gene.log 
    rowid=$(cat $geneids | nl | grep $gene | awk '{print $1 + 1}') #+1 for col header in ge file - row number of the gene in the ge file 
    ge_donors=$(zcat $gefile | head -n $rowid | tail -n 1 | cut -f $colind2)  #gene expression from the individuals of interest
    paste --delimiters=' ' <(cut -d' ' -f1-5 $wd/$gene.fam) <(echo $ge_donors | sed 's/ /\n/g') > $wd/$gene.mod.fam #modifying fam file with ge data 
    mv $wd/$gene.mod.fam $wd/$gene.fam 
    TMP=$tmpdir/${gene}
    OUT=$weights/${gene}    
    Rscript benchmark_no_target_betas.R --gene $gene --bfile $wd/$gene --covar $home_dir/intermediate/Covar_All.txt --hsq_p 1 --tmp $TMP --out $OUT --models lasso --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --gemmaout ${gene}_${tissue} --PATH_gemma /expanse/lustre/projects/ddp412/kakamatsu/gemma-0.98.5-linux-static-AMD64 --verbose 2 --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --datasets $dataset --crossval 10 --ldref ${ldref}
    rm $wd/$gene.* 

fi

done
