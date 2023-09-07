#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 03:00:00
#SBATCH -J MAGEPRO
#SBATCH -A ddp412
#SBATCH -o ../workingerr/MAGEPRO.%j.%N.out
#SBATCH -e ../workingerr/MAGEPRO.%j.%N.err
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

#--- read command line arguments
batch=$1
dataset=$2
gefile=$3
scratch=$4 #make wd and temp, get plink
intermed=$5
weights=$6
plink_exec=$7
plink_exec1=$8
gcta=$9

echo $plink_exec1

#--- create paths
tmpdir=$scratch/tmp
mkdir $tmpdir
wd=$scratch/wd
mkdir $wd
plinkdir=$scratch/plink_gene
batchfile=$intermed/genes_assign_Whole_Blood.txt 
geneids=$intermed/All_Genes_Expressed.txt 
ind=$intermed/All_Individuals.txt #all individuals with both genotype and ge data 

#--- get column numbers (in ge data) of individuals we are interested in
genes=$(awk '$2 == '${batch}' {print $1}' $batchfile) #all genes with the batch number passed in to the script 
#grabbing just the individuals we want
alldonors=$(zcat $gefile | head -n 1)
colind=$(echo $alldonors | sed 's| |\n|g' | nl | grep -f $ind | awk 'BEGIN {ORS=","} {print $1}') 
colind2=${colind%,} #remove last comma - at this point we have a comma delimited list of column numbers of individuals we want to extract from the ge data 

for gene in $genes
do
filecheck=$plinkdir/$gene.bed
if test -f "$filecheck"
then
    $plink_exec --bfile $plinkdir/$gene --make-bed --keep $ind --out $wd/$gene
    rm $wd/$gene.log 
    rowid=$(cat $geneids | nl | grep $gene | awk '{print $1 + 1}') #+1 for col header in ge file - row number of the gene in the ge file 
    ge_donors=$(zcat $gefile | head -n $rowid | tail -n 1 | cut -f $colind2)  #gene expression from the individuals of interest
    paste --delimiters='\t' <(cut -f1-5 $wd/$gene.fam) <(echo $ge_donors | sed 's/ /\n/g') > $wd/$gene.mod.fam #modifying fam file with ge data 
    mv $wd/$gene.mod.fam $wd/$gene.fam 
    TMP=$tmpdir/${gene}
    OUT=$weights/${gene}
    Rscript ComputeMultiPopWeights.R --gene $gene --bfile $wd/$gene --covar $intermed/Covar_All.txt --hsq_p 1 --tmp $TMP --out $OUT --PATH_gcta $gcta --verbose 2 --PATH_plink ${plink_exec1} --datasets $dataset #FUSION USES PLINK1.9 FOR --lasso FLAG
    rm $wd/$gene.* 

fi

done




