#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 05:00:00
#SBATCH -J MAGEPROvalidate
#SBATCH -A ddp412
#SBATCH -o ../../working_err/MAGEPROvalidate.%j.%N.out
#SBATCH -e ../../working_err/MAGEPROvalidate.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"

module purge
module load cpu/0.15.4
module load gcc/9.2.0
module load mvapich2/2.3.6
module load slurm
source ~/.bashrc
conda activate r_env

#--- read command line arguments
batch=$1
gefile=$2
scratch=$3 
intermed=$4
weights=$5
plink_exec=$6
gcta=$7
resid=${8}
hsq_p=${9}
hsq_set=${10}
verbose=${11}
noclean=${12}
save_hsq=${13}
genemodel=${14}

#--- create paths
tmpdir=$scratch/tmp
mkdir $tmpdir
wd=$scratch/wd
mkdir $wd
plinkdir=$scratch/plink_gene
batchfile=$intermed/Genes_Assigned.txt 
geneids=$intermed/All_Genes_Expressed.txt 
ind=$intermed/All_Individuals.txt #all individuals with both genotype and ge data 
samples=$intermed/Sample_IDs.txt 

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
    $plink_exec --bfile $plinkdir/$gene --allow-no-sex --make-bed --keep $ind --out $wd/$gene
    rm $wd/$gene.log 
    rowid=$(cat $geneids | nl | grep $gene | awk '{print $1 + 1}') #+1 for col header in ge file - row number of the gene in the ge file 
    ge_donors=$(zcat $gefile | head -n $rowid | tail -n 1 | cut -f $colind2)  #gene expression from the individuals of interest
    paste --delimiters=' ' <(cut -d' ' -f1-5 $wd/$gene.fam) <(echo $ge_donors | sed 's/ /\n/g') > $wd/$gene.mod.fam #modifying fam file with ge data 
    mv $wd/$gene.mod.fam $wd/$gene.fam 
    TMP=$tmpdir/${gene}
    OUT=$weights/${gene}
    Rscript TEST_GENE_MODELS.R --gene $gene --bfile $wd/$gene --genemodel $genemodel --covar $intermed/Covar_All.txt --tmp $TMP --out $OUT --PATH_gcta $gcta --PATH_plink ${plink_exec} --resid $resid --hsq_p $hsq_p --hsq_set $hsq_set --verbose $verbose --noclean $noclean --save_hsq $save_hsq
    rm $wd/$gene.* 
fi

done




