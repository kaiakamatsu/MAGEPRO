#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 00:15:00
#SBATCH -J MAGEPRO
#SBATCH -A csd832
#SBATCH -o ../working_err/MAGEPRO.%j.%N.out
#SBATCH -e ../working_err/MAGEPRO.%j.%N.err
#SBATCH --export=ALL
#SBATCH --constraint="lustre"
#--- EDIT ABOVE TO SUIT YOUR HPC CLUSTER

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
sumstats_dir=$8
sumstats=$9
models=${10}
ss=${11}
resid=${12}
hsq_p=${13}
lassohsq=${14}
hsq_set=${15}
crossval=${16}
verbose=${17}
noclean=${18}
save_hsq=${19}
ldref_pt=${20}
prune_r2=${21}
threshold_p=${22}
ldref_PRSCSx=${23}
dir_PRSCSx=${24}
phi_shrinkage_PRSCSx=${25}
pops=${26}
susie_pip=${27}
susie_beta=${28}
susie_cs=${29}
impact_path=${30}
ldref_dir=${31}
ldrefs=${32}
cl_thresh=${33}


#--- create directory for temporary/working files 
#--- NOTE: create "tmp" and "wd" where your system can write/delete files efficiently (EDIT TO SUIT YOUR HPC CLUSTER)
#tmpdir=$scratch/tmp
tmpdir=/expanse/lustre/projects/ddp412/sgolzari/job_$SLURM_JOBID/tmp
mkdir -p $tmpdir
#wd=$scratch/wd
wd=/expanse/lustre/projects/ddp412/sgolzari/job_$SLURM_JOBID/wd
mkdir -p $wd

#--- define file paths 
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
    $plink_exec --bfile $plinkdir/$gene --allow-no-sex --make-bed --keep $ind --indiv-sort f $ind --out $wd/$gene
    rm $wd/$gene.log 
    chr=$(cat ${wd}/${gene}.bim | head -n 1 | awk '{print $1}')
    rowid=$(cat $geneids | nl | grep $gene | awk '{print $1 + 1}') #+1 for col header in ge file - row number of the gene in the ge file 
    ge_donors=$(zcat $gefile | head -n $rowid | tail -n 1 | cut -f $colind2)  #gene expression from the individuals of interest
    paste --delimiters=' ' <(cut -d' ' -f1-5 $wd/$gene.fam) <(echo $ge_donors | sed 's/ /\n/g') > $wd/$gene.mod.fam #modifying fam file with ge data 
    mv $wd/$gene.mod.fam $wd/$gene.fam 
    TMP=$tmpdir/${gene}
    OUT=$weights/${gene}
    cmd="Rscript MAGEPRO.R --gene $gene --bfile $wd/$gene --covar $intermed/Covar_All.txt --tmp $TMP --out $OUT --PATH_gcta $gcta --PATH_plink ${plink_exec} --sumstats_dir $sumstats_dir --sumstats $sumstats --models $models --ss $ss --resid $resid --hsq_p $hsq_p --lassohsq $lassohsq --hsq_set $hsq_set --crossval $crossval --verbose $verbose --noclean $noclean --save_hsq $save_hsq --ldref_pt $ldref_pt --prune_r2 $prune_r2 --threshold_p $threshold_p --ldref_PRSCSx $ldref_PRSCSx --dir_PRSCSx $dir_PRSCSx --phi_shrinkage_PRSCSx $phi_shrinkage_PRSCSx --pops $pops --susie_pip $susie_pip --susie_beta $susie_beta --susie_cs $susie_cs --impact_path ${impact_path}${chr}.txt --ldref_dir ${ldref_dir} --ldrefs ${ldrefs} --cl_thresh ${cl_thresh}"
    echo $cmd
    $cmd
    rm $wd/$gene.* 
fi

done





Rscript MAGEPRO.R --gene ENSG00000000419 
--bfile /expanse/lustre/projects/ddp412/sgolzari/job_31624887/wd/ENSG00000000419 
--tmp /expanse/lustre/projects/ddp412/sgolzari/job_31624887/tmp/ENSG00000000419 
--out /expanse/lustre/projects/ddp412/sgolzari/testMAGEPRO/_weights/ENSG00000000419 
--PATH_plink /expanse/lustre/projects/ddp412/sgolzari/plink/plink 
--sumstats_dir /expanse/lustre/projects/ddp412/sgolzari/fine_mapping 
--sumstats genoa --models MAGEPRO --ss 1032 --resid FALSE --hsq_p 1 
--lassohsq 0.05 --hsq_set NA --crossval 5 --verbose 2 --noclean FALSE --save_hsq FALSE 
--ldref_pt /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/1KG_geno/plink_HM3_variants_EUR/maf_hwe_rate_relatedness_HM3_ingenoa_b38_redo/GEUVADIS_EUR 
--prune_r2 NA --threshold_p NA --ldref_PRSCSx NA --dir_PRSCSx PRScsx --phi_shrinkage_PRSCSx 1e-06 
--pops EUR,EUR,AMR,AFR --susie_pip 8 --susie_beta 9 --susie_cs 10 --impact_path NA20.txt 
--ldref_dir /expanse/lustre/projects/ddp412/sgolzari/ldref --ldrefs GENOA_AA --cl_thresh 0.97
