###Begin an interactive session before running the script

#bash 0_v8_AFR.sh <path to genotype plink files + plink file prefix> <path to ge data> <path to covariates file> <path to intermediate directory> <path to scratch directory> <path to output directory of weights file> <path to plink> <comma separated list of datasets to use> <path to gcta>

#<path to genotype plink files + plink file prefix>
	#pathtoplink/plink_prefix_<chr#>.bed/bim/fam -> pathtoplink/plink_prefix_

#--- Set up r environment 
source ~/.bashrc
conda activate r_env

echo "STARTING MAIN SCRIPT"

#--- collect all command line arguments
geno=$1 # PLINK FILES SPLIT BY CHROMOSOME = /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO_gtexEUR/genotype/GTEx_v8_genotype_EUR_HM3_chr
ge=$2 # /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO_gtexEUR/intermediate/Whole_Blood.v8.normalized_expression.bed.gz
covar=$3 # /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO_gtexEUR/intermediate/Whole_Blood.v8.covariates.txt 
intermediate=$4 # /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO_gtexEUR/intermediate
scratch=$5 # /expanse/lustre/scratch/kakamatsu/temp_project/EUR_MAGEPRO
output=$6 # /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO_gtexEUR/weights
plink_executable=$7 # /expanse/lustre/projects/ddp412/kakamatsu/plink2
plink_executable1=$8 # /expanse/lustre/projects/ddp412/kakamatsu/plink #NEED FOR LASSO IN FUSION
datasets=$9 # eqtlgen,ota,his,eur,mesa,genoa
gcta=${10} #/expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust


#--- define paths
gene_beds=${scratch}/gene_beds
plink_pergene=${scratch}/plink_gene

#--- finds individuals with both genotype and gene expression data in GTEx
Rscript 1_v8.NoDownsample.R $geno $ge $intermediate

#Only run once or everytime scratch directory deletes itself
run=2 #run: 1, don't run: 2 
if [ $run -eq 1 ]; 
then
 
 mkdir $gene_beds
 mkdir $plink_pergene

 #--- get genes that have ge data
 Rscript 2_v8.GetMeasuredGenes.R $ge $gene_beds $intermediate
 
 #--- create plink files per gene
 echo "CREATING PLINK FILES PER GENE"
 for i in $gene_beds/* #looping through every bed file created in Rscript 2
 do
 name=`basename $i | sed 's/.bed//g'` #get gene name
 chrom=$(awk '{print $1}' $i) #get chromosome number
 plink_file=${geno}${chrom} 
 $plink_executable --bfile $plink_file --extract bed1 $i --out $plink_pergene/$name --make-bed #extract SNPs in cis-region 
 #NOTE ABOUT PLINK --extract bed1 <file>: file has to be in this specific BED format 
 # chromosome code, start range, end range, set ID, group label 
 rm $plinkdir/$name.log
 done
fi

#--- extract covariate for all people of interest
Rscript 3_v8.covar.AFR.R $covar $intermediate

#--- assign genes in batches to run jobs in "embarassingly" parallel way
# 20 = number of batches -> this runs 20 jobs in parallel
Rscript 4_v8.assignBatch.R 20 $intermediate $plink_pergene

#--- run jobs per batch to compute gene models 
# this part uses the script runmultipopweights.sh -> THIS SHOULD BE CHANGED BASED ON THE HPC CLUSTER
bash multipopweights.sh $datasets $intermediate/genes_assign_Whole_Blood.txt $ge ${scratch} $intermediate $output $plink_executable $plink_executable1 $gcta

