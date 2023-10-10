scratch=/expanse/lustre/projects/ddp412/kakamatsu/testingdir/gene_subset/scratch
tmpdir=$scratch/tmp
wd=$scratch/wd
plinkdir=$scratch/plink_gene
plink_exec=../../plink
intermed=/expanse/lustre/projects/ddp412/kakamatsu/testingdir/gene_subset/intermediate
batchfile=$intermed/Genes_Assigned.txt 
geneids=$intermed/All_Genes_Expressed.txt 
ind=$intermed/All_Individuals.txt #all individuals with both genotype and ge data 
samples=$intermed/Sample_IDs.txt 

#--- get column numbers (in ge data) of individuals we are interested in
gene=$1
#grabbing just the individuals we want
alldonors=$(zcat /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/Whole_Blood/Whole_Blood.v8.normalized_expression.bed.gz | head -n 1)
colind=$(echo $alldonors | sed 's| |\n|g' | nl | grep -f $samples | awk 'BEGIN {ORS=","} {print $1}') 
colind2=${colind%,} #remove last comma - at this point we have a comma delimited list of column numbers of individuals we want to extract from the ge data 

$plink_exec --bfile $plinkdir/$gene --allow-no-sex --make-bed --keep $ind --out $wd/$gene
rm $wd/$gene.log 
rowid=$(cat $geneids | nl | grep $gene | awk '{print $1 + 1}') #+1 for col header in ge file - row number of the gene in the ge file 
ge_donors=$(zcat /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GTEX_ge_covar/Whole_Blood/Whole_Blood.v8.normalized_expression.bed.gz | head -n $rowid | tail -n 1 | cut -f $colind2)  #gene expression from the individuals of interest
paste --delimiters=' ' <(cut -d' ' -f1-5 $wd/$gene.fam) <(echo $ge_donors | sed 's/ /\n/g') > $wd/$gene.mod.fam #modifying fam file with ge data 
mv $wd/$gene.mod.fam $wd/$gene.fam 
TMP=$tmpdir/${gene}
OUT=$weights/${gene}
Rscript ../MAGEPRO.R --gene $gene --bfile $wd/$gene --covar $intermed/Covar_All.txt --tmp $TMP --out $OUT --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --PATH_plink ${plink_exec} --sumstats_dir /expanse/lustre/projects/ddp412/kakamatsu/MAGEPRO_datasets --sumstats eqtlgen,genoa,mesahis,eurgtex,mesaafr,ota_CD16p_Mono,ota_CL_Mono,ota_LDG,ota_mDC,ota_Mem_CD4,ota_Mem_CD8,ota_Naive_B,ota_Naive_CD4,ota_Naive_CD8,ota_Neu,ota_NK,ota_pDC,ota_Plasmablast --models SINGLE,META,MAGEPRO --ss 31684,1031,352,574,233,416,416,416,416,416,416,416,416,416,416,416,416,416 --cell_meta ota --hsq_p 1 --verbose 2
rm $wd/$gene.* 
