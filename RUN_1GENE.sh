#bash RUN_1GENE.sh ENSG00000001630.10 /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_LCL/scratch_filtered_QC /expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_LCL/intermediate_filtered_QC /expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/GE/GEUVADIS_normalized_processed_ge_EUR.bed.gz

gene=$1 
scratch=$2 #/expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_LCL/scratch_filtered_QC
intermed=$3 #/expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GEUVADIS_EUR_LCL/intermediate_filtered_QC
gefile=$4 #/expanse/lustre/projects/ddp412/kakamatsu/GENOTYPES_GE_data/GEUVADIS/GE/GEUVADIS_normalized_processed_ge_EUR.bed.gz

plinkdir=$scratch/plink_gene
batchfile=$intermed/Genes_Assigned.txt
geneids=$intermed/All_Genes_Expressed.txt
ind=$intermed/All_Individuals.txt #all individuals with both genotype and ge data
samples=$intermed/Sample_IDs.txt
wd=/expanse/lustre/projects/ddp412/kakamatsu/scratch_genes/wd
tmpdir=/expanse/lustre/projects/ddp412/kakamatsu/scratch_genes/tmp
weights=/expanse/lustre/projects/ddp412/kakamatsu/scratch_genes/weights

alldonors=$(zcat $gefile | head -n 1)
colind=$(echo $alldonors | sed 's| |\n|g' | nl | grep -f $samples | awk 'BEGIN {ORS=","} {print $1}')
colind2=${colind%,}

/expanse/lustre/projects/ddp412/kakamatsu/plink --bfile $plinkdir/$gene --allow-no-sex --make-bed --keep $ind --out $wd/$gene
    rm $wd/$gene.log
    chr=$(cat ${wd}/${gene}.bim | head -n 1 | awk '{print $1}')
    rowid=$(cat $geneids | nl | grep $gene | awk '{print $1 + 1}') #+1 for col header in ge file - row number of the gene in the ge file
    ge_donors=$(zcat $gefile | head -n $rowid | tail -n 1 | cut -f $colind2)  #gene expression from the individuals of interest
    paste --delimiters=' ' <(cut -d' ' -f1-5 $wd/$gene.fam) <(echo $ge_donors | sed 's/ /\n/g') > $wd/$gene.mod.fam #modifying fam file with ge data
    mv $wd/$gene.mod.fam $wd/$gene.fam
    TMP=$tmpdir/${gene}
    OUT=$weights/${gene}
    Rscript MAGEPRO_IMPACT_PT_resid.R --gene $gene --bfile $wd/$gene --covar $intermed/Covar_All.txt --tmp $TMP --out $OUT --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --PATH_plink /expanse/lustre/projects/ddp412/kakamatsu/plink --sumstats_dir /expanse/lustre/projects/ddp412/kakamatsu/MAGEPRO_datasets_allinfo_PT_IMPACT/ready_LCL --sumstats eqtlgen,mesahis,genoa --models SINGLE,META,MAGEPRO --ss 31684,352,1031 --hsq_p 1 --crossval 5 --verbose 1
    rm $wd/$gene.*
