#bash eurgtex_main.sh /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/GTEx_EUR/Whole_Blood_all /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/GTEx_plink/GTEx_v8_genotype_AFR_HM3_exclude_dups.allchr.bim /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/GTEx_EUR/genes


output_full=$1 #output directory for fullsumstats
ref_bim=$2
output_genes=$3

#download files 
#bash loop_download.sh $output_full

#create a reference file to map chr_bp_a1_a0 format to rsid 
#Rscript create_rsid_ref.R $ref_bim $output_full

#read in sumstats, filter for snps of interest and split by gene
sbatch loop_split_merge.sh $output_full $output_genes
