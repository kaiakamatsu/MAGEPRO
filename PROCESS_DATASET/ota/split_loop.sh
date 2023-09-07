for i in CD16p_Mono Mem_CD4 Naive_B Naive_CD8 Plasmablast mDC CL_Mono LDG Mem_CD8 NK Naive_CD4 Neu pDC
do
	echo $i	
	# Rscript splitGenesFilterSnps.R $i <path to directory where sumstats are stores> <path to snp reference file> <path to directory to store outputs> 
	Rscript splitGenesFilterSnps.R $i /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/OTA_nominal/nominal_data /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/GTEx_plink/GTEx_v8_genotype_AFR_HM3_exclude_dups.allchr.bim /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/OTA_nominal/genes_format
done
