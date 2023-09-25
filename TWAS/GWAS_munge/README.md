GWAS sumstats from: http://www.mhi-humangenetics.org/en/resources/

GTEx lookup table from: https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz

step 1: create a lookup table to convert snp naming convention (build 37) used in BCX LETTRELAB GWAS sumstats to dbSNP151, filter for HM3 snps in gtex
	create_lookup.R
	filter_loop.sh	

step 2: run ldsc/munge_sumstats.py 
	munge_loop.sh	
