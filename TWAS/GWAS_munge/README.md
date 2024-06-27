# BLOOD TRAITS GWAS

GWAS sumstats from: http://www.mhi-humangenetics.org/en/resources/

GTEx lookup table from: https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz

step 1: create a lookup table to convert snp naming convention (build 37) used in BCX LETTRELAB GWAS sumstats to dbSNP151, filter for HM3 snps in gtex

- LETTRELAB_lookup_table
- LETTRELAB_filter

step 2: run ldsc/munge_sumstats.py 

- munge_loop_LETTRE.sh

# GBMI GWAS 

GWAS sumstats from: https://docs.google.com/spreadsheets/d/1sSU_JfPKs6EZLcY9t3gXsSGPZA1LsP_z/edit?gid=1819800283#gid=1819800283

step 1: filter sumstats for HM3 snps 

- GBMI_filter

step 2: run ldsc/munge_sumstats.py

- munge_loop_GBMI.sh
