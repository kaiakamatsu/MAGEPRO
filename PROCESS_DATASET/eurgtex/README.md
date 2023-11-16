# Extracting and preparing data from gtex

see **eurgtex_main.sh** for full pipeline

loop_download.sh
- download raw sumstats file (data posted on google cloud by gtex)

create_rsid_ref.R 
- create a file that maps rsid to chr_#_pos_A1_A0_b38 format 
- NOTE: gtex names snps with "chr_#_pos_A1_A0_b38" in sumstats instead of rsid 

split_merge.R 
- filter full sumstats to only keep snps of interest
- split data by gene 

