# run TWAS using MAGEPRO weights

- gwas summary stats should already be processed - see GWAS_munge directory 

- ld reference file are in-sample plink files

- before running TWAS, run the following pipeline to get a summary of Rdat files 
> 0_getRdatSummary.sh 

- to set up weights and run TWAS 
> 1_TWAS_main.sh

- editted FUSION script to exclude META weights for our analysis 
