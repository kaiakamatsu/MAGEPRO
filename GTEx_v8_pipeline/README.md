---pipeline for running magepro on GTEx data

bash 0_v8_GTEx.sh <path to genotype plink files/plink file prefix> <path to ge data> <path to covariates file> <path to intermediate directory> <path to scratch directory> <path to output directory of weights file> <path to plink2> <path to plink> <comma separated list of datasets to use> <path to gcta>



---extra information on command-line args

<path to genotype plink files/plink file prefix>
	pathtoplink/plink_prefix_<chr#>.bed/bim/fam -> pathtoplink/plink_prefix_

<path to intermediate directory>
	a directory where intermediate files are stored 

<path to scratch directory>
	a directory where scratch files are stored



---format of genotype, ge, and covar

genotypes should be in plink bed/bim/fam format 
compatible ge data is from 
	https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_expression_matrices.tar
compatible covar data is from 
	https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_covariates.tar.gz
	5 PC, 5 inferredCov, pcr, platform, sex are regressed out of the gene expression

---NOTE

runmultipopweights.sh
	may have to be changed to fit the slurm system you are using 
