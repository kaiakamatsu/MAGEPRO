suppressMessages(library("optparse"))
suppressMessages(library("data.table"))

# --- PARSE COMMAND LINE ARGUMENTS 

option_list = list(
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus chr number and bed/bim/fam) [required] \n
	      MAGEPRO will take samples with both GE and GENOTYPE data"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--scratch", action="store", default=NA, type='character',
              help="Path to scratch directory [required]"),
  make_option("--ge", action="store", default=NA, type='character',
              help="Path to individual-level, normalized gene expression data in matrix format [required]"),
  make_option("--covar", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) [optional]"),
  make_option("--num_covar", action="store", default=NA, type='double',
              help="Number of covariate to use (number of rows to extract from --covar file). 'ALL' for all covariates in the file. 
	      Default assumes gtex covariate file (first 5 PC, first 15 inferredcov, pcr, platform, sex: ideal for N < 150) [optional]"),
  make_option("--num_batches", action="store", default=20, type='double',
              help="Number of batch jobs to split the genes into. Default 20 [optional]"),
  make_option("--rerun", action="store_true", default=FALSE,
              help="Are you rerunning the pipeline? TRUE = skip recreating plink files per gene"),
  make_option("--intermed_dir", action="store", default=".", type='character',
              help="Directory to store intermediate files"),
  make_option("--subset_genes", action="store", default=NA, type='character',
              help="Path to file with gene of interests in one column"),
  make_option("--sumstats_dir", action="store", default=NA, type='character',
              help="Path to external sumstats"),
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Comma-separated list of external datasets to include"),
  make_option("--models", action="store", default="SINGLE,META,MAGEPRO",
              help="Comma-separated list of models to use \n 
	      SINGLE = single ancestry FUSION lasso approach \n
	      META = ss-weighted meta-analysis of datasets \n 
	      PT = pruning and thresholding \n
	      SuSiE = sum of single effects regression \n
	      SuSiE_IMPACT = sum of single effects regression with IMPACT scores as priors \n
	      PRSCSx = PRS-CSx multi-ancestry PRS method \n 
	      MAGEPRO_fullsumstats = magepro model, no sparsity \n
	      MAGEPRO = magepro model"),
  make_option("--ss", action="store", default=NA, type='character',
              help="Comma-separated list of sample sizes of sumstats (in the same order)"), 
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gcta", action="store", default="gcta_nr_robust", type='character',
              help="Path to gcta executable [%default]"),
  make_option("--resid", action="store_true", default=FALSE,
              help="Also regress the covariates out of the genotypes [default: %default]"),              
  make_option("--hsq_p", action="store", default=0.01, type='double',
              help="Minimum heritability p-value for which to compute weights [default: %default]"),
  make_option("--lassohsq", action="store", default=0.05, type='double',
              help="Backup heritability value to use in lasso regression if skipping heritability estimate or gcta fails (when ss is very small)"),
  make_option("--hsq_set", action="store", default=NA, type='double',
              help="Skip heritability estimation and set hsq estimate to this value [optional]"),
  make_option("--crossval", action="store", default=5, type='double',
              help="How many folds of cross-validation, 0 to skip [default: %default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"),
  make_option("--save_hsq", action="store_true", default=FALSE,
              help="Save heritability results even if weights are not computed [default: %default]"), 
  make_option("--ldref_pt", action="store", default=NA, type='character',
              help="Path to LD reference file for pruning and thresholding, prefix of plink formatted files (assumed to be split by chr) \n 
	      ex. path/file_chr for path/file_chr1.bed/bim/fam "),
  make_option("--prune_r2", action="store", default=NA, type='numeric',
              help="Pruning threshold to use in P+T. If not provided, it will be tuned via 5-fold cross-validation"),
  make_option("--threshold_p", action="store", default=NA, type='numeric',
              help="p-value threshold to use in P+T. If not provided, it will be tuned via 5-fold cross-validation"),
  make_option("--ldref_PRSCSx", action="store", default=NA, type='character',
              help="Path to LD reference directory for PRS-CSx, made available by PRS-CSx github"),
  make_option("--dir_PRSCSx", action="store", default="PRScsx", type='character',
              help="Path to PRS-CSx directory, containing executable (github repo)"),
  make_option("--phi_shrinkage_PRSCSx", action="store", default=1e-6, type='numeric',
              help="Shrinkage parameter for PRS-CSx"),
  make_option("--pops", action="store", default=NA, type='character',                               
	      help="Comma separated list of ancestries of datasets for PRS-CSx (ex. EUR,EAS,AFR)"),
  make_option("--susie_pip", action="store", default=8, type='numeric',
	      help="Column number in external datasets where susie pips are stored"),
  make_option("--susie_beta", action="store", default=9, type='numeric',
              help="Column number in external datasets where susie coefs are stored"),
  make_option("--susie_cs", action="store", default=10, type='numeric',                                             
	      help="Column number in external datasets where susie credible set groups are stored"),
  make_option("--impact_path", action="store", default=NA, type='character',
              help="path to file with impact scores for each snp")
)

opt = parse_args(OptionParser(option_list=option_list))
if ( opt$verbose >= 1 ) print(opt)

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}

# --- COLLECT PEOPLE WITH BOTH GE AND GENOTYPE DATA
if ( opt$verbose >= 1 ) cat("### EXTRACTING SAMPLES WITH BOTH GE AND GENOTYPE DATA\n")
arg = paste("Rscript MAGEPRO_PIPELINE/1_CollectSamples.R", opt$bfile, opt$ge, opt$intermed_dir, sep=" ")
system( arg , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
if ( opt$verbose >= 1 ) cat("### WROTE RESULTS TO ", opt$intermed_dir,"/All_Individuals.txt\n", sep = "")

# --- WRITE A FILE OF JUST SAMPLE-IDS
arg = paste0("cut -f2 ", opt$intermed_dir,"/All_Individuals.txt > ", opt$intermed_dir,"/Sample_IDs.txt")
system( arg , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )

# --- MAKE DIRECTORY FOR GENE BEDS AND PLINK FILES PER GENE 
gene_beds_dir = paste0(opt$scratch, "/gene_beds")
plink_gene_dir = paste0(opt$scratch, "/plink_gene")

if (! opt$rerun){
system( paste("mkdir", gene_beds_dir, sep = " "), ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
system( paste("mkdir", plink_gene_dir, sep = " "), ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )

# --- CREATE GENE BEDS (skip if --rerun)
if ( opt$verbose >= 1 ) cat("### CREATING GENE BEDS FOR GENES OF INTEREST\n")
arg = paste("Rscript MAGEPRO_PIPELINE/2_GeneBed.R", opt$ge, gene_beds_dir, opt$intermed_dir, opt$subset_genes, sep = " ")
system( arg , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
if ( opt$verbose >= 1 ) cat("### WROTE GENE BEDS IN ", gene_beds_dir, "\n", sep = "")
if ( opt$verbose >= 1 ) cat("### WROTE FILE OF ALL GENES IN ANALYSIS IN ", opt$intermed_dir, "/All_Genes_Expressed.txt\n", sep = "")

# --- CREATE PLINK FILES PER GENE (skip if --rerun)
if ( opt$verbose >= 1 ) cat("### CREATING PLINK FILES PER GENE: MAY TAKE A WHILE\n")

bed_files = list.files(path = gene_beds_dir)
for (b in bed_files){
	data = fread( paste0(gene_beds_dir, "/", b), header = F, sep = "\t" )
	name = data$V4[1]
	chr = data$V1[1]
	start = data$V2[1]
	end = data$V3[1]
	plink_file = paste0(opt$bfile, chr)
	arg = paste(opt$PATH_plink, "--bfile", plink_file, "--chr", chr, "--from-bp", start, "--to-bp", end, "--allow-no-sex", "--make-bed", "--out", paste0(plink_gene_dir, "/", name), sep = " " )
	system( arg , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
	system( paste0("rm -rf ", plink_gene_dir, "/", name, ".log"), ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
}
}else{
if ( opt$verbose >= 1 ) cat("### SKIPPING CREATION OF PLINK FILES PER GENE. MAKE SURE THESE FILE EXISTS IN ", plink_gene_dir, "\n", sep = "")
}

# --- PREPARING COVARIATE FILE, DICTATE NUMBER OF COVARIATES TO EXTRACT 
if (!is.na(opt$covar)){
if ( opt$verbose >= 1 ) cat("### PREPARING COVARIATE FILE\n")
arg = paste("Rscript MAGEPRO_PIPELINE/3_PrepCovar.R", opt$covar, opt$intermed_dir, opt$num_covar, sep = " ")
system( arg , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
if ( opt$verbose >= 1 ) cat("### WROTE COVARIATE FILE IN ", opt$intermed_dir, "/Covar_All.txt\n", sep = "")
}else{
if ( opt$verbose >= 1 ) cat("### SKIPPING COVAR FILE PROCESSING \n")
}

# --- SPLIT UP GENES INTO BATCHES, SUBSET TO GENES OF INTEREST
if ( opt$verbose >= 1 ) cat("### ASSIGNING GENE BATCHES \n")
arg = paste("Rscript MAGEPRO_PIPELINE/4_AssignBatches.R", opt$intermed_dir, opt$num_batches, opt$subset_genes, sep = " ")
system( arg , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
if ( opt$verbose >= 1 ) cat("### WROTE GENE BATCH FILE IN ", opt$intermed_dir, "/Genes_Assigned.txt\n", sep = "")


# --- CREATE DIRECTORY FOR STANDARD OUT AND ERROR FILES FROM BATCH JOBS
if ( opt$verbose >= 1 ) cat("### CREATING DIRECTORY ../working_err FOR BATCH JOB OUT/ERROR FILES \n")
system( "mkdir ../working_err" , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )

# --- RUN SLURM JOBS PER BATCH 
# read batch file and split genes per batch number 
if ( opt$verbose >= 1 ) cat("### RUNNING JOBS \n")
batches <- c(1:opt$num_batches)
for (batch in batches){
arg = paste("sbatch MAGEPRO_PIPELINE/5_RunJobs.sh", batch, opt$ge, opt$scratch, opt$intermed_dir, opt$out, opt$PATH_plink, opt$PATH_gcta, opt$sumstats_dir, opt$sumstats, opt$models, opt$ss, opt$resid, opt$hsq_p, opt$lassohsq, opt$hsq_set , opt$crossval, opt$verbose, opt$noclean, opt$save_hsq, opt$ldref_pt, opt$prune_r2, opt$threshold_p, opt$ldref_PRSCSx, opt$dir_PRSCSx, opt$phi_shrinkage_PRSCSx, opt$pops, opt$susie_pip, opt$susie_beta, opt$susie_cs, opt$impact_path, sep = " ") # you may have to edit this script "5_RunJobs.sh" to suit your HPC cluster
system( arg , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
}

# --- PROCESS RESULTS WHEN JOBS FINISH RUNNING
if ( opt$verbose >= 1 ) cat("### WHEN JOBS FINISH RUNNING, USE THE SCRIPTS IN PROCESS_RESULTS DIRECTORY TO SUMMARIZE GENE MODELS \n")


