suppressMessages(library("optparse"))
suppressMessages(library("data.table"))

# --- PARSE COMMAND LINE ARGUMENTS 

option_list = list(
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) [required] \n
	      MAGEPRO will take samples with both GE and GENOTYPE data"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--genemodel", action="store", default=NA, type='character',
              help="Path to gene model to apply"),
  make_option("--scratch", action="store", default=NA, type='character',
              help="Path to scratch directory [required]"),
  make_option("--ge", action="store", default=NA, type='character',
              help="Path to individual-level, normalized gene expression data in matrix format [required]"),
  make_option("--covar", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) [optional]"),
  make_option("--num_covar", action="store", default=NA, type='double',
              help="Number of covariate to use (number of rows to extract from --covar file). 
	      Default assumes gtex covariate file (first 5 PC, first 5 inferredcov, pcr, platform, sex) [optional]"),
  make_option("--num_batches", action="store", default=20, type='double',
              help="Number of batch jobs to split the genes into. Default 20 [optional]"),
  make_option("--rerun", action="store_true", default=FALSE,
              help="Are you rerunning the pipeline? TRUE = skip recreating plink files per gene"),
  make_option("--intermed_dir", action="store", default="", type='character',
              help="Directory to store intermediate files"),
  make_option("--subset_genes", action="store", default=NA, type='character',
              help="Path to file with gene of interests in one column"),
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gcta", action="store", default="gcta_nr_robust", type='character',
              help="Path to gcta executable [%default]"),
  make_option("--resid", action="store_true", default=FALSE,
              help="Also regress the covariates out of the genotypes [default: %default]"),              
  make_option("--hsq_p", action="store", default=0.01, type='double',
              help="Minimum heritability p-value for which to compute weights [default: %default]"),
  make_option("--hsq_set", action="store", default=NA, type='double',
              help="Skip heritability estimation and set hsq estimate to this value [optional]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"),
  make_option("--save_hsq", action="store_true", default=FALSE,
              help="Save heritability results even if weights are not computed [default: %default]"), 
  make_option("--insample_bim", action="store", default=NA, type = 'character',
              help="Path to gene-specific bim files of in-sample cohort. These files will be used to define cis window snps in this out-of-sample cohort.")
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
arg = paste("Rscript ../MAGEPRO_VALIDATE_PIPELINE/1_CollectSamples.R", opt$bfile, opt$ge, opt$intermed_dir, sep=" ")
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
arg = paste("Rscript ../MAGEPRO_VALIDATE_PIPELINE/2_GeneBed.R", opt$ge, gene_beds_dir, opt$intermed_dir, opt$subset_genes, sep = " ")
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
	base_name = strsplit(name, split = "[.]")[[1]][1]
	plink_file = paste0(opt$bfile, chr)
	arg = paste(opt$PATH_plink, "--bfile", plink_file, "--extract", paste0(opt$insample_bim, "/", base_name, ".bim"), "--allow-no-sex", "--make-bed", "--out", paste0(plink_gene_dir, "/", name), sep = " " )
	system( arg , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
	system( paste0("rm -rf ", plink_gene_dir, "/", name, ".log"), ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
}
}else{
if ( opt$verbose >= 1 ) cat("### SKIPPING CREATION OF PLINK FILES PER GENE. MAKE SURE THESE FILE EXISTS IN ", plink_gene_dir, "\n", sep = "")
}

# --- PREPARING COVARIATE FILE, DICTATE NUMBER OF COVARIATES TO EXTRACT 
if (!is.na(opt$covar)){
if ( opt$verbose >= 1 ) cat("### PREPARING COVARIATE FILE\n")
arg = paste("Rscript ../MAGEPRO_VALIDATE_PIPELINE/3_PrepCovar.R", opt$covar, opt$intermed_dir, opt$num_covar, sep = " ")
system( arg , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
if ( opt$verbose >= 1 ) cat("### WROTE COVARIATE FILE IN ", opt$intermed_dir, "/Covar_All.txt\n", sep = "")
}else{
if ( opt$verbose >= 1 ) cat("### SKIPPING COVAR FILE PROCESSING \n")
}

# --- SPLIT UP GENES INTO BATCHES, SUBSET TO GENES OF INTEREST
if ( opt$verbose >= 1 ) cat("### ASSIGNING GENE BATCHES \n")
arg = paste("Rscript ../MAGEPRO_VALIDATE_PIPELINE/4_AssignBatches.R", opt$intermed_dir, opt$num_batches, opt$subset_genes, sep = " ")
system( arg , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
if ( opt$verbose >= 1 ) cat("### WROTE GENE BATCH FILE IN ", opt$intermed_dir, "/Genes_Assigned.txt\n", sep = "")


# --- CREATE DIRECTORY FOR STANDARD OUT AND ERROR FILES FROM BATCH JOBS
if ( opt$verbose >= 1 ) cat("### CREATING DIRECTORY ../working_err FOR BATCH JOB OUT/ERROR FILES \n")
system( "mkdir ../../working_err" , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )

# --- RUN SLURM JOBS PER BATCH 
# read batch file and split genes per batch number 
if ( opt$verbose >= 1 ) cat("### RUNNING JOBS \n")
batches <- c(1:opt$num_batches)
for (batch in batches){
arg = paste("sbatch ../MAGEPRO_VALIDATE_PIPELINE/5V2_RunTestModels.sh", batch, opt$ge, opt$scratch, opt$intermed_dir, opt$out, opt$PATH_plink, opt$PATH_gcta, opt$resid, opt$hsq_p, opt$hsq_set , opt$verbose, opt$noclean, opt$save_hsq, opt$genemodel, sep = " ") # you may have to edit this script "5_RunJobs.sh" to suit your HPC cluster
system( arg , ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
}

# --- PROCESS RESULTS WHEN JOBS FINISH RUNNING
if ( opt$verbose >= 1 ) cat("### WHEN JOBS FINISH RUNNING, USE THE SCRIPTS IN PROCESS_RESULTS DIRECTORY TO SUMMARIZE GENE MODELS \n")


