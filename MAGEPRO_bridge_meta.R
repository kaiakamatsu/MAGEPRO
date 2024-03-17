# load packages
suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('glmnet'))
suppressMessages(library('data.table'))


option_list = list(
  make_option("--gene", action="store", default=NA, type='character',
              help="ENSG ID, for ex ENSG00000013503.9 [required]"),
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files [required]"),
  make_option("--sumstats_dir", action="store", default=NA, type='character',
              help="Path to external sumstats"),
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Comma-separated list of external datasets to include"),
  make_option("--models", action="store", default="SINGLE,META,MAGEPRO",
              help="Comma-separated list of models to use \n 
	      SINGLE = single ancestry approach \n
	      META = ss-weighted meta-analysis \n 
	      MAGEPRO = magepro model"),
  make_option("--ss", action="store", default=NA, type='character',
              help="Comma-separated list of sample sizes of sumstats (in the same order)"), 
  make_option("--cell_meta", action="store", default=NA, type='character',
              help="Comma-separated list of prefixes of eqtl datasets to ss-meta-analyze (--ss required) \n 
	      for ex. --cell-type-meta ota will meta-analyze [ota_CD16, ota_CD4]. \n 
	      NOTE: this cell-type meta-analysis happens before the --meta flag meta-analysis across all datasets. \n 
	      the average sample size across all cell types will be assigned as the new sample size of the cell-type meta-analyzed eqtl dataset"),
  make_option("--pheno", action="store", default=NA, type='character',
              help="Path to molecular phenotype file (PLINK format) [optional, taken from bfile otherwise]"),
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gcta", action="store", default="gcta_nr_robust", type='character',
              help="Path to gcta executable [%default]"),
  make_option("--covar", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) [optional]"),
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
  make_option("--IMPACT", action="store", default=NA, type='character',
              help="Path to IMPACT txt file")
)

# --- PARSE COMMAND LINE ARGS
opt = parse_args(OptionParser(option_list=option_list))
na <- which(opt == "NA") 
opt[na] <- NA

if ( opt$verbose >= 1 ) print(opt)

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}

# --- PREDICTION MODELS (FROM FUSION FRAMEWORK)

# PLINK: LASSO
weights.lasso = function( input , hsq , snp , out=NA ) {
	# PURPOSE: perform LASSO regression using plink --lasso flag
	# input = plink file (phenotype in fam file)
	# hsq = heritability estimate used for calibrating regression
	# snp = vector of snps used in the gene model
	# RETURN: vector of effect sizes for each snp
	if ( is.na(out) ) out = paste(input,".LASSO",sep='')
	arg = paste( opt$PATH_plink , " --allow-no-sex --bfile " , input , " --lasso " , hsq , " --out " , out , sep='' )
	system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT )
	if ( !file.exists(paste(out,".lasso",sep='')) ) {
	cat( paste(out,".lasso",sep='') , " LASSO output did not exist\n" )
	eff.wgt = rep(NA,length(snp))
	} else {
	eff = read.table( paste(out,".lasso",sep=''),head=T,as.is=T)
	eff.wgt = rep(0,length(snp))
	m = match( snp , eff$SNP )
	m.keep = !is.na(m)	
	m = m[m.keep]
	eff.wgt[m.keep] = eff$EFFECT[m]
	}
	return( eff.wgt )
}

# Marginal Z-scores (used for top1)
weights.marginal = function( genos , pheno , beta=F ) {
	# PURPOSE: compute marginal effect size estimates 
	# genos = plink genotype bed
	# pheno = phenotype as a matrix (1 column)
	# RETURN: vector of effect sizes for each snp
	if ( beta ) eff.wgt = t( genos ) %*% (pheno) / ( nrow(pheno) - 1)
	else eff.wgt = t( genos ) %*% (pheno) / sqrt( nrow(pheno) - 1 )
	return( eff.wgt )
}

if ( opt$verbose == 2 ) cat("Predictive Models Prepared \n")

# --- CLEANUP
cleanup = function() {
	# PURPOSE: clean up temporary directory
	# RETURN: 
	if ( ! opt$noclean ) {
		arg = paste("rm -f " , opt$tmp , "*", sep='')
		system(arg)
	}
}


# --- PROCESS-SUMSTATS
datasets_process <- function(bim, dataset, file, wgts){
	# PURPOSE: process sumstats by flipping signs according to effect allele and matching snps
	# bim = bim file of snps	
	# dataset = name of sumstat dataset
	# file = file path to sumstat (make sure it is properly formatted)
	# wgts = current list of wgts 
	# RETURN: updated list of wgts, including the one that was processed in this function
	table <- fread(file, select = c(2, 3, 4, 5)) 
	for (k in 1:nrow(table)){
		match=which(bim$V2 == table[k, 1])
		if (! identical(match, integer(0))){
			if (!is.na(table[k, 2]) && !is.na(table[k, 3])){
				if (as.character(bim$V5[match]) != as.character(table[k, 2]) || as.character(bim$V6[match]) != as.character(table[k, 3])) { 
					if (as.character(bim$V5[match]) == as.character(table[k, 3]) && as.character(bim$V6[match]) == as.character(table[k, 2])){
						table[k, 4] <- table[k, 4] * -1
					}else{
						table[k, 4] <- 0
					}
				}
			} else {
				table[k, 4] <- 0
			}
		}
	}
	m <- match(bim[,2], table[[1]])
	table <- table[m,]
	w <- which(is.na(table[[1]]))
	if(length(w) > 0){table[w,4] <- 0}
	assign(paste0("pred.wgt.", dataset), table[[4]], envir = parent.frame())
	if(sum(which(eval(parse(text = paste0("pred.wgt.", dataset))) != 0)) > 0){
		wgts <- append(wgts, paste0("pred.wgt.", dataset))
		return (wgts)
	}else {
		return (wgts)
	}
}

# --- SS-WEIGHTED-META-ANALYZE CELL TYPES
meta_cells <- function(dataset, result_wgts, indices){
	# PURPOSE: meta-analyze datasets together (use-case: data from different cell types but same people)
	# PREREQ: a hashmap (list) named "h" has to contain the sample sizes of the datasets you want to meta-analyze
	# dataset = prefix of datasets to combine
	# result_wgts = list to store the name of the meta-analyzed dataset
	# indices = list to store the indices to remove from the original wgts list before combining with result_wgts
	# RETURN: wgts_indices containing both result_wgts and indices
	lookup <- paste0("\\.", dataset)
	i <- grep(lookup, wgts)
	if (!identical(i, integer(0))){
	used <- wgts[i]
	indices <- append(indices, i)
	meta <- 0
	ss <- c()
	for (u in used){
                name <- strsplit(u, split = "[.]")[[1]][3]
                ss <- append(ss, h[[name]])
        }
	sums = sum(ss)
	for (u in used){
		name <- strsplit(u, split = "[.]")[[1]][3]
		meta <- meta + (eval(parse(text=u))*(h[[name]]/sums))
	}
	h[[dataset]] <<- mean(ss, na.rm = T)
	assign(paste0("pred.wgt.", dataset), meta, envir = parent.frame())
	result_wgts <- append(result_wgts, paste0("pred.wgt.", dataset))
	}
	wgts_indices <- list(result_wgts, indices)
	return (wgts_indices)
}

# --- FIND INDICES OF MAGEPRO "SPLIT"
magepro_split <- function(genos_file, h2, ge){
	# PURPOSE: find the indices of potentially predictive snps (used to split sumstats)
	# genos_file = name of plink files 
	# h2 = heritability estimate 
	# ge = gene expression as a matrix with 1 column 
	# RETURN: groupings: vector of strings, "nonzero" stores indices of important SNPs and "zero" stored indices of unimportant ones 
	geno = read_plink(genos_file,impute="avg")
	singleonly = weights.lasso( genos_file , h2 , snp=geno$bim[,2] )	
	if ( sum(is.na(singleonly)) == length(singleonly)) {
		singleonly = weights.marginal( geno$bed , ge , beta=T )
		#distribution of betas, take Z > 1.96 (nominally significant)
		z_scores <- scale(singleonly)
		top <- which(abs(z_scores) > 1.96)
		#if none of the snps are Z > 1.96	
		if (identical(integer(0), top)){
			threshold = max(round(length(singleonly)/4), 1) # split at the top 1/4 of eqtls
			top <- order(singleonly^2, decreasing = TRUE)[1:threshold]
		}
		singleonly[!top] <- 0
	}
	nonzero <<- which(singleonly != 0)
	zero <<- which(singleonly==0)
	groupings <- c()
	if (!identical(integer(0), nonzero)){
		groupings <- append(groupings, "nonzero")
	}
	if (!identical(integer(0), zero)){
		groupings <- append(groupings, "zero")
	}
	return (groupings)
}

# --- BIN SNPs BASED ON IMPACT SCORES
IMPACT_split <- function(IMPACT_file, BIM){
	# PURPOSE: find the indices of potentially predictive snps (used to split sumstats)
	# IMPACT_file = path to IMPACT file 
	# BIM = BIM matrix
	# RETURN: groupings: vector of strings, "nonzero" stores indices of important SNPs and "zero" stored indices of unimportant ones 
	IMPACT_scores <- as.data.frame(fread(IMPACT_file))
	IMPACT_scores_reordered <- IMPACT_scores[match(BIM[,2], IMPACT_scores$SNP),]
	IMPACT_threshold <<- quantile(IMPACT_scores_reordered[,5], 0.90) #CHANGE THIS TO TAKE TOP X%
	nonzero <<- which(IMPACT_scores_reordered[,5] >= IMPACT_threshold)
        zero <<- which(IMPACT_scores_reordered[,5] < IMPACT_threshold)
	groupings <- c()
        if (!identical(integer(0), nonzero)){
                groupings <- append(groupings, "nonzero")
        }
        #if (!identical(integer(0), zero)){
        #        groupings <- append(groupings, "zero")
        #}
        return (groupings)		
}


# --- I/O CHECKS

if (!is.na(opt$gene)){
	name <- strsplit(opt$gene, ".", fixed = TRUE)[[1]][1] # find gene name without version number
}else{
        cat( "ERROR: --gene flag is required\n" , sep='', file=stderr() )
	cleanup()
	q()
}

#models to use
model <- strsplit(opt$models, ",", fixed = TRUE)[[1]]
if ( sum(! model %in% c("SINGLE", "META", "MAGEPRO")) > 0 | length(model) > 3 ){
	cat( "ERROR: Please input valid models: SINGLE or META or MAGEPRO \n" , sep='', file=stderr() )
        cleanup()
        q()
}
if (opt$verbose >= 1) cat("### USING THE FOLLOWING MODELS:", opt$models, "\n")

order <- c("SINGLE","META","MAGEPRO") 
types <- model[match(order, model)]
types <- unique(types[!is.na(types)])

#list of datasets to be used
datasets <- list()

# sumstats to use
if ( ! is.na(opt$sumstats)){
	if ( is.na(opt$sumstats_dir) ){
		cat( "ERROR: --sumstats supplied, but not --sumstats_dir, cannot compute META and MAGEPRO models \n" , sep='', file=stderr() )
                cleanup()
                q()
	}
	sumstats <- strsplit(opt$sumstats, ",", fixed = TRUE)[[1]]
	for (s in sumstats){
		file_dir = paste0(opt$sumstats_dir, "/", s)
		if (dir.exists(file_dir)){
			file_name = paste0(file_dir, "/", name, ".txt")
			if (file.exists(file_name)){
				assign(paste0("file.", s), file_name, envir = parent.frame())
				datasets <- append(datasets, paste0("file.", s))
			}else{
				if ( opt$verbose == 2 ) {
					cat("skipping sumstat", s, "for this gene\n")
				}
			}
		}else{ #if there is no directory of eqtl sumstats in the --sumstats_dir path
			cat( "ERROR: make sure --sumstats_dir has a directory named ", s ," that contains gene-specific files of eqtl data from that sumstat dataset\n" , sep='', file=stderr() )
			cleanup()
			q()
		}
	}
}else{
	if ("META" %in% model | "MAGEPRO" %in% model){
		cat( "ERROR: --sumstats not supplied, cannot compute META and MAGEPRO models \n" , sep='', file=stderr() )
		cleanup()
                q()
	}
}

h <- list() #hashmap for sample sizes per dataset
if ( ("META" %in% model) | ( !is.na(opt$cell_meta) )){
	if (!is.na(opt$ss)){
		sample_sizes <- strsplit(opt$ss, ",", fixed = TRUE)[[1]]
		if (length(sumstats) != length(sample_sizes)){
			cat( "ERROR: --ss flag required an entry for every dataset\n" , sep='', file=stderr() )
                	cleanup()
                	q()
		}
		for (i in 1:length(sumstats)){
			h[[sumstats[i]]] <- as.numeric(sample_sizes[i])
		}
	}else{
		cat( "ERROR: Cannot perform sample-size weighted meta-analysis or cell-type meta-analysis without the --ss flag\n" , sep='', file=stderr() )
                cleanup()
                q()
	}
}

if ( opt$verbose == 2 ) {
	cat("Datasets available for this gene: \n")
	print(datasets)
	if (length(datasets) == 0){
		cat("WARNING: no datasets available for this gene, MAGEPRO and META will have NA weights \n")
	}
}

#if ( "MAGEPRO" %in% model ){
#        if (!is.na(opt$IMPACT)){
#                if ( !file.exists(opt$IMPACT) ){
#                        cat( "ERROR: --IMPACT file does not exist \n" , sep='', file=stderr() )
#                        cleanup()
#                        q()
#                }
#        }else{
#                cat( "ERROR: Cannot perform MAGEPRO without the --IMPACT flag\n" , sep='', file=stderr() )
#                cleanup()
#                q()
#        }
#}

files = paste(opt$bfile,c(".bed",".bim",".fam"),sep='')
if ( !is.na(opt$pheno) ) files = c(files,opt$pheno)
if ( !is.na(opt$covar) ) files = c(files,opt$covar)

for ( f in files ) {
	if ( !file.exists(f) ){
		cat( "ERROR: ", f , " input geno, pheno, or covar file does not exist\n" , sep='', file=stderr() )
		cleanup()
		q()
	}
}

if ( system( paste(opt$PATH_plink,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
	cat( "ERROR: plink executable not found, set with --PATH_plink\n" , sep='', file=stderr() )
	cleanup()
	q()
}

if ( !is.na(opt$hsq_set) && system( opt$PATH_gcta , ignore.stdout=T,ignore.stderr=T ) != 0 ){
	cat( "ERROR: gcta executable not found, set with --PATH_gcta\n" , sep='', file=stderr() )
	cleanup()
	q()
}


if ( opt$verbose == 2 ) cat("I/O checks complete \n")

# --- PREPARE PHENO, LOAD COVAR (FUSION FRAMEWORK)

#fam file read - the modified fam file version (added ge data)  
fam = read.table(paste(opt$bfile,".fam",sep=''),as.is=T, sep = " ")

# Make/fetch the phenotype file
if ( !is.na(opt$pheno) ) {
	pheno.file = opt$pheno
	pheno = read.table(pheno.file,as.is=T)
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(pheno[,1],pheno[,2]) )
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	m = m[m.keep]
	pheno = pheno[m,]
} else { #creating pheno file from fam files
	pheno.file = paste(opt$tmp,".pheno",sep='')
	pheno = fam[,c(1,2,6)]
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}

if ( opt$verbose == 2 ) cat("Made/Fetched phenotype \n")

# Load in the covariates if needed
var_cov = NA
if ( !is.na(opt$covar) ) { 
	covar = ( read.table(opt$covar,as.is=T,head=T) )
	if ( opt$verbose >= 1 ) cat( "### Loaded",ncol(covar)-2,"covariates\n")
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )	
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	pheno = pheno[m.keep,]
	m = m[m.keep]
	covar = covar[m,] #reordering covariates to match fam file
	reg = summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) )) 
        var_cov = reg$r.sq #save to output
	if ( opt$verbose == 2 ) cat( var_cov , "variance in phenotype explained by covariates\n" )
	pheno[,3] = scale(reg$resid) #regressing out covar to single out genetic effect
	raw.pheno.file = pheno.file
	pheno.file = paste(pheno.file,".resid",sep='')
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file) #newly scaled pheno file 
}


if ( opt$verbose == 2 ) cat("Covariates loaded and regressed out of phenotype \n")

geno.file = opt$tmp
# recode to the intersection of samples and ge
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$bfile," --pheno ",pheno.file," --keep ",pheno.file," --make-bed --out ",geno.file,sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)


# --- HERITABILITY ANALYSIS (FUSION FRAMEWORK)
if ( is.na(opt$hsq_set) ) {
	if ( opt$verbose >= 1 ) cat("### ESTIMATING HERITABILITY\n")

	# 1. generate GRM
	arg = paste( opt$PATH_plink," --allow-no-sex --bfile ",opt$tmp," --make-grm-bin --out ",opt$tmp,sep='')
	system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

	# 2. estimate heritability
	if ( !is.na(opt$covar) ) {
		arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",raw.pheno.file," --qcovar ",opt$covar," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
	} else {
		arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",pheno.file," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
	}
	system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

	# 3. evaluate LRT and V(G)/Vp
	#V(G)/Vp = proportion of genotypic to phenotypic variation
	#LRT - likelihood ratio test - compare goodness of fit 
	if ( !file.exists( paste(opt$tmp,".hsq",sep='') ) ) {
		if ( opt$verbose >= 1 ) cat(opt$tmp,"does not exist, GCTA could not converge, forcing h2 to fixed value\n",file=stderr()) 
		hsq = NA
		hsq.pv = NA
	}else{
		hsq.file = read.table(file=paste(opt$tmp,".hsq",sep=''),as.is=T,fill=T)
		hsq = as.numeric(unlist(hsq.file[hsq.file[,1] == "V(G)/Vp",2:3]))
		hsq.pv = as.numeric(unlist(hsq.file[hsq.file[,1] == "Pval",2]))
		if ( hsq[1] > 1 ){hsq = 1} #should not happen but gcta may overestimate at very low sample sizes
		if ( opt$verbose >= 1 ) cat("Heritability (se):",hsq,"LRT P-value:",hsq.pv,'\n')
		if ( opt$save_hsq ) cat( opt$out , hsq , hsq.pv , '\n' , file=paste(opt$out,".hsq",sep='') )

		# 4. stop if insufficient
		if (!is.na(hsq.pv)){
			if ( hsq.pv > opt$hsq_p ) { 
				cat(opt$tmp," : heritability ",hsq[1],"; LRT P-value ",hsq.pv," : skipping gene\n",sep='',file=stderr())
				cleanup()
				q()
			}
		}
	}
} else {
	if ( opt$verbose >= 1 ) cat("### SKIPPING HERITABILITY ESTIMATE\n")
	hsq = opt$hsq_set
	hsq.pv = NA
}

# --- LOAD AND PREP GENOTYPES (FUSION FRAMEWORK)
genos = read_plink(geno.file,impute="avg")
#mafs = apply(genos$bed,2,mean)/2  
#sds = apply(genos$bed,2,sd)  

# important : genotypes and phenotypes are standardized and scaled here:
genos$bed = scale(genos$bed)
pheno = genos$fam[,c(1,2,6)]
pheno[,3] = scale(pheno[,3])

# check if any genotypes are NA
nasnps = apply( is.na(genos$bed) , 2 , sum )
if ( sum(nasnps) != 0 ) {
	cat( "WARNING :",sum(nasnps != 0),"SNPs could not be scaled and were zeroed out, make sure all SNPs are polymorphic\n" , file=stderr())
	genos$bed[,nasnps != 0] = 0
}

# regress covariates out of the genotypes as well (this is more accurate but slower)
if ( !is.na(opt$covar) && opt$resid ) {
	if ( opt$verbose == 2 ) cat("Regressing covariates out of the genotypes\n")
	for ( i in 1:ncol(genos$bed) ) {
		genos$bed[,i] = summary(lm( genos$bed[,i] ~ as.matrix(covar[,3:ncol(covar)]) ))$resid
	}
	genos$bed = scale(genos$bed)
	genos$bed[is.nan(genos$bed)] <- 0 # when there is no variation at a snp, scale() sets genotypes of the snp to NaN. set those to 0
}

N.tot = nrow(genos$bed)
if ( opt$verbose >= 1 ) cat(nrow(pheno),"phenotyped samples, ",nrow(genos$bed),"genotyped samples, ",ncol(genos$bed)," markers\n")

# --- SETUP SUMSTATS FOR META and MAGEPRO 

lasso_h2 <- hsq[1]
if( (lasso_h2 < 0) | (is.na(lasso_h2)) ){
	if ( opt$verbose >= 1 ) cat("forcing lasso heritability to ", opt$lassohsq, " \n")
	lasso_h2 <- opt$lassohsq
		#change to lowest h2 that is considered significant 
			#hsq_afr = 0.064251
			#hsq_eur = 0.008909
}  #when gcta does not converge or yield wild estimates

if ("MAGEPRO" %in% model | "META" %in% model){

if ( opt$verbose >= 1){
	cat("### PROCESSING SUMSTATS \n")
}

wgts <- c() #sumstats weights before splitting (used for meta-analysis)

for (d in datasets){
	name <- strsplit(d, split="[.]")[[1]][2]
	wgts <- datasets_process(genos$bim, name, eval(parse(text = d)), wgts) # Run process dataset function on all datasets
}

if (!is.na(opt$cell_meta)){
if (opt$verbose == 2){
	cat("Meta-analyzing cell types together to create one eqtl profile per --cell_meta prefix \n")
}
new_wgts <- c()
index <- c()
cells_datasets <- strsplit(opt$cell_meta, split = ",")[[1]]
for (c in cells_datasets){
	wgts_index <- meta_cells(c, new_wgts, index)
	new_wgts <<- wgts_index[[1]]
	index <<- wgts_index[[2]]
}
if (length(index) > 0){
	new_wgts <- append(new_wgts, wgts[-index])
	wgts = new_wgts
}
}

ext <- length(wgts)

}

if ("META" %in% model){
# COMPUTE TOTAL SAMPLE SIZE OF SUMSTATS -> USED LATER IN META-ANALYSIS
total_ss_sumstats <- 0 
for (w in wgts){
	dataset <- strsplit(w, split="[.]")[[1]][3]
	total_ss_sumstats = total_ss_sumstats + h[[dataset]]
}
}

if ("MAGEPRO" %in% model){

#wgt2 <- c() #magepro weights

# SPLIT DATASETS 
#groups <- IMPACT_split( opt$IMPACT, genos$bim ) 

#for (w in wgts){
#	for (g in groups){
#		vec <- eval(parse(text = w))
#		vec[-eval(parse(text = g))] <- 0
#		assign(paste0(w, ".", g), vec, envir = parent.frame())
#		wgt2 <- append(wgt2, paste0(w, ".", g))
#	}	
#}	

#ext2 <- length(wgt2)

# --- no split
wgt2 <- wgts
ext2 <- ext
# ---

# PREPARE MATRIX FOR MAGEPRO

dataset_matrix <- matrix(nrow=ncol(genos$bed))
for (w in wgt2){
	dataset_matrix <- cbind(dataset_matrix, eval(parse(text = w)))
}
# remove first column of NA created during matrix initialization
dataset_matrix <- dataset_matrix[,-1]

meta_datasets <- c()

for (i in 1:nrow(dataset_matrix)){
	zeros <- which(dataset_matrix[i,] == 0)
	weighted_beta <- sum( (dataset_matrix[i,]*unlist(h)) )
	if (identical(zeros, integer(0))){
		ss_contributing <- sum(unlist(h))
	}else{
		ss_contributing <- sum(unlist(h[zeros]))
	}
	new_beta <- weighted_beta/ss_contributing
	if (new_beta == 0){
		new_beta = 1
	}
	meta_datasets <- append(meta_datasets, new_beta)
}

magepro_matrix <- diag(meta_datasets)

}


# --- CROSSVALIDATION ANALYSES
set.seed(1)

#default crossval = 5 fold split
if ( opt$crossval <= 1 ) { 
if ( opt$verbose >= 1 ) cat("### SKIPPING CROSS-VALIDATION\n")
avg_training_r2_magepro <- NA
avg_training_r2_single <- NA
avg_training_r2_meta <- NA
cv.performance <- NA
} else {
if ( opt$verbose >= 1 ) cat("### PERFORMING",opt$crossval,"FOLD CROSS-VALIDATION\n")
cv.all = pheno 
n = nrow(cv.all)
cv.sample = sample(n) #sample randomly 
cv.all = cv.all[ cv.sample , ]
folds = cut(seq(1,n),breaks=opt$crossval,labels=FALSE) #5 fold split - split into 5 groups 
cv.calls = matrix(NA,nrow=n,ncol=length(model)) 

# --- for checking cor() of weights in CV 
wgt.cv = matrix(0, nrow=nrow(genos$bim), ncol=opt$crossval)
# ---
# --- keep track of SINGLE_model
SINGLE_top1 <- 0
# ---

r2_training_magepro <- c()
r2_training_meta <- c()
r2_training_single <- c()

if ("META" %in% model){
training_ss <- N.tot * ((opt$crossval - 1)/opt$crossval)
total_ss_cv <- total_ss_sumstats + training_ss # total sample size in cross validation (used as denominator in ss-weighted meta-analysis)
}

# --- Cross-Validation
for ( cv in 1:opt$crossval ) { 		
	colcount = 1
	if ( opt$verbose >= 1 ) cat("- Crossval fold",cv,"\n")
	indx = which(folds==cv,arr.ind=TRUE)
	cv.train = cv.all[-indx,] #training set is the other 4 groups 
	intercept = mean( cv.train[,3] ) 
	cv.train[,3] = scale(cv.train[,3]) 
	cv.file = paste(opt$tmp,".cv",sep='')
	write.table( cv.train , quote=F , row.names=F , col.names=F , file=paste(cv.file,".keep",sep=''))	
	arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$tmp," --keep ",cv.file,".keep --out ",cv.file," --make-bed",sep='')
	system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
	
	# SINGLE ANCESTRY------------------------------------------------------------------
	# lasso_h2 defined when we split datasets
	pred.wgt = weights.lasso( cv.file , lasso_h2 , snp=genos$bim[,2] )	
	if ( sum(is.na(pred.wgt)) == length(pred.wgt)) {
		if ( opt$verbose >= 1 ) cat("LASSO pushed all weights to 0, using top1 as backup \n")
		SINGLE_top1 <- SINGLE_top1 + 1
		pred.wgt = weights.marginal( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3,drop=F]) , beta=T )
		pred.wgt[ - which.max( pred.wgt^2 ) ] = 0
	}
	if (length(pred.wgt) == 1){
		pred.wgt <- t(pred.wgt) # 1 snp in the cis window -> transpose for "matrix" mult
	}
	
	if ("SINGLE" %in% model){
	
	cv.calls[ indx , colcount ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt
	pred_train_single = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], ] %*% pred.wgt)))
	r2_training_single = append(r2_training_single, pred_train_single$adj.r.sq)
	
	colcount = colcount + 1

	}
	#--------------------------------------------------------------------------

	# SS-WEIGHTED META-ANALYSIS-------------------------------------------------------------
	if ("META" %in% model){
	if (ext > 0){
	pred.wgt.meta <- pred.wgt * (training_ss/total_ss_cv)
	for (w in wgts){
		dataset <- strsplit(w, split="[.]")[[1]][3]
		size <- h[[dataset]]
		pred.wgt.meta = pred.wgt.meta + (eval(parse(text=w)) * (size/total_ss_cv))
	}
	if (length(pred.wgt.meta) == 1){
		pred.wgt.meta <- t(pred.wgt.meta) 
	}
	cv.calls[ indx , colcount ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt.meta 
	pred_train_meta = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], ] %*% pred.wgt.meta))) 
        r2_training_meta = append(r2_training_meta, pred_train_meta$adj.r.sq)
	}else{
	cv.calls[ indx , colcount ] = NA
	r2_training_meta = append(r2_training_meta, NA)
	}

	colcount = colcount + 1
	
	}
	#-----------------------------------------------------------------------------

	# MAGEPRO----------------------------------------------------------------------
	
	if ("MAGEPRO" %in% model){
	
	if (ext2 > 0){	

	eq <- genos$bed[ cv.sample[ -indx ] , ] %*% magepro_matrix

	#2. run lasso regression to find optimal coefficients
	y <- cv.glmnet(x = eq , y = cv.train[,3], alpha = 1, nfold = 5, intercept = T, standardize = T)

	cf = coef(y, s = "lambda.min")[2:(ncol(eq)+1)]
	pred.wgt.magepro <- magepro_matrix %*% cf
	if(length(pred.wgt.magepro) == 1){
		pred.wgt.magepro <- t(pred.wgt.magepro)
	}	

	cv.calls[ indx , colcount ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt.magepro

	# --- for checking cor() of weights in CV
	wgt.cv[,cv] = pred.wgt.magepro
	# --- 

	#store the r2 on training set
	pred_train = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], ] %*% pred.wgt.magepro))) #r^2 between predicted and actual 
	r2_training_magepro = append(r2_training_magepro, pred_train$adj.r.sq)	
	}else{
	cv.calls[ indx , colcount ] = NA
	r2_training_magepro = append(r2_training_magepro, NA)
	}
	
	}
	#-------------------------------------------------------------------------------

}

if (opt$verbose >= 1) cat("### COLLECTING CROSS VALIDATION RESULTS \n")

#uses 80% of data to predict into 20%, then uses another 80% to predict into the next 20%, and so on. 
#after everyone has a prediction, compute rsq and p value. 

#compute rsq + P-value for each model
cv.performance = matrix(NA,nrow=2,ncol=length(model))  
rownames(cv.performance) = c("rsq","pval") 
colnames(cv.performance) = types 

for ( mod in 1:ncol(cv.calls) ) { 
	if ( !is.na(sd(cv.calls[,mod])) && sd(cv.calls[,mod]) != 0 ) {   
		reg = summary(lm( cv.all[,3] ~ cv.calls[,mod] )) #r^2 between predicted and actual 
		cv.performance[ 1, mod ] = reg$adj.r.sq
		cv.performance[ 2, mod ] = reg$coef[2,4]
	} else {
		cv.performance[ 1, mod ] = NA
		cv.performance[ 2, mod ] = NA
	}
}
if ( opt$verbose >= 1 ) write.table(cv.performance,quote=F,sep='\t')


#take average of r2 on training set
avg_training_r2_single <- avg_training_r2_meta <- avg_training_r2_magepro <- NA
if ("SINGLE" %in% model) avg_training_r2_single <- mean(r2_training_single)
if ("META" %in% model) avg_training_r2_meta <- mean(r2_training_meta)
if ("MAGEPRO" %in% model) avg_training_r2_magepro <- mean(r2_training_magepro)

}

# ---- Full Analysis 
if (opt$verbose >= 1) cat("### COMPUTING FULL GENE MODELS \n")

wgt.matrix = matrix(0, nrow=nrow(genos$bim), ncol=length(model)) 
colnames(wgt.matrix) = types
rownames(wgt.matrix) = genos$bim[,2]

colcount = 1

# --- SINGLE ANCESTRY
pred.wgtfull = weights.lasso( geno.file , lasso_h2 , snp=genos$bim[,2] )	
if ( sum(is.na(pred.wgtfull)) == length(pred.wgtfull)) {
	pred.wgtfull = weights.marginal( genos$bed , as.matrix(pheno[,3]) , beta=T )
	pred.wgtfull[ - which.max( pred.wgtfull^2 ) ] = 0
}
if(length(pred.wgtfull) == 1){
	pred.wgtfull <- t(pred.wgtfull)
}
if ("SINGLE" %in% model){
wgt.matrix[, colcount] = pred.wgtfull
colcount = colcount + 1
}
# --- SS-WEIGHTED META-ANALYSIS
if ("META" %in% model){
if (ext > 0){
total_ss_full <- total_ss_sumstats + N.tot
pred.wgt.metafull <- pred.wgtfull * (N.tot/total_ss_full)
for (w in wgts){
	dataset <- strsplit(w, split="[.]")[[1]][3]
	size <- h[[dataset]]
	pred.wgt.metafull = pred.wgt.metafull + (eval(parse(text=w)) * (size/total_ss_full))
}
if(length(pred.wgt.metafull) == 1){
        pred.wgt.metafull <- t(pred.wgt.metafull)
}
wgt.matrix[, colcount] = pred.wgt.metafull
}
else{
wgt.matrix[, colcount] = NA
}
colcount = colcount + 1
}
# --- MAGEPRO
cf_total = NA
if ("MAGEPRO" %in% model){
if (ext2 > 0){
eqfull <- genos$bed %*% magepro_matrix	
#run ridge regression to find optimal coefficients and compute multipop weight
yfull <- cv.glmnet(x = eqfull , y = pheno[,3], alpha = 1, nfold = 5, intercept = T, standardize = T)
cf_total = coef(yfull, s = "lambda.min")[2:(ncol(eqfull)+1)]	
pred.wgt.mageprofull <- magepro_matrix %*% cf_total
if(length(pred.wgt.mageprofull) == 1){
	pred.wgt.mageprofull <- t(pred.wgt.mageprofull)
}
wgt.matrix[, colcount] = pred.wgt.mageprofull
}else{
wgt.matrix[, colcount] = NA
}
}

#--- SAVE RESULTS
snps = genos$bim
if ("MAGEPRO" %in% model){
wgtmagepro <- wgt2
coefs <- NA # this version of the model does not have coefs per dataset
}else{
wgtmagepro = NA
coefs = NA
}

# --- for checking cor() of weights in CV
cors_weights <- c()
for (i in 1:(opt$crossval-1)){
	cors_weights <- append(cors_weights, cor(wgt.cv[,i], wgt.cv[,(i+1)]) )
}
avg_cor <- mean(cors_weights)
# ---

save( wgt.matrix, snps, cv.performance, hsq, hsq.pv, N.tot , wgtmagepro, coefs, avg_training_r2_single, avg_training_r2_meta, avg_training_r2_magepro, var_cov, avg_cor, SINGLE_top1, file = paste( opt$out , ".wgt.RDat" , sep='' ) )

# --- CLEAN-UP
if ( opt$verbose >= 1 ) cat("### CLEANING UP\n")
cleanup()
