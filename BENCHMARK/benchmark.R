# load packages
suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('glmnet'))
suppressMessages(library('methods'))
suppressMessages(library('data.table'))
suppressMessages(library('hash'))
suppressMessages(library('tools'))
suppressMessages(library('metafor'))
suppressMessages(library('Rmpfr'))
suppressMessages(library('dplyr'))

option_list = list(
  make_option("--gene", action="store", default=NA, type='character',
              help="ENSG ID, for ex ENSG00000013503.9"),
  make_option("--tissue", action="store", default=NA, type='character',
              help="Tissue name, for ex Whole_Blood"),
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files [required]"),
  make_option("--pheno", action="store", default=NA, type='character',
              help="Path to molecular phenotype file (PLINK format) [optional, taken from bfile otherwise]"),
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gcta", action="store", default="gcta_nr_robust", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gemma", action="store", default="gemma", type='character',
              help="Path to plink executable [%default]"),
  make_option("--covar", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) [optional]"),
  make_option("--resid", action="store_true", default=FALSE,
              help="Also regress the covariates out of the genotypes [default: %default]"),              
  make_option("--hsq_p", action="store", default=0.01, type='double',
              help="Minimum heritability p-value for which to compute weights [default: %default]"),
  make_option("--hsq_set", action="store", default=NA, type='double',
              help="Skip heritability estimation and set hsq estimate to this value [optional]"),
  make_option("--crossval", action="store", default=5, type='double',
              help="How many folds of cross-validation, 0 to skip [default: %default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"),
  make_option("--rn", action="store_true", default=FALSE,
              help="Rank-normalize the phenotype after all QC: [default: %default]"),
  make_option("--save_hsq", action="store_true", default=FALSE,
              help="Save heritability results even if weights are not computed [default: %default]"),			  
  make_option("--gemmaout", action="store", default=NA,type='character',
              help="Gemma must output to output/ directory that is child directory of directory from which fusion is run  [default: %default]"),
  make_option("--models", action="store", default="blup,lasso,top1,enet", type='character',
              help="Comma-separated list of prediction models [default: %default]\n
					Available models:\n
					top1:\tTop eQTL (standard marginal eQTL Z-scores always computed and stored)\n
					blup:\t Best Unbiased Linear Predictor (dual of ridge regression)\n
					bslmm:\t Bayesian Sparse Linear Model (spike/slab MCMC)\n
					lasso:\t LASSO regression (with heritability used as lambda)\n
					enet:\t Elastic-net regression (with mixing parameter of 0.5)\n"),			  
  make_option("--datasets", action="store", default="", type='character', 
	      help="Comma-separated list of external datasets to include"),
  make_option("--ldref", action = "store", default="", type='character',
	      help="Path to ld reference file for pruning and thresholding")

)


# parse command-line arguments
opt = parse_args(OptionParser(option_list=option_list))
print(opt)

# find gene name without version number
name <- strsplit(opt$gene, ".", fixed = TRUE)[[1]][1]
print(name)

# datasets to use
datas <- strsplit(opt$datasets, ",", fixed = TRUE)[[1]]

#PRS-CSx working dir 

PRS_CSx_working_dir = paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/benchmark/working_data/", name, "/")
system(paste0("mkdir ", PRS_CSx_working_dir))

# assign all external dataset file names
file.eur <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/GTEx_EUR/genes/", opt$gene, ".txt")
ota_cells <- list("CD16p_Mono","CL_Mono","LDG","Mem_CD4","Mem_CD8","NK","Naive_B","Naive_CD4","Naive_CD8","Neu","Plasmablast","mDC","pDC")
for (c in ota_cells){
        assign(paste0("file.ota.", c), paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/OTA_nominal/genes/", c, "/", name, ".txt"))
}
file.his <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/MESA_HIS/genes/",name,".txt")
file.genoa <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/genoa/genes_nominal/genes/",name,".txt")
file.mesa <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/MESA/genes/",name,".txt")
file.eqtlgen <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/eQTLGEN/genes/", name, ".txt")


# create a hashmap to store sample size of each datset (used later for sample-size weighted meta-analysis)
h <- hash()
h[["eur"]] <- 574
h[["ota"]] <- 416
h[["genoa"]] <- 1031
h[["mesa"]] <- 233
h[["his"]] <- 352
h[["eqtlgen"]] <- 31684 

pops <- hash()
pops[["eur"]] <- "EUR"
pops[["ota"]] <- "EAS"
pops[["genoa"]] <- "AFR"
pops[["mesa"]] <- "AFR"
pops[["his"]] <- "AMR"
pops[["eqtlgen"]] <- "EUR" 


# all available external datasets to a list
datasets <- list("file.eur","file.ota.CD16p_Mono", "file.ota.CL_Mono", "file.ota.LDG", "file.ota.Mem_CD4", "file.ota.Mem_CD8", "file.ota.NK", "file.ota.Naive_B", "file.ota.Naive_CD4", "file.ota.Naive_CD8", "file.ota.Neu", "file.ota.Plasmablast", "file.ota.mDC", "file.ota.pDC", "file.his", "file.genoa", "file.mesa", "file.eqtlgen")

# remove datasets not being used
datasets <- datasets[grepl(paste(datas, collapse = "|"), datasets)]
print(datasets)

datasets <- lapply(datasets, function(x) eval(parse(text=x)))

# remove datasets that are not available for that gene
for (d in datasets) {
	if (!file.exists(d)){
		datasets <- datasets[datasets != d]
	}
}
print(datasets)


if(length(datasets) >= 0) {
models = unique(unlist(strsplit(opt$models,','))) # a list of prediction models to be used
M = length(models) # number of models 

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}

# --- PREDICTION MODELS

# BSLMM
weights.bslmm = function( input , bv_type , snp , out=NA ) {
	if ( is.na(out) ) out = paste(input,".BSLMM",sep='')
	print(out)
	arg = paste( opt$PATH_gemma , " -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000 -bfile " , input , " -bslmm " , bv_type , " -o " , out , sep='' )
	system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
	eff = read.table( paste(getwd(), "/output/", out,".param.txt",sep=''),head=T,as.is=T)
	eff.wgt = rep(NA,length(snp))
	m = match( snp , eff$rs )
	m.keep = !is.na(m)
	m = m[m.keep]
	eff.wgt[m.keep] = (eff$alpha + eff$beta * eff$gamma)[m]
	return( eff.wgt )
}

# PLINK: LASSO
weights.lasso = function( input , hsq , snp , out=NA ) {
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
	if ( beta ) eff.wgt = t( genos ) %*% (pheno) / ( nrow(pheno) - 1)
	else eff.wgt = t( genos ) %*% (pheno) / sqrt( nrow(pheno) - 1 )
	return( eff.wgt )
}

# Elastic Net
weights.enet = function( genos , pheno , alpha=0.5 ) {
	eff.wgt = matrix( 0 , ncol=1 , nrow=ncol(genos) )
	# remove monomorphics
	sds = apply( genos  , 2 , sd )
	keep = sds != 0 & !is.na(sds)
	enet = cv.glmnet( x=genos[,keep] , y=pheno , alpha=alpha , nfold=5 , intercept=T , standardize=F )
	eff.wgt[ keep ] = coef( enet , s = "lambda.min")[2:(sum(keep)+1)]
	return( eff.wgt )
}

print("predictive models made")

# --- CLEANUP
cleanup = function() {
	if ( ! opt$noclean ) {
		arg = paste("rm -f " , opt$tmp , "*", sep='')
		system(arg)
                arg = paste("rm -f " , opt$gemmaout , "*", sep='')  
		system(arg)
		arg = paste("rm -f " , PRS_CSx_working_dir , "*", sep='')  
		system(arg)
	}
}

# Perform i/o checks here:
files = paste(opt$bfile,c(".bed",".bim",".fam"),sep='')
if ( !is.na(opt$pheno) ) files = c(files,opt$pheno)
if ( !is.na(opt$covar) ) files = c(files,opt$covar)

for ( f in files ) {
	if ( !file.exists(f) ){
		cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
		cleanup()
		q()
	}
}

if ( system( paste(opt$PATH_plink,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
	cat( "ERROR: plink could not be executed, set with --PATH_plink\n" , sep='', file=stderr() )
	cleanup()
	q()
}

if ( !is.na(opt$hsq_set) && system( opt$PATH_gcta , ignore.stdout=T,ignore.stderr=T ) != 0 ){
	cat( "ERROR: gcta could not be executed, set with --PATH_gcta\n" , sep='', file=stderr() )
	cleanup()
	q()
}

if ( sum(models=="bslmm" | models=="blup") != 0 && system( paste(opt$PATH_gemma,"-h") , ignore.stdout=T,ignore.stderr=T ) != 0 ){
	cat( "ERROR: gemma could not be executed, set with --PATH_gemma or remove 'bslmm' and 'blup' from models\n" , sep='', file=stderr() )
	cleanup()
	q()
}

print("io checks complete")

# ---

#fam file read - the modified fam file version (added ge data)  
fam = read.table(paste(opt$bfile,".fam",sep=''),as.is=T)
print("fam file read")

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

print("made/fetched phenotype file")

if ( opt$rn ) {
	library('GenABEL')
	library(preprocessCore)
	pheno[,3] = rntransform( pheno[,3] )
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}

# Load in the covariates if needed
if ( !is.na(opt$covar) ) { #reorder covar file to match pheno file
	covar = ( read.table(opt$covar,as.is=T,head=T) )
	if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )  
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	pheno = pheno[m.keep,]
	m = m[m.keep]
	covar = covar[m,] #reordering covariates to match fam file
        #if there are few people, we may have more covariates, so limit to 5 peer factors maybe 
	reg = summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))  
	if ( opt$verbose >= 1 ) cat( reg$r.sq , "variance in phenotype explained by covariates\n" )
	pheno[,3] = scale(reg$resid) #regressing out covar to single out genetic effect
	raw.pheno.file = pheno.file
	pheno.file = paste(pheno.file,".resid",sep='')
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file) #newly scaled pheno file 
}

print("covariate loaded and regressed out from pheno")

geno.file = opt$tmp
# recode to the intersection of samples and new phenotype
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$bfile," --pheno ",pheno.file," --keep ",pheno.file," --make-bed --out ",geno.file,sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

print("recode samples and new phenotype using plink")

# --- HERITABILITY ANALYSIS
if ( is.na(opt$hsq_set) ) {
if ( opt$verbose >= 1 ) cat("### Estimating heritability\n")

# 1. generate GRM
arg = paste( opt$PATH_plink," --allow-no-sex --bfile ",opt$tmp," --make-grm-bin --out ",opt$tmp,sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
print("grm generated")

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
	cat(opt$tmp,"does not exist, likely GCTA could not converge, forcing h2 = 0.064251\n",file=stderr())
	hsq_afr = 0.064251
	hsq_afr.pv = NA
	#if heritability estimate does not converge, push through with a preset heritability value = 0.064251 - smallest h2 that has p < 0.05
}else{
	hsq.file = read.table(file=paste(opt$tmp,".hsq",sep=''),as.is=T,fill=T)
	hsq_afr = as.numeric(unlist(hsq.file[hsq.file[,1] == "V(G)/Vp",2:3]))
	hsq_afr.pv = as.numeric(unlist(hsq.file[hsq.file[,1] == "Pval",2]))
if ( opt$verbose >= 1 ) cat("Heritability (se):",hsq_afr,"LRT P-value:",hsq_afr.pv,'\n')
if ( opt$save_hsq ) cat( opt$out , hsq_afr , hsq_afr.pv , '\n' , file=paste(opt$out,".hsq",sep='') )

# 4. stop if insufficient
# hsq_p set to 1, compute weights for all genes, regardless of hsq p value
if (!is.na(hsq_afr.pv)){
	if ( hsq_afr.pv > opt$hsq_p ) { 
		cat(opt$tmp," : heritability ",hsq_afr[1],"; LRT P-value ",hsq_afr.pv," : skipping gene\n",sep='',file=stderr())
		cleanup()
		q()
	}
}
}
}else {
if ( opt$verbose >= 1 ) cat("### Skipping heritability estimate\n")
hsq_afr = opt$hsq_set
hsq_afr.pv = NA
}

# read in genotypes
genos = read_plink(geno.file,impute="avg")
mafs = apply(genos$bed,2,mean)/2  
sds = apply(genos$bed,2,sd)  
# important : genotypes are standardized and scaled here:
genos$bed = scale(genos$bed)
pheno = genos$fam[,c(1,2,6)]
pheno[,3] = scale(pheno[,3])

print("read genotype and scaled")

# check if any genotypes are NA
nasnps = apply( is.na(genos$bed) , 2 , sum )
if ( sum(nasnps) != 0 ) {
	cat( "WARNING :",sum(nasnps != 0),"SNPs could not be scaled and were zeroed out, make sure all SNPs are polymorphic\n" , file=stderr())
	genos$bed[,nasnps != 0] = 0
}

# regress covariates out of the genotypes as well (this is more accurate but slower)
if ( !is.na(opt$covar) && opt$resid ) {
	if ( opt$verbose >= 1 ) cat("regressing covariates out of the genotypes\n")
	for ( i in 1:ncol(genos$bed) ) {
		genos$bed[,i] = summary(lm( genos$bed[,i] ~ as.matrix(covar[,3:ncol(covar)]) ))$resid
	}
	genos$bed = scale(genos$bed)
}


N.tot = nrow(genos$bed)
if ( opt$verbose >= 1 ) cat(nrow(pheno),"phenotyped samples, ",nrow(genos$bed),"genotyped samples, ",ncol(genos$bed)," markers\n")

# --- CROSSVALIDATION ANALYSES
set.seed(1)

# 2 fold split for validation and testing for PRS-CSx
if ( opt$crossval <= 1 ) { 
if ( opt$verbose >= 1 ) cat("### skipping cross-validation\n")
} else {
if ( opt$verbose >= 1 ) cat("### performing",opt$crossval,"fold cross-validation\n")
cv.all = pheno #note, cv.all is the gene expression phenotype for all 80 afr individuals for whole blood for ex 
n = nrow(cv.all)
cv.sample = sample(n) #sample randomly 
cv.all = cv.all[ cv.sample , ]
folds = cut(seq(1,n),breaks=opt$crossval,labels=FALSE) #5 fold split - split into 5 groups 
cv.calls = matrix(NA,nrow=n,ncol=2) 
}


#list of all wgt vectors 
wgts <- list()
input <- list()
ss <- list()
pp <- list()

#PRS-CSx ONLY USES SUMMARY STATS, SO EXCLUDE EUR GENE MODEL

#process function to handle ref/alt flips for external datasets---------------------------------
#also prepare the input files for PRS-CSx
datasets_process <- function(dataset, file, cell, snp, A1, A0, B, P){
	#dataset = name of dataset
	#file = file path to data
	#cell = cell type if it exists (NA if data is not split by cell type)
	#wgts = current list of wgts 
	#snp = column number of rsid 
	#A1 = column number of relevant (alt) allele 
	#A0 = column number of ref allele
	#B = column number of effect size Beta
	if (is.na(cell)){
		cell <- ""
	}
	table <- fread(file, select = c(snp, A1, A0, B, P)) 
	for (k in 1:nrow(table)){
		match=which(genos$bim$V2 == table[k, 1])
		if (identical(match, integer(0))){
			print(paste0(dataset, " ", cell, ":", "this snp not in GTEx, skipping ref/alt check iteration"))
		}else{
		if (!is.na(table[k, 2]) && !is.na(table[k, 3])){
			if (as.character(genos$bim$V5[match]) != as.character(table[k, 2]) || as.character(genos$bim$V6[match]) != as.character(table[k, 3])) {  # check if A1 column is effect allele
				if (as.character(genos$bim$V5[match]) == as.character(table[k, 3]) && as.character(genos$bim$V6[match]) == as.character(table[k, 2])){
					table[k, 4] <- table[k, 4] * -1
					table[k, 2] <- genos$bim$V5[match]
					table[k, 3] <- genos$bim$V6[match]
				}else{
					table[k, 4] <- 0
				}
			}
		} else {
			table[k, 4] <- 0
		}
		}
	}
	df <- data.frame(table)
	colnames(df) <- c("SNP", "A1", "A2", "BETA", "P")
	f_name <- paste0(PRS_CSx_working_dir, dataset, cell, ".txt")
	input <<- append(input, f_name)
	ss <<- append(ss, h[[dataset]])
	pp <<- append(pp, pops[[dataset]])
	write.table(df, file = f_name, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}


#function to process ota weights before meta-analyzing 

flip_signs_ota <- function(dataset, file, snp, A1, A0, B, P){
	table <- fread(file, select = c(snp, A1, A0, B, P)) 
	Z <- qnorm( as.numeric( as.bigq(1) - as.bigq(table[[5]]/2) ) )
	Z <- Z * sign(table[[4]])
	table <- cbind(table, Z)
	colnames(table) <- c("snp", "A1", "A0", "B", "P", "Z")
	for (k in 1:nrow(table)){
		match=which(genos$bim$V2 == table[k, 1])
		if (!identical(match, integer(0))){
			if (!is.na(table[k, 2]) && !is.na(table[k, 3])){
				if (as.character(genos$bim$V5[match]) != as.character(table[k, 2]) || as.character(genos$bim$V6[match]) != as.character(table[k, 3])) {  # check if A1 column is effect allele
					if (as.character(genos$bim$V5[match]) == as.character(table[k, 3]) && as.character(genos$bim$V6[match]) == as.character(table[k, 2])){
						table[k, 4] <- table[k, 4] * -1
						temp = table[k,2]
						table[k,2] <- table[k,3]
						table[k,3] <- temp
					}else{
						table[k, 4] <- NA
					}
				}
			} else {
				table[k, 4] <- NA
			}
		}
	}
	SE <- table$B/table$Z
	w <- is.na(SE)
	if (sum(w) > 0){
		SE[w] = NA
	}
	results_return <- list(table$snp, table$A1, table$A0, table$B, SE)
	return (results_return)
}

#function to add missing snps to dataframe - when meta-analyzing ota, snps have to be in common 
add_missing_snps <- function(df, unique) {
  missing_snps <- setdiff(unique, df$snp)  # Find missing SNPs
  # check if missing is undefined
  if (length(missing_snps) == 0){
	  return (df)
  }else{
  	  missing_rows <- data.frame(snp = missing_snps, A1 = 'X', A0 = 'X', weights = NA, se = NA)  # Create missing rows - separate allele info used so 'X' is just a filler 
  	  updated_df <- rbind(df, missing_rows)  # Combine original and missing rows
  	  return(updated_df)
  }
}


#function to take posterior weights and prepare weight vectors 

weights_process <- function(dataset, file, wgts){
	if (file.size(file) != 0L){
		table <- fread(file)
		m <- match(genos$bim[,2], table[[2]])
		table <- table[m,]
		w <- which(is.na(table[[2]]))
		if(length(w) > 0){
			table[w,6] <- 0
		}
		assign(paste0("pred.wgt.", dataset), table[[6]], envir = parent.frame())
		if(sum(which(eval(parse(text = paste0("pred.wgt.", dataset))) != 0)) > 0){
			wgts <- append(wgts, paste0("pred.wgt.", dataset))
			return (wgts)
		}else {
			return (wgts)
		}
	}else{
		return (wgts)
	}
}



#EUR
if (file.eur %in% datasets){
	datasets_process("eur", file.eur, NA, 11, 13, 14, 8, 7)
}


#OTA  - have to meta-analyze cell types before shrinking 
ota_use = F
ota_dfs <- c() # cell-type specific dataframe
ota_cells <- c() # cell types

all_ota <- c(file.ota.CD16p_Mono, file.ota.CL_Mono, file.ota.LDG, file.ota.Mem_CD4, file.ota.Mem_CD8, file.ota.NK, file.ota.Naive_B, file.ota.Naive_CD4, file.ota.Naive_CD8, file.ota.Neu, file.ota.Plasmablast, file.ota.mDC, file.ota.pDC)

for (o in all_ota){
	if(o %in% datasets){
		cell <- strsplit(o, split = "/")[[1]][10]
                ota_cells <- append(ota_cells, cell)
                results_list <- flip_signs_ota("ota", o, 6, 8, 7, 13, 12)
                df <- data.frame(results_list[1], results_list[2], results_list[3], results_list[4], results_list[5])
                colnames(df) <- c("snp", "A1", "A0", "weights", "se")
                ota_dfs <- append(ota_dfs, list(df))
                ota_use = T
	}
}

if (ota_use == T){

# iterate through all dataframes, add to a list and find all unique snps
all_snps <- c()
all_A1 <- c()
all_A0 <- c()
for (c in 1:length(ota_cells)){
        df <- ota_dfs[[c]]
        all_snps <- append(all_snps, df$snp)
        all_A1 <- append(all_A1, df$A1)
        all_A0 <- append(all_A0, df$A0)
}
snp_info <- data.frame(all_snps, all_A1, all_A0)
snp_info <- snp_info[!duplicated(snp_info$all_snps), ]
rownames(snp_info) <- NULL
unique_snps <- snp_info$all_snps
snp_info <- snp_info[order(snp_info$all_snps), ]

# modify each dataframe individually to add missing snps and change order
sorted_order <- snp_info$all_snps
for (c in 1:length(ota_cells)){
        modified <- add_missing_snps(ota_dfs[[c]], unique_snps)
        modified <- modified[order(match(modified$snp, sorted_order)), ]
        rownames(modified) <- NULL
        ota_dfs[[c]] <- modified
}

#do meta-analysis - extract B and P
#iterate per snp
B_ota <- c()
P_ota <- c()
for (s in 1:length(unique_snps)){
	effect_sizes <- c()
        standard_errors <- c()
        sample_sizes <- c()

        for (c in 1:length(ota_cells)){
                effect_sizes <- append(effect_sizes, ota_dfs[[c]][[4]][s])
                standard_errors <- append(standard_errors, ota_dfs[[c]][[5]][s])
        }

	n <- sum(which(!is.na(effect_sizes)))

	if (n == 0){
		B_ota <- append(B_ota, NA)
                P_ota <- append(P_ota, NA)
	}else{
		# calculate B 
		w <- which(!is.na(effect_sizes))
		B_meta <- mean(effect_sizes[w]) #ota sample sizes are some for cell types
		B_ota <- append(B_ota, B_meta)
		# calculate P from Var of linear combination of random variables 
		#https://stats.stackexchange.com/questions/495215/standard-error-standard-deviation-and-variance-confusion
		variances = standard_errors[w]**2
		#https://stats.libretexts.org/Bookshelves/Introductory_Statistics/OpenIntro_Statistics_(Diez_et_al)./02%3A_Probability/2.05%3A_Random_Variables#:~:text=The%20standard%20deviation%20of%20the,square%20root%20of%20the%20variance.
		var_meta = sum(variances)/(n**2)
		se_meta = sqrt(var_meta)
		Z = B_meta/se_meta
		P_ota <- append(P_ota, 2*pnorm(q = abs(Z), lower.tail=FALSE))
	}
}

ota_df <- data.frame(snp_info[,1:3], B_ota, P_ota)
colnames(ota_df) <- c("SNP", "A1", "A2", "BETA", "P")
f_name <- paste0(PRS_CSx_working_dir, "ota.txt")
input <<- append(input, f_name)
ss <<- append(ss, h[["ota"]])
pp <<- append(pp, pops[["ota"]])
ota_df <- na.omit(ota_df)
write.table(ota_df, file = f_name, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

}

#GENOA
if (file.genoa %in% datasets){
	datasets_process("genoa", file.genoa, NA, 3, 5, 6, 8, 10)
}


#MESA
if (file.mesa %in% datasets){
	datasets_process("mesa", file.mesa, NA, 1, 14, 13, 6, 4)
}

#MESA_HIS
if (file.his %in% datasets){
	datasets_process("his", file.his, NA, 1, 14, 13, 6, 4)
}

#eqtlgen
if (file.eqtlgen %in% datasets){
	datasets_process("eqtlgen", file.eqtlgen, NA, 1, 4, 5, 12, 14)
}


#---------------PRS-CSx shrinkage prior

#PRS-CSx needs atleast one summary stat to run 

if (length(datasets) < 1){
	stop("PRS-CSx needs atleast 1 summary stats to run, skipping this gene")
}

#save bim file to PRS-CSx working directory 
write.table(genos$bim, file = paste0(PRS_CSx_working_dir, "snps.bim"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#run PRS-CSx 
ss <- unlist(ss)
pp <- unlist(pp)
input_prs <- paste(input, collapse=',')
ss_prs <- paste(ss, collapse=',')
pp_prs <- paste(pp, collapse=',')

arg = paste0("python PRScsx/PRScsx.py --ref_dir=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/benchmark/LD_ref --bim_prefix=", PRS_CSx_working_dir, "snps", " --sst_file=", input_prs," --n_gwas=", ss_prs, " --pop=", pp_prs, " --chrom=", genos$bim[1,1]," --phi=1e-4 --out_dir=", PRS_CSx_working_dir, " --out_name=results")  
print(arg)
system(arg)

#retrieve posterior snp effect sizes, reassign weights
unique_pp <- unique(pp)
for (n in unique_pp){
	wgts <- weights_process(n, paste0(PRS_CSx_working_dir,"results_", n, "_pst_eff_a1_b0.5_phi1e-04_chr", genos$bim[1,1], ".txt"), wgts)
}


#---------------


ext <- length(wgts)

if(ext <= 0){
	stop("no external weights available")
}

#----------------



for ( i in 1:opt$crossval ) { #for every chunk in crossval
	
	
	if ( opt$verbose >= 1 ) cat("- Crossval fold",i,"\n")
	indx = which(folds==i,arr.ind=TRUE)
	cv.train = cv.all[-indx,] #training set is the other 4 groups 
	intercept = mean( cv.train[,3] ) 
	cv.train[,3] = scale(cv.train[,3]) 
	# hide current fold; writing new plink files for the 80% remaining
	cv.file = paste(opt$tmp,".cv",sep='')
	write.table( cv.train , quote=F , row.names=F , col.names=F , file=paste(cv.file,".keep",sep=''))	
	arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$tmp," --keep ",cv.file,".keep --out ",cv.file," --make-bed",sep='')
	system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

	for ( mod in 1:1 ) {
	
		#PT_SUMSTATS ----------------------------------------------------------------
		# genos$bed has people on rows, snps on column
		#linear regression per snp, extract p value and effect size, adding to dataframe 
		
		Betas <- c()
		Pvals <- c()
		SNPs <- c()
		for (col in 1:ncol(genos$bed)){
			SNPs <- append(SNPs, colnames(genos$bed)[col])
			snp <- genos$bed[ cv.sample[-indx],col]
			model <- lm(cv.train[,3] ~ snp)
			b <- coef(model)[2]
			if (is.na(b)){
				b = 0
			}
			Betas <- append(Betas, b)
			if ( nrow(summary(model)$coefficients) > 1 ){
				Pvals <- append(Pvals, summary(model)$coefficients[2, "Pr(>|t|)"])
			}else{
				Pvals <- append(Pvals, 1)
			}
		}	
		sumstats <- data.frame(SNPs, Betas, Pvals)

		#write dataframe, run plink --clump to P+T, write output to temp directory
		# p1 = 0.5, r2 = 0.2,   see Tiffany's paper here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8049522/
		sumstats_file = paste(opt$tmp,"_sumstats.txt",sep='')
		write.table(sumstats, file = sumstats_file, quote = F, row.names = F, col.names = T, sep = '\t')

		arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$ldref, "/GTEx_v8_genotype_AFR_HM3_exclude_dups.", genos$bim[1,1] ," --clump-p1 0.5 --clump-r2 0.2 --clump-snp-field SNPs --clump-field Pvals --clump ", opt$tmp, "_sumstats.txt"," --out ", sumstats_file,sep='')
        	system(arg)

		#read in result file, create wgt matrix and predict 
		if (file.exists(paste0(sumstats_file, ".clumped"))){
			clumped <- fread(file = paste0(sumstats_file, ".clumped"), header = T)
			clumps <- clumped[[3]]
			sumstats <- filter(sumstats, SNPs %in% clumps)

			m_sumstats <- match(genos$bim[,2], sumstats[[1]])
                	sumstats <- sumstats[m_sumstats,]
                	w_sumstats <- which(is.na(sumstats[[1]]))
                	if(length(w_sumstats) > 0){
                        	sumstats[w_sumstats,2] <- 0
                	}			

			pred.wgt.PT_sumstats = sumstats[[2]]
		}else{
			pred.wgt.PT_sumstats = rep(0, times = nrow(genos$bim))
		}

		if(length(pred.wgt.PT_sumstats) == 1){
			pred.wgt.PT_sumstats <- t(pred.wgt.PT_sumstats)
		}

		cv.calls[ indx , mod*2 ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt.PT_sumstats
		
		system(paste0("rm -rf ", sumstats_file, ".clumped"))

		#----------------------------------------------------------------------------


		#PRS-CSx----------------------------------------------------------------------

		#compute afr gene model weights with 4 fold 
		
		lasso_h2 <- hsq_afr[1]
		if(lasso_h2 < 0){lasso_h2 <- 0.064251}
		pred.wgt = weights.lasso( cv.file , lasso_h2 , snp=genos$bim[,2] )
		
		if ( sum(is.na(pred.wgt)) == length(pred.wgt)) {
			pred.wgt = weights.marginal( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3,drop=F]) , beta=T )
			pred.wgt[ - which.max( pred.wgt^2 ) ] = 0
		}

		if (length(pred.wgt) == 1){
			pred.wgt <- t(pred.wgt)
		}

		#collect adjusted weights from PRS-CSx
		eq = list("genos$bed[cv.sample[ -indx ], ] %*% pred.wgt")
		for (w in wgts){
			
			if (length(eval(parse(text = w))) == 1){
				assign(w, t(eval(parse(text = w))))
			}
			
			eq <- append(eq, paste0("genos$bed[ cv.sample[ -indx ], ] %*% ", w))
		}
		eq <- paste(eq, collapse="+")
		eq <- paste0("cv.train[,3] ~ ", eq)

		#compute optimal linear combination on validation set
		y <- lm(eval(parse(text = eq))) 

		coef <- coef(y)
		coef <- ifelse(is.na(coef), 0, coef)

		w_eq = list("coef[2]*pred.wgt")
		for (nums in 3:(ext+2)){
			w_eq <- append(w_eq, paste0("coef[", nums,"]*",wgts[nums-2]))
		}
		w_eq <- paste(w_eq, collapse="+")

		print(w_eq)

		pred.wgt.prs_csx <- eval(parse(text=w_eq))
	
		#if (length(pred.wgt.prs_csx) == 1){
		#	pred.wgt.prs_csx <- t(pred.wgt.prs_csx)
		#}

		#predict on testing set
		cv.calls[ indx , (mod*2 - 1) ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt.prs_csx
	
		#-------------------------------------------------------------------------------


	}
}

#uses 80% of data to predict into 20%, then uses another 80% to predict into the next 20%, and so on. 
#after everyone has a prediction, compute rsq and p value. 

#compute rsq + P-value for each model
cv.performance = matrix(NA,nrow=2,ncol=2)  
rownames(cv.performance) = c("rsq","pval") 
types <- c("PRS-CSx", "PT_sumstats") 
colnames(cv.performance) = types 

for ( mod in 1:ncol(cv.calls) ) { 
	if ( !is.na(sd(cv.calls[,mod])) && sd(cv.calls[,mod]) != 0 ) { #length 1 or 0 vector -> sd = NA, all values are identical -> sd = 0  
		reg = summary(lm( cv.all[,3] ~ cv.calls[,mod] )) #r^2 between predicted and actual 
		cv.performance[ 1, mod ] = reg$adj.r.sq
		cv.performance[ 2, mod ] = reg$coef[2,4]
	} else {
		cv.performance[ 1, mod ] = NA
		cv.performance[ 2, mod ] = NA
	}
}
if ( opt$verbose >= 1 ) write.table(cv.performance,quote=F,sep='\t')


snps = genos$bim

save( snps, cv.performance , hsq_afr, hsq_afr.pv, N.tot , wgts, models, file = paste( opt$out , ".wgt.RDat" , sep='' ) )

# --- CLEAN-UP
if ( opt$verbose >= 1 ) cat("Cleaning up\n")
cleanup()
}else{
	cat("Gene is not cis heritable in all of the populations")
}

