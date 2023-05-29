# load packages
suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('glmnet'))
suppressMessages(library('methods'))
suppressMessages(library('data.table'))
suppressMessages(library('hash'))

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
					enet:\t Elastic-net regression (with mixing parameter of 0.5)\n")			  
		  
)

# parse command-line arguments
opt = parse_args(OptionParser(option_list=option_list))
print(opt)

# find gene name without version number
name <- strsplit(opt$gene, ".", fixed = TRUE)[[1]][1]
print(name)

# assign all external dataset file names
file.test <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/","v8_allEUR_",opt$tissue,"/",opt$tissue,".",opt$gene,".wgt.RDat")
ota_cells <- list("CD16p_Mono","CL_Mono","LDG","Mem_CD4","Mem_CD8","NK","Naive_B","Naive_CD4","Naive_CD8","Neu","Plasmablast","mDC","pDC")
for (c in ota_cells){
	assign(paste0("file.ota.", c), paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/OTA/genes/", c, "/", name, ".txt"))
}
file.mesahis <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/MESA_HIS/filtered/clumped_sig/",name,".txt") 
file.genoa <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/genoa/significant_clumped/clumped/",name,".txt")
file.mesa <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/MESA/filtered/clumped_sig/",name,".txt") 
file.jung <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/Jung/clumped/",name,".txt")
file.ishi.B <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/Ishigaki/permutation/B_permutation/clumped/",name,".txt")
file.ishi.CD4 <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/Ishigaki/permutation/CD4_permutation/clumped/",name,".txt")
file.ishi.CD8 <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/Ishigaki/permutation/CD8_permutation/clumped/",name,".txt")
file.ishi.Mono <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/Ishigaki/permutation/Mono_permutation/clumped/",name,".txt")
file.ishi.NK <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/Ishigaki/permutation/NK_permutation/clumped/",name,".txt")
file.ishi.PB <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/Ishigaki/permutation/PB_permutation/clumped/",name,".txt")
file.eqtlgen <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/eQTLGEN/clumped/", name, ".txt")
file.peruvian <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/peruvian/genes/", name, ".txt")


# create a hashmap to store sample size of each datset (used later for sample-size weighted meta-analysis)
h <- hash()
h[["eur"]] <- 574
h[["ota"]] <- 416
h[["jung"]] <- 101
h[["genoa"]] <- 1031
h[["ishi"]] <- 105
h[["mesa"]] <- 233
h[["mesahis"]] <- 352
h[["eqtlgen"]] <- 31684 
h[["peruvian"]] <- 259 

# all available external datasets to a list
datasets <- list(file.test, file.ota.CD16p_Mono, file.ota.CL_Mono, file.ota.LDG, file.ota.Mem_CD4, file.ota.Mem_CD8, file.ota.NK, file.ota.Naive_B, file.ota.Naive_CD4, file.ota.Naive_CD8, file.ota.Neu, file.ota.Plasmablast, file.ota.mDC, file.ota.pDC, file.mesahis, file.genoa, file.mesa, file.jung, file.ishi.B, file.ishi.CD4, file.ishi.CD8, file.ishi.Mono, file.ishi.NK, file.ishi.PB, file.peruvian)
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

#default crossval =5 fold split
if ( opt$crossval <= 1 ) { 
if ( opt$verbose >= 1 ) cat("### skipping cross-validation\n")
} else {
if ( opt$verbose >= 1 ) cat("### performing",opt$crossval,"fold cross-validation\n")
cv.all = pheno #note, cv.all is the gene expression phenotype for all 80 afr individuals for whole blood for ex 
n = nrow(cv.all)
cv.sample = sample(n) #sample randomly 
cv.all = cv.all[ cv.sample , ]
folds = cut(seq(1,n),breaks=opt$crossval,labels=FALSE) #5 fold split - split into 5 groups 
cv.calls = matrix(NA,nrow=n,ncol=3) 
}


#list of all wgt vectors 
wgts <- list()


#process EUR gene model here: 
if (file.test %in% datasets){
load(file.test) 
#wgt.RDat includes the following: 
	#snps
	#cv.performance - rsq and pval for all 4 models 
	#wgt.matrix - weights for each snp for each model 
hsq_eur <- hsq
hsq_eur.pv <- hsq.pv

bestmod <- which.min(cv.performance[2,]) #model with the smallest p value 
if(bestmod == 1){
	wgt.matrix[ - which.max( wgt.matrix[,1]^2 ) , 1] = 0
} 
#for top1, identify the snp with all the weight and set everything else to 0  #how do we know which snp carries all the weight? largest weight^2
pred.wgt.eur <- wgt.matrix[,bestmod] #take weights for each snp from the best model 
pred.wgt.eur[is.na(pred.wgt.eur)] <- 0 

for (k in 1:nrow(snps)){
		match=which(genos$bim$V2 == snps$V2[k])
		if (identical(match, integer(0))){
			print(paste0("eur:", "no snp in common, skipping ref/alt check iteration"))
		}else{
		if(!is.na(snps$V5[k]) && !is.na(snps$V6[k])){
		if (genos$bim$V5[match] != snps$V5[k] || genos$bim$V6[match] != snps$V6[k]) {
			if (genos$bim$V5[match] == snps$V6[k] && genos$bim$V6[match] == snps$V5[k]){
				pred.wgt.eur[k] <- pred.wgt.eur[k] * -1
			}else{
				pred.wgt.eur[k] <- 0
			}
		}
		}else {
			pred.wgt.eur[k] <- 0
		}
		}
}

m <- match(genos$bim[,2],names(pred.wgt.eur))
pred.wgt.eur <- pred.wgt.eur[m] #reorder eur if necessary - confirmed that order of the snps are the same 
#if there are more snps in eur than afr, they will be removed in the line above 
#if there are more snps in afr than in eur, they will simply have value of NA from the line above. 
#therefore, weight must be zero for any NA values
w <- which(is.na(pred.wgt.eur))
if(length(w)>0){pred.wgt.eur[w] <- 0}
if(sum(which(pred.wgt.eur != 0)) > 0){
	wgts <- append(wgts, "pred.wgt.eur")
}
}else{
	hsq_eur <- NA
	hsq_eur.pv <- NA
}


#IMPORTANT: for the eQTL stats, use only snps present in both populations above 
#Note: SNPs in the eQTL stats that are not in AFR data are removed, SNPs in the AFR but not in the eQTL stats will have a weight 0 in the eQTL weights 
#To simplify, if the SNP in AFR is in neither EUR or other populations, only the AFR weight is used. 

#function to process ota data for different cell types 
ota_process <- function(file, cell, wgts){
	table.ota <- fread(file, select = c(6, 14, 15, 16))
	for (k in 1:nrow(table.ota)){
		match=which(genos$bim$V2 == table.ota$Variant_ID[k])
		if (identical(match, integer(0))){
			print(paste0("ota,", cell, ":", "no snp in common, skipping ref/alt check iteration"))
		}else{
		if (!is.na(table.ota$A1[k]) && !is.na(table.ota$A0[k])){
		if (genos$bim$V5[match] != table.ota$A1[k] || genos$bim$V6[match] != table.ota$A0[k]) {
			if (genos$bim$V5[match] == table.ota$A0[k] && genos$bim$V6[match] == table.ota$A1[k]){
				table.ota$Backward_slope[k] <- table.ota$Backward_slope[k] * -1
			}else{
				table.ota$Backward_slope[k] <- 0
			}
		}
		}else {
			table.ota$Backward_slope[k] <- 0
		}
		}
	}
	mOTA <- match(genos$bim[,2], table.ota$Variant_ID)
	table.ota <- table.ota[mOTA,]
	wOTA <- which(is.na(table.ota[,1]))
	if(length(wOTA) > 0){table.ota[wOTA,2] <- 0} 
	assign(paste0("pred.wgt.ota.", cell ), table.ota$Backward_slope, envir = parent.frame())
	if(sum(which(eval(parse(text = paste0("pred.wgt.ota.", cell))) != 0)) > 0){
		wgts <- append(wgts, paste0("pred.wgt.ota.", cell))
		return (wgts)
	}else {
		return (wgts)
	}
}


#OTA
all_ota <- c(file.ota.CD16p_Mono, file.ota.CL_Mono, file.ota.LDG, file.ota.Mem_CD4, file.ota.Mem_CD8, file.ota.NK, file.ota.Naive_B, file.ota.Naive_CD4, file.ota.Naive_CD8, file.ota.Neu, file.ota.Plasmablast, file.ota.mDC, file.ota.pDC)

print("processing ota files")
for (o in all_ota){
	if(o %in% datasets){
		splits <- strsplit(o, split = "/")[[1]][10]
		print(splits)
		wgts <- ota_process(o, splits, wgts)
	}
}


#JUNG 
if (file.jung %in% datasets){
table.jung <- fread(file.jung, select = c(3, 13, 14, 15)) 
for (k in 1:nrow(table.jung)){
		match=which(genos$bim$V2 == table.jung$SNP[k])
		if (identical(match, integer(0))){
			print(paste0("JUNG:", "no snp in common, skipping ref/alt check iteration"))
		}else{
		if (!is.na(table.jung$A1[k]) && !is.na(table.jung$A0[k])){
		if (genos$bim$V5[match] != table.jung$A1[k] || genos$bim$V6[match] != table.jung$A0[k]) {
			if (genos$bim$V5[match] == table.jung$A0[k] && genos$bim$V6[match] == table.jung$A1[k]){
				table.jung$SLOPE[k] <- table.jung$SLOPE[k] * -1
			}else{
				table.jung$SLOPE[k] <- 0
			}
		}
		}else {
			table.jung$SLOPE[k] <- 0
		}
		}
}
mJUNG <- match(genos$bim[,2], table.jung$SNP)
table.jung <- table.jung[mJUNG,]
wJUNG <- which(is.na(table.jung[,1]))
if(length(wJUNG) > 0){table.jung[wJUNG,2] <- 0} 
pred.wgt.jung <- table.jung$SLOPE
if (sum(which(pred.wgt.jung != 0 )) > 0){
	wgts <- append(wgts, "pred.wgt.jung")
}
}


#GENOA
if (file.genoa %in% datasets){
table.genoa <- fread(file.genoa, select = c(3, 13, 15, 16)) 
for (k in 1:nrow(table.genoa)){
		match=which(genos$bim$V2 == table.genoa$SNP[k])
		if (identical(match, integer(0))){
			print(paste0("GENOA:", "no snp in common, skipping ref/alt check iteration"))
		}else{
		if (!is.na(table.genoa$A1[k]) && !is.na(table.genoa$A0[k])){
		if (genos$bim$V5[match] != table.genoa$A1[k] || genos$bim$V6[match] != table.genoa$A0[k]) {
			if (genos$bim$V5[match] == table.genoa$A0[k] && genos$bim$V6[match] == table.genoa$A1[k]){
				table.genoa$SLOPE[k] <- table.genoa$SLOPE[k] * -1
			}else{
				table.genoa$SLOPE[k] <- 0
			}
		}
		}else {
			table.genoa$SLOPE[k] <- 0
		}
		}
}
mGENOA <- match(genos$bim[,2], table.genoa$SNP)
table.genoa <- table.genoa[mGENOA,]
wGENOA <- which(is.na(table.genoa[,1]))
if(length(wGENOA) > 0){table.genoa[wGENOA,2] <- 0} 
pred.wgt.genoa <- table.genoa$SLOPE
if (sum(which(pred.wgt.genoa != 0)) > 0){
	wgts <- append(wgts, "pred.wgt.genoa")
}
}

#function to process ishigaki data for different cell types
ishigaki_process <- function(file, cell, wgts){
	table.ishi <- fread(file, select = c(3, 13, 14, 15)) 
	for (k in 1:nrow(table.ishi)){
		match=which(genos$bim$V2 == table.ishi$SNP[k])
		if (identical(match, integer(0))){
			print(paste0("ishi,", cell, ":", "no snp in common, skipping ref/alt check iteration"))
		}else{
		if (!is.na(table.ishi$A1[k]) && !is.na(table.ishi$A0[k])){
			if (genos$bim$V5[match] != table.ishi$A1[k] || genos$bim$V6[match] != table.ishi$A0[k]) {
				if (genos$bim$V5[match] == table.ishi$A0[k] && genos$bim$V6[match] == table.ishi$A1[k]){
					table.ishi$SLOPE[k] <- table.ishi$SLOPE[k] * -1
				}else{
					table.ishi$SLOPE[k] <- 0
				}
			}
		} else {
			table.ishi$SLOPE[k] <- 0
		}
		}
	}
	mIshi <- match(genos$bim[,2], table.ishi$SNP)
	table.ishi <- table.ishi[mIshi,]
	wIshi <- which(is.na(table.ishi[,1]))
	if(length(wIshi) > 0){table.ishi[wIshi,2] <- 0} 
	assign(paste0("pred.wgt.ishi.", cell ), table.ishi$SLOPE, envir = parent.frame())
	if(sum(which(eval(parse(text = paste0("pred.wgt.ishi.", cell))) != 0)) > 0){
		wgts <- append(wgts, paste0("pred.wgt.ishi.", cell))
		return (wgts)
	}else {
		return (wgts)
	}
}


#ISHIGAKI


if(file.ishi.B %in% datasets){
	wgts <- ishigaki_process(file.ishi.B, "B", wgts)
}


if(file.ishi.CD4 %in% datasets){
	wgts <- ishigaki_process(file.ishi.CD4, "CD4", wgts)
}


if(file.ishi.CD8 %in% datasets){
	wgts <- ishigaki_process(file.ishi.CD8, "CD8", wgts)
}


if(file.ishi.Mono %in% datasets){
	wgts <- ishigaki_process(file.ishi.Mono, "Mono", wgts)
}


if(file.ishi.NK %in% datasets){
	wgts <- ishigaki_process(file.ishi.NK, "NK", wgts)
}


if(file.ishi.PB %in% datasets){
	wgts <- ishigaki_process(file.ishi.PB, "PB", wgts)
}


#MESA
if (file.mesa %in% datasets){
table.mesa <- fread(file.mesa, select = c(3, 13, 14, 15)) 
for (k in 1:nrow(table.mesa)){
		match=which(genos$bim$V2 == table.mesa$SNP[k])
		if (identical(match, integer(0))){
			print(paste0("MESA:", "no snp in common, skipping ref/alt check iteration"))
		}else{
		if (!is.na(table.mesa$A1[k]) && !is.na(table.mesa$A0[k]) ){
		if (genos$bim$V5[match] != table.mesa$A1[k] || genos$bim$V6[match] != table.mesa$A0[k]) {
			if (genos$bim$V5[match] == table.mesa$A0[k] && genos$bim$V6[match] == table.mesa$A1[k]){
				table.mesa$SLOPE[k] <- table.mesa$SLOPE[k] * -1
			}else{
				table.mesa$SLOPE[k] <- 0
			}
		}
		} else {
			table.mesa$SLOPE[k] <- 0
		}
		}
}
mMESA <- match(genos$bim[,2], table.mesa$SNP)
table.mesa <- table.mesa[mMESA,]
wMESA <- which(is.na(table.mesa[,1]))
if(length(wMESA) > 0){table.mesa[wMESA,2] <- 0} 
pred.wgt.mesa <- table.mesa$SLOPE
if (sum(which(pred.wgt.mesa != 0)) > 0){
	wgts <- append(wgts, "pred.wgt.mesa")
}
}

#MESA_HIS
if (file.mesahis %in% datasets){
table.mesahis <- fread(file.mesahis, select = c(3, 13, 14, 15)) 
for (k in 1:nrow(table.mesahis)){
		match=which(genos$bim$V2 == table.mesahis$SNP[k])
		if (identical(match, integer(0))){
			print(paste0("MESAHIS:", "no snp in common, skipping ref/alt check iteration"))
		}else{
		if(!is.na(table.mesahis$A1[k]) && !is.na(table.mesahis$A0[k])){
		if (genos$bim$V5[match] != table.mesahis$A1[k] || genos$bim$V6[match] != table.mesahis$A0[k]) {
			if (genos$bim$V5[match] == table.mesahis$A0[k] && genos$bim$V6[match] == table.mesahis$A1[k]){
				table.mesahis$SLOPE[k] <- table.mesahis$SLOPE[k] * -1
			}else{
				table.mesahis$SLOPE[k] <- 0
			}
		}
		}else {
			table.mesahis$SLOPE[k] <- 0
		}
		}
}
mMESAhis <- match(genos$bim[,2], table.mesahis$SNP)
table.mesahis <- table.mesahis[mMESAhis,]
wMESAhis <- which(is.na(table.mesahis[,1]))
if(length(wMESAhis) > 0){table.mesahis[wMESAhis,2] <- 0} 
pred.wgt.mesahis <- table.mesahis$SLOPE
if (sum(which(pred.wgt.mesahis != 0)) > 0) {
	wgts <- append(wgts, "pred.wgt.mesahis")
}
}

#eqtlgen
if (file.eqtlgen %in% datasets){
table.eqtlgen <- fread(file.eqtlgen, select = c(3, 13, 14, 15)) 
for (k in 1:nrow(table.eqtlgen)){
		match=which(genos$bim$V2 == table.eqtlgen$SNP[k])
		if (identical(match, integer(0))){
			print(paste0("eqtlgen:", "no snp in common, skipping ref/alt check iteration"))
		}else{
		if (!is.na(table.eqtlgen$A1[k]) && !is.na(table.eqtlgen$A0[k]) ){
		if (genos$bim$V5[match] != table.eqtlgen$A1[k] || genos$bim$V6[match] != table.eqtlgen$A0[k]) {
			if (genos$bim$V5[match] == table.eqtlgen$A0[k] && genos$bim$V6[match] == table.eqtlgen$A1[k]){
				table.eqtlgen$SLOPE[k] <- table.eqtlgen$SLOPE[k] * -1
			}else{
				table.eqtlgen$SLOPE[k] <- 0
			}
		}
		} else {
			table.eqtlgen$SLOPE[k] <- 0
		}
		}
}
meqtlgen <- match(genos$bim[,2], table.eqtlgen$SNP)
table.eqtlgen <- table.eqtlgen[meqtlgen,]
weqtlgen <- which(is.na(table.eqtlgen[,1]))
if(length(weqtlgen) > 0){table.eqtlgen[weqtlgen,2] <- 0} 
pred.wgt.eqtlgen <- table.eqtlgen$SLOPE
if (sum(which(pred.wgt.eqtlgen != 0)) > 0){
	wgts <- append(wgts, "pred.wgt.eqtlgen")
}
}


#peruvian - flipped signs of beta already 
if (file.peruvian %in% datasets){
table.peruvian <- fread(file.peruvian, select = c(3, 7)) 
mperuvian <- match(genos$bim[,2], table.peruvian$rsid)
table.peruvian <- table.peruvian[mperuvian,]
wperuvian <- which(is.na(table.peruvian[,2]))
if(length(wperuvian) > 0){table.peruvian[wperuvian,1] <- 0} 
pred.wgt.peruvian <- table.peruvian$Beta
if (sum(which(pred.wgt.peruvian != 0)) > 0){
	wgts <- append(wgts, "pred.wgt.peruvian")
}
}

print("edited to only use snps common across all populations")

#metanalyze (by sample size) cell types from same dataset 

new_wgts <- c()
index <- c()
iota <- grep("ota", wgts)
if (!identical(iota, integer(0))){
ota_used <- wgts[iota]
index <- append(index, iota)
len_ota <- length(ota_used)
ota_sum <- 0
for (o in ota_used){
ota_sum <- ota_sum + eval(parse(text=o))
}
pred.wgt.ota <- ota_sum/len_ota
new_wgts <- append(new_wgts, "pred.wgt.ota")
}

iishi <- grep("ishi", wgts)
if (!identical(iishi, integer(0))){
ishigaki_used <- wgts[iishi]
index <- append(index, iishi)
len_ishi <- length(ishigaki_used)
ishi_sum <- 0
for (i in ishigaki_used){
ishi_sum <- ishi_sum + eval(parse(text=i))
}
pred.wgt.ishi <- ishi_sum/len_ishi
new_wgts <- append(new_wgts, "pred.wgt.ishi")
}
if (length(index) > 0){
	new_wgts <- append(new_wgts, wgts[-index])
	wgts = new_wgts
	print(wgts)
	print("cell types combined")
}

ext <- length(wgts)

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

		print("performing LASSO and top1 as backup")
		
		lasso_h2 <- hsq_afr[1]
		if(lasso_h2 < 0){lasso_h2 <- 0.064251}
		pred.wgt = weights.lasso( cv.file , lasso_h2 , snp=genos$bim[,2] )
		
		if ( sum(is.na(pred.wgt)) == length(pred.wgt)) {
			print("all weights pushed to 0, using top1 as backup")
			pred.wgt = weights.marginal( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3,drop=F]) , beta=T )
			pred.wgt[ - which.max( pred.wgt^2 ) ] = 0
		}
		# predict from weights into sample

		print("predicting from weights into sample")
		
		# when there is only 1 snp in the wgt file, we need to reformat	
		if (length(pred.wgt) == 1){
			pred.wgt <- t(pred.wgt)
		}
		
		#AFR ONLY------------------------------------------------------------------
		cv.calls[ indx , mod*3-2 ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt
		#cv.calls is a matrix of predicted value - 1 fold at a time
		#pred.wgt - training on other 4 folds
		#indx = people in 1 fold we removed when we trained using 4 fold
		#--------------------------------------------------------------------------

		#META-ANALYSIS-------------------------------------------------------------
		#sample-size weighted meta-analysis 
		#calculate total sample size 
		total <- 0
		for (w in wgts){
			dataset <- strsplit(w, split="[.]")[[1]][3]
			size <- h[[dataset]]
			total = total + size
		}
		total = total + 80
		pred.wgt.meta <- pred.wgt * (80/total)
		for (w in wgts){
			dataset <- strsplit(w, split="[.]")[[1]][3]
			size <- h[[dataset]]
			pred.wgt.meta = pred.wgt.meta + (eval(parse(text=w)) * (size/total))
		}
		if (length(pred.wgt.meta) == 1){
			pred.wgt.meta <- t(pred.wgt.meta)
		}
		
		cv.calls[ indx , mod*3-1 ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt.meta 
		#-----------------------------------------------------------------------------

		#MAGEPRO----------------------------------------------------------------------

		#create a matrix of weights from all data sets - ridge regression input
		w1 <- genos$bed[ cv.sample[ -indx ] , ] %*% pred.wgt
		eq <- matrix(0, nrow = length (w1), ncol = ext + 1)
		eq[, 1] <- w1
		num = 2
		for (w in wgts){
			eq[,num] <- genos$bed[ cv.sample[ -indx ] , ] %*%  eval(parse(text = w))
			num = num+1
		}	
		
		#cv.glmnet runs if there are two or more columns
		if (ext > 0){
			y <- cv.glmnet(x = eq , y = cv.train[,3], alpha = 0, nfold = 5, intercept = T, standardize = T)
			cf = coef(y, s = "lambda.min")[2:(ext+2)]	
			predtext <- "cf[1]*pred.wgt"
			for (i in 2:(length(cf))){
					predtext <- paste0(predtext, " + cf[", i, "]*", wgts[(i-1)])
			}	
			pred.wgt.magepro <- eval(parse(text = predtext))
		}else{
			pred.wgt.magepro <- pred.wgt
			cf = NA
		}

		if(length(pred.wgt.magepro) == 1){
			pred.wgt.magepro <- t(pred.wgt.magepro)
		}	

		cv.calls[ indx , mod*3 ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt.magepro
		#-------------------------------------------------------------------------------


	}
}


#uses 80% of data to predict into 20%, then uses another 80% to predict into the next 20%, and so on. 
#after everyone has a prediction, compute rsq and p value. 

#compute rsq + P-value for each model
cv.performance = matrix(NA,nrow=2,ncol=3)  
rownames(cv.performance) = c("rsq","pval") 
types <- c("afr","meta","multipop") 
x <- c() 
for (i in models){for (j in types){x <- c(x,paste0(i,",",j))}} #model,type 
colnames(cv.performance) = x  

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

# ---- Full Analysis 
if (opt$verbose >= 1) cat("Computing full-sample weights \n")

wgt.matrix = matrix(0, nrow=nrow(genos$bim), ncol=3) 
colnames(wgt.matrix) = c("AFRONLY","META","MULTIPOP")
rownames(wgt.matrix) = genos$bim[,2]

pred.wgt2 = weights.lasso( geno.file , lasso_h2 , snp=genos$bim[,2] )	
if ( sum(is.na(pred.wgt2)) == length(pred.wgt2)) {
if (opt$verbose >= 1) cat("all weights pushed to 0, using top1 as backup\n")
pred.wgt2 = weights.marginal( genos$bed , as.matrix(pheno[,3]) , beta=T )
pred.wgt2[ - which.max( pred.wgt2^2 ) ] = 0
}

wgt.matrix[, 1] = pred.wgt2

#full meta-analysis
pred.wgt.meta2 <- pred.wgt2 * (80/total)
for (w in wgts){
	dataset <- strsplit(w, split="[.]")[[1]][3]
	size <- h[[dataset]]
	pred.wgt.meta2 = pred.wgt.meta2 + (eval(parse(text=w)) * (size/total))
}
wgt.matrix[, 2] = pred.wgt.meta2

#create matrix for glmnet input
w1 <- genos$bed %*% pred.wgt2
eq2 <- matrix(0, nrow = length (w1), ncol = ext + 1)
eq2[, 1] <- w1
num = 2
for (w in wgts){
	eq2[,num] <- genos$bed %*% eval(parse(text = w))
	num = num+1
}	

#store alpha*beta for all snps per dataset (for downstream analysis of which datasets contribute the most)
alpha_beta <- matrix(0, nrow = length(pred.wgt2), ncol = ext + 1)

#run ridge regression to find optimal coefficients and compute multipop weight
if (ext > 0){
	y2 <- cv.glmnet(x = eq2 , y = pheno[,3], alpha = 0, nfold = 5, intercept = T, standardize = T)
	cf_total = coef(y2, s = "lambda.min")[2:(ext+2)]	
	total_coeff <- cf_total
	predtext2 <- "cf_total[1]*pred.wgt2"
	alpha_beta[, 1] <- eval(parse(text = predtext2))
	for (i in 2:(length(cf_total))){
		predtext2 <- paste0(predtext2, " + cf_total[", i, "]*", wgts[(i-1)])
		alpha_beta[, i] <- cf_total[i]*eval(parse(text = wgts[(i-1)]))
	}	
	pred.wgt.multipop <- eval(parse(text = predtext2))
}else{
	pred.wgt.multipop <- pred.wgt
	total_coeff = c(NA)
	alpha_beta = NA
}		

print(cor(pheno[,3], (genos$bed %*% pred.wgt.multipop))^2 )

wgt.matrix[, 3] = pred.wgt.multipop

#---

snps = genos$bim

#add african weight to the list of weights used
wgts <- append("pred.wgt", wgts)
print(total_coeff)

if (is.matrix(alpha_beta)){
	colnames(alpha_beta) <- wgts
}

#if the alpha is 0 for any dataset, it did not add anything to the model
if (sum(is.na(total_coeff)) == 0){
w <- which(total_coeff == 0)
#remove that dataset from the list of datasets used in the model
if (length(w) > 0 ) {
wgts <- wgts[-w]
total_coeff <- total_coeff[-w]
}
}
print(wgts)


save( wgt.matrix, alpha_beta, snps, cv.performance , hsq_afr, hsq_afr.pv, hsq_eur, hsq_eur.pv, N.tot , wgts, total_coeff, models, file = paste( opt$out , ".wgt.RDat" , sep='' ) )

# --- CLEAN-UP
if ( opt$verbose >= 1 ) cat("Cleaning up\n")
cleanup()
}else{
	cat("Gene is not cis heritable in all of the populations")
}
