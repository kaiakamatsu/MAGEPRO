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
  make_option("--datasets", action="store", default="", type='character', 
	      help="Comma-separated list of external datasets to include")
)


# --- PREPARING EXTERNAL DATASET FILES AND PARSING OPTIONS
# parse command-line arguments
opt = parse_args(OptionParser(option_list=option_list))
print(opt)

# find gene name without version number
name <- strsplit(opt$gene, ".", fixed = TRUE)[[1]][1]

# datasets to use
datas <- strsplit(opt$datasets, ",", fixed = TRUE)[[1]]

# datasets to use
datas <- strsplit(opt$datasets, ",", fixed = TRUE)[[1]]
print(datas)

# assign all external dataset file names
file.test <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/","v8_allEUR_",opt$tissue,"/",opt$tissue,".",opt$gene,".wgt.RDat")
ota_cells <- list("CD16p_Mono","CL_Mono","LDG","Mem_CD4","Mem_CD8","NK","Naive_B","Naive_CD4","Naive_CD8","Neu","Plasmablast","mDC","pDC")
for (c in ota_cells){
	assign(paste0("file.ota.", c), paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/OTA_nominal/genes/", c, "/", name, ".txt"))
}
file.his <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/MESA_HIS/filtered/genes/",name,".txt") 
file.genoa <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/genoa/genes_nominal/genes/",name,".txt")
file.mesa <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/MESA/filtered/genes/",name,".txt") 
file.eqtlgen <- paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/eQTLGEN/genes/", name, ".txt")

# create a hashmap to store sample size of each datset (used later for sample-size weighted meta-analysis)
h <- hash()
h[["eur"]] <- 574
h[["ota"]] <- 416
h[["genoa"]] <- 1031
h[["mesa"]] <- 233
h[["his"]] <- 352
h[["eqtlgen"]] <- 31684 

# all available external datasets to a list
datasets <- list("file.test", "file.ota.CD16p_Mono", "file.ota.CL_Mono", "file.ota.LDG", "file.ota.Mem_CD4", "file.ota.Mem_CD8", "file.ota.NK", "file.ota.Naive_B", "file.ota.Naive_CD4", "file.ota.Naive_CD8", "file.ota.Neu", "file.ota.Plasmablast", "file.ota.mDC", "file.ota.pDC", "file.his", "file.genoa", "file.mesa", "file.eqtlgen")

# filter for datasets you want to provide MAGEPRO
datasets <- datasets[grepl(paste(datas, collapse = "|"), datasets)]
datasets <- lapply(datasets, function(x) eval(parse(text=x)))

# remove datasets that are not available for that gene
for (d in datasets) {
	if (!file.exists(d)){
		datasets <- datasets[datasets != d]
	}
}

if ( opt$verbose == 2 ) {
	cat("Datasets available for this gene: \n")
	print(datasets)
}

if(length(datasets) >= 0) {

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}

# --- PREDICTION MODELS

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

if ( opt$verbose == 2 ) cat("Predictive Models Prepared \n")

# --- CLEANUP
cleanup = function() {
	if ( ! opt$noclean ) {
		arg = paste("rm -f " , opt$tmp , "*", sep='')
		system(arg)
                arg = paste("rm -f " , opt$gemmaout , "*", sep='')  
		system(arg)
	}
}


# --- i/o CHECKS
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


if ( opt$verbose == 2 ) cat("I/O checks complete \n")

# --- PREPARE PHENO, LOAD COVAR

#fam file read - the modified fam file version (added ge data)  
fam = read.table(paste(opt$bfile,".fam",sep=''),as.is=T)

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

if ( opt$rn ) {
	library('GenABEL')
	library(preprocessCore)
	pheno[,3] = rntransform( pheno[,3] )
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}

# Load in the covariates if needed
if ( !is.na(opt$covar) ) { 
	covar = ( read.table(opt$covar,as.is=T,head=T) )
	if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )  
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	pheno = pheno[m.keep,]
	m = m[m.keep]
	covar = covar[m,] #reordering covariates to match fam file
	reg = summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))  
	if ( opt$verbose >= 1 ) cat( reg$r.sq , "variance in phenotype explained by covariates\n" )
	pheno[,3] = scale(reg$resid) #regressing out covar to single out genetic effect
	raw.pheno.file = pheno.file
	pheno.file = paste(pheno.file,".resid",sep='')
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file) #newly scaled pheno file 
}


if ( opt$verbose == 2 ) cat("Covariates loaded and regressed out of phenotype \n")

geno.file = opt$tmp
# recode to the intersection of samples and new phenotype
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$bfile," --pheno ",pheno.file," --keep ",pheno.file," --make-bed --out ",geno.file,sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)


# --- HERITABILITY ANALYSIS
if ( is.na(opt$hsq_set) ) {
	if ( opt$verbose >= 1 ) cat("### Estimating heritability\n")

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
		if ( opt$verbose >= 1 ) cat(opt$tmp,"does not exist, GCTA could not converge, forcing h2 = 0.064251\n",file=stderr())
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
} else {
	if ( opt$verbose >= 1 ) cat("### Skipping heritability estimate\n")
	hsq_afr = opt$hsq_set
	hsq_afr.pv = NA
}


# --- LOAD AND PREP GENOTYPES 
genos = read_plink(geno.file,impute="avg")
mafs = apply(genos$bed,2,mean)/2  
sds = apply(genos$bed,2,sd)  
# important : genotypes are standardized and scaled here:
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
	if ( opt$verbose >= 1 ) cat("Regressing covariates out of the genotypes\n")
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
r2_training_magepro <- c()
r2_training_afr <- c()
}


#list of all wgt vectors 
wgts <- list()


# --- EUR gene model
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

#flip weights based on effect allele
for (k in 1:nrow(snps)){
		
	match=which(genos$bim$V2 == snps$V2[k])
	
	if (identical(match, integer(0))){
		print(paste0("eur:", "no snp in common, skipping ref/alt check iteration"))
	}else{
		if(!is.na(snps$V5[k]) && !is.na(snps$V6[k])){
			if (as.character(genos$bim$V5[match]) != as.character(snps$V5[k]) || as.character(genos$bim$V6[match]) != as.character(snps$V6[k])) {
				if (as.character(genos$bim$V5[match]) == as.character(snps$V6[k]) && as.character(genos$bim$V6[match]) == as.character(snps$V5[k])){
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

if ( opt$verbose == 2) cat("EUR gene models processed \n")

#IMPORTANT: for the eQTL stats, use only snps present in both populations above 
#Note: SNPs in the eQTL stats that are not in AFR data are removed, SNPs in the AFR but not in the eQTL stats will have a weight 0 in the eQTL weights 
#If the SNP in AFR is in neither EUR or other populations, only the AFR weight is used. 

# --- PROCESS OTHER EXTERNAL DATASETS

datasets_process <- function(dataset, file, cell, wgts, snp, A1, A0, B){
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
	table <- fread(file, select = c(snp, A1, A0, B)) 
	for (k in 1:nrow(table)){
		match=which(genos$bim$V2 == table[k, 1])
		if (! identical(match, integer(0))){
			if (!is.na(table[k, 2]) && !is.na(table[k, 3])){
				if (as.character(genos$bim$V5[match]) != as.character(table[k, 2]) || as.character(genos$bim$V6[match]) != as.character(table[k, 3])) {  # check if A1 column is effect allele
					if (as.character(genos$bim$V5[match]) == as.character(table[k, 3]) && as.character(genos$bim$V6[match]) == as.character(table[k, 2])){
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
	m <- match(genos$bim[,2], table[[1]])
	table <- table[m,]
	w <- which(is.na(table[[1]]))
	if(length(w) > 0){table[w,4] <- 0}
	assign(paste0("pred.wgt.", dataset , cell ), table[[4]], envir = parent.frame())
	if(sum(which(eval(parse(text = paste0("pred.wgt.", dataset, cell))) != 0)) > 0){
		wgts <- append(wgts, paste0("pred.wgt.", dataset, cell))
		return (wgts)
	}else {
		return (wgts)
	}
}

#OTA 
all_ota <- c(file.ota.CD16p_Mono, file.ota.CL_Mono, file.ota.LDG, file.ota.Mem_CD4, file.ota.Mem_CD8, file.ota.NK, file.ota.Naive_B, file.ota.Naive_CD4, file.ota.Naive_CD8, file.ota.Neu, file.ota.Plasmablast, file.ota.mDC, file.ota.pDC)

for (o in all_ota){
	if(o %in% datasets){
		cell <- strsplit(o, split = "/")[[1]][10]
		wgts <- datasets_process("ota", o, cell, wgts, 6, 8, 7, 13)
	}
}

#GENOA
if (file.genoa %in% datasets){
	wgts <- datasets_process("genoa", file.genoa, NA, wgts, 3, 5, 6, 8)
}

#MESA
if (file.mesa %in% datasets){
	wgts <- datasets_process("mesa", file.mesa, NA, wgts, 1, 14, 13, 6)
}

#MESA_HIS
if (file.his %in% datasets){
	wgts <- datasets_process("his", file.his, NA, wgts, 1, 14, 13, 6)
}

#eqtlgen
if (file.eqtlgen %in% datasets){
	wgts <- datasets_process("eqtlgen", file.eqtlgen, NA, wgts, 1, 4, 5, 12)
}



# --- META-ANALYZE DIFFERENT CELL TYPES TOGETHER
new_wgts <- c()
index <- c()

#function to meta-analyze cell types - combine into one vector 
#inputs 
	#dataset = name of dataset 
	#new_wgts = list of newly assigned weights 
	#indices = list of indices that will be replaced by one new vector 
meta_cells <- function(dataset, result_wgts, indices){
	i <- grep(dataset, wgts)
	if (!identical(i, integer(0))){
	used <- wgts[i]
	indices <- append(indices, i)
	len <- length(used)
	sum <- 0
	for (u in used){
		sum <- sum + eval(parse(text=u))
	}
	assign(paste0("pred.wgt.", dataset), sum/len, envir = parent.frame())
	result_wgts <- append(result_wgts, paste0("pred.wgt.", dataset))
	}
	assign("new_wgts", result_wgts, envir=parent.frame())
	assign("index", indices, envir = parent.frame())
}

cells_datasets <- c("ota")
for (c in cells_datasets){
	meta <- meta_cells(c, new_wgts, index)
}


if (length(index) > 0){
	new_wgts <- append(new_wgts, wgts[-index])
	wgts = new_wgts
}

ext <- length(wgts)

# --- COMPUTE TOTAL SAMPLE SIZE -> USED LATER IN META-ANALYSIS

total <- 0
for (w in wgts){
	dataset <- strsplit(w, split="[.]")[[1]][3]
	size <- h[[dataset]]
	total = total + size
}
total = total + 80


# --- SPLIT DATASETS 
lasso_h2 <- hsq_afr[1]
if(lasso_h2 < 0){lasso_h2 <- 0.064251}
afronly = weights.lasso( geno.file , lasso_h2 , snp=genos$bim[,2] )	
if ( sum(is.na(afronly)) == length(afronly)) {
	afronly = weights.marginal( genos$bed , as.matrix(pheno[,3]) , beta=T )
	#distribution of betas, take Z > 1.96
	z_scores <- scale(afronly)
	top <- which(abs(z_scores) > 1.96)
	#if none of the snps are Z > 1.96	
	if (identical(integer(0), top)){
		threshold = max(round(length(afronly)/4), 1)
		top <- order(afronly^2, decreasing = TRUE)[1:threshold]
	}
	temp <- rep(0, length(afronly))
	temp[top] <- afronly[top]
	afronly = temp
}
nonzero <- which(afronly != 0)
zero <- which(afronly==0)
groups <- c()
if (!identical(integer(0), nonzero)){
	groups <- append(groups, "nonzero")
}
if (!identical(integer(0), zero)){
	groups <- append(groups, "zero")
}

wgt2 <- c()

for (w in wgts){
	for (g in groups){
		vec <- eval(parse(text = w))
		vec[-eval(parse(text = g))] <- 0
		assign(paste0(w, ".", g), vec, envir = parent.frame())
		wgt2 <- append(wgt2, paste0(w, ".", g))
	}	
}	

ext2 <- length(wgt2)


if ( opt$verbose >= 2){
	cat("External summary stats processed \n", "Printing the datasets being used: \n")
	print(wgts)
}
#------------------------------------------------------------


# --- Cross-Validation
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


	# lasso_h2 defined when we split datasets
	pred.wgt = weights.lasso( cv.file , lasso_h2 , snp=genos$bim[,2] )
	
	if ( sum(is.na(pred.wgt)) == length(pred.wgt)) {
		if ( opt$verbose >= 1 ) cat("LASSO pushed all weights to 0, using top1 as backup \n")
		pred.wgt = weights.marginal( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3,drop=F]) , beta=T )
		pred.wgt[ - which.max( pred.wgt^2 ) ] = 0
	}
	# predict from weights into sample

	# when there is only 1 snp in the wgt file, we need to reformat	
	if (length(pred.wgt) == 1){
		pred.wgt <- t(pred.wgt)
	}
		
	#AFR ONLY------------------------------------------------------------------
	cv.calls[ indx , 1 ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt
	#cv.calls is a matrix of predicted value - 1 fold at a time
	#pred.wgt - training on other 4 folds
	#indx = people in 1 fold we removed when we trained using 4 fold
	
	pred_train_afr = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], ] %*% pred.wgt))) #r^2 between predicted and actual 
	r2_training_afr = append(r2_training_afr, pred_train_afr$adj.r.sq)
		
	#--------------------------------------------------------------------------

	#META-ANALYSIS-------------------------------------------------------------
	#sample-size weighted meta-analysis 
	#calculate total sample size 
	pred.wgt.meta <- pred.wgt * (80/total)
	#add sample size weighted contribution 
	for (w in wgts){
		dataset <- strsplit(w, split="[.]")[[1]][3]
		size <- h[[dataset]]
		pred.wgt.meta = pred.wgt.meta + (eval(parse(text=w)) * (size/total))
	}
	if (length(pred.wgt.meta) == 1){
		pred.wgt.meta <- t(pred.wgt.meta)
	}

	cv.calls[ indx , 2 ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt.meta 
	#-----------------------------------------------------------------------------

	#MAGEPRO----------------------------------------------------------------------
	#create a matrix of weights from all data sets - ridge regression input
	w1 <- genos$bed[ cv.sample[ -indx ] , ] %*% pred.wgt
	eq <- matrix(0, nrow = length (w1), ncol = ext2+1)
	eq[,1] <- w1
	num = 2
	for (w in wgt2){
		eq[,num] <- genos$bed[ cv.sample[ -indx ] , ] %*%  eval(parse(text = w))
		num = num+1
	}	

	#check if any columns are all 0	- situation when there is nonzero Beta but genotype at that snp has no variance 
	remove <- c()
	for (i in 1:ncol(eq)){
		if (length(which(eq[,i] == 0)) == nrow(eq)){
			remove <- append(remove, i)
		}
	}
	if (length(remove) > 0){
		eq <- eq[, -remove]
		ext2 <- ext2 - length(remove)
		wgt2 <- wgt2[-(remove-1)]
	}

	#cv.glmnet runs if there are two or more columns
	if (ext2 > 0){
		y <- cv.glmnet(x = eq , y = cv.train[,3], alpha = 0, nfold = 5, intercept = T, standardize = T)
		cf = coef(y, s = "lambda.min")[2:(ext2+2)]
		predtext <- "cf[1]*pred.wgt"
		for (i in 2:(length(cf))){
			predtext <- paste0(predtext, " + cf[", i, "]*", wgt2[(i-1)])
		}
		pred.wgt.magepro <- eval(parse(text = predtext))
	}else{
		pred.wgt.magepro <- pred.wgt #magepro weight defaults to afronly if no external datasets available
		cf = NA
	}
	
	if(length(pred.wgt.magepro) == 1){
		pred.wgt.magepro <- t(pred.wgt.magepro)
	}	

	cv.calls[ indx , 3 ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt.magepro

	#store the r2 on training set
	pred_train = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], ] %*% pred.wgt.magepro))) #r^2 between predicted and actual 
	r2_training_magepro = append(r2_training_magepro, pred_train$adj.r.sq)
		

	#-------------------------------------------------------------------------------

}

if (opt$verbose >= 1) cat("Cross-validation complete \n")

#uses 80% of data to predict into 20%, then uses another 80% to predict into the next 20%, and so on. 
#after everyone has a prediction, compute rsq and p value. 

#compute rsq + P-value for each model
cv.performance = matrix(NA,nrow=2,ncol=3)  
rownames(cv.performance) = c("rsq","pval") 
types <- c("afr","meta","magepro") 
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

# ---- Full Analysis (all 80 individuals)  
if (opt$verbose >= 1) cat("Computing full-sample weights \n")

wgt.matrix = matrix(0, nrow=nrow(genos$bim), ncol=3) 
colnames(wgt.matrix) = c("AFRONLY","META","MAGEPRO")
rownames(wgt.matrix) = genos$bim[,2]

pred.wgt2 = weights.lasso( geno.file , lasso_h2 , snp=genos$bim[,2] )	
if ( sum(is.na(pred.wgt2)) == length(pred.wgt2)) {
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
if(length(pred.wgt2) == 1){
	pred.wgt2 <- t(pred.wgt2)
}	
w1 <- genos$bed %*% pred.wgt2
eq2 <- matrix(0, nrow = length (w1), ncol = ext2 + 1)
eq2[, 1] <- w1
num = 2
for (w in wgt2){	
	eq2[,num] <- genos$bed %*% eval(parse(text = w))
	num = num+1
}	

#store alpha*beta for all snps per dataset (for downstream analysis of which datasets contribute the most)
alpha_beta <- matrix(0, nrow = length(pred.wgt2), ncol = ext2 + 1)

#run ridge regression to find optimal coefficients and compute multipop weight
if (ext2 > 0){
	y2 <- cv.glmnet(x = eq2 , y = pheno[,3], alpha = 0, nfold = 5, intercept = T, standardize = T)
	cf_total = coef(y2, s = "lambda.min")[2:(ext2+2)]	
	total_coeff <- cf_total #saved for down-stream analysis
	predtext2 <- "cf_total[1]*pred.wgt2"
	alpha_beta[, 1] <- eval(parse(text = predtext2))
	for (i in 2:(length(cf_total))){
		predtext2 <- paste0(predtext2, " + cf_total[", i, "]*", wgt2[(i-1)])
		alpha_beta[, i] <- cf_total[i]*eval(parse(text = wgt2[(i-1)]))
	}	
	pred.wgt.multipop <- eval(parse(text = predtext2))
}else{
	pred.wgt.multipop <- pred.wgt2
	total_coeff = c(NA)
	alpha_beta = NA
}		

wgt.matrix[, 3] = pred.wgt.multipop

#---

snps = genos$bim

#take average of r2 on training set
avg_training_r2_magepro <- mean(r2_training_magepro)
avg_training_r2_afr <- mean(r2_training_afr)

#add african weight to the list of weights used
wgt2 <- append("pred.wgt", wgt2)

if (is.matrix(alpha_beta)){
	colnames(alpha_beta) <- wgt2
}

#if the alpha is 0 for any dataset, it did not add anything to the model
if (sum(is.na(total_coeff)) == 0){
w <- which(total_coeff == 0)
#remove that dataset from the list of datasets used in the model
if (length(w) > 0 ) {
wgt2 <- wgt2[-w]
total_coeff <- total_coeff[-w]
alpha_beta <- alpha_beta[, -w]
}
}

save( wgt.matrix, alpha_beta, snps, cv.performance , hsq_afr, hsq_afr.pv, hsq_eur, hsq_eur.pv, N.tot , wgt2, total_coeff, avg_training_r2_magepro, avg_training_r2_afr, file = paste( opt$out , ".wgt.RDat" , sep='' ) )

# --- CLEAN-UP
if ( opt$verbose >= 1 ) cat("Cleaning up\n")
cleanup()
}else{
	cat("Gene is not cis heritable in all of the populations")
}
