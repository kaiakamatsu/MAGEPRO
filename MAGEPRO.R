# load packages
suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('glmnet'))
suppressMessages(library('data.table'))
suppressMessages(library('dplyr'))
suppressMessages(library('susieR'))

source(('fine_mapping.R'), chdir = TRUE)


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
	      SINGLE = single ancestry FUSION lasso approach \n
	      META = ss-weighted meta-analysis of datasets \n 
	      PT = pruning and thresholding \n
	      SuSiE = sum of single effects regression \n
	      SuSiE_IMPACT = sum of single effect regression with IMPACT scores as prior on SNP selection \n
	      PRSCSx = PRS-CSx multi-ancestry PRS method \n 
	      MAGEPRO_fullsumstats = magepro model, no sparsity \n
	      MAGEPRO = magepro model"),
  make_option("--ss", action="store", default=NA, type='character',
              help="Comma-separated list of sample sizes of sumstats (in the same order)"), 
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
  make_option("--impact_path", action="store", default=NA, type='character',
              help="path to file with impact scores for each snp"), 
  make_option("--ldref_dir", action="store", default=NA, type="character",
  			  help="Directory containing ld reference files used for susie fine mapping"),
  make_option("--ldrefs", action="store", default=NA, type="character",
  			  help="Comma-separated list of ld reference files (plink prefixes) used for susie fine mapping"),
  make_option("--cl_thresh", action="store", default=0.97, type="numeric",
  			  help="Clumping threshold for plink to clump SNPs in summary statistics that are in high LD before running SuSiE [optional]"), 
  make_option("--out_susie", action="store", default=NA, type='character',
              help="Path to susie output directory [required if using MAGEPRO and not skipping susie]"), 
  make_option("--skip_susie", action="store_true", default=FALSE,
              help="Boolean to skip SuSiE preprocessing. This assumes summary statistics in sumstats_dir have columns 8/9/10 with PIP/POSTERIOR/CS from susie")
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

weights.meta = function(target_model, target_ss, wgts, h, total_ss) {
	# PURPOSE: ss-weighted meta-analysis of effect sizes from datasets and target
	# target_model = vector of target population effect sizes
	# target_ss = sample size used to compute the target_model
	# wgts = vector of the names (string) of external dataset vectors to use 
	# h = hashmap that maps the dataset name to sample size 
	# total_ss = total sample size of all datasets and target_ss

	pred.wgt.meta <- target_model * (target_ss/total_ss)
	for (w in wgts){
		dataset <- strsplit(w, split="[.]")[[1]][3]
		size <- h[[dataset]]
		pred.wgt.meta = pred.wgt.meta + (eval(parse(text=w)) * (size/total_ss))
	}
	if (length(pred.wgt.meta) == 1){
		pred.wgt.meta <- t(pred.wgt.meta) 
	}
	return (pred.wgt.meta)
}


tune.pt = function(geno_bed, geno_bim, pheno_df, r2_choices, p_choices, cross, ldref, tmp){
	# PURPOSE: perform 5-fold cross validation on individuals in geno-bed and pheno to fine the best r2 and p parameter out of given choices 
	# geno_bed = genotype matrix (row = people, col = snps)
	# geno_bim = bim file for target population (snps correspond to geno_bed)
	# pheno_df = phenotype data for individuals in geno_bed (dataframe with phenotype in third column)
	# r2_choices = vector of r2 thresholds to search over 
	# p_choices = vector of p value thresholds to search over 
	# cross = number of cross validation folds
	# ldref = ld reference file for pruning and thresholding 
	# tmp = path/name for temporary files	
	# RETURN: c(r2_best, p_best): a vector holding the best r2 and p value threshold 
	# 1. create a matrix of size |r2_choices| x |p_choices|
		# holds the r-squared of the model created using i_th r2 and j_th p 
	performance = matrix(NA, nrow=length(r2_choices), ncol=length(p_choices))
	rownames(performance) <- r2_choices
	colnames(performance) <- p_choices
	# 2. define indices of cross validation 
	all_phenos = pheno_df
	sample_size = nrow(all_phenos)
	index_sample = sample(sample_size) #sample randomly 
	all_phenos = all_phenos[index_sample,]
	folds = cut(seq(1,sample_size),breaks=cross,labels=FALSE)
	# 3. for every fold
	calls = matrix(NA,nrow=sample_size,ncol=(length(r2_choices)*length(p_choices))) # reset calls for every r2/p pair 
	for ( fold in 1:cross ) {
		fold_indx = which(folds==fold,arr.ind=TRUE)
		train_phenos = all_phenos[-fold_indx,]
		train_phenos[,3] = scale(train_phenos[,3]) 
		geno_bed_train = geno_bed[index_sample[-fold_indx], , drop = FALSE]	
		# 4. compute betas in training fold
		sumstats <- eqtl_assoc(geno_bed_train, geno_bim, train_phenos[,3])	
		# 5. for every r2, p pair 
		colnum = 1
		for (i in 1:length(r2_choices)){
			for (j in 1:length(p_choices)){
				r2_current = r2_choices[i]
				p_current = p_choices[j]
				# P+T
	       			current_pred = prune_thresh(sumstats, geno_bim, r2_current, p_current, ldref, tmp)
				# predict 
				if (sum(is.na(current_pred)) == nrow(geno_bim)){
                        		calls[fold_indx, colnum] = NA
                		}else{
                        		calls[fold_indx, colnum] <- geno_bed[index_sample[fold_indx], , drop = FALSE] %*% current_pred
				}
				colnum = colnum + 1
			}
		}
	}
	# 6. compute r-squared of predictions and store in matrix
	colnum = 1
	for (i in 1:length(r2_choices)){
		for (j in 1:length(p_choices)){
			if ( !is.na(sd(calls[,colnum])) && sd(calls[,colnum]) != 0 ) {
				summ = summary(lm( all_phenos[,3] ~ calls[,colnum] )) #r^2 between predicted and actual
				performance[i,j] = summ$adj.r.sq
                        }else{
                                performance[i,j] = NA
                        }
			colnum = colnum + 1
		}
	}
	# 7. return the r2 and p that results in the maximum r-squared
	max_index <- which(performance == max(performance), arr.ind=TRUE)
	if (nrow(max_index) > 0){
		row_max <- max_index[1,1]
		col_max <- max_index[1,2]
		best_r2_p <- c(r2_choices[row_max], p_choices[col_max])
	}else{
		best_r2_p <- c(0.2,0.5) # if none of the tuned parameters result in any clumps, just try 0.2 and 0.5
	}
	return(best_r2_p)
}

eqtl_assoc = function(geno_bed, geno_bim, pheno) {
	# PURPOSE: compute effect sizes of eQTL
	# geno_bed = genotype matrix (row = people, col = snps)
	# geno_bim = bim file for target population
	# pheno = phenotype of interest across people
	
	Betas <- c()
        Pvals <- c()
        SNPs <- c()
	for(col in 1:ncol(geno_bed)){
                SNPs <- append(SNPs, colnames(geno_bed)[col])
                snp <- geno_bed[,col]
                model <- lm(pheno ~ snp)
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
        sumstats_pt <- data.frame(SNPs, Betas, Pvals)
	return (sumstats_pt)
}

prune_thresh = function(sumstats_pt, geno_bim, r2, p, ldref, tmp) {
	# PURPOSE: run pruning+thresholding on dataframe with rsid, betas, pval
	# sumstats_pt = dataframe with SNPs, Betas, Pvals as columns
	# geno_bim = bim file for target population
	# r2 = LD pruning threshold 
	# p = p-value threshold 
	# ldref = path/file to plink ldref file for P+T
	# tmp = path/name for temporary file
	
	#write dataframe, run plink --clump to P+T, write output to temp directory
        sumstats_file = paste(tmp,"_sumstats.txt",sep='')
        write.table(sumstats_pt, file = sumstats_file, quote = F, row.names = F, col.names = T, sep = '\t')
        arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ", ldref, geno_bim[1,1], " --clump-p1 ", p ," --clump-r2 ", r2 ," --clump-snp-field SNPs --clump-field Pvals --clump ", tmp, "_sumstats.txt"," --out ", sumstats_file,sep='')
        system(arg, ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)	

	#read in result file, create wgt matrix and predict
        if (file.exists(paste0(sumstats_file, ".clumped"))){
        	clumped <- fread(file = paste0(sumstats_file, ".clumped"), header = T)
		clumps <- clumped[[3]]
		sumstats_pt <- sumstats_pt[which(sumstats_pt$SNPs %in% clumps),]
		m_sumstats <- match(geno_bim[,2], sumstats_pt[[1]])
                sumstats_pt <- sumstats_pt[m_sumstats,]
		w_sumstats <- which(is.na(sumstats_pt[[1]]))
                if(length(w_sumstats) > 0){
                        sumstats_pt[w_sumstats,2] <- 0
                }
                pred.wgt.PT_sumstats = sumstats_pt[[2]]
        }else{
                pred.wgt.PT_sumstats = rep(NA, times = nrow(geno_bim))
        }
        if(length(pred.wgt.PT_sumstats) == 1){
               pred.wgt.PT_sumstats <- t(pred.wgt.PT_sumstats)
        }
	system(paste0("rm -rf ", sumstats_file, ".clumped"))
	return (pred.wgt.PT_sumstats)

}

weights.pt = function(geno_bed, geno_bim, pheno, r2, p, ldref, tmp) {
	# PURPOSE: compute effect sizes of eQTL and run pruning+thresholding
	# geno_bed = genotype matrix (row = people, col = snps)
	# geno_bim = bim file for target population
	# pheno = phenotype of interest across people
	# r2 = LD pruning threshold 
	# p = p-value threshold 
	# ldref = path/file to plink ldref file for P+T
	# tmp = path/name for temporary file
	sumstats_pt <- eqtl_assoc(geno_bed, geno_bim, pheno)
	pred.wgt.PT_sumstats <- prune_thresh(sumstats_pt, geno_bim, r2, p, ldref, tmp)
	return (pred.wgt.PT_sumstats)
}

weights.susie = function(geno, pheno){
	# PURPOSE: fit the sum of single effects model and return coefficients for prediction
	# geno = nxp genotype matrix 
	# pheno = n-length vector of phenotype 
	# RETURN: p-length vector of effect sizes

	fitted <- susie(geno, pheno)
	coefs <- coef(fitted)[-1]
	pred.wgt.susie <- coefs
	return(pred.wgt.susie)

}


weights.susie_impact = function(geno, pheno, impact_file){
	# PURPOSE: fit the sum of single effects model and return coefficients for prediction
	# geno = nxp genotype matrix 
	# pheno = n-length vector of phenotype 
	# impact_file = path to file with impact scores for every SNP
	# RETURN: p-length vector of effect sizes
	# use IMPACT scores as priors for SuSiE
	
	snps <- colnames(geno)

	# read in impact scores, take snps in geno 
	impact <- fread(impact_file)
	impact_ordered <- impact %>% 
		filter(SNP %in% snps) %>%
		arrange(match(SNP, snps))

	fitted <- susie(geno, pheno, prior_weights = impact_ordered$means)
	coefs <- coef(fitted)[-1]
	pred.wgt.susie <- coefs
	return(pred.wgt.susie)

}


shrinkage.prscsx = function(exec_dir, ldref, working_dir, shrinkage, input, ss, pp, bim) {
	# PURPOSE: run PRS-CSx shrinkage and retrieve list of weights
	# exec_dir: directory where PRScsx.py executable is 
	# ldref: directory where ldref for PRS-CSx is (see PRS-CSx github for preparation)
	# working_dir: directory where external summary statistics files for PRS-CSx and snps list is written 
	# input: comma separated sequence of path/file of summary statistics for PRS-CSx
	# ss: comma separated sequence of sample sizes of summary statistics 
	# pp: comma separated sequence of ancestries of summary statistics
	# bim: bim file of target population 
	
	write.table(bim, file = paste0(working_dir, "snps.bim"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
	arg = paste0("python ", exec_dir, "/PRScsx.py --ref_dir=", ldref, " --bim_prefix=", working_dir, "snps", " --sst_file=", input," --n_gwas=", ss, " --pop=", pp, " --chrom=", bim[1,1]," --phi=", shrinkage," --out_dir=", working_dir, " --out_name=results")
	if ( opt$verbose >= 1 ) print(arg)
	system(arg)
	
	wgts <- c()
	unique_pp <- unique(strsplit(pp, split = ",")[[1]])
	for (n in unique_pp){
		wgts <- weights_process_prscsx(n, paste0(working_dir,"results_", n, "_pst_eff_a1_b0.5_phi", sub(".00e", "e", sprintf("%.2e", shrinkage)) ,"_chr", bim[1,1], ".txt"), wgts, bim)
	}
	return (wgts)

}

weights_process_prscsx <- function(dataset, file, wgts, bim){
        if (file.size(file) != 0L){
                table <- fread(file)
                m <- match(bim[,2], table[[2]])
                table <- table[m,]
                w <- which(is.na(table[[2]]))
                if(length(w) > 0){
                        table[w,6] <- 0
                }
                assign(paste0("pred.wgtprscsx.", dataset), table[[6]], envir = .GlobalEnv)
                if(sum(which(eval(parse(text = paste0("pred.wgtprscsx.", dataset))) != 0)) > 0){
                        wgts <- append(wgts, paste0("pred.wgtprscsx.", dataset))
                        return (wgts)
                }else {
                        return (wgts)
                }
        }else{
                return (wgts)
        }
}

weights.prscsx = function(wgts, geno, pheno, target_model = NA){
	# PURPOSE: compute optimal linear combination of post-shrinkage weights from prscsx
	# wgts: list of strings, name of variables holding post-shrinkage weights 
	# geno: genotype matrix of target cohort 
	# pheno: vector of phenotypes of target cohort, for individuals in geno
	# target_model: name of variable holding gene model for target population (if provided, it is added as one of the features in the linear combination)
	ext <- length(wgts)
	eq = list()
	if (!is.na(target_model)){
		eq <- append(eq, paste0("geno %*% ", target_model))
	}
        for (w in wgts){
                if (length(eval(parse(text = w))) == 1){
                       assign(w, t(eval(parse(text = w))))
                }
                eq <- append(eq, paste0("geno %*% ", w))
        }
	eq <- paste(eq, collapse="+")
        eq <- paste0("pheno ~ ", eq)
	y <- lm(eval(parse(text = eq)))
	coef <- coef(y)
	coef <- ifelse(is.na(coef), 0, coef)
	w_eq = list()
        if (!is.na(target_model)){
		w_eq <- append(w_eq, paste0("coef[", 2,"]*",target_model))
		for (nums in 3:(ext+2)){
                	w_eq <- append(w_eq, paste0("coef[", nums,"]*",wgts[nums-2]))
        	}
	}else{
		for (nums in 2:(ext+1)){
                	w_eq <- append(w_eq, paste0("coef[", nums,"]*",wgts[nums-1]))
        	}
	}
        w_eq <- paste(w_eq, collapse="+")
        pred.wgt.prs_csx <- eval(parse(text=w_eq))
        if (length(pred.wgt.prs_csx) == 1){
                pred.wgt.prs_csx <- t(pred.wgt.prs_csx)
        }
	return(pred.wgt.prs_csx)
}

weights.magepro_marquezluna = function(basemodel, wgts, geno, pheno_df, geno_bim, tmp, lasso_h2, save_alphas) {
	# PURPOSE: compute MAGEPRO weights, use MarquezLuna approach to compute Betas and alphas in the same training fold
	# basemodel = target population prediction model, a vector of effect sizes (this base model will be one of the features in the regression)
	# wgts = list of name of variables holding magepro_processed weights 
	# geno = genotype matrix for training mixing weights 
	# pheno_df = data frame with gene expression values for individuals in geno in column 3 (also for training mixing weights)
	# geno_bim = bim file of target population
	# tmp = directory/name for temporary files
	# lasso_h2 = heritability of gene, used for plink --lasso
	# save_alphas = T/F, save alpha coefficients to 'cf_total'?
	ext <- length(wgts)
	#1a. format glmnet input
	eq <- matrix(0, nrow = nrow(geno), ncol = ext+1)
	# 10-fold cv to compute betas and "predict" on all people
	all_phenos <- pheno_df
	sample_size <- nrow(all_phenos)
	index_sample <- sample(sample_size) #sample randomly 
	all_phenos <- all_phenos[index_sample,]
	folds <- cut(seq(1,sample_size),breaks=10,labels=FALSE)
	for ( fold in 1:10 ) {
		fold_indx = which(folds==fold,arr.ind=TRUE)
		train_phenos = all_phenos[-fold_indx,]
		train_phenos[,3] = scale(train_phenos[,3]) 
		geno_bed_train = geno[index_sample[-fold_indx], , drop = FALSE]	
		train.file = paste(tmp,".cv_nested",sep='')
		write.table( train_phenos , quote=F , row.names=F , col.names=F , file=paste(train.file,".keep",sep=''))
		arg = paste(opt$PATH_plink ," --allow-no-sex --bfile ",tmp," --keep ",train.file,".keep --out ",train.file," --make-bed",sep='')
		system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
		pred.wgt_current = weights.lasso( train.file , lasso_h2 , snp=geno_bim[,2] )
		if ( sum(is.na(pred.wgt_current)) == nrow(geno_bim)) {
			pred.wgt_current = weights.marginal( geno_bed_train , as.matrix(train_phenos[,3,drop=F]) , beta=T )
			pred.wgt_current[ - which.max( pred.wgt_current^2 ) ] = 0
		}
		if (length(pred.wgt_current) == 1){
			pred.wgt_current <- t(pred.wgt_current) # 1 snp in the cis window -> transpose for "matrix" mult
		}
		if (sum(is.na(pred.wgt_current)) == nrow(geno_bim)){
                        eq[index_sample[fold_indx], 1] = NA
                }else{
			eq[index_sample[fold_indx], 1] <- geno[index_sample[fold_indx], , drop = FALSE] %*% pred.wgt_current
		}
	}
	for (c in 1:length(wgts)){
		eq[,(c+1)] <- geno %*%  eval(parse(text = wgts[c])) 
	}	
	#1b. when geno %*% wgt -> all 0 (snps at nonzero wgt have no variaion among people) -> remove sumstat
	zero_cols <- colSums(eq == 0) == nrow(eq) 
	if (any(zero_cols)) {
  		eq <- eq[, !zero_cols]
  		ext <- ext - sum(zero_cols)
		wgts <- wgts[-(which(zero_cols) - 1)]
	}
	#2. run ridge regression to find optimal coefficients for each dataset
	y <- cv.glmnet(x = eq , y = pheno_df[,3], alpha = 0, nfold = 5, intercept = T, standardize = T)
	cf <- coef(y, s = "lambda.min")[2:(ext+2)]
	predtext <- "cf[1]*basemodel"
	for (i in 2:(length(cf))){
		predtext <- paste0(predtext, " + cf[", i, "]*", wgts[(i-1)])
	}
	pred.wgt.magepro <- eval(parse(text = predtext))
	if(length(pred.wgt.magepro) == 1){
		pred.wgt.magepro <- t(pred.wgt.magepro)
	}
	if(save_alphas) cf_total <<- cf
	return(pred.wgt.magepro)
}

weights.magepro = function(basemodel, wgts, geno, pheno, save_alphas) {
	# PURPOSE: compute MAGEPRO weights
	# basemodel = target population prediction model, a vector of effect sizes (this base model will be one of the features in the regression)
	# wgts = list of name of variables holding magepro_processed weights 
	# geno = genotype matrix for training mixing weights 
	# pheno = vector of gene expression values for individuals in geno (also for training mixing weights)
	# save_alphas = T/F, save alpha coefficients to 'cf_total'?
	ext <- length(wgts)
	#1. format glmnet input
	eq <- matrix(0, nrow = nrow(geno), ncol = ext+1)
	eq[,1] <- geno %*% basemodel
	
	for (c in 1:length(wgts)){
		eq[,(c+1)] <- geno %*%  eval(parse(text = wgts[c])) 
	}	
	#2. run ridge regression to find optimal coefficients for each dataset
	y <- cv.glmnet(x = eq , y = pheno, alpha = 0, nfold = 5, intercept = T, standardize = T)
	cf = coef(y, s = "lambda.min")[2:(ext+2)]
	predtext <- "cf[1]*basemodel"
	for (i in 2:(length(cf))){
		predtext <- paste0(predtext, " + cf[", i, "]*", wgts[(i-1)])
	}
	pred.wgt.magepro <- eval(parse(text = predtext))
	if(length(pred.wgt.magepro) == 1){
		pred.wgt.magepro <- t(pred.wgt.magepro)
	}
	if(save_alphas) cf_total <<- cf
	return(pred.wgt.magepro)
}


if ( opt$verbose == 2 ) cat("Predictive Models Prepared \n")

# --- CLEANUP
cleanup = function() {
	# PURPOSE: clean up temporary directory
	# RETURN: 
	if ( ! opt$noclean ) {
		arg = paste("rm -rf " , opt$tmp , "*", sep='')
		system(arg)
	}
}

# --- LOAD SUMSTATS AND FIX FLIPPED ALLELES
load_flip_dataset <- function(bim, dataset, file, datasets_loaded, select, susie = FALSE){
	# PURPOSE: read summary statistics, flip signs according to effect allele 
	# bim = bim file of snps
	# dataset = name of sumstat dataset
	# file = file path to sumstat (make sure it is properly formatted)
	# datasets_loaded = current list of loaded datasets
	# select = columns numbers to extract from file 
	# susie = true/false: using susie posterior weights? 
	# RETURN: updated list of loaded datasets, including the one processed here
	table <- fread(file, select = select) 
	for (k in 1:nrow(table)){
		match=which(bim$V2 == table[k, 1])
		if (! identical(match, integer(0))){
			if (!is.na(table[k, 2]) && !is.na(table[k, 3])){
				if (as.character(bim$V5[match]) != as.character(table[k, 2]) || as.character(bim$V6[match]) != as.character(table[k, 3])) { 
					if (as.character(bim$V5[match]) == as.character(table[k, 3]) && as.character(bim$V6[match]) == as.character(table[k, 2])){
						table[k, 4] <- table[k, 4] * -1
						if (susie){
							table[k, 7] <- table[k, 7] * -1
						}
					}else{
						table[k, 4] <- 0
						if (susie){
							table[k, 7] <- 0
						}
					}
				}
			} else {
				table[k, 4] <- 0
				if (susie){
                                	table[k, 7] <- 0
                                }
			}
		}
	}
	m <- match(bim[,2], table[[1]])
	table <- table[m,]
	w <- which(is.na(table[[1]]))
	if(length(w) > 0){
		table[w,4] <- 0
		if (susie){			
			table[w,7] <- 0
		}
	}
	assign(paste0("loaded.", dataset), table, envir = .GlobalEnv)
	datasets_loaded <- append(datasets_loaded, paste0("loaded.", dataset))
	return(datasets_loaded)
}

# --- PROCESS-SUMSTATS
datasets_process_fullsumstats <- function(dataset, loaded, wgts){
	# PURPOSE: process effect sizes from full summary statistics
	# dataset = name of sumstat dataset
	# loaded = dataframe containing summary statistics from this dataset 
		# SNP(rsid) A1(allele1) A2(allele2) b(effect size) p(p-value) pip(pip from SuSiE) cs(credible set membersip from SuSiE) 
	# wgts = current list of the names of the variables containing these wgts
	# RETURN: updated list of wgts, including the one that was processed in this function
	assign(paste0("pred.wgt.", dataset), loaded[[4]], envir = .GlobalEnv)
	if(sum(which(eval(parse(text = paste0("pred.wgt.", dataset))) != 0)) > 0){
		wgts <- append(wgts, paste0("pred.wgt.", dataset))
		return (wgts)
	}else {
		return (wgts)
	}
}

# --- PROCESS-SUMSTATS-USE-SUSIE
datasets_process_susie <- function(dataset, loaded, wgts){
	# PURPOSE: process sumstats using SuSiE results: ONLY keep nonzero effect sizes for high PIP SNPs. 
	# dataset = name of sumstat dataset
	# loaded = dataframe containing summary statistics from this dataset
		# SNP(rsid) A1(allele1) A2(allele2) b(effect size) p(p-value) pip(pip from SuSiE) b(coef from SuSiE) cs(credible set membersip from SuSiE)
	# wgts = current list of the names of the variables containing these wgts
	# RETURN: updated list of wgts including the one that was processed in this function

	# STRATEGY:
		# SuSiE coefficients for all snps

	assign(paste0("pred.wgtsusie.", dataset), loaded[[7]], envir = .GlobalEnv)
	if(sum(which(eval(parse(text = paste0("pred.wgtsusie.", dataset))) != 0)) > 0){
		wgts <- append(wgts, paste0("pred.wgtsusie.", dataset))
		return (wgts)
	}else {
		return (wgts)
	}
}

datasets_process_prscsx <- function(name, loaded, working_dir){
	# PURPOSE = write summary statistics files to prepare to run PRS-CSx
	# name = name of the dataset
	# loaded = dataframe containing summary statistics from this dataset
		# SNP(rsid) A1(allele1) A2(allele2) b(effect size) p(p-value) pip(pip from SuSiE) cs(credible set membersip from SuSiE)
	# working_dir = temporary directory to write intermediate files for running PRS-CSx
	df <- data.frame(loaded[,c(1,2,3,4,5)])
	colnames(df) <- c("SNP", "A1", "A2", "BETA", "P")
	f_name <- paste0(working_dir, name, ".txt")
	write.table(df, file = f_name, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
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
if ( sum(! model %in% c("SINGLE", "META", "PT", "SuSiE", "SuSiE_IMPACT","PRSCSx", "MAGEPRO_fullsumstats", "MAGEPRO")) > 0 | length(model) > 8 ){
	cat( "ERROR: Please input valid models \n" , sep='', file=stderr() )
        cleanup()
        q()
}

if (opt$verbose >= 1) cat("### USING THE FOLLOWING MODELS:", opt$models, "\n")

order <- c("SINGLE","META", "PT", "SuSiE", "SuSiE_IMPACT", "PRSCSx", "MAGEPRO_fullsumstats", "MAGEPRO")
types <- model[match(order, model)]
types <- unique(types[!is.na(types)])

#list of datasets to be used
datasets <- list()

# sumstats to use
if ( ! is.na(opt$sumstats)){
	if ( is.na(opt$sumstats_dir) ){
		cat( "ERROR: --sumstats supplied, but not --sumstats_dir \n" , sep='', file=stderr() )
                cleanup()
                q()
	}
	sumstats <- strsplit(opt$sumstats, ",", fixed = TRUE)[[1]]
	for (s in sumstats){
		file_dir = paste0(opt$sumstats_dir, "/", s)
		if (dir.exists(file_dir)){
			file_name = paste0(file_dir, "/", name, ".txt")
			if (file.exists(file_name)){
				assign(paste0("file.", s), file_name, envir = .GlobalEnv)
				datasets <- append(datasets, paste0("file.", s))
			}else{
				if ( opt$verbose == 2 ) {
					cat("skipping sumstat", s, "for this gene\n")
				}
			}
		}else{ #if there is no directory of eqtl sumstats in the --sumstats_dir path
			cat( "ERROR: make sure --sumstats_dir has a directory named ", s ," that contains gene-specific files of eqtl data\n" , sep='', file=stderr() )
			cleanup()
			q()
		}
	}
}else{
	if ("META" %in% model | "PRSCSx" %in% model | "MAGEPRO_fullsumstats" %in% model | "MAGEPRO" %in% model){
		cat( "ERROR: --sumstats not supplied, cannot compute META, PRS-CSx, MAGEPRO models \n" , sep='', file=stderr() )
		cleanup()
                q()
	}
}

hashmap_ss <- list() #hashmap for sample sizes per dataset
if ( ("META" %in% model | "PRSCSx" %in% model |  ("MAGEPRO" %in% model & !opt$skip_susie) ) ){
	if (!is.na(opt$ss)){
		sample_sizes <- strsplit(opt$ss, ",", fixed = TRUE)[[1]]
		if (length(sumstats) != length(sample_sizes)){
			cat( "ERROR: --ss flag required an entry for every dataset\n" , sep='', file=stderr() )
                	cleanup()
                	q()
		}
		hashmap_ss <<- setNames(as.numeric(sample_sizes), sumstats)
	}else{
		cat( "ERROR: Cannot perform sample-size weighted meta-analysis or PRS-CSx without the --ss flag\n" , sep='', file=stderr() )
                cleanup()
                q()
	}
}

cohort_map <- list()
if ("MAGEPRO" %in% model & (!opt$skip_susie) ) {
	if ( !is.na(opt$ldref_dir) ) {
		if (!file.exists(opt$ldref_dir)) {
			cat( "ERROR: --ldref_dir directory does not exist, please check the path\n", sep='', file=stderr())
			cleanup()
			q()
		}
	} else {
		cat( "ERROR: --ldref_dir not supplied, cannot perform fine mapping of sumstats for MAGEPRO model\n", sep='', file=stderr() )
				cleanup()
				q()
	}

	if (!is.na(opt$ldrefs)) {
		ldrefs_list <- strsplit(opt$ldrefs, ",")[[1]]
		if (length(sumstats) != length(ldrefs_list)) {
			cat ("ERROR: --ldrefs flag requires an entry for every dataset\n", sep='', file=stderr())
				cleanup()
				q()
		}
		# create map with cohort as keys 
		cohort_map <- setNames(
		lapply(seq_along(sumstats), function(i) {
			list(sample_size = sample_sizes[i], ldref = ldrefs_list[i])
		}),
		sumstats
		)
	} else {
		cat( "ERROR: --ldrefs not supplied, cannot perform fine mapping for MAGEPRO model\n", sep='', file=stderr() )
				cleanup()
				q()
	}

	if (!is.na(opt$out_susie)) {
		if (!file.exists(opt$out_susie)) {
			cat( "ERROR: --out_susie directory does not exist, please check the path\n", sep='', file=stderr())
			cleanup()
			q()
		}
	} else {
		cat( "ERROR: --out_susie not supplied, necessary to store susie outputs from MAGEPRO\n", sep='', file=stderr() )
		cleanup()
		q()
	}
	
}


if ( "PT" %in% model ){

	if (is.na(opt$ldref_pt)){
		cat( "ERROR: Cannot perform PT without ld reference file (--ldref_pt) \n" , sep='', file=stderr() )
                cleanup()
                q()
	}

	if ( !is.na(opt$prune_r2) ){
		if ( opt$prune_r2 < 0 | opt$prune_r2 > 1 ){
                	cat( "ERROR: --prune_r2 has to be between 0 and 1 \n" , sep='', file=stderr() )
                	cleanup()
                	q()
        	}
	}

	if ( !is.na(opt$threshold_p) ){
		if ( opt$threshold_p < 0 | opt$threshold_p > 1 ){
			cat( "ERROR: --threshold_p has to be between 0 and 1 \n" , sep='', file=stderr() )
                	cleanup()
                	q()
		}	
	}

}

pops <- list()
if ( "PRSCSx" %in% model ){

	if (opt$verbose >= 1) cat( "Using PRS-CSx model. Please make sure all dependencies for PRS-CSx are installed (see our Github or PRS-CSx Github) \n " )

	if (is.na(opt$ldref_PRSCSx)){
        	cat( "ERROR: input ldref directory for PRS-CSx required (check PRS-CSx github to learn how to download and prepare the ldref files) \n" , sep='', file=stderr() )
        	cleanup()
        	q()
	}

	if (!is.na(opt$pops)){
		populations <- strsplit(opt$pops, ",", fixed = TRUE)[[1]]
		if (length(sumstats) != length(populations)){
			cat( "ERROR: --pops flag required an entry for every dataset to run PRSCSx\n" , sep='', file=stderr() )
			cleanup()
                        q()
                }
                for (i in 1:length(sumstats)){
                        pops[[sumstats[i]]] <- populations[i]
                }
        }else{
                cat( "ERROR: Cannot perform PRS-CSx without the --pops flag\n" , sep='', file=stderr() )
                cleanup()
                q()
        }

	#PRS-CSx working dir
	PRS_CSx_working_dir = paste0(opt$tmp, "PRSCSx/")
	system(paste0("mkdir ", PRS_CSx_working_dir))

}


if ( "SuSiE" %in% model | "SuSiE_IMPACT" %in% model ){
	if ( "SuSiE_IMPACT" %in% model ){
		if(is.na(opt$impact_path)){
			cat( "ERROR: Cannot perform SuSiE_IMPACT without the --impact_path flag\n" , sep='', file=stderr() )
                	cleanup()
                	q()
		}
	}
}


if ( opt$verbose == 2 ) {
	cat("Datasets available for this gene: \n")
	print(datasets)
	if (length(datasets) == 0){
		cat("WARNING: no datasets available for this gene, META, PRS-CSx, and MAGEPRO will have NA weights \n")
	}
}

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
	if ( opt$verbose >= 1 ) cat("### Loaded",ncol(covar)-2,"covariates\n")
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
}

N.tot = nrow(genos$bed)
if ( opt$verbose >= 1 ) cat(nrow(pheno),"phenotyped samples, ",nrow(genos$bed),"genotyped samples, ",ncol(genos$bed)," markers\n")

# --- prepare heritability value to use for LASSO in plink
lasso_h2 <- hsq[1]
if( (lasso_h2 < 0) | (is.na(lasso_h2)) ){
	if ( opt$verbose >= 1 ) cat("forcing lasso heritability to ", opt$lassohsq, " \n")
	lasso_h2 <- opt$lassohsq
}  #when gcta does not converge or yield wild estimates
# --- SETUP SUMSTATS
ext <- length(datasets)

# --- READ SUMSTATS AND FLIP ALLELES AS NECESSARY
loaded_datasets <- c()
if (ext > 0){
	select_cols <- c(2,3,4,5,7) # CAN EDIT THIS LINE WITH CUSTOMIZED COL NUMBERS
	susie <- FALSE
	if ( "MAGEPRO" %in% model ){
		select_cols <- append(select_cols, c(8, 9, 10))
		susie <- TRUE
		if (!opt$skip_susie){
			cohort_fine_mapping(cohort_map, opt$sumstats_dir, opt$tmp, opt$ldref_dir, opt$out_susie, opt$gene, opt$PATH_plink, opt$cl_thresh, opt$verbose)
		}
	}
	for (d in datasets) {
		name <- strsplit(d, split="[.]")[[1]][2]
		loaded_datasets <- load_flip_dataset(genos$bim, name, eval(parse(text = d)), loaded_datasets, select_cols, susie)
	}
}

# --- PREPARE SUMMARY STATISTICS FOR META AND MAGEPRO_fullsumstats
if ( ("META" %in% model | "MAGEPRO_fullsumstats" %in% model)  & (ext > 0) ){

	if ( opt$verbose >= 1){
		cat("### PROCESSING SUMSTATS FOR META AND MAGEPRO_fullsumstats \n")
	}

	wgt_fullsumstats <- c() #sumstats weights before splitting (used for meta-analysis)

	for (loaded in loaded_datasets){
		name <- strsplit(loaded, split="[.]")[[1]][2]
		wgt_fullsumstats <- datasets_process_fullsumstats(name, eval(parse(text = loaded)), wgt_fullsumstats)
	}

	ext_fullsumstats <- length(wgt_fullsumstats)

	# COMPUTE TOTAL SAMPLE SIZE OF SUMSTATS -> USED LATER IN META-ANALYSIS
	total_ss_sumstats <- 0 
	for (w in wgt_fullsumstats){
		dataset <- strsplit(w, split="[.]")[[1]][3]
		total_ss_sumstats = total_ss_sumstats + hashmap_ss[[dataset]]
	}

}else{
	ext_fullsumstats <- 0
}

# --- PREPARE SUMMARY STATISTICS FOR PRS-CSx
if ( ("PRSCSx" %in% model) & (ext > 0) ){

	if ( opt$verbose >= 1){
			cat("### PROCESSING SUMSTATS FOR PRS-CSx \n")
	}

	# process datasets 
	input <- c()
	ss <- c()
	pp <- c()

	for (loaded in loaded_datasets){
		name <- strsplit(loaded, split="[.]")[[1]][2]
		datasets_process_prscsx( name, eval(parse(text = loaded)), PRS_CSx_working_dir)
		input <- append(input, paste0(PRS_CSx_working_dir, name, ".txt"))
		ss <- append(ss, hashmap_ss[[name]])
		pp <- append(pp, pops[[name]])
	}

	input_prs <- paste(input, collapse=',')
	ss_prs <- paste(ss, collapse=',')
	pp_prs <- paste(pp, collapse=',')

	# run the shrinkage 
	wgt_prscsx <- shrinkage.prscsx(opt$dir_PRSCSx, opt$ldref_PRSCSx, PRS_CSx_working_dir, opt$phi_shrinkage_PRSCSx, input_prs, ss_prs, pp_prs, genos$bim)

	ext_prscsx <- length(wgt_prscsx)

}else{
	ext_prscsx <- 0
}

# --- PREPARE SUMMARY STATISTICS FOR MAGEPRO
if ( ("MAGEPRO" %in% model) & (ext > 0) ){

	if ( opt$verbose >= 1){
			cat("### PROCESSING SUMSTATS FOR MAGEPRO \n")
	}

	wgt_magepro <- c() #magepro weights
	for (loaded in loaded_datasets){
        name <- strsplit(loaded, split="[.]")[[1]][2]
        wgt_magepro <- datasets_process_susie(name, eval(parse(text = loaded)), wgt_magepro)
	}

	ext_magepro <- length(wgt_magepro)

}else{
	ext_magepro <- 0
}


# --- CROSSVALIDATION ANALYSES
set.seed(1)
avg_training_r2_single <- avg_training_r2_meta <- avg_training_r2_pt <- avg_training_r2_susie <- avg_training_r2_susie_impact <- avg_training_r2_prscsx <- avg_training_r2_magepro_fullsumstats <- avg_training_r2_magepro <- NA
#default crossval = 5 fold split
if ( opt$crossval <= 1 ) { 
if ( opt$verbose >= 1 ) cat("### SKIPPING CROSS-VALIDATION\n")
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
wgt.cv = matrix(NA, nrow=nrow(genos$bim), ncol=opt$crossval)
# ---

# --- keep track of SINGLE_model
SINGLE_top1 <- 0
# ---

r2_training_single <- c()
r2_training_meta <- c()
r2_training_pt <- c()
r2_training_susie <- c()
r2_training_susie_impact <- c()
r2_training_prscsx <- c()
r2_training_magepro_fullsumstats <- c()
r2_training_magepro <- c()

if ( ("META" %in% model) & (ext_fullsumstats > 0) ){
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
	if ( sum(is.na(pred.wgt)) == nrow(genos$bim)) {
		if ( opt$verbose >= 1 ) cat("LASSO pushed all weights to 0, using top1 as backup \n")
		SINGLE_top1 <- SINGLE_top1 + 1
		pred.wgt = weights.marginal( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3,drop=F]) , beta=T )
		pred.wgt[ - which.max( pred.wgt^2 ) ] = 0
	}
	if (length(pred.wgt) == 1){
		pred.wgt <- t(pred.wgt) # 1 snp in the cis window -> transpose for "matrix" mult
	}
	assign("pred.wgt", pred.wgt, envir = .GlobalEnv)
	if ("SINGLE" %in% model){	
	cv.calls[ indx , colcount ] = genos$bed[ cv.sample[ indx ] , , drop = FALSE] %*% pred.wgt
	pred_train_single = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], , drop = FALSE] %*% pred.wgt)))
	r2_training_single = append(r2_training_single, pred_train_single$adj.r.sq)
	
	colcount = colcount + 1

	}
	#--------------------------------------------------------------------------

	# SS-WEIGHTED META-ANALYSIS-------------------------------------------------------------
	if ("META" %in% model){
	if (ext_fullsumstats > 0){
	pred.wgt.meta = weights.meta(pred.wgt, training_ss, wgt_fullsumstats, hashmap_ss, total_ss_cv)	
	cv.calls[ indx , colcount ] = genos$bed[ cv.sample[ indx ] , , drop = FALSE] %*% pred.wgt.meta 
	pred_train_meta = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], , drop = FALSE] %*% pred.wgt.meta))) 
        r2_training_meta = append(r2_training_meta, pred_train_meta$adj.r.sq)
	}else{
	cv.calls[ indx , colcount ] = NA
	r2_training_meta = append(r2_training_meta, NA)
	}

	colcount = colcount + 1
	
	}
	#-----------------------------------------------------------------------------

	# P+T--------------------------------------------------------------------------
	if ("PT" %in% model){
	# if --prune_r2 and --threshold_p are not provided, cross validation on the training fold again to tune r2 and p values 
	if ( is.na(opt$prune_r2) & is.na(opt$threshold_p) ){
	tuned_pt <- tune.pt(genos$bed[cv.sample[-indx], , drop = FALSE], genos$bim, cv.train, c(0.2,0.5,0.8), c(0.001, 0.01, 0.1, 0.5), 2, opt$ldref_pt, opt$tmp)
	}else if ( is.na(opt$prune_r2) & !is.na(opt$threshold_p) ){
	tuned_pt <- tune.pt(genos$bed[cv.sample[-indx], , drop = FALSE], genos$bim, cv.train, c(0.2,0.5,0.8), c(opt$threshold_p), 2, opt$ldref_pt, opt$tmp)
	}else if ( !is.na(opt$prune_r2) & is.na(opt$threshold_p) ){
	tuned_pt <- tune.pt(genos$bed[cv.sample[-indx], , drop = FALSE], genos$bim, cv.train, c(opt$prune_r2), c(0.001, 0.01, 0.1, 0.5), 2, opt$ldref_pt, opt$tmp)
	}else{
	tuned_pt <- c(opt$prune_r2, opt$threshold_p)
	}
	pred.wgt.PT_sumstats <- weights.pt(genos$bed[cv.sample[-indx], , drop = FALSE], genos$bim, cv.train[,3], tuned_pt[1], tuned_pt[2], opt$ldref_pt, opt$tmp)
	if (sum(is.na(pred.wgt.PT_sumstats)) == nrow(genos$bim)){
	if ( opt$verbose >= 1 ) cat("No clumps remaining after P+T, NA results \n")
	cv.calls[ indx , colcount ] = NA
	r2_training_pt = append(r2_training_pt, NA)	
	}else{
	cv.calls[ indx , colcount ] = genos$bed[cv.sample[indx], , drop = FALSE] %*% pred.wgt.PT_sumstats
	pred_train_pt = summary(lm( cv.all[-indx,3] ~ (genos$bed[cv.sample[-indx], , drop = FALSE] %*% pred.wgt.PT_sumstats)))
        r2_training_pt = append(r2_training_pt, pred_train_pt$adj.r.sq)
	}
	
	colcount = colcount + 1

	}
	#------------------------------------------------------------------------------

	# SuSiE------------------------------------------------------------------------

	if ("SuSiE" %in% model){
	pred.wgt.susie <- weights.susie(genos$bed[cv.sample[-indx], , drop = FALSE], cv.train[,3])
	if ( sum(pred.wgt.susie == 0) == nrow(genos$bim) ) {
		pred.wgt.susie = weights.marginal( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3,drop=F]) , beta=T )
		pred.wgt.susie[ - which.max( pred.wgt.susie^2 ) ] = 0
	}
	if (length(pred.wgt.susie) == 1){
		pred.wgt.susie <- t(pred.wgt.susie) # 1 snp in the cis window -> transpose for "matrix" mult
	}
	cv.calls[ indx , colcount ] = genos$bed[ cv.sample[ indx ] , , drop = FALSE] %*% pred.wgt.susie
	pred_train_susie = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], , drop = FALSE] %*% pred.wgt.susie)))
        r2_training_susie = append(r2_training_susie, pred_train_susie$adj.r.sq)

	colcount = colcount + 1
	
	}

	#------------------------------------------------------------------------------

	# SuSiE_IMPACT-----------------------------------------------------------------

	if ("SuSiE_IMPACT" %in% model){
	pred.wgt.susieimpact <- weights.susie_impact(genos$bed[cv.sample[-indx], , drop = FALSE], cv.train[,3], opt$impact_path)
	cv.calls[ indx , colcount ] = genos$bed[ cv.sample[ indx ] , , drop = FALSE] %*% pred.wgt.susieimpact
	pred_train_susieimpact = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], , drop = FALSE] %*% pred.wgt.susieimpact)))
        r2_training_susie_impact = append(r2_training_susie_impact, pred_train_susieimpact$adj.r.sq)

	colcount = colcount + 1
	

	}

	#------------------------------------------------------------------------------

	# PRS-CSx----------------------------------------------------------------------

	if ("PRSCSx" %in% model){
	if (ext_prscsx > 0){
	pred.wgt.prs_csx <- weights.prscsx(wgt_prscsx, genos$bed[cv.sample[-indx], , drop = FALSE], cv.train[,3], "pred.wgt")
	cv.calls[ indx , colcount ] = genos$bed[ cv.sample[ indx ] , , drop = FALSE] %*% pred.wgt.prs_csx
	pred_train_prscsx = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], , drop = FALSE] %*% pred.wgt.prs_csx)))
        r2_training_prscsx = append(r2_training_prscsx, pred_train_prscsx$adj.r.sq)
	}else{
	cv.calls[ indx , colcount ] = NA
        r2_training_prscsx = append(r2_training_prscsx, NA)
	}

	colcount = colcount + 1

	}	

	#------------------------------------------------------------------------------

	# MAGEPRO_fullsumstats---------------------------------------------------------

	if ("MAGEPRO_fullsumstats" %in% model){
	if (ext_fullsumstats > 0){
	##pred.wgt.magepro_fullsumstats <- weights.magepro_marquezluna(pred.wgt, wgt_fullsumstats, genos$bed[cv.sample[-indx], , drop = FALSE], cv.train, genos$bim, opt$tmp, lasso_h2, FALSE)
	pred.wgt.magepro_fullsumstats <- weights.magepro(pred.wgt, wgt_fullsumstats, genos$bed[cv.sample[-indx], , drop = FALSE], cv.train[,3], FALSE)
	cv.calls[ indx , colcount ] = genos$bed[ cv.sample[ indx ], , drop = FALSE] %*% pred.wgt.magepro_fullsumstats
	pred_train_magepro_fullsumstats = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], , drop = FALSE] %*% pred.wgt.magepro_fullsumstats))) #r^2 between predicted and actual
        r2_training_magepro_fullsumstats = append(r2_training_magepro_fullsumstats, pred_train_magepro_fullsumstats$adj.r.sq)
	}else{
	cv.calls[ indx , colcount ] = NA
        r2_training_magepro_fullsumstats = append(r2_training_magepro_fullsumstats, NA)
	}
	
	colcount = colcount + 1

	}

	#------------------------------------------------------------------------------


	# MAGEPRO----------------------------------------------------------------------
	
	if ("MAGEPRO" %in% model){	
	if (ext_magepro > 0){	
	#pred.wgt.magepro <- weights.magepro_marquezluna(pred.wgt, wgt_magepro, genos$bed[cv.sample[-indx], , drop = FALSE], cv.train, genos$bim, opt$tmp, lasso_h2, FALSE)
	pred.wgt.magepro <- weights.magepro(pred.wgt, wgt_magepro, genos$bed[cv.sample[-indx], , drop = FALSE], cv.train[,3], FALSE)
	cv.calls[ indx , colcount ] = genos$bed[ cv.sample[ indx ], , drop = FALSE] %*% pred.wgt.magepro
	# --- for checking cor() of weights in CV
	wgt.cv[,cv] = pred.wgt.magepro
	# --- 
	#store the r2 on training set
	pred_train_magepro = summary(lm( cv.all[-indx,3] ~ (genos$bed[ cv.sample[-indx], , drop = FALSE] %*% pred.wgt.magepro))) #r^2 between predicted and actual 
	r2_training_magepro = append(r2_training_magepro, pred_train_magepro$adj.r.sq)	
	}else{
	cv.calls[ indx , colcount ] = NA
	r2_training_magepro = append(r2_training_magepro, NA)
	}
	
	}
	#-------------------------------------------------------------------------------

}

if (opt$verbose >= 1) cat("### COLLECTING CROSS VALIDATION RESULTS \n")

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
if ("SINGLE" %in% model) avg_training_r2_single <- mean(r2_training_single)
if ("META" %in% model) avg_training_r2_meta <- mean(r2_training_meta)
if ("PT" %in% model) avg_training_r2_pt <- mean(r2_training_pt)
if ("SuSiE" %in% model) avg_training_r2_susie <- mean(r2_training_susie)
if ("SuSiE_IMPACT" %in% model) avg_training_r2_susie_impact <- mean(r2_training_susie_impact)
if ("PRSCSx" %in% model) avg_training_r2_prscsx <- mean(r2_training_prscsx)
if ("MAGEPRO_fullsumstats" %in% model) avg_training_r2_magepro_fullsumstats <- mean(r2_training_magepro_fullsumstats)
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
	if ( opt$verbose >= 1 ) cat("LASSO pushed all weights to 0, using top1 as backup \n")
	pred.wgtfull = weights.marginal( genos$bed , as.matrix(pheno[,3]) , beta=T )
	pred.wgtfull[ - which.max( pred.wgtfull^2 )] = 0
}
if(length(pred.wgtfull) == 1){
	pred.wgtfull <- t(pred.wgtfull)
}
assign("pred.wgtfull", pred.wgtfull, envir = .GlobalEnv)
if ("SINGLE" %in% model){
wgt.matrix[, colcount] = pred.wgtfull
colcount = colcount + 1
}
# --- SS-WEIGHTED META-ANALYSIS
if ("META" %in% model){
if (ext_fullsumstats > 0){
total_ss_full <- total_ss_sumstats + N.tot
pred.wgt.metafull <- weights.meta(pred.wgtfull, N.tot, wgt_fullsumstats, hashmap_ss, total_ss_full)
wgt.matrix[, colcount] = pred.wgt.metafull
}else{
wgt.matrix[, colcount] = NA
}
colcount = colcount + 1
}
# --- P+T 
if ("PT" %in% model){
if ( is.na(opt$prune_r2) & is.na(opt$threshold_p) ){
tuned_pt <- tune.pt(genos$bed, genos$bim, pheno, c(0.2,0.5,0.8), c(0.001, 0.01, 0.1, 0.5), 5, opt$ldref_pt, opt$tmp)
}else if ( is.na(opt$prune_r2) & !is.na(opt$threshold_p) ){
tuned_pt <- tune.pt(genos$bed, genos$bim, pheno, c(0.2,0.5,0.8), c(opt$threshold_p), 5, opt$ldref_pt, opt$tmp)
}else if ( !is.na(opt$prune_r2) & is.na(opt$threshold_p) ){
tuned_pt <- tune.pt(genos$bed, genos$bim, pheno, c(opt$prune_r2), c(0.001, 0.01, 0.1, 0.5), 5, opt$ldref_pt, opt$tmp)
}else{
tuned_pt <- c(opt$prune_r2, opt$threshold_p)
}
pred.wgt.PT_sumstatsfull <- weights.pt(genos$bed, genos$bim, pheno[,3], tuned_pt[1], tuned_pt[2], opt$ldref_pt, opt$tmp)
wgt.matrix[, colcount] = pred.wgt.PT_sumstatsfull
colcount = colcount + 1
}
# --- SuSiE
if ("SuSiE" %in% model){
pred.wgt.susiefull <- weights.susie(genos$bed, pheno[,3])
if ( sum(pred.wgt.susiefull == 0) == nrow(genos$bim) ) {
pred.wgt.susiefull = weights.marginal( genos$bed , as.matrix(pheno[,3]) , beta=T ) # use marginal weights for susie if NA
pred.wgt.susiefull[ - which.max( pred.wgt.susiefull^2 ) ] = 0
}
if (length(pred.wgt.susiefull) == 1){
pred.wgt.susiefull <- t(pred.wgt.susiefull) # 1 snp in the cis window -> transpose for "matrix" mult
}
wgt.matrix[, colcount] = pred.wgt.susiefull
colcount = colcount + 1
}
# --- SuSiE_IMPACT
if ("SuSiE_IMPACT" %in% model){
pred.wgt.susieimpactfull <- weights.susie_impact(genos$bed, pheno[,3], opt$impact_path)
wgt.matrix[, colcount] = pred.wgt.susieimpactfull
colcount = colcount + 1
}
# --- PRS-CSx
if ("PRSCSx" %in% model){
if (ext_prscsx > 0){
pred.wgt.prs_csxfull <- weights.prscsx(wgt_prscsx, genos$bed, pheno[,3], "pred.wgtfull")
wgt.matrix[, colcount] = pred.wgt.prs_csxfull
}else{
wgt.matrix[, colcount] = NA
}
colcount = colcount + 1
}
# --- MAGEPROfullsumstats
if ("MAGEPRO_fullsumstats" %in% model){
if (ext_fullsumstats > 0){
#pred.wgt.magepro_fullsumstatsfull <- weights.magepro_marquezluna(pred.wgtfull, wgt_fullsumstats, genos$bed, pheno, genos$bim, opt$tmp, lasso_h2, FALSE)
pred.wgt.magepro_fullsumstatsfull <- weights.magepro(pred.wgtfull, wgt_fullsumstats, genos$bed, pheno[,3], FALSE)
wgt.matrix[, colcount] = pred.wgt.magepro_fullsumstatsfull
}else{
wgt.matrix[, colcount] = NA
}
colcount = colcount + 1
}
# --- MAGEPRO
cf_total = NA
if ("MAGEPRO" %in% model){
if (ext_magepro > 0){
#pred.wgt.mageprofull <- weights.magepro_marquezluna(pred.wgtfull, wgt_magepro, genos$bed, pheno, genos$bim, opt$tmp, lasso_h2, TRUE)
pred.wgt.mageprofull <- weights.magepro(pred.wgtfull, wgt_magepro, genos$bed, pheno[,3], TRUE)
wgt.matrix[, colcount] = pred.wgt.mageprofull
}else{
wgt.matrix[, colcount] = NA
}
}

#--- SAVE RESULTS
snps = genos$bim
if ( ("MAGEPRO" %in% model) & (ext_magepro > 0) ){
wgtmagepro <- append("pred.wgt", wgt_magepro)
}else{
wgtmagepro = NA
}

# --- for checking cor() of weights in CV
cors_weights <- c()
for (i in 1:(opt$crossval-1)){
	cors_weights <- append(cors_weights, cor(wgt.cv[,i], wgt.cv[,(i+1)]) )
}
avg_cor <- mean(cors_weights)
# ---

save( wgt.matrix, snps, cv.performance, hsq, hsq.pv, N.tot , wgtmagepro, cf_total, avg_training_r2_single, avg_training_r2_meta, avg_training_r2_pt, avg_training_r2_susie, avg_training_r2_susie_impact, avg_training_r2_prscsx, avg_training_r2_magepro_fullsumstats, avg_training_r2_magepro, var_cov, avg_cor, SINGLE_top1, file = paste( opt$out , ".wgt.RDat" , sep='' ) )

# --- CLEAN-UP
if ( opt$verbose >= 1 ) cat("### CLEANING UP\n")
cleanup()
