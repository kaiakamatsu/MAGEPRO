# load packages
suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('glmnet'))
suppressMessages(library('data.table'))
suppressMessages(library('dplyr'))


# --- UNUSED FUNCTIONS BELOW

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
	assign(paste0("pred.wgt.", dataset), meta, envir = .GlobalEnv)
	result_wgts <- append(result_wgts, paste0("pred.wgt.", dataset))
	}
	wgts_indices <- list(result_wgts, indices)
	return (wgts_indices)
}

# --- BIN SNPs BASED ON IMPACT SCORES
IMPACT_split <- function(IMPACT_file, BIM){
	# PURPOSE: find the indices of potentially predictive snps (used to split sumstats)
	# IMPACT_file = path to IMPACT file 
	# BIM = BIM matrix
	# RETURN: groupings: vector of strings, "nonzero" stores indices of important SNPs and "zero" stored indices of unimportant ones 
	IMPACT_scores <- as.data.frame(fread(IMPACT_file))
	IMPACT_scores_reordered <- IMPACT_scores[match(BIM[,2], IMPACT_scores$SNP),]
	IMPACT_threshold <<- quantile(IMPACT_scores_reordered[,5], 0.95) #CHANGE THIS TO TAKE TOP X%
	nonzero <<- which(IMPACT_scores_reordered[,5] >= IMPACT_threshold)
	groupings <- c()
        if (!identical(integer(0), nonzero)){
                groupings <- append(groupings, "nonzero")
        }
        return (groupings)		
}


# --- UNUSED FUNCTIONS END