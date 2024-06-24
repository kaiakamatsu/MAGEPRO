# load packages
suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('glmnet'))
suppressMessages(library('data.table'))
suppressMessages(library('dplyr'))
suppressMessages(library('tools'))


option_list = list(
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--scratch", action="store", default=NA, type='character',
              help="Path to temporary files [required]"),
  make_option("--ldref_dir", action="store", default=NA, type="character",
  			  help="Directory containing ld reference files [required]"),
  make_option("--ldrefs", action="store", default=NA, type="character",
  			  help="Comma-separated list of ld reference files [required]"),
  make_option("--sumstats_dir", action="store", default=NA, type='character',
              help="Path to external sumstats"),
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Comma-separated list of external datasets to include"),
  make_option("--ss", action="store", default=NA, type='character',
              help="Comma-separated list of sample sizes of sumstats (in the same order)"), 
  make_option("--n_threads", action="store", default=1, type="numeric",
  			  help="Number of threads used in parallelization of running SuSiE [optional]"),
  make_option("--cl_thresh", action="store", default=0.97, type="numeric",
  			  help="Clumping threshold for plink to clump SNPs that have high R2 [optional]"),
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"),        
  make_option("--susie_pip", action="store", default=8, type='numeric',
              help="Column number in external datasets where susie pips are stored"),
  make_option("--susie_beta", action="store", default=9, type='numeric',
              help="Column number in external datasets where susie coefs are stored"),
  make_option("--susie_cs", action="store", default=10, type='numeric',                                             
              help="Column number in external datasets where susie credible set groups are stored")
  make_option("--resume_partial", action="store_true", default=FALSE,
            help="Skip completed jobs if the code was interrupted, to avoid rerunning the entire dataset [optional]")
)

# --- PARSE COMMAND LINE ARGS
opt = parse_args(OptionParser(option_list=option_list))


# only run SuSiE on specified folders
sumstats_list <- strsplit(opt$sumstats, ",")[[1]]
cohort_sizes <- strsplit(opt$ss, ",")[[1]]
ldrefs_list <- strsplit(opt$ldrefs, ",")[[1]]

# create map with cohort as keys and 
cohort_map <- setNames(
  lapply(seq_along(sumstats_list), function(i) {
    list(cohort_size = cohort_sizes[i], ldref = ldrefs_list[i])
  }),
  sumstats_list
)

expected_header <- "Gene SNP A1 A2 b SE P"

for (cohort in names(cohort_map)) {
	cohort_data <- cohort_map[[cohort]]
	cohort_path <- file.path(opt$sumstats_dir, cohort)
	cohort_ld_directory <- file.path(opt$scratch, paste0(cohort, "ld"))
	cohort_ldref_path <- file.path(opt$ldref_dir, cohort_data$ldref)
	genes <- list.files(cohort_path)

	# create directories for output and ld_matrix_path
	make_ld_matrix_path <- paste0(
		"mkdir -p ", opt$scratch, "/", cohort, "ld/"
	)
	make_output_path <- paste0(
		"mkdir -p ", opt$out, "/", cohort, "/"
	)

	system(make_ld_matrix_path, wait = TRUE)
	system(make_output_path, wait = TRUE)

	for (gene_txt in genes) {
		
		# TODO this checker might slow down the runtime (not significantly, but still)
		# TODO could make a boolean logic for the for loop by checking for the flag outside of for loop.
		if (opt$resume_partial) {
			output_file <- file.path(opt$out, cohort, gene_txt)
			if (file.exists(output_file)) {
				cat(paste0("Output for ", cohort, "/", gene_txt, " already exists. Skipping...\n"))
				next
			}
		}

		# TODO this checker might slow down the runtime (not significantly, but still)
		if (readLines(gene_txt, n = 1) != "Gene SNP A1 A2 b SE P") {
			stop(paste0("Header in", file.path(cohort_path, gene_txt), "does not match the expected", expected_header, "format.\nMake sure all files have the expected header.\n"))
		}

		gene <- file_path_sans_ext(gene_txt)
		
		clump_cmd <- paste0(
			opt$PATH_plink, 
			" --bfile ", cohort_ldref_path,
			" --r2 inter-chr --ld-window-r2 0",
			" --clump ", file.path(cohort_path, gene_txt),
			" --clump-p1 1 --clump-r2 ", opt$cl_thresh,
			" --extract ", file.path(cohort_path, gene_txt),
			" --out ", file.path(cohort_ld_directory, gene),
			" > /dev/null 2> ", gene, "error.log"
		)

		tryCatch({
			exit_status <- system(clump_cmd, wait = TRUE)
			# Check if the command failed and if error.log exists and is not empty
			if (exit_status != 0 && file.exists(paste0(gene, "error.log")) && file.info(paste0(gene, "error.log"))$size > 0) {
				cat(readLines(paste0(gene, "error.log")), sep = "\n")
			}
			file.remove(paste0(gene, "error.log"))
		}, error = function(e) {
			cat("Error in clump command: ", e$message, " Moving to other genes...\n")
			next
		})

		get_ld_cmd <- paste0(
			opt$PATH_plink,
			" --bfile ", cohort_ldref_path,
			" --r2 inter-chr --ld-window-r2 0",
			" --extract ", file.path(cohort_ld_directory, paste0(gene, ".clumped")), 
			" --out ", file.path(cohort_ld_directory, gene),
			" > /dev/null 2> ", gene, "error.log"
		)

		tryCatch({
			exit_status <- system(get_ld_cmd, wait = TRUE)
			# Check if the command failed and if error.log exists and is not empty
			if (exit_status != 0 && file.exists(paste0(gene, "error.log")) && file.info(paste0(gene, "error.log"))$size > 0) {
				cat(readLines(paste0(gene, "error.log")), sep = "\n")
			}
			file.remove(paste0(gene, "error.log"))
		}, error = function(e) {
			cat("Error in get_ld command: ", e$message, " Moving to other genes...\n")
			next
		})

		rcommand <- paste0("Rscript SUSIE/get_pips.R -g ", cohort_path, "/",
		" -c ", gene_txt, 
		" -l ", file.path(cohort_ld_directory, paste0(gene, ".ld")),
		" -n ", cohort_data$cohort_size,
		" -o ", file.path(opt$out, cohort, gene_txt)
		)

		tryCatch({
			system(rcommand, wait = TRUE)
		}, error = function(e) {
			cat("Error in R script command:", e$message, "\n")
			next
		})
	}
}
