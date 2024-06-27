expected_header <- c("Gene", "SNP","A1", "A2", "BETA", "SE", "P")
suppressMessages(library('tools'))
suppressMessages(library('data.table'))

cohort_fine_mapping <- function(cohort_map, sumstats_dir, tmp, ldref_dir, out, gene, PATH_plink, cl_thresh, verbose) {
	for (cohort in names(cohort_map)) {
		cohort_data <- cohort_map[[cohort]]
		cohort_path <- file.path(sumstats_dir, cohort)
		cohort_ld_directory <- file.path(tmp, paste0(cohort, "ld"))
		cohort_ldref_path <- file.path(ldref_dir, cohort_data$ldref)
		out_cohort_path <- file.path(out, cohort)

		# create directories for output and ld_matrix_path
		system(paste0("mkdir -p ", cohort_ld_directory), wait = TRUE)
		system(paste0("mkdir -p ", out_cohort_path), wait = TRUE)

		gene_txt <- paste0(gene, '.txt')
		check <- file.path(cohort_path, gene_txt)
		if (!file.exists(check)){
			if (verbose == 2) {
				cat("skipping sumstat", cohort, "for this gene\n")
			}
			next
		}
		path_to_fine_mapping_output <- gene_fine_mapping(gene_txt, cohort, cohort_data, cohort_path, cohort_ld_directory, cohort_ldref_path, PATH_plink, cl_thresh, out_cohort_path)
		if (path_to_fine_mapping_output == "Error") {
			next
		}
		assign(paste0("file.", cohort), path_to_fine_mapping_output, envir = .GlobalEnv) # reassign file path to the path to susie results
		cat("successfully fine mapped ", gene, " for ", cohort, "\n")
	}
}


gene_fine_mapping <- function(gene_txt, cohort, cohort_data, cohort_path, cohort_ld_directory, cohort_ldref_path, plink, cl_thresh, out) {
	file_path <- file.path(cohort_path, gene_txt)
	if ( ! identical( colnames(fread(file_path, nrows = 0)), expected_header) ) {
		cat( "WARNING: rewriting summary statistics file with correct headers. Dropping extra columns. Make sure statistics are in the correct order specified in README \n" , sep='', file=stderr() )
		temp_df <- fread(file_path, select = ( c(1:length(expected_header)) ) )
		colnames(temp_df) <- expected_header
		write.table(temp_df, file = file_path, sep = ' ', quote = F, col.names = T, row.names = F)
	}

	gene <- file_path_sans_ext(gene_txt)
	
	clump_cmd <- paste0(
		plink, 
		" --bfile ", cohort_ldref_path,
		" --clump ", file_path,
		" --clump-p1 1 --clump-r2 ", cl_thresh,
		" --extract ", file_path,
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
		cat("Error in clump command: ", e$message, " Moving to next cohort...\n")
		return("Error")
	})

	get_ld_cmd <- paste0(
		plink,
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
		cat("Error in get_ld command: ", e$message, " Moving to next cohort...\n")
		return("Error")
	})

	rcommand <- paste0("Rscript SUSIE/get_pips.R -g ", cohort_path, "/",
	" -c ", gene_txt, 
	" -l ", file.path(cohort_ld_directory, paste0(gene, ".ld")),
	" -n ", cohort_data$sample_size,
	" -o ", file.path(out, gene_txt)
	)

	tryCatch({
		system(rcommand, wait = TRUE)
	}, error = function(e) {
		cat("Error in R script command:", e$message, "\n")
		return("Error")
	})

	to_return <- file.path(out, gene_txt)

	return(to_return)
}