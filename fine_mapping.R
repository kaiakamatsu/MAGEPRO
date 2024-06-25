expected_header <- "Gene SNP A1 A2 b SE P"
suppressMessages(library('tools'))

gene_fine_mapping <- function(gene_txt, cohort, cohort_data, cohort_path, cohort_ld_directory, cohort_ldref_path, plink, cl_thresh, out) {
	file_path <- file.path(cohort_path, gene_txt)
	if (trimws(readLines(file_path, n = 1)) != expected_header) {
		stop(paste0("Error: Header in", file.path(cohort_path, gene_txt), "does not match the expected", expected_header, "format.\nMake sure all sumstat files have the expected header.\n"))
		cleanup()
		q()
	}

	gene <- file_path_sans_ext(gene_txt)
	
	clump_cmd <- paste0(
		plink, 
		" --bfile ", cohort_ldref_path,
		" --r2 inter-chr --ld-window-r2 0",
		" --clump ", file.path(cohort_path, gene_txt),
		" --clump-p1 1 --clump-r2 ", cl_thresh,
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

	print(paste0("rcommand: ", rcommand))

	tryCatch({
		system(rcommand, wait = TRUE)
	}, error = function(e) {
		cat("Error in R script command:", e$message, "\n")
		return("Error")
	})

	cat("finished running gene_fine_mapping\n")

	to_return <- file.path(out, gene_txt)

	return(to_return)
}