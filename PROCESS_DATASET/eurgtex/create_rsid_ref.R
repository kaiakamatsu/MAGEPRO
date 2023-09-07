library(data.table)
library(stringr)
library(dplyr)


args = commandArgs(trailingOnly=TRUE)
print("arguments read")

if (length(args)==0){ #make sure there are input files 
	print("no input files")
	stop()
}else{ #assign the paths given in the command line arguments to variables 
	SNPS_file = args[1]
	OUT = args[2]
	print("arguments assigned to variables")
}


if (file.exists(SNPS_file)){
	SNPS = fread(file = SNPS_file, header = F)
	print("snps reference read") 
}else{
	print("snps reference file does not exist")
	stop()
}


chroms <- split(SNPS, SNPS$V1)

for (c in chroms){
	# preallocate memory
	n_rows <- nrow(c)
	variant_id <- character(n_rows)
	id <- character(n_rows)
	chr <- character(n_rows)
	A_effect <- character(n_rows)
	A_alt <- character(n_rows)

	# Loop through the data using vectorized operations
	for (i in 1:n_rows) {
    		n <- paste0('chr', c[i, 1], '_', c[i, 4], '_', c[i, 6], '_', c[i, 5], '_b38')
    		print(n)
    		variant_id[i] <- n
    		id[i] <- c[i, 2]
    		chr[i] <- c[i, 1]
    		A_effect[i] <- c[i, 5]
    		A_alt[i] <- c[i, 6]
	}

	# Create the data frame
	id_name <- matrix(c(id, variant_id, chr, A_effect, A_alt), nrow = n_rows, ncol = 5)
	#id_name <- data.frame(id, variant_id, chr, A_effect, A_alt)
	print("dataframe created")

	write.table(id_name, file = paste0(OUT, "/", id_name[1,3],"snpMAP.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)

}
