library(data.table)
library(stringr)
library(dplyr)

base <- "/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/weights_MAGEPRO_fullsumstats/Whole_Blood/"
files <- list.files(base)
len <- length(files)
print(paste0("number of wgt files: ", len))

output <- matrix(0, nrow = len, ncol = 1)

for (f in 1:len){
	path <- paste0(base, files[f])
	print(path)
	output[f,1] <- path 
}

write.table(output, file = "magepro_weights_paths.txt", row.names = F, col.names = F, sep = "\t", quote = F)
