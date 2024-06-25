library(data.table)
library(stringr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
base <- args[1]
weights_name <- args[2]
files <- list.files(base)
len <- length(files)

output <- matrix(0, nrow = len, ncol = 1)
for (f in 1:len){
	path <- paste0(base, files[f])
	output[f,1] <- path 
}
write.table(output, file = paste0("magepro_weights_paths_", weights_name,".txt"), row.names = F, col.names = F, sep = "\t", quote = F)
