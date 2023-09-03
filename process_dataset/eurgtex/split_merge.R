library(data.table)
library(stringr)
library(dplyr)
library(arrow)

args = commandArgs(trailingOnly=TRUE)
print("arguments read")

if (length(args)==0){ #make sure there are input files 
	print("no input files")
	stop()
}else{ #assign the paths given in the command line arguments to variables 
	EQTL = args[1]
	SNPS_file = args[2]
	OUT = args[3]
	print("arguments assigned to variables")
}


if (file.exists(SNPS_file)){
	SNPS = fread(file = SNPS_file, header = T)
	print("snps reference read") 
}else{
	print("snps reference file does not exist")
	stop()
}
colnames(SNPS) <- c("id", "variant_id", "chr", "A_effect", "A_alt")

if (file.exists(EQTL)){
	df_sumstats <- read_parquet(EQTL)
	print("sumstats read") 
}else{
	print("eqtl file does not exist")
	stop()
}

df_sumstats <- merge(df_sumstats, SNPS, by = 'variant_id')
# only variants present in both are kept

df_sumstats <- df_sumstats[, c(2, 11, 13, 14, 8)]

print(head(df_sumstats))

groups <- split(df_sumstats, df_sumstats$phenotype_id)
print("file split")

for (x in groups){
	gene <- strsplit(x[1,1], split = "[.]")[[1]][1]
	print(gene)
	write.table(x, file = paste(OUT, paste0(gene, ".txt"), sep = "/"), quote = FALSE, row.names = F)
}
print("split by genes")
