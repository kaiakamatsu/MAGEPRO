library(data.table)
library(stringr)
library(dplyr)


getSNPs = function( eQTL , SNPs) {
  eQTL2 <- filter(eQTL, match(Variant_ID, SNPs$V2, nomatch = -1) != -1) #filter the eQTL file so that only rows with rs values that are in the SNPs file are left
  return (eQTL2)
}
print("snp filter function created")


args = commandArgs(trailingOnly=TRUE)
cell_type = args[1]
eQTL = args[2]
SNP = args[3]
OUT = args[4]


if (file.exists(SNP)){
	SNPS = fread(SNP, header = F)
	print("read SNPs")
}else{
	print("no SNPs reference file")
	stop()
}


file = paste0(eQTL, "/", cell_type, "_nominal.txt")

if (file.exists(file)){
	file = fread(file, header = T)
	print("file read") 
}else{
	print(paste0(cell_type, " file does not exist"))
}


file <- getSNPs(file, SNPS)
print("snps extracted")

file <- file[, c(1, 6, 8, 7, 13)]
print(head(file))

groups <- split(file, file$Gene_id)
print("file split by genes")

system(paste0("mkdir ", OUT, "/", cell_type))

for (x in groups){
	basename <- strsplit(as.character(x[1,1]), ".", fixed = TRUE)[[1]][1] 
	write.table(x, file = paste(OUT, cell_type, paste0(basename, ".txt"), sep = "/"), quote = FALSE, row.names = F)
}
print("split genes for one cell type")

