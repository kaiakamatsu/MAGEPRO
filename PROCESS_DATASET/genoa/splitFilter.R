#--- read in dplyr package for filtering 
library(dplyr)
library(data.table)
print("packages is ready")


#--- function to filter for SNPs we want - CHANGE COLUMN NAME OF eQTL corresponding to the SNPs 
getSNPs = function( eQTL , SNPs) {
  eQTL2 <- filter(eQTL, match(rs, SNPs$V2, nomatch = -1) != -1) 
  return (eQTL2)
}
print("getsnp function created")


#--- split PER GENE - change name of gene column 
splitGene = function( eQTL ) {
	genes <- split(eQTL, eQTL$GENE)
	return (genes)
}
print("splitGene function created")


#--- read in command line arguments as a vector 
args = commandArgs(trailingOnly=TRUE)
print("arguments read")

if (length(args)==0){
	print("no input files")
	stop()
}else{ 
	EQTL = args[1]
	SNPS = args[2]
	OUT = args[3]
	print("arguments assigned to variables")
}


#--- read in eqtl files 
if (file.exists(EQTL)){
	eQTL = fread(EQTL, header = T) 
	print("eQTL table read")
}else{
	print("no eqtl file")
	stop()
}

#--- read snps to keep
if (file.exists(SNPS)){ #check if SNPS file exists 
	SNPs = fread(SNPS, header = F) #read in SNPs we want 
	print("SNP table read")
}else{
	print("no SNP reference file")
	stop()
}

#--- extract snps of interest
extracted <- getSNPs(eQTL, SNPs)
print("extraction completed")

extracted <- extracted[, c(2, 3, 5, 6, 8)]
print(head(extracted))

#--- split the eQTL file by gene into multiple data frames using the splitGene function 
split <- splitGene(extracted)
print("file split by gene")


#--- write a new txt file into the output directory (OUT) for each gene table, change column number of gene depending on eQTL file  
for (x in split){
	write.table(x, file = paste(OUT, paste0(x[1,1], ".txt"), sep = "/"), quote = FALSE, row.names = F)
}
