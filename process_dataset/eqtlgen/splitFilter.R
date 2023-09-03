#read in dplyr package for filtering 
library(dplyr)
library(data.table)
print("packages is ready")

#function to filter and keep SNP we want - CHANGE COLUMN NAME OF eQTL corresponding to the SNPs 
getSNPs = function( eQTL , SNPs) {
  eQTL2 <- filter(eQTL, match(SNP, SNPs$V2, nomatch = -1) != -1) #filter the eQTL file so that only rows with rs values that are in the SNPs file are left
  return (eQTL2)
}
print("getsnp function created")


#split PER GENE and write a new txt file for each in the OUT directory, change the column name depending on the eQTL file 
splitGene = function( eQTL ) {
	genes <- split(eQTL, eQTL$Gene)
	return (genes)
}
print("splitGene function created")


#read in command line arguments as a vector 
args = commandArgs(trailingOnly=TRUE)
print("arguments read")

if (length(args)==0){ #make sure there are input files 
	print("no input files")
	stop()
}else{ #assign the paths given in the command line arguments to variables 
	EQTL = args[1]
	SNPS = args[2]
	OUT = args[3]
	print("arguments assigned to variables")
}


if (file.exists(EQTL)){ #check if eQTL file exists
	eQTL = fread(EQTL, header = T) 
	print("eQTL table read")
}else{
	print("no eqtl file")
	stop()
}


if (file.exists(SNPS)){ #check if SNPS file exists 
	SNPs = read.table(SNPS, header = F, as.is = F) #read in SNPs we want 
	print("SNP table read")
}else{
	print("no snp reference file")
	stop()
}


#run get SNP function 
extracted <- getSNPs(eQTL, SNPs)
print("extraction completed")

extracted <- extracted[, c(10, 1, 4, 5, 12)]
print(head(extracted))

#split the eQTL file by gene into multiple data frames using the splitGene function 
split <- splitGene(extracted)
print("file split by gene")


#write a new txt file into the scratch directory (OUT) per gene table, change column number of gene depending on eQTL file  
for (x in split){
	write.table(x, file = paste(OUT, paste0(x[1,1], ".txt"), sep = "/"), quote = FALSE, row.names = F)
}

