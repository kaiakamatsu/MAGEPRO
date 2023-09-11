library(data.table)

#--- read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
intermediate_dir <- args[1]
batches <- as.numeric(args[2]) #number of batches

#--- get all of genes names from the scratch directory and split them up into jobs 
allgenenames <- fread(file = paste0(intermediate_dir, "/All_Genes_Expressed.txt"), header = F)
#goal: first col = gene name, second col = batch number 
batchnum <- cut(seq(1,nrow(allgenenames)),breaks=batches,labels=FALSE) #break up genes in batches 
x <- cbind(allgenenames,batchnum) #each gene name has a corresponding batch number ranging from 1-20
write.table(x, file = paste0(intermediate_dir,"/Genes_Assigned.txt"), row.names = F, col.names = F, sep = "\t", quote= F) 
