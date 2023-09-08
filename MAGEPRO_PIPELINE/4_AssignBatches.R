library(data.table)

#--- read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
batches <- as.numeric(args[1]) #number of batches
intermediate_dir <- args[2]
scratchplink <- args[3]

#--- get all of genes names from the scratch directory and split them up into jobs 
allgenes <- list.files(path = paste0(scratchplink,"/"), pattern = ".fam")
allgenenames <- as.data.frame(gsub(patter = ".fam", replacement = "", x=allgenes))
#goal: first col = gene name, second col = batch number 
batchnum <- cut(seq(1,nrow(allgenenames)),breaks=batches,labels=FALSE) #break up genes in batches 
x <- cbind(allgenenames,batchnum) #each gene name has a corresponding batch number ranging from 1-20
write.table(x, file = paste0(intermediate_dir,"/genes_assign_Whole_Blood.txt"), row.names = F, col.names = F, sep = "\t", quote= F) 
