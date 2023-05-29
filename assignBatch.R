library(data.table)
tissues <- "Whole_Blood"
batches <- 20 #20 jobs 

#splitting to submit jobs 
for (i in 1:length(tissues)){ 
print(tissues[i])

#get all of genes names from the scratch plink directory
allgenes <- list.files(path = "/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/plink_cissnps_AFR/", pattern = ".fam")
allgenenames <- gsub(patter = ".fam", replacement = "", x=allgenes)
#goal: first col = gene name, second col = batch number 
print(nrow(allgenenames))
batchnum <- cut(seq(1,nrow(allgenenames)),breaks=batches,labels=FALSE) #breaks up gene names in 20 groups 
print(batchnum)
x <- cbind(allgenenames,batchnum) #each gene name has a corresponding batch number ranging from 1-20
write.table(x, file = paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO/intermedfiles/genes_assign_",tissues[i],"2.txt"), row.names = F, col.names = F, sep = "\t", quote= F) 
}

