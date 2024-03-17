library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
output <- args[1] #path to output directory
weightsdir <- args[2]
genes_assign <- args[3] #path to genes assigned file
models <- strsplit(args[4], split = ",")[[1]] #comma-separated list of models used (SINGLE,META,MAGEPRO)

num_models <- length(models)

genes <- read.table(file = genes_assign, header = F)$V1 #all genes in analysis

r2_h2 <- matrix(0, length(genes), 2*num_models+4)

for (j in 1:length(genes)){

file <- paste0(weightsdir,"/",genes[j],".wgt.RDat") 
print(genes[j])

if(file.exists(file)){
if(file.info(file)$size > 0){
load(file)
r2_h2[j,1] <- genes[j]
r2_h2[j,2:(num_models+1)] <- cv.performance[1,]
r2_h2[j,((num_models+2): (2*num_models+1) )] <- cv.performance[2,]
r2_h2[j, 2*num_models+2] <- hsq[1]
r2_h2[j, 2*num_models+3] <- hsq[2]
r2_h2[j, 2*num_models+4] <- hsq.pv

}
}
}

h <- which(r2_h2[,1] == 0)
if(length(h) > 0){r2_h2 <- r2_h2[-h,]}

newcolnamesr2h2 <- c("gene", paste0(models, "_r2"), paste0(models, "_pv"), "hsq", "hsq_se", "hsq.pv")
colnames(r2_h2) <- newcolnamesr2h2

write.table(r2_h2, file = paste0(output, "/PRSCSX_results.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
