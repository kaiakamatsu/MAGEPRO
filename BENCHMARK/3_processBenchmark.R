library(data.table)
library(dplyr)


tissues <- "Whole_Blood"

for (i in 1:length(tissues)){
genes <- read.table(file = paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/genes_assign_Whole_Blood3",".txt"), header = F)$V1 #all genes in analysis

r2_h2 <- matrix(0, length(genes), 9)

for (j in 1:length(genes)){

file <- paste0("/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/weights_benchmark_confirm/",tissues[i],"/",tissues[i],".",genes[j],".wgt.RDat") 

print(genes[j])

if(file.exists(file)){
load(file)
r2_h2[j,1] <- genes[j]
r2_h2[j,2:3] <- cv.performance[1,]
r2_h2[j,4:5] <- cv.performance[2,]
r2_h2[j, 6] <- paste(wgts, collapse = ",")
r2_h2[j, 7] <- hsq_afr[1]
r2_h2[j, 8] <- hsq_afr[2]
r2_h2[j, 9] <- hsq_afr.pv

}
}


newcolnamesr2h2 <- c("gene", "PRS_CSx_r2", "PT_r2","PRS_CSx_pv", "PT_pv","datasets", "hsq_afr", "hsq_afr_se", "hsq_afr_pv")
colnames(r2_h2) <- newcolnamesr2h2

h <- which(r2_h2[,2] == 0)
if(length(h) > 0){r2_h2 <- r2_h2[-h,]}


write.table(r2_h2, file = "benchmark_results_confirm.txt", row.names = F, col.names = T, sep = "\t", quote = F)
}

