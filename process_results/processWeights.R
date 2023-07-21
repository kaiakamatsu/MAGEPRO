library(data.table)
library(dplyr)


tissues <- "Whole_Blood"

for (i in 1:length(tissues)){
genes <- read.table(file = paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/genes_assign_Whole_Blood2",".txt"), header = F)$V1 #all genes in analysis

r2_h2 <- matrix(0, length(genes), 14)

for (j in 1:length(genes)){

file <- paste0("/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/weights_MAGEPRO_final/",tissues[i],"/",tissues[i],".",genes[j],".wgt.RDat") 

print(genes[j])

if(file.exists(file)){
if(file.info(file)$size > 0){
load(file)
r2_h2[j,1] <- genes[j]
r2_h2[j,2:(ncol(cv.performance)+1)] <- cv.performance[1,]
r2_h2[j,((ncol(cv.performance)+2): 7)] <- cv.performance[2,]
r2_h2[j, 8] <- hsq_afr[1]
r2_h2[j, 9] <- hsq_afr[2]
r2_h2[j, 10] <- hsq_afr.pv
if (length(wgt2) > 0){
	print(paste(wgt2, collapse = ","))
	r2_h2[j, 11] <- paste(wgt2, collapse = ",")
}else{
	r2_h2[j, 11] <- NA
}

nonzeroAFR <- length(which(wgt.matrix[, 1] != 0))
nonzeroMAGEPRO <- length(which(wgt.matrix[, 3] != 0))

r2_h2[j, 12] <- nonzeroAFR
r2_h2[j, 13] <- nonzeroMAGEPRO

ab <- c()

if (is.matrix(alpha_beta)){
	for (c in 1:ncol(alpha_beta)){
		ab <- append(ab, sum(alpha_beta[,c]^2))
	}
	ab <- paste(ab, collapse = ",")
}else{
	ab <- NA
}

print(ab)

r2_h2[j, 14] <- ab
}
}
}


newcolnamesr2h2 <- c("gene", "lasso.top1_afr", "lasso.top1_meta","lasso.top1_magepro", "lasso.top1_afr_pv", "lasso.top1_meta_pv","lasso.top1_magepro_pv", "hsq_afr", "hsq_afr_se", "hsq_afr.pv", "datasets", "nonzeroAFR", "nonzeroMAGEPRO", "alpha_beta")
colnames(r2_h2) <- newcolnamesr2h2

h <- which(r2_h2[,2] == 0)
if(length(h) > 0){r2_h2 <- r2_h2[-h,]}


write.table(r2_h2, file = "/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO/process_results/MAGEPRO_r2_h2.txt", row.names = F, col.names = T, sep = "\t", quote = F)
}

