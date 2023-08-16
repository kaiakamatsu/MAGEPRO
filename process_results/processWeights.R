library(data.table)
library(dplyr)


tissues <- "Whole_Blood"
args <- commandArgs(trailingOnly = TRUE)
output <- args[1] #path to output directory
genes_assign <- args[2] #path to genes assigned file

for (i in 1:length(tissues)){
genes <- read.table(file = genes_assign, header = F)$V1 #all genes in analysis

r2_h2 <- matrix(0, length(genes), 15)

for (j in 1:length(genes)){

file <- paste0(output,"/",genes[j],".wgt.RDat") 

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

r2_h2[j, 15] <- paste(total_coeff, collapse = ",")

}
}
}


newcolnamesr2h2 <- c("gene", "lasso.top1_eur", "lasso.top1_meta","lasso.top1_magepro", "lasso.top1_eur_pv", "lasso.top1_meta_pv","lasso.top1_magepro_pv", "hsq_eur", "hsq_eur_se", "hsq_eur.pv", "datasets", "nonzeroEUR", "nonzeroMAGEPRO", "alpha_beta", "alphas")
colnames(r2_h2) <- newcolnamesr2h2

h <- which(r2_h2[,2] == 0)
if(length(h) > 0){r2_h2 <- r2_h2[-h,]}


write.table(r2_h2, file = "MAGEPRO_results.txt", row.names = F, col.names = T, sep = "\t", quote = F)
}

