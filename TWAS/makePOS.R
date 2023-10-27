library(data.table)
library(stringr)
library(dplyr)

genebed <- "/expanse/lustre/projects/ddp412/kakamatsu/GENE_MODELS/GTExAFR_WHOLEBLOOD_MAGEPRO/scratch/gene_beds/"
home <- "/expanse/lustre/projects/ddp412/kakamatsu/fusion_multipop"

paths <- fread(paste0(home, "/magepro_weights_paths.txt"), header = FALSE)

output <- matrix(0, nrow = nrow(paths), ncol = 5)
colnames(output) = c("WGT", "ID", "CHR", "P0", "P1")

for (n in 1:length(paths$V1)){
	wgt <- strsplit(paths$V1[n], "/")[[1]][10]
	genebase <- strsplit(wgt, "[.]")[[1]][1]
	genever <- strsplit(wgt, "[.]")[[1]][2]
	gene <- paste0(genebase, ".", genever)
	print(gene)
	bed <- fread(paste0(genebed, gene, ".bed"), header = FALSE)
	output[n, 1] <- paste0(gene, ".wgt.RDat")
	output[n, 2] <- bed$V4[1]
	output[n, 3] <- bed$V1[1]
	output[n, 4] <- bed$V2[1]
	output[n, 5] <- bed$V3[1]
	print(output[n, ])
}


write.table(output, file = paste0(home,"/magepro_weight.pos"), row.names = F, col.names = T, sep = "\t", quote = F)
