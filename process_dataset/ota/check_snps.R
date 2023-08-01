#script to check if the same snp-gene pairs are tested in every cell type - helps with meta-analysis later
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
gene = args[1]

dir = "/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/OTA_nominal/genes"

cell_types = c("CD16p_Mono", "CL_Mono", "LDG", "mDC", "Mem_CD4", "Mem_CD8", "Naive_B", "Naive_CD4", "Naive_CD8", "Neu", "NK", "pDC", "Plasmablast")
cell_avail = c()

for (c in cell_types){
	file = paste0(dir, "/", c, "/", gene, ".txt")
	if (file.exists(file)){
		cell_avail <- append(cell_avail, c)
	}
}

for (i in 1:length(cell_avail)){
	filename = paste0(dir, "/", cell_avail[i], "/", gene, ".txt")
	df <- fread(filename, header = T)
	print(nrow(df)) #same num rows
	if (i < length(cell_avail)){
		filename2 = paste0(dir, "/", cell_avail[(i+1)], "/", gene, ".txt")
		df2 <- fread(filename2, header = T)
		print(identical(df$Variant_ID, df2$Variant_ID)) #true 
	}
}
