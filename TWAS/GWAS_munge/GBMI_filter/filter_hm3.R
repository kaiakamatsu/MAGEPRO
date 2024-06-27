library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pop <- args[1]

traits <- c("Asthma","COPD","Gout","HF","IPF","Stroke","VTE")

hm3 <- fread("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.bim", header = F)

for(trait in traits){

	sumstats <- fread( file = paste0(pop, "/", trait, "_Bothsex_", pop,"_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"), header = T)
	sumstats_filtered <- filter(sumstats, rsid %in% hm3$V2)
	out_file <- paste0(pop, "/filtered_hm3/", trait ,"_", pop, "_GBMI.txt")
	print(paste0(trait, ", number of snps: ", nrow(sumstats_filtered)))
	write.table(sumstats_filtered, file = out_file, quote = F, row.names = F, col.names = T, sep = '\t')

}
