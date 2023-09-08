library(data.table)
library(dplyr)

#--- read in command-line args
args <- commandArgs(trailingOnly = TRUE)
ge_file <- args[1]
gene_beds <- args[2]
intermed <- args[3]
subset <- args[4]

#--- parsing through ge file to grab TSS information for each gene
ensg <- data.frame() 
ge <- fread(ge_file, header = T)
chrs <- sapply(1:nrow(ge), function(x) strsplit(ge$`#chr`[x], split = "chr")[[1]][2]) #modify the #chr column: "1" instead of "chr1" 
ge$`#chr` <- chrs
wremove <- which(chrs == "X") #remove X chr genes, Y is not included  
if(length(wremove)>0){ge <- ge[-wremove,]}
dump <- ge[,1:4] #grabbing just the first 4 columns (#chr, start, end, gene_id)
ensg <- rbind(ensg,dump) #add it to the ensg
ensg <- unique(ensg) #delete duplicates 
if (subset != 'NA'){
subset_genes <- fread(subset, header = F)
ensg <- filter(ensg, gene_id %in% subset_genes$V1)
}
print( paste0( "number of genes with ge data in analysis: ", nrow(ensg))) #19696 for whole_blood 
write.table(ensg[,4], paste0(intermed, "/All_Genes_Expressed.txt"), row.names = F, col.names = F, sep ="\t", quote = F)

#--- write files detailing the cis-region bp of every gene
for (i in 1:nrow(ensg)){  
bounds_l <- max(0,ensg$start[i] - 500000) #left bound = whichever is bigger 0 or start - 500kb
bounds_r <- ensg$end[i] + 500000 #it's not the gene end, it's just the transcription start site (TSS), so we have done +/- 500kb of gene TSS, because GTEx doesn't give us the end windows to look for SNPs based on genes 
chr <- ensg$`#chr`[i] #store chr #
bed <- c(chr,bounds_l,bounds_r,ensg$gene_id[i],0,"+")
bed <- matrix(bed, nrow = 1, ncol = 6) #make a 1 row x 6 column matrix 
#make bed file 
write.table(bed, paste0(gene_beds, "/",ensg$gene_id[i],".bed"), row.names = F, col.names = F, sep ="\t", quote = F)
}
#note: some genes do not have variants. 
