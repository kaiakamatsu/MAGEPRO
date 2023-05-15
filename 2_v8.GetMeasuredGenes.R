library(data.table)

tissues <- c("Whole_Blood")

ensg <- data.frame() 
for (i in 1:length(tissues)){
print(i)
ge <- fread(paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/",tissues[i],".v8.normalized_expression.bed.gz"), header = T)
chrs <- sapply(1:nrow(ge), function(x) strsplit(ge$`#chr`[x], split = "chr")[[1]][2]) #modify the #chr column: "1" instead of "chr1" 
ge$`#chr` <- chrs
wremove <- which(chrs == "X") #remove X chr genes, I checked with grep and Y is not included  
if(length(wremove)>0){ge <- ge[-wremove,]}
dump <- ge[,1:4] #grabbing just the first 4 columns (#chr, start, end, gene_id)
ensg <- rbind(ensg,dump) #add it to the ensg
}

ensg <- unique(ensg) #delete duplicates 
print(nrow(ensg)) #19696 for whole_blood 

#comment out if starting from scratch 
#already_existing_genes <- list.files("/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/gene_beds/") #list of files in gene_beds
#already_existing_genes <- gsub(already_existing_genes, pattern  = ".bed", replacement = "") #remove the .bed part
#w <- which(is.na(match(ensg$gene_id,already_existing_genes))) #match it with the gene_id just read #1492
#ensg <- ensg[w,] 

#ensg <- ensg[sample(nrow(ensg), size = 1000, replace = F),] #downsampling
#don't downsample, run whole pipeline. 
for (i in 1:nrow(ensg)){ #make this quick so just do 1000 random ones - 
print(i)
bounds_l <- max(0,ensg$start[i] - 500000) #left bound = whichever is bigger 0 or start - 500kb
bounds_r <- ensg$end[i] + 500000 #it's not the gene end, it's just the transcription start site (TSS), so we have done +/- 500kb of gene TSS, because GTEx doesn't give us the end windows to look for SNPs based on genes 
chr <- ensg$`#chr`[i] #store chr#
bed <- c(chr,bounds_l,bounds_r,ensg$gene_id[i],0,"+")
bed <- matrix(bed, nrow = 1, ncol = 6) #make a 1 row x 6 column matrix 
#make bed file 
write.table(bed, paste0("/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/gene_beds/",ensg$gene_id[i],".bed"), row.names = F, col.names = F, sep ="\t", quote = F)
}
#note: some genes do not have variants. 



