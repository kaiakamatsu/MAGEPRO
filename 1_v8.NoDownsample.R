library(data.table)

tissues <- "Whole_Blood"  #just whole blood

pop_people <- fread("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO/data/GTEx_plink/GTEx_v8_genotype_AFR_HM3_exclude_dups.22.fam", header = F)$V2 #read in just the second column of the fam file (people id) 
for (i in 1:length(tissues)){ 
ge <- fread(paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO/data/",tissues[i],".v8.normalized_expression.bed.gz"), header = T, nrows=1) #first row has all people ids 
all_people_with_geneexp <- colnames(ge)[-c(1:4)] #a vector of column names EXCEPT 1-4 (exclude #chr, start, end, gene_id columns) 
pop_people_with_geneexp <- intersect(all_people_with_geneexp,pop_people) #finds people common in both ge and geno data 
write.table(pop_people_with_geneexp,file=paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/MAGEPRO/intermedfiles/All_Individuals_",tissues[i],".txt"), row.names = F, col.names = F, quote = F)
}

