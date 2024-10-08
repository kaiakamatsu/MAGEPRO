suppressMessages(library(data.table))
suppressMessages(library(dplyr))

#--- read in command line arguments 
args <- commandArgs(trailingOnly = TRUE)
fam_file <- paste0(args[1], "1.fam") #path to plink files containing all individuals of interest 
ge_file <- args[2] #path to normalized gene expression file with people ids on first row
intermed_dir <- args[3] #directory to store intermediate files - will store all people ids 

#--- extract all people with both ge and genotype data
fam <- fread(fam_file, header = F) #read in just the second column of the fam file (people id) 
pop_people <- fam$V2
ge <- fread(ge_file, header = T, nrows=1) #first row has all people ids 
all_people_with_geneexp <- colnames(ge)[-c(1:4)] #a vector of column names EXCEPT 1-4 (exclude #chr, start, end, gene_id columns) 
people_interest <- fam %>%
  filter(V2 %in% all_people_with_geneexp) %>%
  arrange(match(V2, all_people_with_geneexp))
print(paste0( "number of individuals with both genotype and ge data: " , nrow(people_interest)) )
pop_people_with_geneexp <- people_interest[,c(1,2)]
write.table(pop_people_with_geneexp,file=paste0(intermed_dir, "/All_Individuals.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
