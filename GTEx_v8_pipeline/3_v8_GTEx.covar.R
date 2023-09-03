library(data.table)

#--- read in command-line args
args <- commandArgs(trailingOnly = TRUE)
covar_file <- args[1]
intermediate_dir <- args[2]

#--- matching covar to people
print("MATCHING COVARS DATA TO PEOPLE")
covar <- fread(covar_file, header = T)
#format needs people on the rows, col1 = fam col 1, col2 = fam col 2, rest = covar
covar_t <- t(covar) #transpose matrix to get people on row 1
covar_t <- cbind(0,rownames(covar_t),covar_t)  
rownames(covar_t) <- NULL
colnames(covar_t) <- covar_t[1,] #add first row as column names 
covar_t <- covar_t[-1,] #make sure first row is deleted - it is now the column name
ind <- fread(paste0(intermediate_dir, "/All_Individuals.txt"), header = F)$V1 #read in people ID
covar_ind <- covar_t[,2] #column 2 has people IDS #colnames(covar)[-1]
w <- which(covar_ind %in% ind) #get indices where people ID match
covar_DS <- covar_t[w,] #extract the indices where people ID match  #cbind(covar[,1],covar_mat[,w])
if(nrow(covar_DS) < 2*(ncol(covar)-3)){covar_DS <- covar_DS[,-grep("InferredCov",colnames(covar_DS))[-c(1:5)]]}

#--- write output
write.table(covar_DS, file=paste0( intermediate_dir, "/Covar_All.txt"),row.names = F, col.names = T, sep = "\t", quote = F)
