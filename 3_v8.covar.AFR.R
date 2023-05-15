library(data.table)


#tissues <- c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_EBV-transformed_lymphocytes","Cells_Cultured_fibroblasts","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood")

tissues <- "Whole_Blood"

#matching covar to people
for (i in 1:length(tissues)){
covar <- fread(paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/",tissues[i],".v8.covariates.txt"), header = T)
#Sasha's format needs people on the rows, col1 = fam col 1, col2 = fam col 2, rest = covar
covar_t <- t(covar) #transpose matrix to get people on row 1
covar_t <- cbind(0,rownames(covar_t),covar_t) #add a column of 0 and a column of rownames 
rownames(covar_t) <- NULL
colnames(covar_t) <- covar_t[1,] #add first row as column names 
covar_t <- covar_t[-1,] #make sure first row is deleted - it is now the column name
ind <- fread(paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/All_Individuals_",tissues[i],".txt"), header = F)$V1 #read in people ID
covar_ind <- covar_t[,2] #column 2 has people IDS #colnames(covar)[-1]
w <- which(covar_ind %in% ind) #get indices where people ID match
covar_DS <- covar_t[w,] #extract the indices where people ID match  #cbind(covar[,1],covar_mat[,w])
if(nrow(covar_DS) < 2*(ncol(covar)-3)){covar_DS <- covar_DS[,-grep("InferredCov",colnames(covar_DS))[-c(1:5)]]} #just InferredCov 1-5
write.table(covar_DS, file=paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/Covar_All_",tissues[i],".txt"),row.names = F, col.names = T, sep = "\t", quote = F)
}

