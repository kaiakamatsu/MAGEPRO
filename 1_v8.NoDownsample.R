library(data.table)

#tissues <- c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_EBV-transformed_lymphocytes","Cells_Cultured_fibroblasts","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood")

tissues <- "Whole_Blood"  #just whole blood

pop_people <- fread("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/GTEx_plink/GTEx_v8_genotype_AFR_HM3_exclude_dups.22.fam", header = F)$V2 #read in just the second column of the fam file (people id) 
for (i in 1:length(tissues)){ #just whole blood for this analysis
ge <- fread(paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/",tissues[i],".v8.normalized_expression.bed.gz"), header = T, nrows=1) #first row has all people ids 
all_people_with_geneexp <- colnames(ge)[-c(1:4)] #a vector of column names EXCEPT 1-4 (exclude #chr, start, end, gene_id columns) 
pop_people_with_geneexp <- intersect(all_people_with_geneexp,pop_people) #finds people common in both ge and geno data 
write.table(pop_people_with_geneexp,file=paste0("/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/All_Individuals_",tissues[i],".txt"), row.names = F, col.names = F, quote = F)
}

#repeat for EAS data?

