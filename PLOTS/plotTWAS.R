install.packages("devtools")
install.packages("R.utils")
install.packages("ggplot2")
install.packages("data.table")
install.packages("dplyr", verbose=T)
library(ggplot2)
library(dplyr)
library(data.table)
library(R.utils)
library(devtools)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

install.packages("VennDiagram")
library(VennDiagram)

# --- analysis for tiffany 

# create a set of genes either significantly heritable in EUR or AFR 

stats_afr <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_fullsumstats_r2_h2.txt', as.is = T, header = T)
stats_eur <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_results_eur.txt', as.is = T, header = T)

colnames(stats_afr)
stats_afr_sigh2 <- filter(stats_afr, hsq_afr.pv < 0.01 & hsq_afr > 0)

colnames(stats_eur)
stats_eur_sigh2 <- filter(stats_eur, hsq_eur.pv < 0.01 & hsq_eur > 0)

afr_genes <- stats_afr_sigh2$gene
eur_genes <- stats_eur_sigh2$gene
length(afr_genes)
length(eur_genes)

result <- union(afr_genes, eur_genes)
bonfer <- length(result)
thresh <- 0.05/bonfer

genes <- as.data.frame(result)
write.table(genes, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/genes_sigh2_either.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

#https://www.biotools.fr/human/ensembl_symbol_converter
#genes_names <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/sigh2_AFR_EA_name.txt', as.is = T, header = F, sep = "\t")
#colnames(genes_names) <- c("ID", "NAME")

#try biomaRt
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensemble <- c()
for (r in genes$result){
  ensemble <- append(ensemble, strsplit(r, split="[.]")[[1]][1])
}
genes_list = ensemble
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes_list,mart= mart)
genes <- cbind(genes, ensemble)
colnames(genes) <- c("ID", "ensembl_gene_id")
genes_table <- merge(genes, G_list, by = "ensembl_gene_id", all = TRUE)
nrow(genes_table)

# find ancestry specific TWAS hits (bonferroni)

phenos <- c("BAS", "HGB", "MCHC", "MPV", "RBC", "EOS", "LYM", "MCV", "NEU", "RDW", "HCT", "MCH", "MON", "PLT", "WBC")
cols <- NA
sig_magepro <- matrix(nrow = 0, ncol = 21)
for (pheno in phenos){
  twasdata <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputMAGEPRO_BCX_LETTER/", pheno,".combined.txt"), header= T)
  twasdata_sig_genes <- filter(twasdata, ID %in% genes_table$ID)
  twasdata_sig <- filter(twasdata_sig_genes, TWAS.P < thresh)
  PHENO <- rep(pheno, times = nrow(twasdata_sig))
  sig <- cbind(twasdata_sig, PHENO)
  if (nrow(sig) > 0){
    sig_magepro <- rbind(sig_magepro, as.matrix(sig))
    cols <<- colnames(sig)
  }
}

magepro <- as.data.frame(sig_magepro)
colnames(magepro) <- cols
nrow(magepro) # 144 which is more than 116 reported by MA-FOCUS

magepro_names <- filter(genes_table, ID %in% magepro$ID)
combined_magepro <- merge(magepro, magepro_names, by = "ID", all = TRUE)
nrow(combined_magepro)

write.table(combined_magepro, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/MAGEPRO_bonferroni_sigh2_either.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)

# AFR 

phenos <- c("BAS", "HGB", "MCHC", "MPV", "RBC", "EOS", "LYM", "MCV", "NEU", "RDW", "HCT", "MCH", "MON", "PLT", "WBC")
cols <- NA
sig_afr <- matrix(nrow = 0, ncol = 21)
for (pheno in phenos){
  twasdata <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputAFR_BCX_LETTER/", pheno,".combined.txt"), header= T)
  twasdata_sig_genes <- filter(twasdata, ID %in% genes_table$ID)
  twasdata_sig <- filter(twasdata_sig_genes, TWAS.P < thresh)
  PHENO <- rep(pheno, times = nrow(twasdata_sig))
  sig <- cbind(twasdata_sig, PHENO)
  if (nrow(sig) > 0){
    sig_afr <- rbind(sig_afr, as.matrix(sig))
    cols <<- colnames(sig)
  }
}

afr <- as.data.frame(sig_afr)
colnames(afr) <- cols
nrow(afr) #131

afr_names <- filter(genes_table, ID %in% afr$ID)
combined_afr <- merge(afr, afr_names, by = "ID", all = TRUE)

write.table(combined_afr, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/AFR_bonferroni_sigh2_either.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)


# EUR model 

phenos <- c("BAS", "HGB", "MCHC", "MPV", "RBC", "EOS", "LYM", "MCV", "NEU", "RDW", "HCT", "MCH", "MON", "PLT", "WBC")
cols <- NA
sig_eur <- matrix(nrow = 0, ncol = 21)
for (pheno in phenos){
  twasdata <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputEUR_BCX_LETTER/", pheno,".combined.txt"), header= T)
  twasdata_sig_genes <- filter(twasdata, ID %in% genes_table$ID)
  twasdata_sig <- filter(twasdata_sig_genes, TWAS.P < thresh)
  PHENO <- rep(pheno, times = nrow(twasdata_sig))
  sig <- cbind(twasdata_sig, PHENO)
  if (nrow(sig) > 0){
    sig_eur <- rbind(sig_eur, as.matrix(sig))
    cols <<- colnames(sig)
  }
}

eur <- as.data.frame(sig_eur)
colnames(eur) <- cols
nrow(eur) #6827 compared to 6,236 in MA-FOCUS

eur_names <- filter(genes_table, ID %in% eur$ID)
combined_eur <- merge(eur, eur_names, by = "ID", all = TRUE)
nrow(combined_eur)

write.table(combined_eur, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/EUR_bonferroni_sigh2_either.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)

# compare our TWAS hits to MA-FOCUS 

load("C:/Users/kaiak/Downloads/twas.RData") #https://github.com/mancusolab/MA-FOCUS-data-code/blob/main/real-data/data/twas.RData
df_mafocus <- as.data.frame(twas_all)

unique(df_mafocus$POP)
df_eur_mafocus <- filter(df_mafocus, POP == "EA")
df_afr_mafocus <- filter(df_mafocus, POP == "AA")

eur_mafocus <- split(df_eur_mafocus, df_eur_mafocus$PHEN)
afr_mafocus <- split(df_afr_mafocus, df_afr_mafocus$PHEN)

# what proportion of our hits are also hits in MA-FOCUS? (p < 0.05/4579)
phenos <- c("BAS", "HGB", "MCHC", "MPV", "RBC", "EOS", "LYM", "MCV", "NEU", "RDW", "HCT", "MCH", "MON", "PLT", "WBC")
#AFR(MAGEPRO)
percent_common <- c() #0.00000000 0.00000000 0.00000000 0.20000000 0.33333333 0.07692308 0.00000000 0.37500000 0.00000000 0.33333333 0.06382979
for (pheno in phenos){
  twasdata <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputMAGEPRO_BCX_LETTER/", pheno,".combined.txt"), header= T)
  twasdata_sig_genes <- filter(twasdata, ID %in% genes_table$ID)
  twasdata_sig <- filter(twasdata_sig_genes, TWAS.P < thresh)
  if (nrow(twasdata_sig) > 0){
    
    names <- filter(genes_table, ID %in% twasdata_sig$ID)
    combined <- merge(twasdata_sig, names, by = "ID", all = TRUE)
    
    comparison <- afr_mafocus[[pheno]]
    comparison_sig <- filter(comparison, TWAS.P < (0.05/4579))
    
    total <- nrow(combined)
    common <- length(which(combined$hgnc_symbol %in% comparison_sig$ID))
    
    percent_common <- append(percent_common, common/total)
  }
}
#EUR
percent_common <- c() #0.1592357 0.2094241 0.1709845 0.1876209 0.1970534 0.1876380 0.2004008 0.1953642 0.1760870 0.1926782 0.2170330 0.1780000 0.2054507 0.2116317 0.2014787
for (pheno in phenos){
  twasdata <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputEUR_BCX_LETTER/", pheno,".combined.txt"), header= T)
  twasdata_sig_genes <- filter(twasdata, ID %in% genes_table$ID)
  twasdata_sig <- filter(twasdata_sig_genes, TWAS.P < thresh)
  if (nrow(twasdata_sig) > 0){
    
    names <- filter(genes_table, ID %in% twasdata_sig$ID)
    combined <- merge(twasdata_sig, names, by = "ID", all = TRUE)
    
    comparison <- eur_mafocus[[pheno]]
    comparison_sig <- filter(comparison, TWAS.P < (0.05/4579))
    
    total <- nrow(combined)
    common <- length(which(combined$hgnc_symbol %in% comparison_sig$ID))
    
    percent_common <- append(percent_common, common/total)
  }
}

# subset to hits in common and check correlation? 

# across all results, what is the correlation? 
#Afr
phenos <- c("BAS", "HGB", "MCHC", "MPV", "RBC", "EOS", "LYM", "MCV", "NEU", "RDW", "HCT", "MCH", "MON", "PLT", "WBC")
cors <- c() #0.1673255 0.1809370 0.2019416 0.1706386 0.1638782 0.1609890 0.2665420 0.1873377 0.2079294 0.1773886 0.1660443 0.1822299 0.2042523 0.1986272 0.2338586
for (p in phenos){
  twasdata <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputMAGEPRO_BCX_LETTER/", p,".combined.txt"), header= T)
  names <- filter(genes_table, ID %in% twasdata$ID)
  combined <- merge(twasdata, names, by = "ID", all = TRUE)
  comparison <- afr_mafocus[[p]]
  # filter for genes in common and reorder both to be in the same order
  comparison = filter(comparison, ID %in% combined$hgnc_symbol)
  combined = filter(combined, hgnc_symbol %in% comparison$ID)
  m1 <- match(comparison$ID, combined$hgnc_symbol)
  combined <- combined[m1, ]
  df <- data.frame(combined$TWAS.Z, comparison$TWAS.Z)
  df <- na.omit(df)
  colnames(df) <- c("MAGEPRO","MAFOCUS")
  cors <- append(cors, cor(df$MAGEPRO, df$MAFOCUS))
}

#EUR
phenos <- c("BAS", "HGB", "MCHC", "MPV", "RBC", "EOS", "LYM", "MCV", "NEU", "RDW", "HCT", "MCH", "MON", "PLT", "WBC")
cors <- c() #0.11659153 0.17198161 0.25387615 0.19224023 0.18605941 0.05944016 0.14975165 0.31742256 0.10850904 0.19631324 0.19116171 0.27288473 0.05538389 0.16896935 0.15632017
for (p in phenos){
  twasdata <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputEUR_BCX_LETTER/", p,".combined.txt"), header= T)
  names <- filter(genes_table, ID %in% twasdata$ID)
  combined <- merge(twasdata, names, by = "ID", all = TRUE)
  comparison <- eur_mafocus[[p]]
  # filter for genes in common and reorder both to be in the same order
  comparison = filter(comparison, ID %in% combined$hgnc_symbol)
  combined = filter(combined, hgnc_symbol %in% comparison$ID)
  m1 <- match(comparison$ID, combined$hgnc_symbol)
  combined <- combined[m1, ]
  df <- data.frame(combined$TWAS.Z, comparison$TWAS.Z)
  df <- na.omit(df)
  colnames(df) <- c("EURour","MAFOCUS")
  cors <- append(cors, cor(df$EURour, df$MAFOCUS))
}

#--- ancestry specific hits 

#all TWAS eur
eur <- matrix(nrow = 0, ncol = 21)
cols <- NA
for (pheno in phenos){
  twasdata <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputEUR_BCX_LETTER/", pheno,".combined.txt"), header= T)
  twasdata_sig_genes <- filter(twasdata, ID %in% genes_table$ID)
  PHENO <- rep(pheno, times = nrow(twasdata_sig_genes))
  all <- cbind(twasdata_sig_genes, PHENO)
  if (nrow(all) > 0){
    eur <- rbind(eur, as.matrix(all))
    cols <<- colnames(all)
  }
}
eur <- as.data.frame(eur)
colnames(eur) <- cols
eursplit <- split(eur, eur$PHENO)

aa_specific <- combined_magepro[which(! combined_magepro$ID %in% combined_eur$ID),] #in magepro afr but not in eur
aa_specific <- combined_magepro[which( (! combined_magepro$ID %in% combined_eur$ID) & (! combined_magepro$ID %in% combined_afr$ID) ),] #in magepro afr, not in afronly and eur
TWASeur_p <- c()
TWASeur_z <- c()
GWASeur_id <- c()
GWASeur_z <- c()
EQTLeur.ID <- c()
EQTLeur.R2 <- c()
EQTLeur.Z <- c()

for(a in 1:nrow(aa_specific)){
  pheno <- aa_specific$PHENO[a]
  eur_data <- eursplit[[pheno]]
  w <- which(eur_data$ID == aa_specific$ID[a])
  TWASeur_p <- append(TWASeur_p, eur_data$TWAS.P[w])
  TWASeur_z <- append(TWASeur_z, eur_data$TWAS.Z[w])
  GWASeur_id <- append(GWASeur_id, eur_data$BEST.GWAS.ID[w])
  GWASeur_z <- append(GWASeur_z, eur_data$BEST.GWAS.Z[w])
  EQTLeur.ID <- append(EQTLeur.ID, eur_data$EQTL.ID[w])
  EQTLeur.R2 <- append(EQTLeur.R2, eur_data$EQTL.R2[w])
  EQTLeur.Z <- append(EQTLeur.Z, eur_data$EQTL.Z[w])
}

aa_specific_with_eur <- cbind(aa_specific, EQTLeur.ID, EQTLeur.R2, EQTLeur.Z, GWASeur_id, GWASeur_z, TWASeur_z, TWASeur_p)
nrow(aa_specific_with_eur)
rownames(aa_specific_with_eur) <- NULL

write.table(aa_specific_with_eur, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPROaa_specific_gene_trait_assoc.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)


# --- MANHATTAN 
twasdata <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputMAGEPRO/BMI.combined.txt", header= T)
colnames(twasdata)

genes <- twasdata$ID
chr <- twasdata$CHR
pos <- unlist(lapply(twasdata$P0, function(x) x + 500000))
model <- twasdata$MODEL
r2 <- twasdata$MODELCV.R2
z <- twasdata$TWAS.Z
p <- twasdata$TWAS.P
cat <- c(rep("", nrow(twasdata)))

df <- data.frame(genes, chr, pos, model, r2, z, p, cat)

for (i in 1:nrow(df)){
  if ((df$chr[i] %% 2) == 0){
    df$cat[i] = "even"
  }else{
    df$cat[i] = "odd"
  }
}

#---compute Z threshold 
df2 <- filter(df, r2 > 0)
p <- p.adjust(df$p, method = "fdr")
w <- which(p < 0.05)
df2[w, ]
threshold <- max(df2$p[w])
z <- abs(qnorm(threshold))
print(z)
nrow(df2)
df2$genes[w]
df2$z[w]


#POSITIVE AND NEGATIVE Z?
ggplot(df2, aes(x = pos, y = z, color = cat)) + 
  geom_point(alpha = 1, size = 1.5, position = position_jitter(width = 1000000)) +
  facet_wrap( ~ chr, nrow = 1, strip.position="bottom" , scales = "free_x") + 
  theme_bw() + 
  theme(strip.background=element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust=-0.5), panel.spacing = unit(0, "lines"), panel.border = element_blank(), axis.text.x=element_blank(), legend.position="none") + 
  labs(x = "Gene Position", y = "TWAS Z Score", title = "FLEX AFR White Blood Cell Count") + 
  scale_colour_manual(values = c("black", "#0B3295")) + 
  geom_hline(yintercept=0, linetype="solid") +
  geom_hline(yintercept=z, linetype="dashed", color = "#5B5E65", linewidth = 1.5) + 
  geom_hline(yintercept=-z, linetype="dashed", color = "#5B5E65", linewidth = 1.5) +
  scale_y_continuous(breaks = seq(-8, 8, by = 1), limits = c(-8, 8))+
  coord_cartesian(clip = "off")
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/TWASAFR_WhiteCount.png", width = 10, height = 8)



#--- compute number of associations found in each TWAS 

#MAGEPRO

#phenos <- c("Diabetes", "Hypertension", "MonocyteCount", "RedBloodCellCount", "BMI", "EosinophilCount", "LymphocyteCount", "PlateletCount", "RedBloodCellDistributionWidth", "WhiteBloodCellCount")
phenos <- c("HGB", "MCHC", "MPV", "RBC", "EOS", "LYM", "MCV", "NEU", "RDW", "HCT", "MCH", "MON", "PLT", "WBC")
count = 0

#run first pheno, set up dataframe 
twasdata <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputMAGEPRO_BCX_LETTER/BAS.combined.txt", header= T)
genes <- twasdata$ID
chr <- twasdata$CHR
pos <- unlist(lapply(twasdata$P0, function(x) x + 500000))
model <- twasdata$MODEL
r2 <- twasdata$MODELCV.R2
z <- twasdata$TWAS.Z
p <- twasdata$TWAS.P
gwasid <- twasdata$BEST.GWAS.ID
gwasz <- twasdata$BEST.GWAS.Z
df <- data.frame(genes, chr, pos, model, r2, z, p, gwasid, gwasz)
df <- na.omit(df)
df2 <- filter(df, r2 > 0)
p <- p.adjust(df2$p, method = "fdr")
w <- which(p < 0.05)
threshold <- max(df2$p[w])
z <- abs(qnorm(threshold))
df3 <- df2[w, ]
addpheno <- rep(c("BAS"), each=nrow(df3))
addzthresh <- rep(c(z), each = nrow(df3))
result <- data.frame(addpheno, df3, addzthresh)
colnames(result) <- c("pheno", "gene", "chr", "pos", "model", "r2", "z", "p", "bestgwasid", "bestgwasz","z_threshold")


for (pheno in phenos){
  twasdata <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputMAGEPRO_BCX_LETTER/", pheno,".combined.txt"), header= T)
  genes <- twasdata$ID
  chr <- twasdata$CHR
  pos <- unlist(lapply(twasdata$P0, function(x) x + 500000))
  model <- twasdata$MODEL
  r2 <- twasdata$MODELCV.R2
  z <- twasdata$TWAS.Z
  p <- twasdata$TWAS.P
  gwasid <- twasdata$BEST.GWAS.ID
  gwasz <- twasdata$BEST.GWAS.Z
  df <- data.frame(genes, chr, pos, model, r2, z, p, gwasid, gwasz)
  #filtering and count significantly associated genes
  df2 <- filter(df, r2 > 0)
  p <- p.adjust(df2$p, method = "fdr")
  w <- which(p < 0.05)
  if (length(w) != 0){
    count = count+length(w)
    #record significant genes
    threshold <- max(df2$p[w])
    z <- abs(qnorm(threshold))
    df3 <- df2[w, ]
    addpheno <- rep(c(pheno), each=nrow(df3))
    addzthresh <- rep(c(z), each = nrow(df3))
    df4 <- data.frame(addpheno, df3, addzthresh)
    colnames(df4) <- c("pheno", "gene", "chr", "pos", "model", "r2", "z", "p", "bestgwasid", "bestgwasz","z_threshold")
    result <- rbind(result, df4)
  }
}

nrow(result) #57

write.table(result, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/MAGEPRO_significant.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)

#AFR

phenos <- c("HGB", "MCHC", "MPV", "RBC", "EOS", "LYM", "MCV", "NEU", "RDW", "HCT", "MCH", "MON", "PLT", "WBC")
count = 0

#run first pheno, set up dataframe 
twasdata <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputAFR_BCX_LETTER/BAS.combined.txt", header= T)
genes <- twasdata$ID
chr <- twasdata$CHR
pos <- unlist(lapply(twasdata$P0, function(x) x + 500000))
model <- twasdata$MODEL
r2 <- twasdata$MODELCV.R2
z <- twasdata$TWAS.Z
p <- twasdata$TWAS.P
df <- data.frame(genes, chr, pos, model, r2, z, p)
df2 <- filter(df, r2 > 0)
p <- p.adjust(df2$p, method = "fdr")
w <- which(p < 0.05)
count = count + length(w)
threshold <- max(df2$p[w])
z <- abs(qnorm(threshold))
df3 <- df2[w, ]
addpheno <- rep(c("BAS"), each=nrow(df3))
addzthresh <- rep(c(z), each = nrow(df3))
result2 <- data.frame(addpheno, df3, addzthresh)
colnames(result2) <- c("pheno", "gene", "chr", "pos", "model", "r2", "z", "p", "z_threshold")


for (pheno in phenos){
  twasdata <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputAFR_BCX_LETTER/", pheno,".combined.txt"), header= T)
  genes <- twasdata$ID
  chr <- twasdata$CHR
  pos <- unlist(lapply(twasdata$P0, function(x) x + 500000))
  model <- twasdata$MODEL
  r2 <- twasdata$MODELCV.R2
  z <- twasdata$TWAS.Z
  p <- twasdata$TWAS.P
  df <- data.frame(genes, chr, pos, model, r2, z, p)
  #filtering and count significantly associated genes
  df2 <- filter(df, r2 > 0)
  p <- p.adjust(df2$p, method = "fdr")
  w <- which(p < 0.05)
  if (length(w) != 0){
    count = count+length(w)
    #record significant genes
    threshold <- max(df2$p[w])
    z <- abs(qnorm(threshold))
    df3 <- df2[w, ]
    addpheno <- rep(c(pheno), each=nrow(df3))
    addzthresh <- rep(c(z), each = nrow(df3))
    df4 <- data.frame(addpheno, df3, addzthresh)
    colnames(df4) <- c("pheno", "gene", "chr", "pos", "model", "r2", "z", "p", "z_threshold")
    result2 <- rbind(result2, df4)
  }
}


nrow(result2) #40


write.table(result2, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/AFR_significant.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)

length(unique(result$pheno)) #found associations for 13 phenotypes
length(unique(result2$pheno)) #found associations for 4 phenotypes 

# --- venn diagram between FLEX and AFR

FLEX <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/MAGEPRO_significant.txt", header= T)
AFR <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/AFR_significant.txt", header= T)

flexgroup <- split(FLEX, FLEX$pheno)
afrgroup <- split(AFR, AFR$pheno)

nrow(AFR)
nrow(FLEX)

numflex = 0 #26
numafr = 0 #9
numcommon = 0 #32

phenos <- c("Asthma", "Diabetes", "Hypertension", "MonocyteCount", "RedBloodCellCount", "BMI", "EosinophilCount", "LymphocyteCount", "PlateletCount", "RedBloodCellDistributionWidth", "WhiteBloodCellCount")

for (p in phenos){
  flex <- flexgroup[[p]]
  afr <- afrgroup[[p]]
  if (!is.null(afr)){
    w <- which(!(flex$gene %in% afr$gene))
    flex <- flex[w, ]
  }
  if ( !is.null(flex)){
    if (nrow(flex) > 0){
      numflex <<-sum(numflex, nrow(flex))
    }
  }
}

for (p in phenos){
  flex <- flexgroup[[p]]
  afr <- afrgroup[[p]]
  if (!is.null(afr)){
    w <- which(!(afr$gene %in% flex$gene))
    afr <- afr[w, ]
  }
  if (!is.null(afr)){
    if (nrow(afr) > 0){
        numafr <<-sum(numafr, nrow(afr))
    }
  }
}

for (p in phenos){
  flex <- flexgroup[[p]]
  afr <- afrgroup[[p]]
  if (!is.null(afr)){
    w <- which((flex$gene %in% afr$gene))
    flex <- flex[w, ]
  }
  if (!is.null(flex)){
    if (nrow(flex) > 0){
      numcommon <<-sum(numcommon, nrow(flex))
    }
  }
}

#EUR

phenos <- c("Diabetes", "Hypertension", "MonocyteCount", "RedBloodCellCount", "BMI", "EosinophilCount", "LymphocyteCount", "PlateletCount", "RedBloodCellDistributionWidth", "WhiteBloodCellCount")
count = 0

#run first pheno, set up dataframe 
twasdata <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputEUR/Asthma.combined.txt", header= T)
colnames(twasdata)

twasdata <- filter(twasdata, r2 > 0)
zscores <- twasdata$z
zscorespositive <- abs(zscores)
zscoresnegative <- zscorespositive*-1
pvalues <- pnorm(zscoresnegative)
p <- p.adjust(pvalues, method = "fdr")
w <- which(p < 0.05)
df <- twasdata[w, ]
addpheno <- rep(c("Asthma"), each=nrow(df))
result3 <- data.frame(addpheno, df)
colnames(result3) <- c("pheno", "gene", "z", "r2")


for (pheno in phenos){
  twasdata <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputEUR/", pheno,".combined.txt"), header= T)
  twasdata <- filter(twasdata, r2 > 0)
  zscores <- twasdata$z
  zscorespositive <- abs(zscores)
  zscoresnegative <- zscorespositive*-1
  pvalues <- pnorm(zscoresnegative)
  p <- p.adjust(pvalues, method = "fdr")
  w <- which(p < 0.05)
  df <- twasdata[w, ]
  addpheno <- rep(c(pheno), each=nrow(df))
  addframe <- data.frame(addpheno, df)
  colnames(addframe) <- c("pheno", "gene", "z", "r2")
  result3 <- rbind(result3, addframe)
}

print(nrow(result3)) #34917


write.table(result3, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/EUR_significant.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)

#---

#--- load in data and split
FLEX <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/MAGEPRO_significant.txt", header= T)
AFR <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/AFR_significant.txt", header= T)
EUR <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/EUR_significant.txt", header= T)

nrow(FLEX)
nrow(AFR)
nrow(EUR)

flexgroup <- split(FLEX, FLEX$pheno)
afrgroup <- split(AFR, AFR$pheno)
eurgroup <- split(EUR, EUR$pheno)

#---

#--- afr specific split by afronly vs magepro

numflex = 0 #19
numafr = 0 #0
numcommon = 0 #12

phenos <- c("Asthma", "Diabetes", "Hypertension", "MonocyteCount", "RedBloodCellCount", "BMI", "EosinophilCount", "LymphocyteCount", "PlateletCount", "RedBloodCellDistributionWidth", "WhiteBloodCellCount")

df_flex_unique <- data.frame()

for (p in phenos){
  flex <- flexgroup[[p]]
  afr <- afrgroup[[p]]
  eur <- eurgroup[[p]]
  if (!is.null(afr)){
    w <- which(!(flex$gene %in% afr$gene))
    flex <- flex[w, ]
  }
  if (!is.null(eur)){
    w <- which(!(flex$gene %in% eur$gene))
    flex <- flex[w, ]
  }
  if (!is.null(flex)){
    df_flex_unique <<- rbind(df_flex_unique, flex)
    numflex <<-sum(numflex + nrow(flex))
  }
}

write.table(df_flex_unique, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/MAGEPRO_notfoundAFRonly_ancestry_specific_significant.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)


#--- venn diagram across EUR, AFRonly, MAGEPRO

FLEX <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/MAGEPRO_significant.txt", header= T)
AFR <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/AFR_significant.txt", header= T)
EUR <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/TWAS/EUR_significant.txt", header= T)

nrow(FLEX)
nrow(AFR)
nrow(EUR)

flexgroup <- split(FLEX, FLEX$pheno)
afrgroup <- split(AFR, AFR$pheno)
eurgroup <- split(EUR, EUR$pheno)

numflex = 0 #19
numafr = 0 #7
numeur = 0 #34892
num_c_flexafr = 0 #15
num_c_flexeur = 0 #7
num_c_afreur = 0 #2
num_c_all = 0 #16

phenos <- c("Asthma", "Diabetes", "Hypertension", "MonocyteCount", "RedBloodCellCount", "BMI", "EosinophilCount", "LymphocyteCount", "PlateletCount", "RedBloodCellDistributionWidth", "WhiteBloodCellCount")

df_flex_unique <- data.frame()

for (p in phenos){
  flex <- flexgroup[[p]]
  afr <- afrgroup[[p]]
  eur <- eurgroup[[p]]
  if (!is.null(flex)){
    w <- which(!(afr$gene %in% flex$gene))
    afr <- afr[w, ]
  }
  if (!is.null(eur)){
    w <- which(!(afr$gene %in% eur$gene))
    afr <- afr[w, ]
  }
  if (!is.null(afr)){
    numafr <<-sum(numafr, nrow(afr))
  }
}


for (p in phenos){
  flex <- flexgroup[[p]]
  afr <- afrgroup[[p]]
  eur <- eurgroup[[p]]
  if (!is.null(flex)){
    w <- which(!(afr$gene %in% flex$gene))
    afr <- afr[w, ]
    if (!is.null(eur)){
      w <- which((afr$gene %in% eur$gene))
      afr <- afr[w, ]
      if (!is.null(afr)){
        num_c_afreur <-sum(num_c_afreur, nrow(afr))
      }
    }
  }
}


#make venn diagram 
# Set up the Venn diagram data
numflex <- 19
numafr <- 7
numeur <- 34892
num_c_flexafr <- 15
num_c_flexeur <- 7
num_c_afreur <- 2
num_c_all <- 16





#------


#---


#--- how many of the newly found ancestry specific flex genes don't have a GWAS hit in its cis window? 


#---check correlation between ancestries
FLEX <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputFLEX/WhiteBloodCellCount.combined.txt"), header= T)
EUR <- fread(paste0("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/outputEUR/WhiteBloodCellCount.combined.txt"), header= T)
FLEX <- filter(FLEX, MODELCV.R2 > 0)
EUR <- filter(EUR, r2 > 0)

flex <- FLEX[, c(3, 19)]
eur <- EUR[, c(1, 2)]

flex <- na.omit(flex)
eur <- na.omit(eur)

flex <- filter(flex, ID %in% eur$gene)
eur <- filter(eur, gene %in% flex$ID)

flexorder <- flex[order(flex$ID), ]
eurorder <- eur[order(eur$gene), ]

identical(flexorder$ID, eurorder$gene)

cor.test(flexorder$TWAS.Z, eurorder$z)$p.val

#Diabetes p = 0.2215715
#BMI p = 0.519237
#Hypertension p = 0.4594502
#RedBloodCellCount p = 0.7507574
#Asthma p = 0.07152377
#MonocyteCount p = 0.001277968
#EosinophilCount p = 0.1521763
#LymphocyteCount p = 0.1772624
#PlateletCount p = 0.206132
#RedBloodCellDistributionWidth p = 0.001446149
#WhiteBloodCellCount p = 0.9871662
