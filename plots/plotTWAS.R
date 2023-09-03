install.packages("devtools")
install.packages("R.utils")
install.packages('ggplot2')
install.packages('dplyr')
library(ggplot2)
library(dplyr)
library(data.table)
library(R.utils)
library(devtools)

install.packages("VennDiagram")
library(VennDiagram)

# --- analysis for tiffany 




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
