#load packages 
install.packages('ggplot2')
library(ggplot2)
install.packages('dplyr')
library(dplyr)
install.packages('data.table')
library(data.table)
library(RColorBrewer)
library(tidyr)

#below is the code I used to create every plot 

#r2 plot ------------------------------------------------------------------------------------------------------------
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_fullsumstats_r2_h2.txt', as.is = T, header = T)

#other cell types - whole blood external datasets
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_AFR_Muscle_results.txt', as.is = T, header = T)
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_LCLafr_results.txt', as.is = T, header = T)
#---

colnames(stats)
nrow(stats)#19602 genes 
#min(stats$hsq_afr[which(stats$hsq_afr.pv < 0.05 & stats$hsq_afr > 0)]) #0.064251
#filter out genes with...
    #negative heritability 
    #invalid gcta output (for invalid gcta genes, h2 set to 0.064251 -> max h2 with p < 0.05)
    #not present in any external datasets

#w <- which((stats$datasets != "pred.wgt") & (stats$hsq_afr.pv < 0.01))
#w <- which((stats$datasets != "pred.wgt"))
w <- which((stats$datasets != "pred.wgt") & (stats$hsq_afr > 0) & (stats$hsq_afr != 0.064251))
stats <- stats[w, ]
nrow(stats) #19567 genes in this comparison 
mean_h2 = mean(stats$hsq_afr) #0.06321698
stats_r2 <- stats[, 2:4] #extract columns of r2 values for afronly, meta, and magepro
#median(stats_r2$lasso.top1_magepro)

#quick t-test - does r2 increase with our model?
t.test(stats_r2$lasso.top1_magepro, stats_r2$lasso.top1_afr, paired = T, alternative = "greater")$p.value # p = 0 
t.test(stats_r2$lasso.top1_magepro, stats_r2$lasso.top1_meta, paired = T, alternative = "greater")$p.value # p = 0

#populate dataframe for plotting
g <- grep("lasso",colnames(stats_r2))
average_r2 <- sapply(g, function(x) mean(stats_r2[, x], na.rm = T))
sem_r2 <- sapply(g, function(x) sd(stats_r2[,x], na.rm = T)/sqrt(length(which(!is.na(stats_r2[,x])))))
models=c("AFRONLY", "META","MAGEPRO")
df <- data.frame(average_r2,sem_r2, models)
write.table(df, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/magepro_AFR_r2.txt", quote = F, row.names = F, col.names = T)

#check percent increase
(df$average_r2[3] - df$average_r2[1])/df$average_r2[1] #184 %
(df$average_r2[3] - df$average_r2[2])/df$average_r2[2] #466 %

#plot results
cc <- c("#F5E453","#808080","#0666A1")
ggplot(df, aes(x= factor(models, levels =c("META","AFRONLY","MAGEPRO")), y=average_r2, fill=models)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=average_r2-sem_r2, ymax=average_r2+sem_r2), width=.2, position=position_dodge(.9)) +
  labs(x="Model", y = "Gene Expression Prediction R2 - LCL") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=cc, breaks=c("META", "AFRONLY", "MAGEPRO")) + 
  geom_hline(yintercept=0, linetype="solid", color="black") + 
  scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  theme_bw() +
  theme_classic(base_size = 15) + 
  geom_hline(aes(yintercept = mean_h2), linetype = "dotted", color = "red", linewidth = 1)
  

#save plot
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/MAGEPRO_AFR_wholeblood_LCL.png", width = 8, height = 8)

#--------------------------------------------------------------------------------------------------------------------

#--- how many genes have a raw r2 increase of > 0.05? 
nrow(stats) #10745
delta <- stats$lasso.top1_magepro - stats$lasso.top1_afr
length(which(delta > 0.05)) #4956
#--------------------------------------------------

#---num genes increase
length(which(stats$lasso.top1_magepro > stats$lasso.top1_afr))
#---

#---if we subset to significantly heritable & present in at least 1 dataset (h2p < 0.05 and positive)
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_fullsumstats_r2_h2.txt', as.is = T, header = T)
w <- which((stats$hsq_afr > 0) & (stats$hsq_afr.pv < 0.05) & (stats$hsq_afr != 0.064251) & (stats$datasets != "pred.wgt"))
stats <- stats[w, ]
nrow(stats) #2263 genes in this comparison 
mean_h2 = mean(stats$hsq_afr) #0.4062262
stats_r2 <- stats[, 2:4] #extract columns of r2 values for afronly, meta, and magepro

g <- grep("lasso",colnames(stats_r2))
average_r2 <- sapply(g, function(x) mean(stats_r2[, x], na.rm = T))
sem_r2 <- sapply(g, function(x) sd(stats_r2[,x], na.rm = T)/sqrt(length(which(!is.na(stats_r2[,x])))))
models=c("AFRONLY", "META","MAGEPRO")
df <- data.frame(average_r2,sem_r2, models)
write.table(df, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/magepro_r2_h2p0.05_AFR.txt", quote = F, row.names = F, col.names = T)

#---

# benchmark against PRS-CSx -----------------------------------------------------------------------------------------
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_fullsumstats_r2_h2.txt', as.is = T, header = T)
w <- which((stats$datasets != "pred.wgt"))
stats <- stats[w, ]
nrow(stats)

PRSCSx <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/benchmark_results.txt', as.is = T, header = T)
colnames(PRSCSx)
nrow(PRSCSx)

stats = filter(stats, gene %in% PRSCSx$gene)
PRSCSx = filter(PRSCSx, gene %in% stats$gene)
 
#nrow(stats) #19566
#nrow(PRSCSx) #19566

all(stats$gene == PRSCSx$gene)
mean(stats$hsq_afr)

r2 <- c(mean(stats_filter$lasso.top1_meta, na.rm = T),mean(stats_filter$lasso.top1_afr, na.rm = T), mean(stats_filter$lasso.top1_magepro, na.rm = T), mean(PRSCSx_filter$PRS_CSx_r2, na.rm = T), mean(PRSCSx_filter$PT_r2, na.rm = T))
sem <- c(sd(stats_filter$lasso.top1_meta, na.rm = T)/sqrt(length(which(!is.na(stats_filter$lasso.top1_meta)))),sd(stats_filter$lasso.top1_afr, na.rm = T)/sqrt(length(which(!is.na(stats_filter$lasso.top1_afr)))), sd(stats_filter$lasso.top1_magepro, na.rm = T)/sqrt(length(which(!is.na(stats_filter$lasso.top1_magepro)))), sd(PRSCSx_filter$PRS_CSx_r2, na.rm = T)/sqrt(length(which(!is.na(PRSCSx_filter$PRS_CSx_r2)))), sd(PRSCSx_filter$PT_r2, na.rm = T)/sqrt(length(which(!is.na(PRSCSx_filter$PT_r2)))))
Model <- c("META", "AFRONLY", "MAGEPRO", "PRS-CSx", "P+T")
df <- data.frame(Model, r2, sem)

#check percent increase
(df$r2[3] - df$r2[2])/df$r2[2] # MAGEPRO vs AFRONLY %184
(df$r2[3] - df$r2[1])/df$r2[1] # MAGEPRO vs META %466
(df$r2[3] - df$r2[5])/df$r2[5] # MAGEPRO vs PT %291
(df$r2[3] - df$r2[4])/df$r2[4] # MAGEPRO vs PRS-CSx %159


#cc <- c("#00b894", "#17becf","#ff6f31", "#e74c3c","#d3d3d9")
#cc <- c("#E7E4CA", "#40C0AE","#E57F4F", "#DE5D4C", "#202F6A")
cc <- c(
        "#687cae",
        "#465c6d",
        "#b6a752",
        "#aab5c1",
        "#188085"
        )
ggplot(df, aes(x= factor(Model, levels =c("META", "P+T", "AFRONLY","PRS-CSx","MAGEPRO")), y=r2, fill=Model)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=r2-sem, ymax=r2+sem), width=.2, position=position_dodge(.9)) +
  labs(x="Model", y = "Gene Expression Prediction R2")  +
  scale_fill_manual(values=cc, breaks=c("META", "P+T", "AFRONLY","PRS-CSx","MAGEPRO")) + 
  geom_hline(yintercept=0, linetype="solid", color="black") + 
  scale_y_continuous(breaks = seq(0, 0.07, by = 0.01), limits = c(0, 0.07)) + 
  theme_bw() +
  theme_classic(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6)
        )
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/benchmark.png", width = 3, height = 4)


#--------------------------------------------------------------------------------------------------------------------

#r2 vs h2 line plot -------------------------------------------------------------------------------------------------
#IMPORTANT NOTE: I used AFR h2 
r2h2 <- stats #continue working with dataframe created above

r2h2$lasso.top1_afr[r2h2$lasso.top1_afr < 0] <- 0
r2h2$lasso.top1_meta[r2h2$lasso.top1_meta < 0] <- 0
r2h2$lasso.top1_magepro[r2h2$lasso.top1_magepro < 0] <- 0
r2h2$hsq_afr[r2h2$hsq_afr < 0] <- 0
r2h2$hsq_afr[r2h2$hsq_afr == 0.064251] <- 0
groups <- rep(c("META", "AFRonly", "MAGEPRO"), each=nrow(r2h2))
combined <- c(r2h2$lasso.top1_meta, r2h2$lasso.top1_afr, r2h2$lasso.top1_magepro)
h2combined <- c(r2h2$hsq_afr, r2h2$hsq_afr, r2h2$hsq_afr)
df <- data.frame(groups, combined, h2combined)

cc <- c("#F5E453","#808080","#0666A1")
df$groups <- factor(df$groups,levels=unique(as.character(df$groups)))
ggplot(df, aes(x= h2combined, y=combined, color = groups)) + 
  geom_point(alpha = .05) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method='lm', se = TRUE) + 
  labs(x="gcta Heritability h2", y = "Gene Expression Prediction R2") +
  theme_bw() +
  theme_classic(base_size = 15) + 
  scale_color_manual(values=cc)
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/MAGEPRO_line.png", width = 8, height = 8)


#------------------------------------------------------------------------

# --- BENCHMARK LINE R2 H2 PLOT

stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_fullsumstats_r2_h2.txt', as.is = T, header = T)
colnames(stats)
nrow(stats)
w <- which((stats$datasets != "pred.wgt"))
stats <- stats[w, ]
nrow(stats)

PRSCSx <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/benchmark_results.txt', as.is = T, header = T)
colnames(PRSCSx)
nrow(PRSCSx)


# match up genes 
stats = filter(stats, gene %in% PRSCSx$gene)
PRSCSx = filter(PRSCSx, gene %in% stats$gene)
all(stats$gene == PRSCSx$gene)
nrow(PRSCSx)
nrow(stats)

#merge two dataframes 
PRSCSx = PRSCSx[, -c(6:9)]
merged <- merge(stats, PRSCSx, by = 'gene')
nrow(merged)#9767
colnames(merged)

#split into r2 dataframe and h2 dataframe
r2 <- as.matrix(merged[,c(2,3,4,16,17)])
h2 <- merged$hsq_afr
r2[r2 < 0] <- 0 
h2[h2 < 0] <- 0 #if h2 is negative, set to 0 (for cleaner plot)
h2[h2 == 0.064251] <- 0

#populate dataframe
vafr <- r2[,1]
vmeta <- r2[,2]
vmagepro <- r2[,3]
vPRSCSx <- r2[,4]
vPT <- r2[,5]
groups <- rep(c("PT", "META", "AFRonly", "PRSCSx", "MAGEPRO"), each=length(vafr))
h2combined <- c(h2, h2, h2, h2, h2) #append 3 vectors of h2 values
combined <- c(vPT, vmeta, vafr, vPRSCSx, vmagepro) #append vectors 
df <- data.frame(groups, combined, h2combined)

cc <- c("#687cae",
        "#465c6d",
        "#b6a752",
        "#aab5c1",
        "#188085")
#cc <- c("#800080", "#008080","#808000", "#800000", "#000080")
df$groups <- factor(df$groups,levels=c("META", "PT", "AFRonly", "PRSCSx", "MAGEPRO"))
ggplot(df, aes(x= h2combined, y=combined, color = groups)) + 
  geom_point(alpha = .03, size = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method='lm', se = TRUE, size = 0.5) + #intercept?
  labs(x="GCTA h2", y = "Gene Expression Prediction R2") +
  theme_bw() + 
  scale_color_manual(values=cc) +
  theme_classic(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/MAGEPRO_line_benchmark.png", width = 3, height = 3)


#--------------------------------------------------------------------------------------------------------------------

# --- EUR gene models made by MAGEPRO 
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_results_eur.txt', as.is = T, header = T)
colnames(stats)
nrow(stats)
w <- which((stats$datasets != "pred.wgt"))
length(which((stats$hsq_eur.pv < 0.01)))
stats <- stats[w, ]
nrow(stats) #19567 genes - same 19567 genes as AFR
mean_h2 = mean(stats$hsq_eur) #0.04211786

#subset r2 for eur average
stats_r2 <- stats[, 2:4] #not filtered for genes in afr analysis
nrow(stats_r2) #19567 -> all genes in eur analysis with h2 > 0, at least 1 external dataset 

#populate dataframe for eur average
g <- grep("lasso",colnames(stats_r2))
average_r2 <- sapply(g, function(x) mean(stats_r2[, x], na.rm = T))
sem_r2 <- sapply(g, function(x) sd(stats_r2[,x], na.rm = T)/sqrt(length(which(!is.na(stats_r2[,x])))))
models=c("EURONLY", "META","MAGEPRO")
df <- data.frame(average_r2,sem_r2, models)
write.table(df, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/magepro_EUR_r2.txt", quote = F, row.names = F, col.names = T)


#check percent increase
(df$average_r2[3] - df$average_r2[1])/df$average_r2[1] #53%
(df$average_r2[3] - df$average_r2[2])/df$average_r2[2] #74%

#plot results
cc <- c("#F5E453","#808080","#0666A1")
ggplot(df, aes(x= factor(models, levels =c("META","EURONLY","MAGEPRO")), y=average_r2, fill=models)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=average_r2-sem_r2, ymax=average_r2+sem_r2), width=.2, position=position_dodge(.9)) +
  labs(x="Model", y = "Gene Expression Prediction R2") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=cc, breaks=c("META", "EURONLY", "MAGEPRO")) + 
  geom_hline(yintercept=0, linetype="solid", color="black") + 
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01), limits = c(0, 0.05)) + 
  theme_bw() +
  theme_classic(base_size = 15) + 
  geom_hline(aes(yintercept = mean_h2), linetype = "dotted", color = "red", linewidth = 1)

#save plot
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/MAGEPRO_EUR.png", width = 8, height = 8)


# --- 

#--- subset to significantly heritable in EUR
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_results_eur.txt', as.is = T, header = T)
w <- which((stats$hsq_eur > 0) & (stats$hsq_eur.pv < 0.01) & (stats$hsq_eur != 0.008909) & (stats$datasets != "pred.wgt"))
stats <- stats[w, ]
nrow(stats) #5674 genes in this comparison 
mean_h2 = mean(stats$hsq_eur) #0.1259003
stats_r2 <- stats[, 2:4] #extract columns of r2 values for afronly, meta, and magepro

g <- grep("lasso",colnames(stats_r2))
average_r2 <- sapply(g, function(x) mean(stats_r2[, x], na.rm = T))
sem_r2 <- sapply(g, function(x) sd(stats_r2[,x], na.rm = T)/sqrt(length(which(!is.na(stats_r2[,x])))))
models=c("EURONLY", "META","MAGEPRO")
df <- data.frame(average_r2,sem_r2, models)
write.table(df, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/magepro_r2_h2p0.01_EUR.txt", quote = F, row.names = F, col.names = T)
#---


# --- applying EUR MAGEPRO gene models to AFR cohort 
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/EUR_cross_ancestry_AFR.txt', as.is = T, header = T)
nrow(stats)
colnames(stats)
#rows_with_na <- apply(stats, 1, function(row) any(is.na(row)))
#stats[rows_with_na, ] #some na rows - why? eur model could have created a sparse model with a few nonzero snps that do not exist in AFR or where there is no genotype variability 
stats <- na.omit(stats)

#get genes from afr analysis 
stats2 <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_fullsumstats_r2_h2.txt', as.is = T, header = T)
w <- which((stats2$datasets != "pred.wgt"))
stats2 <- stats2[w, ]

stats <- filter(stats, gene %in% stats2$gene)
stats2 <- filter(stats2, gene %in% stats$gene)
all(stats$gene == stats2$gene)
nrow(stats) #19257
nrow(stats2) #19257

mean_h2 = mean(stats2$hsq_afr)

stats2_r2 <- stats2[,2:4]
g <- grep("lasso",colnames(stats2_r2))
average_r2 <- sapply(g, function(x) mean(stats2_r2[, x], na.rm = T))
sem_r2 <- sapply(g, function(x) sd(stats2_r2[,x], na.rm = T)/sqrt(length(which(!is.na(stats2_r2[,x])))))
models=c("AFRONLY", "META","MAGEPRO")
df2 <- data.frame(average_r2,sem_r2, models)


stats_r2 <- stats[, 2:3]
g <- grep("r2",colnames(stats_r2))
average_r2 <- sapply(g, function(x) mean(stats_r2[,x], na.rm = T))
sem_r2 <- sapply(g, function(x) sd(stats_r2[,x], na.rm = T)/sqrt(length(which(!is.na(stats_r2[,x])))))
models=c("EURONLY", "EURMAGEPRO")
df <- data.frame(average_r2,sem_r2, models)


df_plot <- rbind(df2, df)

cc <- c(
  "#465c6d",
  "#83c4b5",
  "#d1625b",
  "#b6a752",
  "#188085"
)
ggplot(df_plot, aes(x= factor(models, levels =c("META", "EURMAGEPRO", "EURONLY","AFRONLY","MAGEPRO")), y=average_r2, fill=models)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=average_r2-sem_r2, ymax=average_r2+sem_r2), width=.2, position=position_dodge(.9)) +
  labs(x="Model", y = "Gene Expression Prediction R2")  +
  scale_fill_manual(values=cc, breaks=c("META", "EURMAGEPRO", "EURONLY","AFRONLY","MAGEPRO")) + 
  geom_hline(yintercept=0, linetype="solid", color="black") + 
  scale_y_continuous(breaks = seq(0, 0.07, by = 0.01), limits = c(0, 0.07)) + 
  theme_bw() +
  theme_classic(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6)
  ) + 
  geom_hline(aes(yintercept = mean_h2), linetype = "dotted", color = "red", linewidth = 1)

ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/Cross_Ancestry_EURtoAFR.png", width = 8, height = 8)

# ---


#delta r2 vs # external datasets ------------------------------------------------------------------------------------
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_fullsumstats_r2_h2.txt', as.is = T, header = T)
w <- which((stats$hsq_afr > 0) & (stats$hsq_afr != 0.064251) & (stats$datasets != "pred.wgt"))
stats <- stats[w, ]
#analyze all genes with atleast 1 external dataset
colnames(stats)
nrow(stats) #10745
numdatasets <- c(rep(0, nrow(stats)))
df <- data.frame(stats$gene, stats$lasso.top1_magepro - stats$lasso.top1_afr, stats$datasets, numdatasets)
colnames(df) <- c("gene", "delta", "datasets", "nums")

# "pred.wgt,pred.wgt.ota.nonzero,pred.wgt.ota.zero" -> 2
for (i in 1:nrow(df)){
  datas <- strsplit(df$datasets[i], split=",")[[1]]
  sets <- c()
  for (d in datas){
    name <- strsplit(d, split = "[.]")[[1]][3]
    if (!is.na(name)){
      sets <- append(sets, name)
    }
  }
  sets <- unique(sets)
  num <- length(sets)
  df[i,4] <- num
}

cor.test(df$delta, df$nums)$p.value #p < 4.033168e-51

x <- split(df, df$nums)
max(df$nums) #6 datasets

averages <- c()
semeans<- c()


#jump by 2
for (num in seq(from=1, to=6, by=2)){
  print(x[[num]][,4])
  print(x[[num+1]][,4])
  l <- x[[num]][, 2]
  l <- append(l, x[[(num+1)]][, 2])
  avg <- mean(l, na.rm = T)
  print(avg)
  averages <- append(averages, avg)
  semean <- sd(l, na.rm = T)/sqrt(length(which(!is.na(l))))
  print(semean)
  semeans <- append(semeans, semean)
}

numberExternalDataSet <- c('1-2', '3-4', '5-6')

df2 <- data.frame(numberExternalDataSet, averages, semeans)
colnames(df2) <- c("External_Datasets", "Average_delta_r2", "SEM")

#df2[1, 2] = 0
#df2[1, 3] = 0

#cc <- c("#d7f1fb", "#87b8e2", "#4b88d6")
cc <- c('#0d74af', '#064a96', '#002754')


#plot 
ggplot(df2, aes(x = factor(External_Datasets, levels =numberExternalDataSet), y=Average_delta_r2, fill=External_Datasets)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=Average_delta_r2 - SEM, ymax=Average_delta_r2 + SEM), width=.3, position=position_dodge(.9)) +
  labs(x="Number of External Datasets", y = "Average Δr2") +     
  scale_fill_manual(breaks=numberExternalDataSet,values=cc) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10)) +
  geom_hline(yintercept=0, linetype="solid", color="black") + 
  scale_y_continuous(breaks = seq(0, 0.065, by = 0.015), limits = c(0, 0.065)) + 
  theme_bw() +
  theme_classic(base_size = 16)

#save 
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/MAGEPRO_datasets.png", width = 8, height = 8)

#--------------------------------------------------------------------------------------------------------------------

#heritability across number of datasets-----------------------------------------------------------------------------

stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_fullsumstats_r2_h2.txt', as.is = T, header = T)
w <- which((stats$datasets != "pred.wgt"))
stats <- stats[w, ]
#analyze all genes with atleast 1 external dataset
colnames(stats)
nrow(stats) #17451
numdatasets <- c(rep(0, nrow(stats)))
df <- data.frame(stats$gene, stats$hsq_afr , stats$datasets, numdatasets)
colnames(df) <- c("gene", "h2", "datasets", "nums")

# "pred.wgt,pred.wgt.ota.nonzero,pred.wgt.ota.zero" -> 2
for (i in 1:nrow(df)){
  datas <- strsplit(df$datasets[i], split=",")[[1]]
  sets <- c()
  for (d in datas){
    name <- strsplit(d, split = "[.]")[[1]][3]
    if (!is.na(name)){
      sets <- append(sets, name)
    }
  }
  sets <- unique(sets)
  num <- length(sets)
  df[i,4] <- num
}

cor.test(df$h2, df$nums)$p.value #p < 

x <- split(df, df$nums)
max(df$nums) #10 datasets

averages <- c()
semeans<- c()

#jump by 2
for (num in seq(from=1, to=10, by=2)){
  print(x[[num]][,4])
  print(x[[num+1]][,4])
  l <- x[[num]][, 2]
  l <- append(l, x[[(num+1)]][, 2])
  avg <- mean(l, na.rm = T)
  print(avg)
  averages <- append(averages, avg)
  semean <- sd(l, na.rm = T)/sqrt(length(which(!is.na(l))))
  print(semean)
  semeans <- append(semeans, semean)
}

numberExternalDataSet <- c('1-2', '3-4', '5-6', '7-8', '9-10')

df2 <- data.frame(numberExternalDataSet, averages, semeans)
colnames(df2) <- c("External_Datasets", "Average_h2", "SEM")

#df2[1, 2] = 0
#df2[1, 3] = 0

cc <- c("#CCCCCC", "#999999", "#666666", "#333333", "#000000")

#plot 
ggplot(df2, aes(x = factor(External_Datasets, levels =numberExternalDataSet), y=Average_h2, fill=External_Datasets)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=Average_h2 - SEM, ymax=Average_h2 + SEM), width=.3, position=position_dodge(.9)) +
  labs(x="Number of External Datasets", y = "Average h2") +     
  scale_fill_manual(breaks=numberExternalDataSet,values=cc) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10)) +
  geom_hline(yintercept=0, linetype="solid", color="black") + 
  scale_y_continuous(breaks = seq(0, 0.015, by = 0.005), limits = c(0, 0.015)) + 
  theme_bw() +
  theme_classic(base_size = 16)

#save 
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/MAGEPRO_h2_datasets.png", width = 8, height = 8)

#-----------------------------------------------------------------------------------------------------------------------

#delta r2 vs quantiles of h2 ----------------------------------------------------------------------------------------

#reload and filter differently 
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/MAGEPRO_fullsumstats_r2_h2.txt', as.is = T, header = T)

#take genes present in atleast 1 external dataset, use AFR h2 values 
w <- which((stats$hsq_afr > 0) & (stats$hsq_afr != 0.064251) & (stats$datasets != "pred.wgt"))
stats <- stats[w, ]
nrow(stats) #10745 genes

#make dataframe with delta r2 and h2 afr
deltar2_h2 <- data.frame(stats$gene, (stats$lasso.top1_magepro - stats$lasso.top1_afr), stats$hsq_afr)
colnames(deltar2_h2) <- c("gene", "delta_r2", "hsq_afr")

cor.test(deltar2_h2$hsq_afr, deltar2_h2$delta_r2)$p.value # p < 5.294431e-154, r = 0.2509672

#quantile h2
q_h2 <- quantile(deltar2_h2$hsq_afr, probs = seq(0,1,0.2))

#populate plotting dataframe
df <- data.frame()
for (x in 1:(length(q_h2)-1)){
  w <- which(deltar2_h2$hsq_afr > q_h2[x] & deltar2_h2$hsq_afr <= q_h2[x+1])
  val <- deltar2_h2$delta_r2[w]
  dump <- cbind(mean(val, na.rm=T), x, sd(val, na.rm = T)/sqrt(length(which(!is.na(val)))))
  df <- rbind(df, dump)
}
quantile_val <- c("0 ~ 0.048", "0.048 ~ 0.11", "0.11 ~ 0.18", "0.18 ~ 0.29", "0.29 ~ 1")
colnames(df) <- c("avg", "quantile", "sem")
df$quantile <- quantile_val

cc <- c("#ffe49e", "#deb54b", "#b47d2d","#a25c24", "#633913")
quantiles <- factor(df$quantile)
ggplot(df, aes(x = quantile, y=avg, fill=quantile_val)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=avg - sem, ymax=avg + sem), width=.3, position=position_dodge(.9)) +
  labs(x="Quantiles of h2", y = "Average Δr2") +    
  scale_fill_manual(breaks=quantile_val,values=cc) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=0, linetype="solid", color="black") + 
  theme_bw() +
  theme_classic(base_size = 16) + 
  scale_y_continuous(breaks = seq(0, 0.08, by = 0.01), limits = c(0, 0.08))
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/MAGEPRO_deltar2_h2.png", width = 8, height = 8)

#--------------------------------------------------------------------------------------------------------------------

#delta r2 vs additional nonzero SNPs -------------------------------------------------------------------------------
nrow(stats) #17380
deltar2 <- stats$lasso.top1_flex - stats$lasso.top1_afr
additionalsnps <- stats$nonzeroFLEX - stats$nonzeroAFR
cor.test(deltar2, additionalsnps)$p.value # p < 

df <- data.frame(deltar2, additionalsnps)

#try number of additional nonzero SNPs as buckets 
q_num <- quantile(additionalsnps, probs = seq(0,1,0.25))
df2 <- data.frame()
for (x in 1:(length(q_num)-1)){
  w <- which(df$additionalsnps > q_num[x] & df$additionalsnps <= q_num[x+1])
  val <- df$deltar2[w]
  dump <- cbind(mean(val, na.rm=T), x, sd(val, na.rm = T)/sqrt(length(which(!is.na(val)))))
  df2 <- rbind(df2, dump)
}


quantile <- c('0-26', '27-64', '65-116', '117-601')
df2 <- cbind(df2, quantile)
colnames(df2) <- c("delta", "quantiles", "sem", "quantile")

cc <- c("#AD4FC4", "#9D3DB5", "#84259B", "#640F78")

ggplot(df2, aes(x = factor(quantile, levels =quantile), y=delta, fill=quantile)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin= pmax(delta - sem, 0), ymax=delta + sem), width=.3, position=position_dodge(.9)) +
  labs(x="Additional Nonzero SNPs", y = "Average Δr2") +     
  scale_fill_manual(breaks=quantile,values=cc) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=0, linetype="solid", color="black") + 
  scale_y_continuous(breaks = seq(0, 0.025, by = 0.005), limits = c(0, 0.025)) + 
  theme_bw() +
  theme_classic(base_size = 16)
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/MAGEPRO_nonzeroSNP.png", width = 8, height = 8)

#--------------------------------------------------------------------------------------------------------------------

#ratio between alphas for zero and nonzero vectors, averaged across all genes --------------------------------------
nrow(stats) #10745
colnames(stats)

# for every gene, take the average ratio (absolute value) between nonzero and zero vectors, add them to a list 
ratios <- c()
for (i in 1:nrow(stats)){
  ratios_pergene <- c()
  #hash map to store dataset and corresponding alphas 
  h_done <- list()
  h <- list()
  datasets <- strsplit(stats$datasets[i], split = ",")[[1]]
  alphas <- strsplit(stats$alphas[i], split = ",")[[1]]
  for (k in 1:length(datasets)){
    h[[datasets[k]]] <- abs(as.numeric(alphas[k]))
    h_done[[datasets[k]]] <- FALSE
  }
  
  #iterate through hashmap and compute ratios, if either the nonzero or zero vector is missing, do not add to average ratio
  for (key in names(h)){
    name <- strsplit(key, split="[.]")[[1]][3]
    if (! is.na(name)){
      zero_name <- paste0("pred.wgt.",name, ".zero")
      nonzero_name <- paste0("pred.wgt.",name, ".nonzero")
      if ( ((zero_name %in% names(h)) & (nonzero_name %in% names(h))) ){
        if (( (h_done[[nonzero_name]] == FALSE) &  (h_done[[zero_name]] == FALSE) )){
          r <- h[[nonzero_name]]/h[[zero_name]]
          ratios <- append(ratios, r)
          h_done[[nonzero_name]] <- TRUE
          h_done[[zero_name]] <- TRUE
        }
      }
    }
  }
}


# take the average of the list
mean(ratios) #1257.41


#--------------------------------------------------------------------------------------------------------------------

#number of genes and percent increase buckets ----------------------------------------------------------------------
nrow(stats)
colnames(stats)
stats_percents <- stats # copy to new dataframe
#zero out negative r2 to prevent sign changes
stats_percents[which(stats_percents$lasso.top1_afr < 0), 2] <- 0 
stats_percents[which(stats_percents$lasso.top1_magepro < 0), 4] <- 0 
percent_deltar2 <- (stats_percents$lasso.top1_magepro - stats_percents$lasso.top1_afr)/stats_percents$lasso.top1_afr
percent_deltar2 <- na.omit(percent_deltar2)
#percent_deltar2 <- percent_deltar2[is.finite(percent_deltar2)]


#populate dataframe
percent_buckets <- c(length(which(percent_deltar2 <= 0)), length(which(percent_deltar2 > 0)), length(which(percent_deltar2 > 0.5)), length(which(percent_deltar2 > 1)), length(which(percent_deltar2 > 2)), length(which(percent_deltar2 > 3)))
percent_increase <- c("<=0%", ">0%", ">50%", ">100%", ">200%", ">300%")
df <- data.frame(percent_increase, percent_buckets)

ggplot(df, aes(x = factor(percent_increase, levels = c("<=0%", ">0%", ">50%", ">100%", ">200%", ">300%")), y = percent_buckets, fill = percent_increase)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  labs(x = "Percent Increase in R2", y = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_y_continuous(limits = c(0, 9000), breaks = seq(0, 9000, by = 1000)) +
  theme_bw() +
  theme_classic(base_size = 15) +
  scale_fill_manual(values = c(
    "<=0%" = "#F7FCF5",
    ">0%" = "#C7E9C0",
    ">50%" = "#74C476",
    ">100%" = "#31A354",
    ">200%" = "#006D2C",
    ">300%" = "#00441B"
  ))

#save plot
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/MAGEPRO_percent_increase.png", width = 8, height = 8)
#--------------------------------------------------------------------------------------------------------------------

#afr vs nonafr datasets -------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------

#for top quantile of delta r2 genes, what is the average change in Beta vs Allele Frequencies -------------------------
#extract the high r2 genes 
delta <- stats$lasso.top1_flex - stats$lasso.top1_afr
max(delta)
q_delta <- quantile(delta, probs = seq(0,1,0.05))
w <- which(delta > q_delta[20] & delta <= q_delta[21])
df2 <- data.frame(stats[w, 1])
colnames(df2) <- c("gene")
write.table(df2, file = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/high_deltar2_genes.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)

#load the snps from those genes
betamaf <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/beta_maf_high.txt', as.is = T, header = T)
betamaf <- betamaf[- which((betamaf$mageprobeta - betamaf$afrbeta) == 0), ]
for (i in 1:nrow(betamaf)){
  if (betamaf$mafafr[i] < 0.01){
    betamaf[i,7] <- "rare"
  }
  if (betamaf$mafafr[i] > 0.01 & betamaf$mafafr[i] < 0.05){
    betamaf[i,7] <- "lowfreq"
  }
  if (betamaf$mafafr[i] > 0.05){
    betamaf[i,7] <- "common"
  }
}
for (i in 1:nrow(betamaf)){
  if (betamaf$mafeur[i] < 0.01){
    betamaf[i,8] <- "rare"
  }
  if (betamaf$mafeur[i] > 0.01 & betamaf$mafeur[i] < 0.05){
    betamaf[i,8] <- "lowfreq"
  }
  if (betamaf$mafeur[i] > 0.05){
    betamaf[i,8] <- "common"
  }
}

common_afr_common_eur <- which(betamaf$V7 == "common" & betamaf$V8 == "common")
common_afr_low_eur <- which(betamaf$V7 == "common" & betamaf$V8 == "lowfreq")
common_afr_rare_eur <- which(betamaf$V7 == "common" & betamaf$V8 == "rare")
low_afr_common_eur <- which(betamaf$V7 == "lowfreq" & betamaf$V8 == "common")
low_afr_low_eur <- which(betamaf$V7 == "lowfreq" & betamaf$V8 == "lowfreq")
low_afr_rare_eur <- which(betamaf$V7 == "lowfreq" & betamaf$V8 == "rare")
rare_afr_common_eur <- which(betamaf$V7 == "rare" & betamaf$V8 == "common")
rare_afr_low_eur <- which(betamaf$V7 == "rare" & betamaf$V8 == "lowfreq")
rare_afr_rare_eur <- which(betamaf$V7 == "rare" & betamaf$V8 == "rare")
AF <- c("common_afr_common_eur","common_afr_low_eur","common_afr_rare_eur", "low_afr_common_eur","low_afr_low_eur", "low_afr_rare_eur", "rare_afr_common_eur", "rare_afr_low_eur", "rare_afr_rare_eur")

avg <- c()
sem <- c()

for (a in AF){
  avg <- append(avg, mean(abs(betamaf$mageprobeta[eval(parse(text=a))] - betamaf$afrbeta[eval(parse(text=a))]), na.rm= T))
  sem <- append(sem, sd(abs(betamaf$mageprobeta[eval(parse(text=a))] - betamaf$afrbeta[eval(parse(text=a))]), na.rm = T)/sqrt(length(which(!is.na(abs(betamaf$mageprobeta[eval(parse(text=a))] - betamaf$afrbeta[eval(parse(text=a))]))))))
}

afr_maf <- rep(c("common_afr", "lowfreq_afr", "rare_afr"), each = 3)

df <- data.frame(AF, avg, sem, afr_maf)


ggplot(df, aes(x = AF, y=avg, fill=afr_maf)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=avg - sem, ymax=avg + sem), width=.3, position=position_dodge(.9)) +
  labs(x="Allele Frequency", y = "Average abs(ΔB)") +
  geom_hline(yintercept=0, linetype="solid", color="black") + 
  theme_bw() +
  theme_classic(base_size = 16) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_y_continuous(breaks = seq(0, 0.02, by = 0.005), limits = c(0, 0.02))
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/MAGEPRO_high_mafbeta.png", width = 14, height = 8)

#-----------------------------------------------------------------------------------------------------------------------

#rank datasets (run one at a time)--------------------------------------------------------------------------------------
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/ranked_datasets.txt', as.is = T, header = T)
stats[,2] <- stats$percent*100
#run each dataset individually, varying number of genes, filtered for h2 > 0 and 1 external dataset present
order <- stats$datasets
stats$datasets <- factor(stats$datasets, levels = order)
ggplot(stats, aes(x = datasets, y=percent, fill=datasets)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  labs(x="Datasets", y = "Percent increase in average r2") +
  geom_hline(yintercept=0, linetype="solid", color="black") + 
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits = c(0, 100)) + 
  theme_classic(base_size = 16) +
  scale_fill_manual(values = c("#004c6d", "#0f5e80", "#2a85a5", 
                                 "#46aeca", "#55c3dc", "#64d8ee")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/MAGEPRO_datasets_ranked_percentIncrease.png", width = 8, height = 8)

#-----------------------------------------------------------------------------------------------------------------------



#rank datasets average (aB)^2

#datasets added in rank order, r2 for train and testing set -------------------------------------------------------------

#normal r2 
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/normal_r2.txt', as.is = T, header = T)
stats[which(stats$r2_test == max(stats$r2_test)),] #eqtlgen,ota,his,test,ishi,mesa,genoa
number_datasets <- c()
for (i in 1:nrow(stats)){
  num <- length(strsplit(stats$datasets[i], split = ",")[[1]])
  number_datasets <- append(number_datasets, num)
}
stats <- data.frame(stats, number_datasets)
stats$number_datasets[7] = 0
colnames(stats)

ggplot(stats, aes(x = number_datasets)) +
  geom_line(aes(y = r2_test, color = "R-squared (Testing)"), linewidth = 0.5) +
  geom_line(aes(y = r2_train, color = "R-squared (Training)"), linewidth = 0.5, linetype = "dashed") +
  geom_errorbar(aes(ymin = r2_test - r2_test_sem, ymax = r2_test + r2_test_sem), width = 0.2) +
  geom_errorbar(aes(ymin = r2_train - r2_train_sem, ymax = r2_train + r2_train_sem), width = 0.2) +
  geom_point(aes(y = r2_test, color = "R-squared (Testing)"), size = 1) +
  geom_point(aes(y = r2_train, color = "R-squared (Training)"), size = 1, shape = 21, fill = "white") +
  scale_color_manual(
    values = c("R-squared (Testing)" = "blue", "R-squared (Training)" = "red"),
    guide = guide_legend(title = "R-squared", override.aes = list(linetype = c("solid", "dashed")))
  ) +
  labs(x = "Number of Datasets", y = "R-squared") +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  scale_y_continuous(breaks = seq(0, 0.31, 0.05), limits = c(0, 0.31))
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/r2_datasets_test_train.png", width = 8, height = 10)


ggplot(stats, aes(x = number_datasets)) +
  geom_line(aes(y = r2_test), color = "darkblue", linewidth = 0.5) +
  geom_errorbar(aes(ymin = r2_test - r2_test_sem, ymax = r2_test + r2_test_sem), width = 0.2) +
  geom_point(aes(y = r2_test), color = "darkblue", size = 1.5) +
  labs(x = "Number of Datasets", y = "R-squared") +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  scale_y_continuous(breaks = seq(0.02, 0.09, 0.01), limits = c(0.02, 0.09))
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/r2_datasets.png", width = 8, height = 10)



#adjusted r2
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/adj_r2.txt', as.is = T, header = T)
number_datasets <- c()
for (i in 1:nrow(stats)){
  num <- length(strsplit(stats$datasets[i], split = ",")[[1]])
  number_datasets <- append(number_datasets, num)
}
stats <- data.frame(stats, number_datasets)
stats$number_datasets[7] = 0
colnames(stats)

ggplot(stats, aes(x = number_datasets)) +
  geom_line(aes(y = r2_adjusted_test, color = "adj R-squared (Testing)"), linewidth = 0.5) +
  geom_line(aes(y = r2_adjusted_train , color = "adj R-squared (Training)"), linewidth = 0.5, linetype = "dashed") +
  geom_errorbar(aes(ymin = r2_adjusted_test - r2_adjusted_test_sem , ymax = r2_adjusted_test + r2_adjusted_test_sem ), width = 0.3) +
  geom_errorbar(aes(ymin = r2_adjusted_train - r2_adjusted_train_sem , ymax = r2_adjusted_train + r2_adjusted_train_sem ), width = 0.3) +
  geom_point(aes(y = r2_adjusted_test, color = "adj R-squared (Testing)"), size = 1) +
  geom_point(aes(y = r2_adjusted_train, color = "adj R-squared (Training)"), size = 1, shape = 21, fill = "white") +
  scale_color_manual(
    values = c("adj R-squared (Testing)" = "blue", "adj R-squared (Training)" = "red"),
    guide = guide_legend(title = "adj R-squared", override.aes = list(linetype = c("solid", "dashed")))
  ) +
  labs(x = "Number of Datasets", y = "adj R-squared") +
  theme_bw() +
  theme(text = element_text(size = 15))

ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/adjr2_datasets_test_train.png", width = 8, height = 10)

ggplot(stats, aes(x = number_datasets)) +
  geom_line(aes(y = r2_adjusted_test), color = "darkblue", linewidth = 0.5) +
  geom_errorbar(aes(ymin = r2_adjusted_test - r2_adjusted_test_sem, ymax = r2_adjusted_test + r2_adjusted_test_sem), width = 0.2) +
  geom_point(aes(y = r2_adjusted_test), color = "darkblue", size = 1.5) +
  labs(x = "Number of Datasets", y = "adj R-squared") +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  scale_y_continuous(breaks = seq(0.02, 0.05, 0.01), limits = c(0.02, 0.05))
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/adjr2_datasets.png", width = 8, height = 10)



#percent increase
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/percent_normal_r2.txt', as.is = T, header = T)
stats$percent <- stats$percent * 100
number_datasets <- c()
for (i in 1:nrow(stats)){
  num <- length(strsplit(stats$datasets[i], split = ",")[[1]])
  number_datasets <- append(number_datasets, num)
}
stats <- data.frame(stats, number_datasets)
stats$number_datasets[7] = 0
colnames(stats)
stats <- arrange(stats, datasets)

cc <- c("#FF7F00", "#4C78A8")

ggplot(stats, aes(x = number_datasets, y = percent)) +
  geom_line(color = cc[1]) +
  geom_point(color = cc[2]) +
  labs(x = "Number of Datasets", y = "Percent Increase") +
  theme_bw(base_size = 15) +
  scale_color_manual(values = cc) +
  scale_y_continuous(breaks = seq(0, 200, 50), limits = c(0, 200))
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/percent_increase_datasets.png", width = 9, height = 8)


# number of genes increase 

#percent increase
stats <- read.table(file='C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/num_genes_increase.txt', as.is = T, header = T)
number_datasets <- c()
for (i in 1:nrow(stats)){
  num <- length(strsplit(stats$datasets[i], split = ",")[[1]])
  number_datasets <- append(number_datasets, num)
}
stats <- data.frame(stats, number_datasets)
stats$number_datasets[7] = 0
colnames(stats)
stats <- arrange(stats, datasets)

cc <- c("#1BBCEA", "#277F99")

ggplot(stats, aes(x = number_datasets, y = num_genes_increase)) +
  geom_line(color = cc[1]) +
  geom_point(color = cc[2]) +
  labs(x = "Number of Datasets", y = "Number of Genes with Improved r2") +
  theme_bw(base_size = 15) +
  scale_color_manual(values = cc) +
  scale_y_continuous(breaks = seq(0, 8000, 1000), limits = c(0, 8000))
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/num_genes_improved_datasets.png", width = 9, height = 8)

#-----------------------------------------------------------------------------------------------------------------------


#rerun TWAS and analyze 

#r2 increase difference between tissues 

#simulations -----------------------------------------------------------------------------------------------------------

# --- same causal snp simulations in various settings 

simresults <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/simulation_results.txt", header = T)
simresults <- fread("C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/simulation_results_LD.txt", header = T)
nrow(simresults)
colnames(simresults)

df_plot_r2 <- simresults[, c(1, 2, 9, 10, 11, 12)]
colnames(df_plot_r2)

h2_set = c()
ss_set = c()
r2 = c()
sem = c()
model = c()

for (i in 1:nrow(df_plot_r2)){
  h2_set <- append(h2_set, df_plot_r2$h2_set[i])
  h2_set <- append(h2_set, df_plot_r2$h2_set[i])
  ss_set <- as.character(append(ss_set, df_plot_r2$ss_set[i]))
  ss_set <- as.character(append(ss_set, df_plot_r2$ss_set[i]))
  
  r2 <- append(r2, df_plot_r2$average_r2_afronly[i])
  sem <- append(sem, df_plot_r2$sem_r2_afronly[i])
  model <- append(model, "AFRONLY")
  r2 <- append(r2, df_plot_r2$average_r2_magepro[i])
  sem <- append(sem, df_plot_r2$sem_r2_magepro[i])
  model <- append(model, "MAGEPRO")
}

df_plot_r2 <- data.frame(h2_set, ss_set, r2, sem, model)

facetlables <- function(string) {
  string <- paste0("h2=", string)
  return (string)
}


df_plot_r2$ss_set <- factor(df_plot_r2$ss_set, levels = c("80", "160", "240", "400", "500"))
plot <- ggplot(df_plot_r2, aes(x = ss_set, y = r2, color = model, group = model)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = r2 - sem, ymax = r2 + sem),
                width = 0.2) +
  theme_bw() +
  labs(x = "AFR Sample Sizes", y = "Average R2", color = "Method") +
  scale_color_manual(values = c("AFRONLY" = "#b6a752", "MAGEPRO" = "#188085")) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), strip.text = element_text(size = 12)) + 
  facet_wrap(~h2_set, nrow = 1, labeller = labeller(
    h2_set = facetlables
  )) +
  scale_y_continuous(limits = c(0, 0.21))+  # Set the global y-axis limits
  geom_hline(data = df_plot_r2, aes(yintercept = h2_set), linetype = "dotted", color = "red", linewidth = 0.5)
plot
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/simulation_r2_samplesizes_h2.png", width = 6, height = 5)
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/simulation_r2_samplesizes_h2_LD.png", width = 6, height = 5)


# --- try violin plots 

df_plot_r2 <- fread(file= "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/r2_violin.txt")
df_plot_r2 <- fread(file= "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/r2_violin_LD.txt")
colnames(df_plot_r2)

df_plot_r2$ss_set_violin <- factor(df_plot_r2$ss_set_violin, levels = c("80", "160", "240", "400", "500"))
plot <- ggplot(df_plot_r2, aes(x = ss_set_violin, y = r2_violin, fill = model)) +
  geom_violin(trim = FALSE, position = position_dodge(0.9)) + 
  facet_wrap(~h2_set_violin, nrow = 1, labeller = labeller(h2_set_violin = facetlables)) +
  geom_hline(data = df_plot_r2, aes(yintercept = h2_set_violin), linetype = "dotted", color = "red", linewidth = 0.5) +
  stat_summary(aes(x = ss_set_violin, shape = "Mean"), fun = mean, geom = "point", position = position_dodge(0.9), size = 1.5, color = "black") +
  stat_summary(aes(x = ss_set_violin, shape = "Median"), fun = median, geom = "point", position = position_dodge(0.9), size = 1.5, color = "white", alpha = 0.7) +
  theme_bw() +
  labs(x = "AFR Sample Sizes", y = "R2 Values", fill = "Method", shape = "Statistic") +
  scale_shape_manual(values = c("Mean" = 19, "Median" = 24), labels = c("Mean", "Median")) +
  scale_fill_manual(values = c("AFRONLY" = "#b6a752", "MAGEPRO" = "#188085")) + 
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), strip.text = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(shape = c(NA, NA))))
plot

ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/simulation_r2_samplesizes_h2_violin.png", width = 10, height = 8)
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/simulation_r2_samplesizes_h2_violin_LD.png", width = 10, height = 8)


# ---


#plot Beta^2 of causal

df_plot_B2 <- simresults[, c(1, 2, 15, 16, 17, 18)]
colnames(df_plot_B2)

h2_set = c()
ss_set = c()
B2 = c()
sem = c()
model = c()

for (i in 1:nrow(df_plot_B2)){
  h2_set <- append(h2_set, df_plot_B2$h2_set[i])
  h2_set <- append(h2_set, df_plot_B2$h2_set[i])
  ss_set <- as.character(append(ss_set, df_plot_B2$ss_set[i]))
  ss_set <- as.character(append(ss_set, df_plot_B2$ss_set[i]))
  
  B2 <- append(B2, df_plot_B2$average_B2_afronly[i])
  sem <- append(sem, df_plot_B2$sem_B2_afronly[i])
  model <- append(model, "AFRONLY")
  B2 <- append(B2, df_plot_B2$average_B2_magepro[i])
  sem <- append(sem, df_plot_B2$sem_B2_magepro[i])
  model <- append(model, "MAGEPRO")
}

df_plot_B2 <- data.frame(h2_set, ss_set, B2, sem, model)

facetlables <- function(string) {
  string <- paste0("h2=", string)
  return (string)
}


df_plot_B2$ss_set <- factor(df_plot_B2$ss_set, levels = c("80", "160", "240", "400", "500"))
plot <- ggplot(df_plot_B2, aes(x = ss_set, y = B2, color = model, group = model)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = B2 - sem, ymax = B2 + sem),
                width = 0.2) +
  theme_bw() +
  labs(x = "AFR Sample Sizes", y = "Average B^2 Causal", color = "Method") +
  scale_color_manual(values = c("AFRONLY" = "#b6a752", "MAGEPRO" = "#188085")) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), strip.text = element_text(size = 12)) + 
  facet_wrap(~h2_set, nrow = 1, labeller = labeller(
    h2_set = facetlables
  )) +
  scale_y_continuous(limits = c(0, 0.25))  # Set the global y-axis limits
plot
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/simulation_B2_samplesizes_h2.png", width = 6, height = 5)
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/simulation_B2_samplesizes_h2_LD.png", width = 6, height = 5)


# plot diff (actual - predicted)^2 effect size causal

df_plot_B2 <- simresults[, c(1, 2, 19, 20, 21, 22)]
colnames(df_plot_B2)

h2_set = c()
ss_set = c()
diff_causal = c()
sem = c()
model = c()

for (i in 1:nrow(df_plot_B2)){
  h2_set <- append(h2_set, df_plot_B2$h2_set[i])
  h2_set <- append(h2_set, df_plot_B2$h2_set[i])
  ss_set <- as.character(append(ss_set, df_plot_B2$ss_set[i]))
  ss_set <- as.character(append(ss_set, df_plot_B2$ss_set[i]))
  
  diff_causal <- append(diff_causal, df_plot_B2$avg_causal_effect_diff_afr[i])
  sem <- append(sem, df_plot_B2$sem_causal_effect_diff_afr[i])
  model <- append(model, "AFRONLY")
  diff_causal <- append(diff_causal, df_plot_B2$avg_causal_effect_diff_magepro[i])
  sem <- append(sem, df_plot_B2$sem_causal_effect_diff_magepro[i])
  model <- append(model, "MAGEPRO")
}

df_plot_B2 <- data.frame(h2_set, ss_set, diff_causal, sem, model)

facetlables <- function(string) {
  string <- paste0("h2=", string)
  return (string)
}


df_plot_B2$ss_set <- factor(df_plot_B2$ss_set, levels = c("80", "160", "240", "400", "500"))
plot <- ggplot(df_plot_B2, aes(x = ss_set, y = diff_causal, color = model, group = model)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = diff_causal - sem, ymax = diff_causal + sem),
                width = 0.2) +
  theme_bw() +
  labs(x = "AFR Sample Sizes", y = "Average (Actual-Pred)^2 Causal", color = "Method") +
  scale_color_manual(values = c("AFRONLY" = "#b6a752", "MAGEPRO" = "#188085")) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), strip.text = element_text(size = 12)) + 
  facet_wrap(~h2_set, nrow = 1, labeller = labeller(
    h2_set = facetlables
  )) +
  scale_y_continuous(limits = c(0, 0.25))  # Set the global y-axis limits
plot
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/causalSNP_effectsize_diff.png", width = 6, height = 5)
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/causalSNP_effectsize_diff_LD.png", width = 6, height = 5)

# plot h2 of all simulation settings 

df_plot_h2 <- simresults[, c(1, 2, 3, 4, 5, 6)]
colnames(df_plot_h2)

h2_set = c()
ss_set = c()
h2 = c()
sem = c()
model = c()

for (i in 1:nrow(df_plot_h2)){
  h2_set <- append(h2_set, df_plot_h2$h2_set[i])
  h2_set <- append(h2_set, df_plot_h2$h2_set[i])
  ss_set <- as.character(append(ss_set, df_plot_h2$ss_set[i]))
  ss_set <- as.character(append(ss_set, df_plot_h2$ss_set[i]))
  
  h2 <- append(h2, df_plot_h2$avg_gctah2_afr[i])
  sem <- append(sem, df_plot_h2$sem_gctah2_eur[i])
  model <- append(model, "AFR")
  h2 <- append(h2, df_plot_h2$avg_gctah2_eur[i])
  sem <- append(sem, df_plot_h2$sem_gctah2_eur[i])
  model <- append(model, "EUR")
}

df_plot_h2 <- data.frame(h2_set, ss_set, h2, sem, model)

facetlables <- function(string) {
  string <- paste0("h2=", string)
  return (string)
}


df_plot_h2$ss_set <- factor(df_plot_h2$ss_set, levels = c("80", "160", "240", "400", "500"))
plot <- ggplot(df_plot_h2, aes(x = ss_set, y = h2, color = model, group = model)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = h2 - sem, ymax = h2 + sem),
                width = 0.2) +
  theme_bw() +
  labs(x = "AFR Sample Sizes", y = "Average GCTA h2", color = "POP") +
  scale_color_manual(values = c("AFR" = "#b6a752", "EUR" = "#188085")) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), strip.text = element_text(size = 12)) + 
  facet_wrap(~h2_set, nrow = 1, labeller = labeller(
    h2_set = facetlables
  )) +
  scale_y_continuous(limits = c(0, 0.24))  +  # Set the global y-axis limits
  geom_hline(data = df_plot_h2, aes(yintercept = h2_set), linetype = "dotted", color = "red", linewidth = 0.5)
plot
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/simulation_h2_samplesizes_h2.png", width = 6, height = 5)
ggsave(filename = "C:/Users/kaiak/OneDrive/Desktop/AMARIUTALAB/final_plots/simulation_h2_samplesizes_h2_LD.png", width = 6, height = 5)



#-----------------------------------------------------------------------------------------------------------------------

