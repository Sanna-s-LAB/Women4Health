# Rsript

# Script to Choose 10 samples (analyzed with 16S metod) with different microbial profile to repeate the analysis wit WGS method
# Note: to obtain samples with different compositions we took the samples with greater distance of beta diversity between them

# Author: Elena Vinerbi (elenavinerbi@cnr.it) and Fabio Chillotti
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(tidyverse)
library(vegan)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

# path
setwd("/home/")
getwd()

# Import table (output previous script '8a.Relative_abundance.R')
# Table of Relative abundance of Species (16S data); row = taxa, column = samples
df_sp <- read.table("Tot_Species.txt", sep = "\t", header= T)
View(df_sp)

# traspone the table
df_sp_t <- t(df_sp)
dim(df_sp_t)
View(df_sp_t)
 
# Calculate beta diversity
beta <- as.matrix(vegdist(df_sp_t,method = 'bray'))
View(beta)
beta <- as.data.frame(beta)

# Detect cluster (three)
hc_beta <- hclust(as.dist(beta))
plot(hc_beta)
labels_beta <- cutree(hc_beta, k = 3)

beta$hc <- factor(labels_beta)
View(beta)

df_three <- df_sp[1:3,]
View(df_three)
df_three <- t(df_three)
df_three <- as.data.frame(df_three)
view(df_three)

# Detect samples with max value of first trhee taxa most abundant
max_iners <- max(df_three$`Lactobacillus Lactobacillus iners`)
max_iners

df_iners <- subset(df_three, df_three$`Lactobacillus Lactobacillus iners` == max_iners)
df_iners

max_crispatus <- max(df_three$`Lactobacillus Lactobacillus crispatus`)
max_crispatus

df_crispatus <- subset(df_three, df_three$`Lactobacillus Lactobacillus crispatus` == max_crispatus)
df_crispatus

max_gasseri <- max(df_three$`Lactobacillus Lactobacillus gasseri`)
max_gasseri

df_gasseri <- subset(df_three, df_three$`Lactobacillus Lactobacillus gasseri` == max_gasseri)
df_gasseri 

# Merge table
sample_bacteria <- rbind(df_iners, df_crispatus, df_gasseri)
View(sample_bacteria)
class(sample_bacteria)

#beta <- as.data.frame(beta)
View(beta)

# Start with L.iners cluster
df_1 <- subset(beta, hc == 1)
df_1

start_1<- c ("Code_women")
start_1 

for (i in 1:4){ 
  vec <- list()
  a <- select(df_1,all_of(start_1))
  
  for (col in start_1){
    vec[[col]] <- a[,col]
  }
  
  a$min <- do.call(pmin, vec)
  max_value <- max(a[,"min"])
  max <- subset (a, min == max_value)
  
  start_1[i+1] <- rownames(max)
  
}

class(start_1)
start_1

# Second cluster L.crispatus
df_2 <- subset(beta, hc == 2)
df_2

start_2<- c ("Code_women")
start_2 

for (i in 1:4){ 
  vec <- list()
  a <- select(df_2,all_of(start_2))
  
  for (col in start_2){
    vec[[col]] <- a[,col]
  }
  
  a$min <- do.call(pmin, vec)
  max_value <- max(a[,"min"])
  max <- subset (a, min == max_value)
  
  start_2[i+1] <- rownames(max)
  
}

class(start_2)
start_2

# Third cluster of L.gasseri
df_3 <- subset(beta, hc == 3)
df_3

start_3 <- c (" Code_women")
start_3 

for (i in 1:4){ 
  vec <- list()
  a <- select(df_3,all_of(start_3))
  
  for (col in start_3){
    vec[[col]] <- a[,col]
  }
  
  a$min <- do.call(pmin, vec)
  max_value <- max(a[,"min"])
  max <- subset (a, min == max_value)
  
  start_3[i+1] <- rownames(max)
  
}

class(start_3)
start_3
length(start_3)

# Calculate PCoA value (for 5 PCoA axis)
pcoa_result <- cmdscale(as.dist(beta))
pcoa_df <- data.frame(PCoA1 = pcoa_result[, 1], PCoA2 = pcoa_result[, 2])
pcoa_df <- data.frame(rownames(pcoa_df), pcoa_df)
colnames(pcoa_df)[1] <- "Sample_ID"
pcoa_df$WGS <- ifelse(pcoa_df$Sample_ID %in% total, pcoa_df$Sample_ID, "")
View(pcoa_df)

### Plot PCoA
P1 <-ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2 , color = WGS )) +
  geom_point()+
  geom_text(aes(label = WGS), nudge_x = 0.02, nudge_y = 0.02)+
  labs(x = "PCoA1 (55.60%)", y = "PCoA2 (16.52%)", color = "Visit" ,title = "PCoA of Genera")
P1 

total <- union(start_1, start_2)
total <- union(total, start_3)
length (total)

# PLot relative abbundance 
View (df_sp)

sample <-  select(df_sp, all_of(total)) [1:3,]
View(sample)

B_1 <- rowMeans(sample[,1:5])
B_1

B_2 <- rowMeans(sample[,6:10])
B_2

B_3 <- rowMeans(sample[,11:15])
B_3

sample_top <-  select(df_sp, all_of(total)) [1:20,]
View(sample_top)

df_top <- colSums(df_sp[21:158, total ])
length(df_top)
df_top <- t(df_top)
View(df_top)

#sample_tot <- data.frame(rownames(sample_tot), sample_tot)
#colnames(sample_tot)[1] <- "Species"
#View(sample_tot)

sample_tot <- rbind(sample_top, df_top)
rownames(sample_tot)[21] <- "Others"
sample_tot <- data.frame(rownames(sample_tot), sample_tot)
colnames(sample_tot)[1] <- "Species"
View(sample_tot)

palette1 <- brewer.pal(8, "Paired")    
palette2 <- brewer.pal(8, "Accent")  
palette3 <- brewer.pal(8, "Set3")
combined_palette <- c(palette1, palette2, palette3)
combined_palette <- combined_palette[1:21]

## Plot for top taxa
Top_tab_melted = melt(sample_tot , value = sample_tot$Species)
Top_tab_melted <- transform(Top_tab_melted, Species = reorder(Species, -value))
colnames(Top_tab_melted)[1:2] = c("Species", "Sample")

P2 <- ggplot(Top_tab_melted ,aes(x =Sample, y = value, fill = Species))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = combined_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  ylab("Relative abundance")+
  xlab("Samples")+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 9),  
        legend.title = element_text(size = 9))
P2

ggsave("Samples_for_WGS.png", P2, width = 12, height = 7)

write.table(sample_tot, "Samples_WGS.csv", sep = "\t")
write.table(sample_tot, "samples_WGS.txt", sep = "\t")



  
