# Rscript

# Script to compare relative abundance of taxa between 16S data and WGS data

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(scales)
library(stringr)
library(gridExtra)

# path
setwd("/home/")
getwd()

# Import table (relative abundance) WGS (output previous script '11b.Preparetable_of_WGS.R')
df <- read.table("WGS_species_20240715.txt", sep = "\t", header=T)
dim(df)
View(df)

#Rel_ab_taxa <- Genus
#View(Rel_ab_taxa)

# Check
# colSums(df[, -1])  


# Calculate means of relative abundance to order the taxa 
rownames(df) <- df$Species

mean <- rowMeans (df[,-1])
mean

tot <- names(mean [order(mean, decreasing=TRUE)])
tot

a <- match(tot, rownames(df))
a
Tot_df <- df[a,]
head(Tot_df)
dim(Tot_df)

#View(Tot_df)

Tot_df[, 2:11] <- lapply(Tot_df[, 2:11], as.numeric)

str(Tot_df)

## Obtain first 15 most abundance taxa
top <- names(mean[order(mean, decreasing=TRUE)]) [1:10]
top

b <- match(top, rownames(df))
b
Top <- df[b,]
head(Top)
dim(Top)
View(Top)

Top [, 2:11] <- lapply(Top[, 2:11], as.numeric)


str(Top)

## Create a columx "Others" 
M <- colSums(Tot_df [11:34, -1]) 
M_t <- t(M)
rownames(M_t)<- "Others"
head(M_t)
dim(M_t)
M_t <- data.frame(rownames(M_t), M_t)
colnames(M_t) [1] <- "Species"
dim(M_t)

## Create a table with top taxa

Top_WGS <- rbind (Top, M_t)
View(Top_WGS)

colSums(Top_WGS[,-1])

## Save output
write.table(Top_WGS, "top10_others_Species_WGS.csv", sep = "\t")
write.table(Top_WGS, "top10_others_Species_WGS.txt", sep = "\t")

# Import table of 16S (output previuos script '7.Change_name_ASV_ambiguos.R')
# Table with counts
df <- read.table("ASV_&_taxa_Subgenus_correctASV.txt", sep = "\t", header=T)
#View(Tab)
dim(df)
head(df)[2:10]

selected_columns <- names(df)[10:ncol(df)]
Tab <- aggregate(. ~ Species, data = df[, c("Species", selected_columns)], FUN = sum)
head(Tab) [1:5]
dim(Tab)

rownames(Tab) <- Tab$Species

Tab <- select(Tab, Code_women1, Code_women2, Code_women3, Code_women4, Code_women5, Code_women6, Code_women7, Code_women8, Code_women9, Code_women10)
head(Tab)
str(Tab)
dim(Tab)

#Taxa <- rownames(Top_WGS)
#View(Tab)
Tab[, 2:10] <- lapply(Tab[, 2:10], as.numeric)

# Take taxa with more 0 counts
Tab$Counts <- rowSums(Tab[,-1])
S <- subset(Tab, Counts > 0)[,-11]
View(S)
dim(S)

S <- data.frame(rownames(S),S)
colnames(S)[1] <- "Species"
head(S)

## Remove words after "_" at the end of secondo word
S$Species <- gsub("^((\\w+\\s+\\w+))_.*", "\\1", S$Species)

## Remove words after "_" at the end of third word
S$Species <- gsub("^((\\w+\\s+\\w+\\s+\\w+))_.*", "\\1", S$Species)
View(S)

selected_columns <- names(S)[2:ncol(S)]
M <- aggregate(. ~ Species, data = S[, c("Species", selected_columns)], FUN = sum)
head(M) [1:5]
dim(M)

M$Species <- gsub(" ", "_", M$Species)
head(M)

M$Species

U <- M[60,]
U

M <- M [-60,]
M$Species
rownames(M) <-M$Species
View(M)

## Calculate relative abbundance of total taxa
Ab_taxa <- colSums(M[, -1])  
head(Ab_taxa)

Rel_ab_taxa <- matrix (nrow=dim(M)[1], ncol=dim(M)[2])
for (i in 1:(ncol(M)-1) ) {  Rel_ab_taxa[, i+1] <- as.numeric(M[, i+1] / Ab_taxa[i])} 
head(colSums(Rel_ab_taxa[,-1]))

Rel_ab_taxa = Rel_ab_taxa[,-1]
Rel_ab_taxa = Rel_ab_taxa * 100
rownames(Rel_ab_taxa) <- M$Species
colnames(Rel_ab_taxa) <- colnames(M)[-1]
head(Rel_ab_taxa)[,5]
View(Rel_ab_taxa)

# Order taxa
mean <- rowMeans (Rel_ab_taxa[,-1])
mean

tot <- names(mean [order(mean, decreasing=TRUE)])
tot

a <- match(tot, rownames(M))
a
Tot_df <- Rel_ab_taxa[a,]
head(Tot_df)
Tot_df <- data.frame(rownames(Tot_df), Tot_df)
colnames(Tot_df)[1] <- "Species"
View(Tot_df)

dim(Tot_df)

#View(Tot_df)

#Tot_df[, 1:10] <- lapply(Tot_df[, 1:10], as.numeric)

str(Tot_df)

## Obtain first 15 most abundance taxa
top <- names(mean[order(mean, decreasing=TRUE)]) [1:10]
top

b <- match(top, rownames(M))
b
Top <- Rel_ab_taxa[b,]
head(Top)
dim(Top)
Top <- data.frame(rownames(Top), Top)
colnames(Top)[1] <- "Species"
View(Top)

Top [, 2:11] <- lapply(Top[, 2:11], as.numeric)
str(Top)

## Create a columx "Others" 
D <- colSums(Tot_df [11:64, -1]) 
D_t <- t(D)
rownames(D_t)<- "Others"
head(D_t)
dim(D_t)
D_t <- data.frame(rownames(D_t), D_t)
colnames(D_t) [1] <- "Species"
dim(D_t)
View(D_t)

Top_16S <- rbind (Top, D_t)
View(Top_df)
colSums(Top_16S[,-1])

write.table(Top_16S, "top10_others_Species_16S.csv", sep = "\t")
write.table(Top_16S, "top10_others_Species_16S.txt", sep = "\t")

## Make a plot
# obtai a species name
Taxa_WGS <- Top_WGS$Species
Taxa_WGS
length(Taxa_WGS)

taxa_16S <- Top_16S$Species
taxa_16S
length(taxa_16S)

# Take a common taxa
common_taxa <- intersect(Taxa_WGS, taxa_16S)
common_taxa
length(common_taxa)

# No common_taxa
no_common_taxa_16S <- setdiff(taxa_16S, Taxa_WGS)
no_common_taxa_16S 
length(no_common_taxa_16S)

no_common_taxa_WGS <- setdiff(Taxa_WGS, taxa_16S)
no_common_taxa_WGS 
length(no_common_taxa_WGS)

# Specific palette for common taxa
common_taxa <- intersect(Taxa_WGS, taxa_16S)
col_common_taxa <- brewer.pal(8, "Paired") 
names(col_common_taxa) <- common_taxa

# Specific palette for NO common taxa
no_common_taxa_WGS <- setdiff(Taxa_WGS, taxa_16S)
col_no_common_taxa_WGS <- c("#d99853", "#bcbd22", "#ff9896")
names(col_no_common_taxa_WGS) <- no_common_taxa_WGS

no_common_taxa_16S <- setdiff(taxa_16S, Taxa_WGS)
col_no_common_taxa_16s <- c("#778953", "#7e693b", "#62a980")
names(col_no_common_taxa_16s) <- no_common_taxa_16S

## Specific palette for specific taxa
col_no_common_taxa_WGS["Gardnerella_vaginalis"] <- "aquamarine"
col_no_common_taxa_16s["Bifidobacterium_vaginale"] <- "aquamarine"

# Combine palette for both WGS and 16S
palette_WGS <- c(col_common_taxa, col_no_common_taxa_WGS)
palette_16S <- c(col_common_taxa, col_no_common_taxa_16s)

# Melt WGS table
Top_tab_melted_WGS <- melt(Top_WGS, value = Top_WGS$Species)
Top_tab_melted_WGS <- transform(Top_tab_melted_WGS, Species = reorder(Species, -value))
colnames(Top_tab_melted_WGS)[1:2] <- c("Species", "Sample")

# plot WGS
P_WGS <- ggplot(Top_tab_melted_WGS, aes(x = Sample, y = value, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette_WGS) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18)) +
  labs(title = "WGS") +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
P_WGS

# Save output
ggsave("WGS_Top10_Species.png", P_WGS, width = 12, height = 7)

# Melt 16S table
Top_tab_melted_16S <- melt(Top_16S, value = Top_16S$Species)
Top_tab_melted_16S <- transform(Top_tab_melted_16S, Species = reorder(Species, -value))
colnames(Top_tab_melted_16S)[1:2] <- c("Species", "Sample")

# plot  16S
P_16S <- ggplot(Top_tab_melted_16S, aes(x = Sample, y = value, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette_16S) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18)) +
  labs(title = "16S") +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
P_16S

# Save output
ggsave("16S_Top10_Species.png", P_WGS, width = 12, height = 7)


# Create a panell
A <- grid.arrange(P_WGS, P_16S)
A

# Save output
ggsave("Plot_comparison_16S_WGS_Top10.png", A, width = 12, height = 7)
