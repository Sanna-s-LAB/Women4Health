#!/bin/bash Rscript

# Script to make a plot of relative abundance following the same order of the sample of relative abundance of V3V4V7V9 regions 
# (see script '8b.Relative_abundance_women')
# Note: for this script a complete metadata file, with all the sample IDs even when no data were present was used.

# Author: Elena Vinerbi
# Last update: 03/02/2025

# R version: v4.4.1,

# Upload library
library(tidyverse)
library(ggplot2)

# path
setwd('/home/')
getwd()

# Import counts table (output previous script '1.Dada2_ITS.R')
df <- read.table('ASV_&_taxa_order_ITS.txt', sep = '\t', header = T)
#View(df)
head(df)[2:9]

# aggregate a column of species
selected_columns <- names(df)[9:ncol(df)]
Tab <- aggregate(. ~ Species, data = df[, c("Species", selected_columns)], FUN = sum)
colnames(Tab)[1] <- "Taxa"
head(Tab) [1:5]
dim(Tab)
rownames(Tab) <- Tab$Taxa
head(Tab)
View(Tab)

# transpose table
df1 <- t(Tab)
df1 <- df1[-1,]
df1 <- data.frame(rownames(df1), df1)
colnames(df1)[1] <- 'ID'
View(df1)
dim(df1)
str(df1)

df1[, 2:59] <- lapply(df1[, 2:59], as.numeric)
str(df1)

# calculate a counts for each samples
df1$Counts <- rowSums(df1[,-1])
View(df1)

# Take only samples with more 200 reads
M <- subset(df1, Counts >= 200)
dim(M)
View(M)

# remove 'counts' coloumn
M <- M[,-60]

# calculate a relative abundance
Tab <- t(M)
Tab <- Tab [-1,]
Tab <- data.frame(rownames(Tab), Tab)
colnames(Tab)[1] <- 'Taxa'
View(Tab)
str(Tab)
dim(Tab)

Tab[, 2:26] <- lapply(Tab[, 2:26], as.numeric)
str(Tab)

Ab_taxa <- colSums(Tab[, -1])  
head(Ab_taxa)

Rel_ab_taxa <- matrix (nrow=dim(Tab)[1], ncol=dim(Tab)[2])
for (i in 1:(ncol(Tab)-1) ) {  Rel_ab_taxa[, i+1] <- as.numeric(Tab[, i+1] / Ab_taxa[i])} 
head(colSums(Rel_ab_taxa[,-1]))

Rel_ab_taxa = Rel_ab_taxa[,-1]
Rel_ab_taxa = Rel_ab_taxa * 100
rownames(Rel_ab_taxa) <- Tab[,1]
colnames(Rel_ab_taxa) <- colnames(Tab)[-1]
colSums(Rel_ab_taxa)
head(Rel_ab_taxa)[,5]
#View(Rel_ab_taxa)

Rel_ab_taxa[is.na(Rel_ab_taxa)] <- 0


# Remove unclassified taxa 
rownames(Tab)

U <- Rel_ab_taxa[rownames(Rel_ab_taxa) == 'Unclassified.Unclassified', ]
head(U) [1:5]

G <- Rel_ab_taxa[rownames(Rel_ab_taxa) != 'Unclassified.Unclassified', ]
rownames(G)
dim(G)

# calculate an average to order the plot
mean_rel_ab <- rowMeans (G) # Use a mean to order taxa
mean_rel_ab

tot_ab.rel <- names(mean_rel_ab[order(mean_rel_ab, decreasing=TRUE)])
#tot_ab.rel

## Obtain first 15 most abundance taxa
top_ab.rel <- names(mean_rel_ab[order(mean_rel_ab, decreasing=TRUE)])[1:10]
top_ab.rel

## Check the nemes of samples with specific order
a <- match(tot_ab.rel, rownames(G))
a
Tot_ab.rel<- G[a,]
head(Tot_ab.rel)
dim(Tot_ab.rel)

## Create a columx "Others" 
df <- colSums(Tot_ab.rel [11:57, ]) # change if the number is not the same
df_t <- t(df)
rownames(df_t)<- "Others"
head(df_t)[1:5]

# Create a table with top taxa
b <- match(top_ab.rel , rownames(Rel_ab_taxa))
b
Top_ab.rel<- Rel_ab_taxa[b,]
head(Top_ab.rel) [1:5]

Top_tab <- rbind(Top_ab.rel, df_t)
Top_tab <- data.frame(rownames(Top_tab), Top_tab)
colnames(Top_tab)[1] <- "Taxa"
head(Top_tab)[1:5]
dim(Top_tab)

## Sum unclassified and Other --> Now in Others levels we also find the unclassified
# View(U)
# View(df_t)
Others <- rbind (U,df_t)
Others <- colSums(Others)
Others <- data.frame(Others)
Others_t <- t(Others)
Others <- data.frame(rownames(Others_t), Others_t)
colnames(Others)[1] <- "Taxa"
head(Others)

## Remove Others (without unclassified to Top table)
Top <- Top_tab [-11,] 

# Merge Top table with Others (with unclassified value)
Top <- rbind(Top, Others)
colSums(Top[,-1])
head(Top)[1:5]
dim(Top)

# Add Unclassified to totat taxa table
Tot <- rbind (Tot_ab.rel, U)
Tot <- data.frame(rownames(Tot), Tot)
colnames(Tot)[1] <-"Taxa"
Tot$Taxa[Tot$Taxa == "U"] <- "Unclassified Unclassified"
rownames(Tot) <- Tot$Taxa
head(Tot)[1:5]
colSums(Tot[,-1])
dim(Tot)

# make a plot with all samples and colored only samples with more 200 reads

# Metadata with all IDs
# example of metadata columns: Plate; ID; Code_women; Visits; ID_paper; Var1; Var2
Meta<- read.table (paste(FOLDER_1, "Metadata_complete.csv", sep =""), sep = ";", header= TRUE)
dim(Meta)
head(Meta)

# Import index create in previus script (see '8b.Relative_abundance_women')
Index <- read.table ("Bacteria_first_visit_INDEX.txt", sep = "\t", header= TRUE)
View(Index)

# Select specific columns
IDx <- select(Index, Index, Women)
View(IDx)
str(IDx)
dim(IDx)

colSums(Top[,-1])  
head(Top)[1:5]

# Create a table with total samples with all visits
Samples <- Meta$Code_complete
missing_cols <- setdiff(Samples, colnames(Top))
length(missing_cols)

A <- as.data.frame(matrix(0, nrow = nrow(Top), ncol = length(missing_cols)))
rownames(A) <- rownames(Top)
colnames(A) <- missing_cols
dim(A)

# Merge table
Top_df <- cbind(Top, A)
dim(Top_df)
View(Top_df)
dim(Top_df)

# Order samples
df_1 <- Top_df$Taxa
head(df_1)


df_2 <- select(Top_df, 2:245)
head(df_2)[1:4]

col_names <- names(df_2)
split_names <- strcapture("X(\\d+)_(\\d+)", col_names, proto = list(first = numeric(), second = numeric()))
order_col <- order(split_names$first, split_names$second)
df_order <- df_2[, order_col]
colnames(df_order)

Top_df <- cbind(df_1, df_order)
colnames(Top_df)[1] <- 'Taxa'
#View(Top_df)
head(Top_df)
dim(Top_df)
length(colnames(Top_df)[2:213])


## Create a specific palette
palette1 <- brewer.pal(8, "Set3")    
palette2<- brewer.pal(3, "Set2")  

combined_palette <- c(palette1, palette2)
combined_palette <- combined_palette[1:11]

## Melt of Top table
D = melt(Top_df)
head(D)
View(D)

D$Women <- sub("_.*", "", D$variable)
head(D)
View(D)

Top_tab_melted <- merge(D, IDx, by = 'Women')
head(Top_tab_melted)
View(Top_tab_melted)
str(Top_tab_melted)

## Reorder table
Top_tab_melted <- transform(Top_tab_melted, Taxa = reorder(Taxa, -value))
#Top_tab_melted$Taxa <- factor(Top_tab_melted$Taxa, levels = taxa)

P <- ggplot(Top_tab_melted, aes(x = reorder(variable, Index), y = value, fill = Taxa)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = function(x) x*100)+
  scale_fill_manual(values = combined_palette) + 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) + 
  labs (title = "") + # cHANGE NAME
  ylab("Relative abundance") +
  xlab("Samples")+
  theme(axis.text.x = element_blank(),       # Remove samples name       
        axis.title = element_text(size = 7, family = 'sans'),
        legend.text = element_text(size = 7, family = 'sans'))
P

# Save output
ggsave("Top10_fungi_All_samples_order.png", P, width = 12, height = 7)

