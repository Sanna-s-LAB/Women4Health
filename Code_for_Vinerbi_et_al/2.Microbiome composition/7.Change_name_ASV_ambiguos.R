# Rscript

# Script to change name of ASVs with ambiguos name recognize after reclassification of lactobacillus genus (see script '5.Reclassification_subgenus')

# Note: before changing the name, the sequences of ASVs with blast were checked

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(tidyverse)

# files path
setwd("/home/")
getwd()

# Import table (Output of previus script 6.Change_subgenus_group_species.R)
# Table with ASVs (row), taxon (row) and samples (column) after change name of species
df <- read.table ("ASV_&_taxa_subgenus.txt", sep = "\t", header= TRUE)
dim(df)
head(df)[2:10]
View(df)

# Detect ASVs ambiguos, remove them from the main df and use a "unclassified" for species of this specific genus
Amylolactobacillus <- subset(df, Genus == "Amylolactobacillus")
dim(Amylolactobacillus)

df <- subset(df, Genus != "Amylolactobacillus")
dim(df)


Amylolactobacillus$Species <- " Unclassified"

Agrilactobacillus <- subset(df, Genus == "Agrilactobacillus")
dim(Agrilactobacillus)

df <- subset(df, Genus != "Agrilactobacillus")
dim(df)

Agrilactobacillus$Species <- "Unclassified"

Paralactobacillus <- subset(df, Genus == "Paralactobacillus")
dim(Paralactobacillus)

df <- subset(df, Genus != "Paralactobacillus")
dim(df)

Paralactobacillus$Species <- "Unclassified"

# Merge df
df <- rbind(df, Amylolactobacillus)
dim(df)

df <- rbind(df, Agrilactobacillus)
dim(df)

df <- rbind(df, Paralactobacillus)
dim(df)
View(df)

# If in the species is 'Unclassified' merge genus name
df$Species <- ifelse(df$Genus == 'Amylolactobacillus', paste(df$Genus, df$Species), df$Species)
df$Species <- ifelse(df$Genus == 'Agrilactobacillus', paste(df$Genus, df$Species), df$Species)
df$Species <- ifelse(df$Genus == 'Paralactobacillus', paste(df$Genus, df$Species), df$Species)

# specific case with ambiguos name
df$Species <- ifelse(df$Genus == 'Lactobacillaceae Unclassified','Lactobacillaceae Unclassified Unclassified',  df$Species)
head(df)
dim(df)

# Order samples
df_1 <- select(df, 1:9)
head(df_1)

df_2 <- select(df, 10:221)
head(df_2)[1:4]

col_names <- names(df_2)
split_names <- strcapture("X(\\d+)_(\\d+)", col_names, proto = list(first = numeric(), second = numeric()))
order_col <- order(split_names$first, split_names$second)
df_order <- df_2[, order_col]

Tab <- cbind(df_1, df_order)
dim(Tab)
colnames(Tab)[10:221]
View(Tab)

#Save output
write.table(Tab,  "ASV_&_taxa_Subgenus_correctASV.csv" ,sep="\t", row.names = T) # table with row = taxon; column = samples
write.table(Tab,  "ASV_&_taxa_Subgenus_correctASV.txt" ,sep="\t", row.names = T)


