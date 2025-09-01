# Rscript

# Script to calculate an alpha diversiry for both Shannon and Simpson index 

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

#Upload library
library(readxl)
library(tidyverse)
library(vegan)
library(rstatix)
library(ggpubr)
library(psych)
library(xtable)
library(tidyverse)
library(patchwork)

# path 
setwd("/home/")
getwd()

# Import table (output previus script 'Change_name_ASV_ambiguos.R')
# table of counts (row = taxon; column = samples)
df <- read.table("ASV_&_taxa_Subgenus_correctASV.txt", sep="\t", header=T)
head(df)[2:9]
dim(df)

# Obtain name of Sample_ID
ID <- colnames(select(df, 10:221)) 
length(ID)
ID

# Aggregate taxa a species level
selected_columns <- names(df)[10:ncol(df)]
Tab <- aggregate(. ~ Species, data = df[, c("Species", selected_columns)], FUN = sum)
colnames(Tab)[1] <- "Taxa"
head(Tab) [1:5]
dim(Tab)

# Manipulation taxa
Taxa_t <- t(Tab)
tax <-Tab$Taxa
Taxa_t <- as.data.frame(t(Tab [,-1]))
colnames(Taxa_t)<- tax

# Shannon_index 
Shannon <-diversity(Taxa_t, index="shannon")
head(Shannon)[1:5]

# Simpson_index
Simpson <- diversity(Taxa_t, index="simpson")
head(Simpson)[1:5]

# Create the table with all information
tab_alpha <- cbind(as.matrix (ID), as.matrix(Richness), as.matrix(Shannon), as.matrix(Simpson))
colnames(tab_alpha)<- c("ID","Shannon_index(H)","Simpson_index (D)")
dim(tab_alpha)
head(tab_alpha)

# Save the output
write.table(tab_alpha, "Tab_alpha_Species.csv", sep="\t")
write.table(tab_alpha, "Tab_alpha_Species.txt", sep="\t")

