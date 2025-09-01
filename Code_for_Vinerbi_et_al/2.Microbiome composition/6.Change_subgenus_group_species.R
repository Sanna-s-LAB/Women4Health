# Rscript

# Script to reclassify species name of lactobacillus with Subgenus group identified (Script above '5.Reclassification_subgenus.R')
# Note: this script requires a metadata file to change name of samples (optional)

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(stringr)
library(gridExtra)

# Files path
setwd('/home/')
getwd()

# Import table (output of Script above '5.Reclassification_subgenus.R')
# table with ASVs (row), taxon (row) and samples (column) with counts for each taxon after reclassification
df <- read.table("ASV_&_taxa_reclassification.txt", sep="\t", header=T)
dim(df)
head(df)[2:10]

# Metadata file
# example of metadata columns: Plate; ID; Code_women; Visits; ID_paper; Var1; Var2
Metadata <- read.table ("Metadata.txt", sep = "\t", header= TRUE)
head(Metadata)
dim(Metadata)

# Change species name with subgenus group name of lactobacilllus
# 1.Detect the groups in dataframe
Gasseri_df <- subset(df, Genus_subgroups == "Lactobacillus gasseri group")
dim(Gasseri_df)

iners_df <- subset(df, Genus_subgroups == "Lactobacillus iners group")
dim(iners_df)

crispatus_df <- subset(df, Genus_subgroups == "Lactobacillus crispatus group")
dim(crispatus_df)

jensenii_df <- subset(df, Genus_subgroups == "Lactobacillus jensenii group")
dim(jensenii_df)

delbrueckii_df <- subset(df, Genus_subgroups == "Lactobacillus delbrueckii group")
dim(delbrueckii_df)

# 2.Remove Lactobacillus groups in dataframe
df <- subset(df, Genus_subgroups != "Lactobacillus gasseri group")
dim(df)

df <- subset(df, Genus_subgroups != "Lactobacillus iners group")
dim(df)

df <- subset(df, Genus_subgroups != "Lactobacillus crispatus group")
dim(df)

df <- subset(df, Genus_subgroups != "Lactobacillus jensenii group")
dim(df)

df <- subset(df, Genus_subgroups != "Lactobacillus delbrueckii group")
dim(df)

# 3.Change species name and merge new names a main dataframe
Gasseri_df$Species <- "Lactobacillus gasseri"
dim(Gasseri_df)

df <- rbind(df, Gasseri_df)
dim(df)

iners_df$Species <- "Lactobacillus iners"
dim(iners_df)

df <- rbind(df, iners_df)
dim(df)

crispatus_df$Species <- "Lactobacillus crispatus"
dim(crispatus_df)

df <- rbind(df, crispatus_df)
dim(df)

jensenii_df$Species <- "Lactobacillus jensenii"
dim(jensenii_df)

df <- rbind(df, jensenii_df)
dim(df)

delbrueckii_df$Species <- "Lactobacillus delbrueckii"
dim(delbrueckii_df)

df <- rbind(df, delbrueckii_df)
dim(df)
#View(df)

# Change name of samples 
colnames(df)[10:221]

colnames(df)[10:221] <- Metadati$Code_women
colnames(df)[10:221]

colnames(df)[10:221] == Metadati$Code_women

# Save output
write.table(df,  "ASV_&_taxa_subgenus.csv" ,sep="\t", row.names = T) # row = taxon and ASVs; column = taxonomic level and samples
write.table(df,  "ASV_&_taxa_subgenus.txt" ,sep="\t", row.names = T)



