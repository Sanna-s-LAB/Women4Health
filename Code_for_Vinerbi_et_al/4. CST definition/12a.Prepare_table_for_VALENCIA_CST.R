# Rscript

# Script to prepare a table for VALENCIA tool
# Note: The information and file necessary is available here https://github.com/ravel-lab/VALENCIA
# We need to:
# file ASV_table
# file taxa_table --> modify taxonomy (e.g: g_Lactobacillus) -> aggregate an genus level except: Lactobacillus spp., Gardnerella spp., 
# Prevotella spp., Atopobium spp., Sneathia spp, -->(e.g., Lactobacillus_crispatus) and I have removed the unclassified. 
# Input table --> col1 = sampleID; col2 = read-count; remaining columns = taxa; rows = samples. 

# For this script part we have recalled some ASVs of bifidobacterium as Gardnerella vaginalis. 
# We used previous data analysed with the Silva database. 

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: 4.4.1; Python version: v3.12.7
# Dependence numpy v1.26.4; pandas v2.2.2 

# Upload library
library(tidyverse)
library(ggplot2)
library(vegan)
library(reshape2)
library(ggpubr)
library(gridExtra)

# Path
setwd("/home/")
getwd()

# Silva's file path
Folder <- "/home/Data_SilvaDB/" 

# Import table (outout previous script: '7.change_name_ASV_ambiguos.R')
# table with counts
tab <- read.table("ASV_&_taxa_Subgenus_correctASV.txt", sep = "\t", header= T)
head(tab)[4:10]
dim(tab)

# Remove Genus_subgroups
tab$Genus_subgroups <- NULL 
head(tab)[4:10]
dim(tab)

# Remove unclassified termns
remove_genus <- function(x) {
  words <- strsplit(x, " ")[[1]]
  if (length(words) == 3) {
    return(words[3])
  } else {
    return(x)
  }
}

tab$Species <- sapply(tab$Species,remove_genus)
head(tab)[4:10]
View(tab)

# Remove genus termns
remove_genus <- function(x) {
  words <- strsplit(x, " ")[[1]]
  if (length(words) == 2) {
    return(words[2])
  } else {
    return(x)
  }
}

tab$Species <- sapply(tab$Species,remove_genus)

# Col "genus" remove a family name 
remove_family <- function(x) {
  words <- strsplit(x, " ")[[1]]
  if (length(words) == 2) {
    return(words[2])
  } else {
    return(x)
  }
}

tab$Genus <- sapply(tab$Genus,remove_family)
head(tab)[4:10]
# View(tab)

# Remove name after "_" both species and genera value
tab$Genus <- sub("_.*", "", tab$Genus)
#View(tab)

tab$Species <- sub("_.*", "", tab$Species)
#View(tab)
dim(tab) 


# Check wich ASVs "Bifidobacterium" are "Gardnerella" 
# Import table with prevous data classify with Silva DB
df_SilvaDB <- read.table(paste(Folder,"ASV_&_taxa_SilvaDB.txt", sep=""), sep="\t", header=T)
dim(df_SilvaDB)

name_filter <- "Gardnerella vaginalis"
A <- subset(df_SilvaDB, Species == name_filter)
dim(A)
head(A)[3:10]

ASV_gen <- A$Row.names 
#ASV_gen

gard <- subset(tab, ASV %in% ASV_gen)
head(gard)[3:10]
dim(gard) 

ASV_remove <- gard$ASV
length(ASV_remove)

gard$Genus <- "Gardnerella"
head(gard)[5:10]

gard$Species <- "vaginalis"
head(gard)[5:10]

# Remove these ASV from main table
tab <- subset(tab, !(ASV %in% ASV_remove))
head(tab)[3:10]
dim(tab) #1129

# merge main table with table of gardnerella vaginalis
tab_df <- rbind(tab, gard)
dim(tab_df) 

# Merge genus and specie only for the genera required by VALENCIA
gen <- c("Lactobacillus", "Gardnerella", "Prevotella", "Atopobium", "Sneathia")

tab_df$Species <- ifelse(
  tab_df$Genus %in% gen & tab_df$Species != "Unclassified", 
  paste(tab_df$Genus, tab_df$Species, sep = "_"),
  tab_df$Species
)

#View(tab_df)

taxa <- tab_df
View(taxa)

## Add labels for taxonomic rank (for each taxa) --> except for a species level
taxa$Kingdom <- ifelse(taxa$Kingdom != "Unclassified", paste0("k_", taxa$Kingdom), taxa$Kingdom)
taxa$Phylum <- ifelse(taxa$Phylum != "Unclassified",  paste0("p_",  taxa$Phylum), taxa$Phylum)
taxa$Class <-  ifelse(taxa$Class != "Unclassified",  paste0("c_",  taxa$Class), taxa$Class)
taxa$Order <- ifelse(taxa$Order != "Unclassified",  paste0("o_",  taxa$Order), taxa$Order)
taxa$Family <-ifelse (taxa$Family !="Unclassified",  paste0("f_",  taxa$Family), taxa$Family)
taxa$Genus <- ifelse(taxa$Genus  !="Unclassified", paste0("g_",  taxa$Genus), taxa$Genus)
head(taxa) [2:9]
dim(taxa)

## save output
write.table(taxa, "ASV_&_taxa_Valencia_table.csv", sep="\t", row.names = F)
write.table(taxa, "ASV_&_taxa_Valencia_table.txt", sep="\t", row.names = F)

# Take genus, species and samples colums
df <- select(taxa, Genus, Species, 9:220)
head(df)[1:5]
dim(df)
#View(df)

# Remove Unclassified value
df <- df %>% filter(Genus != 'Unclassified')
df <- df %>% filter(Species != "Unclassified")
dim(df) 
#View(df)
head(df) [1:5]

# Aggregate genus level
selected_columns <- names(df)[3:ncol(df)]
tab_gen <- aggregate(. ~ Genus, data = df[, c("Genus", selected_columns)], FUN = sum)
head(tab_gen)[1:5] 
dim(tab_gen)

# Aggregate species level
selected_columns <- names(df)[3:ncol(df)]
tab_sp <- aggregate(. ~ Species, data = df[, c("Species", selected_columns)], FUN = sum)
head(tab_sp)[1:5]
dim(tab_sp) 

# take a species with complete name
df_sp <- tab_sp[grepl("_", tab_sp$Species), ]
#View(df_sp)
dim(df_sp) 

# Remove genera with a complete species name
sp <- df_sp$Species
sp <- sub("_.*", "", df_sp$Species)
sp <- unique(sp)
sp <- paste0("g_", sp)
sp

tab_gen <- tab_gen[!tab_gen$Genus %in% sp, ]
View(tab_gen)

# Create a table of genus
rownames(tab_gen) <- tab_gen$Genus
tab_gen <- tab_gen [-1]
tab_gen <- t(tab_gen)
tab_gen <- data.frame(rownames(tab_gen), tab_gen)
colnames(tab_gen)[1] <- "sampleID"
dim(tab_gen)
head(tab_gen) [1:5]

# Create a table of species
rownames(df_sp) <- df_sp$Species
head(df_sp) [1:5]
df_sp <- df_sp[,-1]
df_sp <- t(df_sp)
#View(df_sp)
df_sp <- data.frame(rownames(df_sp), df_sp)
colnames(df_sp)[1] <- "sampleID"
dim(df_sp)
head(df_sp) [1:5]

df_sp$sampleID == tab_gen$sampleID

df_sp <- df_sp[,-1]
head(df_sp) [1:5]

# Union table
df_valencia <- cbind(tab_gen, df_sp)
View(df_valencia)


# Calculate a total counts for each samples
df_valencia$read_count <- rowSums(df_valencia [,-1])
#View(df_valencia)
dim(df_valencia)
View(df_valencia)

df_valencia <- df_valencia %>%
  select(sampleID, read_count, 2:124)
#View(df_valencia)
dim(df_valencia)

## salvo output
write.table(df_valencia, "Tab_for_valencia.txt", sep="\t", row.names = F)
write.table(df_valencia, "Tab_for_valencia.csv", sep=",", row.names = F)

#####################################################################################################################
# To running VALENCIA

#-ref, --reference : path to the reference centroids file (provided)
#-i, --input csv  --> our table
#-o, --output 

# pass to bash enviroment
# python3 Valencia.py -r CST_centroids_012920.csv -i Tab_for_valencia.csv -o Output_Valencia_CST

