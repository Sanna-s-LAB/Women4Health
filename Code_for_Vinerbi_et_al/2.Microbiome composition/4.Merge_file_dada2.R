# Rscript

# Script to remove rare ASVs and merge the outputs od Dada2 (ASV table and taxa table)
# Note: this script requires a metadata file to sort samples (optional)

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(tidyverse)

# Path of file
setwd("/home/")
getwd()

# Import table (output of Dada2, previous script '3.Dada2.R')
# ASV table (row = samples; col = ASVs)
ASV <- read.table ("ASV_table.txt", sep = "\t", header= TRUE)
dim(ASV)
#head(ASV)[1:5]
rownames(ASV)

# Taxa table (row = taxon; column = taxonomic levels)
Taxa <- read.table ('Table_taxa.txt', sep = "\t", header= TRUE)
dim(Taxa)
head(Taxa)[1:5]

# Metadata file
# example of metadata columns: Plate; ID; Code_women; Visits; ID_paper; Var1; Var2
Metadata <- read.table ("Metadata.txt", sep = "\t", header= TRUE)
head(Metadata)
dim(Metadata)

## Order samples of ASV_tab like samples of metadata 
sample <- Metadata$ID
sample
length(sample)

ASV <- ASV[sample, ] # order the samples
dim(ASV)
rownames(ASV)

rownames(ASV) == Metadati$ID # Check if the order is the same

# Remove ASVs rare -> keep only ASV have at least 5 reads at least in 2 samples
selected <- c()    
for ( i in 1:dim(ASV)[2]){
  L = length(which(ASV [, i] >=5))     
  if (L >= 2){selected = c(selected, i)}
}
#selected  

ASV_filt <- ASV [, selected]
dim(ASV_filt)

# merge ASV table with taxa 
ASV_taxa <- merge(Taxa, as.data.frame(t(ASV_filt)), by="row.names")
colnames(ASV_taxa)[1] <- 'ASV'
dim(ASV_taxa)
head(ASV_taxa) [1:5]

# Remove any Mitochondria and Chloroplast (with GTDB is not necessary this step)
#name_filter <- "Mitochondria"
#data_filter <- subset(ASV_taxa, Family == name_filter)
#dim(data_filter)
# ASV_taxa<- ASV_taxa[-which(ASV_taxa$Family == name_filter), ]
# dim(ASV_taxa)

#name_filter <- "Chloroplast"
#data_filter <- subset(ASV_taxa, Order == name_filter)
#dim(data_filter)
# ASV_taxa_filt <- ASV_taxa[-which(ASV_taxa$Order == name_filter), ]
# dim(ASV_taxa_filt)

# Change names of NA with Unclassified 
columns <- colnames(ASV_taxa[2:8])
columns

for (col in columns) {
  ASV_taxa[[col]] <- ifelse(is.na(ASV_taxa[[col]]), 'Unclassified', ASV_taxa[[col]])
}
#View(ASV_taxa)

# Merge name of species level with genus level
ASV_taxa$Species <- ifelse(ASV_taxa$Species != 'Unclassified', paste(ASV_taxa$Genus, ASV_taxa$Species), ASV_taxa$Species)
head(ASV_taxa) [2:5]

# Merge unclassified from family to species
Tab <- ASV_taxa

Tab$Genus <- ifelse(Tab$Genus == 'Unclassified', paste(Tab$Family, Tab$Genus), Tab$Genus)
Tab$Species <- ifelse(Tab$Species == 'Unclassified', paste(Tab$Genus, Tab$Species), Tab$Species)
dim(Tab)
head(Tab)[1:12] # row = taxon and ASVs; column = taxonomic levels and samples with Unclassified level from Family to Species

# Save output
write.table(ASV_filt, "ASVtab_filt.txt", sep="\t", row.names = T) # table with ASVs after remove rare ASVs
write.table(ASV_filt, "ASVtab_filt.csv", sep="\t", row.names = T)

write.table(Tab, "ASV_&_taxa.txt", sep="\t", row.names = T) # table with ASVs, taxa (row) and samples (columns)
write.table(Tab, "ASV_&_taxa.csv", sep="\t", row.names = T)




