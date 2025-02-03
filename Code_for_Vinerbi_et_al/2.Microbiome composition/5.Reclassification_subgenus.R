# Rscript

# Script to reclassify genus and species of Lactobacillus in agreement with Isala protocol

# Note: the information and files about Isala protocol are aviilable here : https://github.com/LebeerLab/Citizen-science-map-of-the-vaginal-microbiome/
# Majority of the code is available here: https://github.com/LebeerLab/Citizen-science-map-of-the-vaginal-microbiome/
# and for this script requires a metadata file to sort samples (necessary for the script to work)

# Author: Elena Vinerbi (elenavinerbi@cnr.it) and Fabio Chillotti 
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(tidyverse)
library(tidyamplicons)
library(dplyr)

# Path
setwd('/home/')
getwd()

fin_refdb_lactobacillaceae <- "DATABASE/Isala_DB/SSUrRNA_GTDB05-lactobacillaceae-all_DADA2.fna.gz"
fin_refdb_lactobacillus_subgenera <- "DATABASE/Isala_DB/SSUrRNA_GTDB05-lactobacillus-subgenera-all_DADA2.fna.gz"

## we need to
# sample -> metadata 
# taxa -> taxa_tab 
# abbundance -> ASV_tab

# Import ASV table (after remove rare ASVs,  output of previuos script '4.Merge_file_dada2.R')
abundances <- read.table('ASVtab_filt.txt', header=T)
dim(abundances)
#View(abundances)
rownames(abundances)

## Import taxa table (Output of Dada2)
taxa <- read.table('ASV_taxa.csv',sep = '\t',header=T)
dim(taxa)
colnames(taxa)
#class(taxa)

# Import metadata to sort the samples 
# # example of metadata columns: Plate; ID; Code_women; Visits; ID_paper; Var1; Var2 
meta <- read.table ("Metadata.txt", sep = "\t", header= TRUE)
head(meta)
meta <- select(meta, 1:4)
head(meta)
# change a column name and use a one column like a description (require below)
colnames(meta) <- c("Plate","ID","sample", "description") 

sample <- meta$sample
sample

# Order the samples
abundances <- abundances [sample, ]
rownames(abundances)
#View(abundances)

## Create a tydy table (Abundance, taxa and samples)
abundances_tidy <- 
  abundances %>% 
  as_tibble() %>%
  mutate(sample = rownames(abundances)) %>%
  gather(key = "taxon", value = "abundance", - sample)
class(abundances_tidy)
#View(abundances_tidy)

samples_tidy <- 
  meta %>% 
  select(sample, description) 
#View(samples_tidy)
#class(samples_tidy)

## Check that sample name of abundance is the same in samples file (order and name)
A <- abundances_tidy$sample
B <- samples_tidy$sample
all(A == B) # check 

taxa_tidy <-  
  taxa %>% ##Taxa_table
  as_tibble() %>%
  rename(species = Species) %>%
  mutate(taxon = rownames(taxa))
#View(taxa_tidy)
#class(taxa_tidy)

## create the tidyamplicons object
ta <- make_tidyamplicons(samples_tidy, taxa_tidy, abundances_tidy)
#class(ta)

#ta %>% write_tidyamplicons(dout_ta)
#ta %>% saveRDS(file = fout_ta)

## Change name of ta is cross file
cross <- ta

# merge family Leuconostocaceae into Lactobacillaceae
cross$taxa <-
  cross$taxa %>%
  {.$Family[.$Family == "Leuconostocaceae"] <- "Lactobacillaceae"; .}

# reclassify Lactobacillaceae ASVs to the new taxonomy
cross <- 
  cross %>% 
  mutate_taxa(genus_oldtaxonomy = Genus) %>% 
  classify_taxa(
    fin_refdb_lactobacillaceae, Family == "Lactobacillaceae", 
    ranks = c("Genus", "species"), sequence_var = "taxon"
  )

# reclassify Lactobacillus ASVs to custom-defined subgenera
cross <- 
  cross %>% 
  mutate_taxa(genus_nosubgroups = Genus) %>%
  classify_taxa(
    fin_refdb_lactobacillus_subgenera, Genus == "Lactobacillus", 
    ranks = c("genus"), sequence_var = "taxon"
  )

#class(cross)
cross

# save integrated dataset as rds file
#saveRDS(cross, file = paste0(dout, "cross_reclassified.rds"))
#cross_reclassified$taxa

df_reclassify <- as.data.frame(cross$taxa)
#View(df_reclassify)

## Create a new table of taxa 
# taxon = Row.names (for match with ASV table)
# From Kingdom to Family (reclassify)
# Genus_old_taxonomy (is the same of Genus and genus_no_subgroups
# species (reclassifay)

df_taxa_new <- select(df_reclassify, taxon, Kingdom, Phylum, Class, Order, Family, genus_oldtaxonomy, genus, species)
colnames(df_taxa_new) <- c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Genus_subgroups", "Species")
#View(df_taxa_new)

# Save output
# table without samples, contains only taxonomic level (column) and taxon (row) with new classification
write.table(df_taxa_new, "Taxa_newclass.txt", sep="\t", row.names = T) 
write.table(df_taxa_new, "Taxa_newclass.csv", sep="\t", row.names = T)

# Merge ASVtable (after filter) and taxa table with new name
df <- abundances
#head(df)
df <- t(df)
df <- data.frame(rownames(df), df)
colnames(df)[1] <- 'ASV'
#View(df)
dim(df)
head(df)[1:5]

taxa <- df_taxa_new
head(taxa)
dim(taxa)

data <- merge(taxa, df, by='ASV')
dim(data)
#View(data)

# change NA with Unclassifed
columns <- colnames(data[2:9])
columns

for (col in columns) {
  data[[col]] <- ifelse(is.na(data[[col]]), 'Unclassified', data[[col]])
}
View(data)

# Taxa belonging to the family 'Lactobacilacea' that are reclassified , be careful when joining genus names with species names otherwise the genus name is repeated 2 times
A <- c('Lactiplantibacillus', 'Pediococcus', 'Limosilactobacillus', 'Lacticaseibacillus', 'Lentilactobacillus', 'Lactobacillus')
A

# Merge genus name with species name
data$Species <- ifelse(data$Species != 'Unclassified' & !(data$Genus %in% A), paste(data$Genus, data$Species), data$Species)

# Merge unclassified from Family to species
data$Genus <- ifelse(data$Genus == 'Unclassified', paste(data$Family, data$Genus), data$Genus)
data$Species <- ifelse(data$Species == 'Unclassified', paste(data$Genus, data$Species), data$Species)
#View(data)

data$Species <- ifelse(data$Species == 'Unclassified Lactobacillus iners', ' Lactobacillaceae Unclassified unclassified', data$Species)
View(data)

# Save output --> table with collapse name 'unclassified' from family to species
write.table(data, "ASV_&_taxa_reclassification.txt", sep="\t", row.names = T) # row = taxon and ASV, colum = taxonomic levels and samples with 
write.table(data, "ASV_&_taxa_reclassification.csv", sep="\t", row.names = T) # unclassified name from family to species 

