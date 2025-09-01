# Rscript

# Script to check the ASvs counts of Bifidobacterium spp. find on GTDB what to correspond on Silva database

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 27/08/2025

# R version: R 4.4.1 

# upload library
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)

# path 
setwd('/home/')
getwd()

# Import table GTDB (output previuos script '7.Change_name_ASV_ambiguos.R')
GTDB <- read.table('ASV_&_taxa_Subgenus_correctASV.txt', sep = '\t', header = T)
head(GTDB)
dim(GTDB)

# Selected specific columns
df <- select(GTDB, ASV,Genus, Species)
colnames(df) <- c('ASV_GTDB', 'Genus_GTDB', 'Species_GTDB') # change name of columns
head(df)

# check counts of Bifidobacterium spp
df1 <- subset(df, Genus_GTDB == 'Bifidobacterium')
dim(df1)

# vector with ASV of bifidobacterium
ASV_GTDB <- df1$ASV_GTDB

# Import table SilvaDB 
SilvaDB <- read.table('ASV_&_taxa_Subgenus_correctASV_SilvaDB.txt', sep = '\t', header = T)
head(SilvaDB)
dim(SilvaDB) 

# Selected specific columns
df2 <- select(SilvaDB, Row.names, Genus, Species)
colnames(df2) <- c('ASV_Silva', 'Genus_Silva', 'Species_Silva')
head(df2)

# take ASVs of Bifidobacterium in results about silva database
Tab <- subset(df2, ASV_Silva %in% ASV_GTDB)
head(Tab)
dim(Tab)

# Merge two tabele
H <- cbind(df1, Tab)
head(H)

# How many vaginal Bifidobacterium are there in GTDB? 
# create a vector with all the species that contain Bifidobacterium spp.
vett <- H$Species_GTDB[grepl("Bifidobacterium vaginale", H$Species_GTDB)]
length(vett)

Bif <- subset(H, Species_GTDB %in% vett) 
head(Bif)

# How many Gardnerella vaginals are there in SIlvaDB?
Gard <- subset(H, Species_Silva == 'Gardnerella vaginalis')
head(Gard)
dim(Gard)  

## Avere un idea di quate specie sono presenti nei due database 
# GTDB --> TOT ASV per Bifidobacterium (genus) = 46 ASV
G1 <- subset(H, Species_GTDB == 'Bifidobacterium Unclassified')
dim(G1)

G2 <- subset(H, Species_GTDB == 'Bifidobacterium bifidum')
dim(G2) 

G3 <- subset(H, Species_GTDB == 'Bifidobacterium infantis')
dim(G3)# 2


# Plot 
head(H)

# Summary data for GTDB
gtdb_counts <- H %>%
  count(Species_GTDB) %>%
  rename(Taxa = Species_GTDB, Count = n) %>%
  mutate(Database = "GTDB")

# Summary data for SilvaDB
silva_counts <- H %>%
  count(Species_Silva) %>%
  rename(Taxa = Species_Silva, Count = n) %>%
  mutate(Database = "SILVA")

# Merge data
combined_counts <- bind_rows(gtdb_counts, silva_counts)
head(combined_counts)
dim(combined_counts)

# Specific palette
taxa_colors <- c(
  "Bifidobacterium Unclassified"     = "#1B9E77",  # Dark2
  "Bifidobacterium bifidum"          = "#D95F02",  # Dark2
  "Bifidobacterium breve"            = "#7570B3",  # Dark2
  "Bifidobacterium infantis"         = "#E7298A",  # Set1
  "Bifidobacterium kashiwanohense_A" = "#66A61E",  # Dark2
  "Bifidobacterium vaginale"         = "#E6AB02",  # Dark2
  "Bifidobacterium vaginale_A"       = "#A6761D",  # Dark2
  "Bifidobacterium vaginale_B"       = "#A6CEE3",  # Set1
  "Bifidobacterium vaginale_F"       = "#1F78B4",  # Set1
  "Bifidobacterium vaginale_H"       = "#B2DF8A",  # Set1
  "Bifidobacterium dentium"          = "#33A02C",  # Set1
  "Bifidobacterium kashiwanohense"   = "#FB9A99",  # Set1
  "Bifidobacterium longum"           = "#FDBF6F",  # Set1
  "Gardnerella vaginalis"            = "#CAB2D6",  # Set1
  "Unclassified"                     = "#6A3D9A"   # Set1
)

# bar plot
P1 <- ggplot(combined_counts, aes(x = Taxa, y = Count, fill = Taxa)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
  facet_wrap(~ Database, scales = "free_x") +
  scale_fill_manual(values = taxa_colors) +
  ylim(0,50) +
  theme_minimal() +
  labs(title = "", x = "", y = "Counts of ASV") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = 'italic', family = 'serif') ,
        axis.title.y = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = 'none')

P1


# Heatmap 
# prepare table

counts <- H %>%
  count(Species_GTDB, Species_Silva) %>%
  rename(Freq = n)
head(counts)

# Heatmap plot
P2 <- ggplot(conteggi, aes(x = Species_GTDB, y = Species_Silva, fill = Freq)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Freq), color = "black", size = 3) +
  scale_fill_gradient(low = "aquamarine3", high = "darkorange") +
  labs(title = "", x = "GTDB", y = "Silva") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = 'italic', family = 'serif'),panel.grid = element_blank(),
        axis.text.y = element_text(size = 10, face = 'italic', family = 'serif'),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = 'none')
P2

ggsave("Comparison_Bifidobacterium_spp_Garnerella_spp.png", P, width = 12, height = 7, dpi = 300)


# merge plots
combined_plot <- ggarrange (P1,P2, nrow = 2, labels = c('A', 'B'))
combined_plot

ggsave("panel_comparison_Bifidobacterium_spp_Garnerella_spp.png", combined_plot, width = 15, height = 10, dpi = 300, bg = 'white')