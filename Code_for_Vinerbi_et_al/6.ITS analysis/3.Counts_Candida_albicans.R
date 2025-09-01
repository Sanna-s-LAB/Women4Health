# Rscript

# Script to check counts of Candida albicans of ITS region 
# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 27/08/2025

# R version: R 4.4.1 

# upload library
library(tidyverse)

# path
setwd('/home/')
getwd()

# Import table
# Metadata file with information about plate, code name, input read and survived reads after QC procces
df <- read.table('Metadata_ITS.csv', sep = ';', header = T)
head(df)

# Import table with ASV and Taxonomic classification of ITS (output of '1.Dada2_ITS.R')
Taxa <- read.table('ASV_&_taxa_ITS.txt', sep = '\t', header = T)
head(Taxa)

# Select only Candida albicans species
Can_alb <- subset(Taxa, Species == 'Candida_albicans')
head(Can_alb)

# Select specific columns
Can_alb <- select(Can_alb, Species, 9:173)
dim(Can_alb)

# Aggregate value (counts)
selected_columns <- names(Can_alb)[2:ncol(Can_alb)]
G <- aggregate(. ~ Species, data = Can_alb[, c("Species", selected_columns)], FUN = sum) # CHANGE NAME
head(G)
rownames(G) <- G$Species

# Traspone table
H <- t(G)
H <- data.frame(rownames(H), H)
H <- H[-1,]
colnames(H) <- c('Code', 'Counts_CA')
dim(H)
head(H)

# Check that code of metadata is the same of candida albicans table
H$Code == B$Code

# Merge tables
M <- merge(B,H, by = 'Code')
head(M)
dim(M)

# Select specific columns
M <- select(M, ID, Counts_CA)
head(M)

# Merge metadata table with candida albicans table
S <- merge(df, M, by = 'ID')
dim(S)
S$Counts_CA <- as.numeric(S$Counts_CA)
str(S)

head(S)
colnames(S)[4] <- 'Counts Candida_albicans'

# Take samples with susrvived reads more to 200 reads
W <- subset(S, Survived.reads >= 200)
dim(W)
head(W)

# Check how many women have candida albicans 
head(H)

# Obtain a women code 
H$Women_code <- sub("_.*", "", H$Code)
head(H)
A <- unique(H$Women_code)
length(A) 

H$Counts_CA <- as.numeric(H$Counts_CA)
str(H$Counts_CA)

# Select specific columns
H1 <- select(H, Women_code, Counts_CA)
head(H1)

# aggregate information about women
selected_columns <- names(H1)[2:ncol(H1)]
H2 <- aggregate(. ~ Women_code, data = H1[, c("Women_code", selected_columns)], FUN = sum) # CHANGE NAME
head(H2)
dim(H2)

# Take samples/women with more 200 reads
H3 <- subset(H1, Counts_CA >= 200)


A <- unique(H3$Women_code)
length(A)
A
# Calculate the mean
mean(H3$Counts_CA)

# Compare the women with Candida albicans with the same women with counts of Lactobacillus
# Import table with taxonomic classification and ASVs of bacteria (output by '7.Change_name_ASV_ambiguos.R')
I <- read.table('ASV_&_taxa.txt', sep = '\t', header = T)
head(I)
dim(I)

# Select all samples and genus column
I1 <- select(I, Genus, 10:221)
head(I1)

# Aggregate counts of genus
selected_columns <- names(I1)[2:ncol(I1)]
I2 <- aggregate(. ~ Genus, data = I1[, c("Genus", selected_columns)], FUN = sum) # CHANGE NAME
View(I2)
rownames(I2) <- I2$Genus

# Take only lactobacillus genus
I2 <- subset(I2, Genus == 'Lactobacillus')
dim(I2)
head(I2) [1:5]
rowMeans(I2[,-1])

# traspone table
I3 <- t(I2)
I3 <- data.frame(rownames(I3), I3)
I3 <- I3[-1,]
colnames(I3)[1] <- 'Code'
head(I3)

# Obtain a code of women
I3$Women_code <- sub("_.*", "", I3$Code)
head(I3)

# select specific columns
I3 <- select(I3, Women_code, Lactobacillus)
head(I3)
I3$Lactobacillus <- as.numeric(I3$Lactobacillus)
str(I3)

# aggregate column
selected_columns <- names(I3)[2:ncol(I3)]
I4 <- aggregate(. ~ Women_code, data = I3[, c("Women_code", selected_columns)], FUN = mean) # CHANGE NAME
head(I4)
str(I4)

# create a table with women with candida and counts more to 200 reads and same women with lactobacillus 
# Merge table
Tab <- merge(H3, I4, by = 'Women_code')
View(Tab) 
str(Tab)

# Mean of candida
meanCand <- mean(Tab$Counts_CA)
meanCand

# Mean of lactobacillus
meanLacto <- mean(Tab$Lactobacillus)
meanLacto 

# Save output
write.table(Tab,  "Counts_Candida_and_lactobacillus.txt" ,sep="\t", row.names = T)
write.table(Tab,  "Counts_Candida_and_lactobacillus.csv" ,sep="\t", row.names = T)