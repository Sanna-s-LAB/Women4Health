# Rscript

# Script to prepare the output of metaphlan table

# Author: Elena Vinerbi (elenavinerbi@cnr.it) and Fabio Chillotti
# Last update: 03/02/2025

# R version: R 4.4.1

# upload library
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

FOLDER_1 ="/home/data16S"

# Import table
# WGS_data --> it is necessary to manipulate the table
df <- read.table("mpa-4.1.1_vJun23_CHOCOPhlAnSGB.tsv", sep = "\t", header=T)
View(df)

df <- df[-1,]
View(df)

# Divided value after "|"
split_values <- function(x) {
  parts <- strsplit(as.character(x), "\n|\\|")[[1]]
  return(parts)
}

split_columns <- lapply(df$clade_name, split_values)
split_columns

## Create a col with taxonomy ranks
df$Kingdom <- NA
df$Phylum <- NA
df$Class <- NA
df$Order <- NA
df$Family <- NA
df$Genus <- NA
df$Species <- NA
df$Type <- NA
View(df)

## Insert name of taxa
for(row in 1:121){
  df$Kingdom[row] <-  sub(".*__", "", ifelse(length(split_columns[[row]]) >=1, split_columns[[row]][1], NA))
  df$Phylum[row] <- sub(".*__", "", ifelse(length(split_columns[[row]]) >=2, split_columns[[row]][2], NA))
  df$Class[row] <- sub(".*__", "",ifelse(length(split_columns[[row]]) >=3, split_columns[[row]][3], NA))
  df$Order[row] <- sub(".*__", "", ifelse(length(split_columns[[row]]) >=4, split_columns[[row]][4], NA))
  df$Family[row] <- sub(".*__", "", ifelse(length(split_columns[[row]]) >=5, split_columns[[row]][5], NA))
  df$Genus[row] <- sub(".*__", "", ifelse(length(split_columns[[row]]) >=6, split_columns[[row]][6], NA))
  df$Species[row] <- sub(".*__", "", ifelse(length(split_columns[[row]]) >=7, split_columns[[row]][7], NA))
  df$Type[row] <- sub(".*__", "", ifelse(length(split_columns[[row]]) >=8, split_columns[[row]][8], NA))
  
  
}
View(df)

## Select specific columns
Tab <- df[,c(12:19,2:11)]
View(Tab)

## Change name of Samples 
colnames(Tab)[9:18] <- substring(colnames(Tab)[9:18], 12, 16)
colnames(Tab)[9:18] <- sapply(colnames(Tab)[9:18], function(stringa){paste("X", stringa, sep = "")})
View(Tab)

colnames(Tab)[9:18]

## Merge Species and Type name
Tab$Type <- ifelse(!is.na(Tab$Type),
                   paste(Tab$Species, Tab$Type, sep = "_"), NA)
View(Tab)

## Select tabele
#M <- select (Tab, Type, 9:18)
#M<- na.omit(M)
#View(M)

#colSums(M[,-1])

## Select from Genus to type
DF <- select(Tab, 6:18)
View(DF)

# Take only Genus name with NA in Species (the same for Species level and type level)
Genus <- matrix(ncol = 11)
Genus <- as.data.frame(Genus)
Species <- matrix(ncol = 11)
Species <- as.data.frame(Species)

k <- 1
x <- 1
for (i in  unique(na.omit(DF$Genus))){
  A <- subset(DF, Genus == i)
  for(j in 1:nrow(A)){
    if (is.na(A$Species[j])){
      ab <- select(A, Genus, 4:13)[j,]
      Genus[k,] <- as.vector(ab)
      k <- k + 1
    } else {
      ab_sp <- select(A, Species, 4:13)[j,]
      if (is.na(A$Type[j])){
        Species[x,] <- ab_sp
        x <- x +1 
      }
      
    }
    
  }
}
View(Genus) 
View(Species) 

# save output
write.table(Genus, "WGS_genus.csv", sep = "\t")
write.table(Genus, "WGS_genus.txt", sep = "\t")

write.table(Species, "WGS_species.csv", sep = "\t")
write.table(Species, "WGS_species.txt", sep = "\t")