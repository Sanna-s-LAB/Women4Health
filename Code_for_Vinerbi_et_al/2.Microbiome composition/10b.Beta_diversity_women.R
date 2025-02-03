# Rscript

# Script to calculate an averge and median of  beta diversity between same women and random women

# Note: This script is designed to analyze both pre-clr (bray curtis index) and post clr (euclidean distance) data
# just change the name of the method when requested.

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(tidyverse)
library(ggplot2)
library(vegan)
library(reshape2)
library(ggpubr)
library(gridExtra)

# path
setwd("/home/")
getwd()


# Import table (output of previuos script '8a.Relative abundance' or '7.Statistical_analyses_pipeline/4.Bacterial_residual.py')
# Table of data after CLR trasformation and adjustment with covariates of laboratory (euclidean distance)
# or table with relative abundance (bray curtis index)
df_rel <- read.table( "data.csv", sep = ",", header= T)
dim(df_rel)

# Metadata
# example of metadata columns: Plate; ID; Code_women; Visits; ID_paper; Var1; Var2
Metadati<- read.table("Metadata.csv", sep = ";", header= TRUE)
dim(Metadati)
head(Metadati)

# for data after CLR trasformation and adjustment with covariates of laboratory
df_rel <- select(df_rel, Code_women, 2:422) # select code of samples and other column
colnames(df_rel)[1] <- "ID"
dim(df_rel)
head(df_rel)[1:5]
rownames(df_rel) <- df_rel$ID
head(df_rel) [1:5]
dim(df_rel)


# Obtain column with code for each women (take code before '_')
Metadati$Women <- sub("_.*", "", Metadati$ID)
head(Metadati) [1:7]

A <- select(Metadati, Women, ID)
head(A)

df_women <- merge(A, df_rel, by = "ID")
head(df_women)[1:5]
dim(df_women)

# Obtain column with visit for each women (take code after '_')
df_women$Visit_number <- sapply(df_women$ID, function(string)substring (string, 6,6)) # obtain a value 1,2 ecc
dim(df_women)
head(df_women)[1:5]
dim(df_women)

# create a specific table
df_women <- select(df_women, 1:2, Visit_number, 3:423)
head(df_women)[1:6]
dim(df_women)

## == Same women  == ##

## Create a vector with unique code for women
G <- df_women$Women
G <- unique(G)
length(G)
G

## empty dataframe for same women results
results <- data.frame(
  Women = character(),
  Mean  = numeric(),
  Median = numeric(),
  stringsAsFactors = FALSE
)

## Calculate Beta diversity for each women inside herself across weeks and after calculate an average and median (for each women)
for (g in G) {
  
  A <- subset(df_women, Women == g)
  A 
  beta_div <- as.matrix(vegdist(A[,-c(1:3)], method = "euclidean")) # may need to vary the method parameter ('bray' for untrasform data)
  
  S <- melt(beta_div)
  S <- unique(S$value)
  S <- setdiff(S, 0.00) ## Remove 0 value (diagonal)
  
  mean_value <- mean(S)
  median_value <- median(S)
  
  #Save results
  results <- rbind(results, data.frame(Women = g, mean  = mean_value, Median = median_value))
}

colnames(results)[1:3] <- c("Women", "Average", "Median")
head(results) ## Obtain a single value (1 for average and 1 for median) for each woman

## save output
write.table(results,"Average_median_Beta_same_women.csv",sep = "\t")
write.table(results,"Average_median_Beta_same_women.txt",sep = "\t")

## Remove NA value
results <- na.omit(results)
dim(results)


## == Random women between visits == ##

# Create 4 vector for each visit groups
V1 <- subset(df_women, Visit_number == "1")
head(V1) [1:5]
dim(V1) #59

V2 <- subset(df_women, Visit_number == "2")
head(V2) [1:5]
dim(V2) #58

V3 <- subset(df_women, Visit_number == "3")
head(V3)[1:5]
dim(V3) #52

V4 <- subset(df_women, Visit_number == "4")
head(V4)[1:5]
dim(V4) #43

#Loop to beta diversity with ramndom women filtering for different visits (repeat 10 times)
boxplot_list_average <- list()
boxplot_list_median <- list()

for (i in 1:10) {
  results_rm<- data.frame(
    Women_rm = character(),
    Average = numeric(),
    Median = numeric(), # change name for data pre CLR-trasformation
    stringsAsFactors = FALSE)
  
  for (r in V1$ID) {
    set.seed(i + 10) ## For repeat the analysis
    selected_women_2 <- sample(subset(V2 , ID != r)$ID,1) ## Select in random way a group of 4 women
    selected_women_3 <- sample(subset(V3 , !(ID %in% c(r, selected_women_2)))$ID,1)
    selected_women_4 <- sample(subset(V4 , !(ID %in% c(r, selected_women_2, selected_women_3)))$ID,1)
    
    B <- subset(df_women, ID %in% c(r, selected_women_2, selected_women_3, selected_women_4))
    rownames(B) <- B$ID
    
    beta_div_rm <- as.matrix(vegdist(B[,-c(1:3)], method = "euclidean")) # change parameter for untrasform data (method = 'bray')
    
    Z <- melt(beta_div_rm)
    Z <- unique(Z$value)
    Z <- setdiff(Z, 0.00) 
    
    mean_value <- mean(Z)
    median_value <- median(Z)
    
    #Salvo i risultati  
    results_rm <- rbind(results_rm, data.frame(Women_rm = r, Average = mean_value, Median = median_value))
  }
  print(dim(results_rm))
  
  ## save results 
  write.table(results_rm, paste("FC/Beta_diversity_random_women_rm_", i, ".csv", sep = ""))
  write.table(results_rm, paste("FC/Beta_diversity_random_women_rm_", i, ".txt", sep = ""))
  
}

