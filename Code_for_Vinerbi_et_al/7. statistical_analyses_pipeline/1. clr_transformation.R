# Rscript
# This script aims to make a function to compute the CLR
# Author: Fabio Chillotti (fabiochillotti@cnr.it)
# Last change: 03/02/2025
# R version: R v4.4.1

library(vegan)
library(ggplot2)
library(parallel)
library(MASS)


#get the data
rel_abb_genus <- read.csv('data/Relative_abundances_genus.csv',sep = '\t') 
rel_abb_genus <- as.data.frame(t(rel_abb_genus))

rel_abb_species <- read.csv('data/Relative_abundances_species.csv',sep = '\t')
rel_abb_species <- as.data.frame(t(rel_abb_species))

#function to compute geometric mean
geom_mean <- function(x) {
  return(exp(mean(log(x))))
}

#function to compute the CLR transformation
clr_transformation <- function(data, eps) {
  
  num_cores <- detectCores() - 1 
  
  # Function to process each row
  process_row <- function(row) { #manage the zero problems
    row <- as.numeric(row)
    row <- row + eps  #translate the abundances by adding a pseudocount
    row <- row / sum(row) #normalize the abundances to sum up 1
    g <- geom_mean(row) #compute the geometric mean of the new row
    res <- log(row/g) #take the logarithm og the row divided by g
    return(res)
  }
  
  res_data <- mclapply(1:nrow(data), function(i) process_row(data[i, ]), mc.cores = num_cores) #apply process_row to our abundances
  
  #make the dataframe and return it
  res_data <- do.call(rbind, res_data)
  res_data <- as.data.frame(res_data)
  colnames(res_data) <- colnames(data)
  rownames(res_data) <- rownames(data)
  
  return(res_data)
}

#compute the clr on our data
clr_genus <- clr_transformation(rel_abb_genus,1e-5)

clr_species <- clr_transformation(rel_abb_species,1e-5)

#save the results in a csv files
write.csv(clr_genus,'~/complete_pipeline/data/clr_genus.csv')
write.csv(clr_species,'~/complete_pipeline/data/clr_species.csv')