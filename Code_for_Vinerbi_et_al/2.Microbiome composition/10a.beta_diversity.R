# Rscript

# Script to calculate a beta diversity

# Note: This script is designed to analyze both pre-clr (bray curtis index) and post clr (euclidean distance) data
# just change the name of the method when requested.

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(vegan)
library(ggplot2)
library(tidyverse)
library(reshape2)

# path
setwd("/home/")
getwd()

# Import table (output previous script '8a.Relative_abundance.R' or '7.Statistical_analyses_pipeline/4.Bacterial_residual.py' )
# Relative abundance for bray curtis index or data post CLR and adjusted for tech variable
df_rel <- read.table( "data.csv", sep = ",", header= T)
dim(df_rel)

# Metadata file
# example of metadata columns: Plate; ID; Code_women; Visits; ID_paper; Var1; Var2
Metadata <- read.table ("Metadata.txt", sep = "\t", header = TRUE)
head(Metadata) 
dim(Metadata)

# for data after CLR trasformation and adjustment with covariates of laboratory
df <- select(df_rel, Code_women, 2:422) # select code of samples (for data output from CLR trasformation)
head(df) [1:5]
rownames(df) <- df$Code
head(df) [1:5]
dim(df)
#View(df)

# Create a column with a number of visits
get_visit <- function(stringa){
  substring(stringa,6,6)
}

df$visit <- sapply(rownames(df),get_visit)
head(df) [420:423]
dim(df)
str(df)

df <- df[,-1]
head(df)[1:5]
dim(df)
View(df)
str(df)


# Calculate Beta diversity 
# change parmeters of 'method' with 'euclidean' for data postCLR or 'bray' untrasform data
BC_matrix <- as.matrix(vegdist(df[,-422],method = 'euclidean')) # remove column with number of visits
dim(BC_matrix)

# Creta a melt table with information for each weeks
A <- melt(BC_matrix)
hist(A$value, 50)

temp1 <- grep("_1",A[,1])
temp2 <- grep("_1",A[temp1,2])
A_visit1 <- A[temp1,][temp2,]
A_visit1 <- subset(A_visit1, Var1 != Var2)
head(A_visit1)

temp1 <- grep("_4",A[,1])
temp2 <- grep("_4",A[temp1,2])
A_visit4 <- A[temp1,][temp2,]
A_visit4 <- subset(A_visit4, Var1 != Var2)
head(A_visit4)

temp1 <- grep("_2",A[,1])
temp2 <- grep("_2",A[temp1,2])
A_visit2 <- A[temp1,][temp2,]
A_visit2 <- subset(A_visit2, Var1 != Var2)
head(A_visit2)

temp1 <- grep("_3",A[,1])
temp2 <- grep("_3",A[temp1,2])
A_visit3 <- A[temp1,][temp2,]
A_visit3 <- subset(A_visit3, Var1 != Var2)
head(A_visit3)

summary(A_visit1$value)
summary(A_visit2$value)
summary(A_visit3$value)
summary(A_visit4$value)

# change name of column
label <- c(rep(1,length(A_visit1$value)), rep(2,length(A_visit2$value)), rep(3,length(A_visit3$value)),rep(4,length(A_visit4$value)))
matrix_to_plot <- data.frame(label, c(A_visit1$value,A_visit2$value,A_visit3$value,A_visit4$value))
names(matrix_to_plot)[2]="beta"
tab <- matrix_to_plot
colnames(tab)[1] <- 'Weeks'
head(tab)

# Violin plot of distribution of beta diversity
# specific palette
col <- c( "1" = "#FFE066","2" = "#FFF5CC", "3" = "#FFC54C", "4" = "#FFA32B") 
col

P_violin <- tab %>% 
  ggplot() +
  theme_bw() +
  geom_violin(aes(x = factor(Weeks), y = beta, fill=factor(Weeks))) + 
  geom_boxplot(aes(x = factor(Weeks), y = beta), width=0.1, alpha=0.6) +
  scale_x_discrete(labels = c("F", "O", "EL", "LL")) +
  labs(title = "", x = "Weeks", y = "Beta diversity") +
  scale_fill_manual(name = "Weeks", values = col, labels = c("F", "O", "EL", "LL"))+
  theme(
    axis.title.x = element_text(size = 7),                       
    axis.title.y = element_text(size = 7),                      
    axis.text.x = element_text(size = 7),                        
    axis.text.y = element_text(size = 7),
    legend.position = 'none')

P_violin

# Save output
ggsave("Violin_beta_PostCLR.png", P_violin, width = 15, height = 7, dpi = 300)


# Calculate Median and average of Beta diversity on each visits
dim(df)
head(df) 

Median_Average_no_sample  <- function(vis) {
  df_i <- subset(df, visit == vis)
  BC_matrix_i <- as.matrix(vegdist(df_i[,-422
  ],method = 'euclidean')) # change method
  num_sample <- ncol(BC_matrix_i)
  results <- matrix(0, nrow = num_sample, ncol = 2, 
                    dimnames = list(colnames(BC_matrix_i), c("Average", "Median")))
  
  for (i in 1:num_sample) {
    no_sample <- BC_matrix_i[i, -i, drop = FALSE]
    results [i, "Average"] <- mean(no_sample, na.rm = TRUE)
    results [i, "Median"] <- median(no_sample, na.rm = TRUE)
  }
  results <- as.data.frame(results)
  return(results)
}

First_beta <- Median_Average_no_sample (1)
head(First_beta)

Second_beta <- Median_Average_no_sample (2)
head(Second_beta)

Third_beta <- Median_Average_no_sample (3)
head(Third_beta)

Fourth_beta <- Median_Average_no_sample (4)
head(Fourth_beta)

# Create table with average/median for all samples (all visits)
BC_tot <- rbind(First_beta,Second_beta,Third_beta,Fourth_beta)
head(BC_tot)
dim(BC_tot)

head(tab)

## Satab## Save output
write.table(tab,  "Distribution_Species_Beta.csv", sep = "\t") # Distribution of beta divesrity for each weeks
write.table(tab,  "Distribution_Species_Beta.txt", sep = "\t")

write.table(BC_tot, "Average_median_Beta_Species.csv", sep = "\t") # Median and average of beta diversity across weeks
write.table(BC_tot, "Average_median_Beta_Species.txt", sep = "\t")
