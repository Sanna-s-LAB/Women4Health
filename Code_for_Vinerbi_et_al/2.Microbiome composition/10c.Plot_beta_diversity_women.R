# Rscript

# Script to make a plot with average of beta diversity across weeks 
# and a plot of beta diversity across women (see above script "10b.Beta_diversity_women") 

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

library(tidyverse) 
library(ggplot2)
library(vegan) 
library(reshape2)
library(ggpubr)
library(gridExtra)


setwd("/home/")
getwd()

# Import table 
# Metada file
# example of metadata columns: Plate; ID; Code_women; Visits; ID_paper; Var1; Var2
Metadata <- read.table ("Metadata.csv", sep = ";", header= TRUE)
head(Metadata) 
colnames(Metadata)[3] <- 'samples'

# Import table (output previous script '10b.Beta_diversity_women.R')
# Average and median of beta diversity
Same_w <- read.table("Average_median_Beta_same_women.txt", sep = "\t", header= T)
colnames(Same_w) <- c("Women_itself", "Average_itself", "Median_itself")
head(Same_w)

rm_1 <- read.table("Beta_diversity_random_women_rm_1.txt", sep = "", header= T)
colnames(rm_1) <- c("Women_rm1", "Average_rm1", "Median_rm1")
dim(rm_1)
head(rm_1)

rm_2 <- read.table("Beta_diversity_random_women_rm_2.txt", sep = "", header= T)
colnames(rm_2) <- c("Women_rm2", "Average_rm2", "Median_rm2")
head(rm_2)

rm_3 <- read.table("Beta_diversity_random_women_rm_3.txt", sep = "", header= T)
colnames(rm_3) <- c("Women_rm3", "Average_rm3", "Median_rm3")
head(rm_3)

rm_4 <- read.table("Beta_diversity_random_women_rm_4.txt", sep = "", header= T)
colnames(rm_4) <- c("Women_rm4", "Average_rm4", "Median_rm4")
dim(rm_4)
head(rm_4)

rm_5 <- read.table("Beta_diversity_random_women_rm_5.txt", sep = "", header= T)
colnames(rm_5) <- c("Women_rm5", "Average_rm5", "Median_rm5")
head(rm_5)

rm_6 <- read.table("Beta_diversity_random_women_rm_6.txt", sep = "", header= T)
colnames(rm_6) <- c("Women_rm6", "Average_rm6", "Median_rm6")
head(rm_6)

rm_7 <- read.table("Beta_diversity_random_women_rm_7.txt", sep = "", header= T)
colnames(rm_7) <- c("Women_rm7", "Average_rm7", "Median_rm7")
head(rm_7)

rm_8 <- read.table("Beta_diversity_random_women_rm_8.txt", sep = "", header= T)
colnames(rm_8) <- c("Women_rm8", "Average_rm8", "Median_rm8")
head(rm_8)

rm_9 <- read.table("Beta_diversity_random_women_rm_9.txt", sep = "", header= T)
colnames(rm_9) <- c("Women_rm9", "Average_rm9", "Median_rm9")
head(rm_9)

rm_10 <- read.table("Beta_diversity_random_women_rm_10.txt", sep = "", header= T)
colnames(rm_10) <- c("Women_rm10", "Average_rm10", "Median_rm10")
head(rm_10)

#box_plot for Average
combined_data <- data.frame(Value = c(Same_w$Average_itself, rm_1$Average_rm1, rm_2$Average_rm2, rm_3$Average_rm3, 
                                      rm_4$Average_rm4, rm_5$Average_rm5, rm_6$Average_rm6, rm_7$Average_rm7, 
                                      rm_8$Average_rm8, rm_9$Average_rm9, rm_10$Average_rm10),
                            Group = c(rep("Same_w", length(Same_w$Average_itself)),
                                      rep("rm_1", length(rm_1$Average_rm1)),
                                      rep("rm_2", length(rm_2$Average_rm2)),
                                      rep("rm_3", length(rm_3$Average_rm3)),
                                      rep("rm_4", length(rm_4$Average_rm4)),
                                      rep("rm_5", length(rm_5$Average_rm5)),
                                      rep("rm_6", length(rm_6$Average_rm6)),
                                      rep("rm_7", length(rm_7$Average_rm7)),
                                      rep("rm_8", length(rm_8$Average_rm8)),
                                      rep("rm_9", length(rm_9$Average_rm9)),
                                      rep("rm_10", length(rm_10$Average_rm10))))

# Change name of column
combined_data$Group[combined_data$Group == "Same_w"] <- "Self"
combined_data$Group[combined_data$Group == "rm_1"] <- "Random1"
combined_data$Group[combined_data$Group == "rm_2"] <- "Random2"
combined_data$Group[combined_data$Group == "rm_3"] <- "Random3"
combined_data$Group[combined_data$Group == "rm_4"] <- "Random4"
combined_data$Group[combined_data$Group == "rm_5"] <- "Random5"
combined_data$Group[combined_data$Group == "rm_6"] <- "Random6"
combined_data$Group[combined_data$Group == "rm_7"] <- "Random7"
combined_data$Group[combined_data$Group == "rm_8"] <- "Random8"
combined_data$Group[combined_data$Group == "rm_9"] <- "Random9"
combined_data$Group[combined_data$Group == "rm_10"] <- "Random10"

combined_data$Group <- factor(combined_data$Group, levels = c("Self", "Random1", "Random2", "Random3",
                                                              "Random4", "Random5", "Random6", "Random7",
                                                              "Random8", "Random9", "Random10"))

# Palette
colors <- c("Self" = "#FFCC99",  
            "Random1" = "#FFE6B3", "Random2" = "#FFE6B3", "Random3" = "#FFE6B3",
            "Random4" = "#FFE6B3", "Random5" = "#FFE6B3", "Random6" = "#FFE6B3",
            "Random7" = "#FFE6B3", "Random8" = "#FFE6B3", "Random9" = "#FFE6B3", "Random10" = "#FFE6B3")


plot1 <- ggplot(combined_data, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE, color = "black") +  
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) +
  scale_fill_manual(values = colors) +  
  theme_bw() +
  labs(title = "",
       x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, family = 'sans'),
        axis.title = element_text(size = 7, family = 'sans'),
        axis.text.y = element_text(size = 7, family = 'sans'),
        legend.position = "none")

plot1
ggsave("Plot_average_women_beta_diversity.png", plot1, width = 15, height = 7, dpi = 300)


# Create a tabele for plot of average of beta diversity

# beta diversity(Average and median) table
tab <- read.table("Average_median_Beta_Species.txt", header=TRUE, sep = "\t", stringsAsFactors = TRUE)
tab <- data.frame(rownames(tab), tab)
colnames(tab)[1] <- "samples"
head(tab)

# Merge metadata file and beta diversity table
metadata_BC <- merge (Metadati, tab, by ="samples", all = TRUE)
head(metadata_BC)

# Create a table
Beta <- data.frame(metadata_BC$Visit, metadata_BC$Average, metadata_BC$Group) # Change median/average
colnames(Beta)<-c("Visit","Average","Group")
head(Beta)
ord <- c("First","Second","Third","Fourth")
Beta$Visit <- factor(Beta$Visit, levels = ord)
Beta$Average <- as.numeric(Beta$Average)
head(Beta)

# Palette
col <- c( "First" = "#FFE066","Second" = "#FFF5CC", "Third" = "#FFC54C", "Fourth" = "#FFA32B") 
col

# Violin plot
# for the values of the segments see the results of the statistical analyses (see script '10d.test_stats_beta_diversity')

Plot2 <- Beta %>%
  ggplot() +
  theme_bw() +
  geom_violin(aes(x = factor(Visit), y = Average, fill=factor(Visit))) +
  geom_boxplot(aes(x = factor(Visit), y = Average), width=0.1, alpha=0.6) +
  scale_x_discrete(labels = c("F", "O", "EL", "LL")) +
  labs(title = "", x = "", y = "Average (Beta diversity)") +
  scale_fill_manual(name = "Weeks", values = col, labels = c("F", "O", "EL", "LL"))+
  theme(axis.title = element_text(size = 7, family = 'sans'),
        axis.text.x = element_text(size = 7, family = 'sans'),
        axis.text.y = element_text(size = 7, family = 'sans'),
        legend.position = 'none') +
  geom_segment(aes(x = 1, xend = 2, y = 80, yend = 80), color = "black", size = 0.5) +
  geom_text(aes(x = 1.5, y = 82, label = "*"), size = 3) +
  geom_segment(aes(x = 1, xend = 3, y = 75, yend = 75), color = "black", size = 0.5) +
  geom_text(aes(x = 2.5, y = 77, label = "*"), size = 3)
Plot2

# save plot
ggsave("Plot_average_beta_diversity.png", plot2, width = 15, height = 7, dpi = 300)
