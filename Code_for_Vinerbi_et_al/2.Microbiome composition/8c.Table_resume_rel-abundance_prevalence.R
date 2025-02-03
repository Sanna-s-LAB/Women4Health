# Rscript

# Script to prepare the table with all information about relative abundance and prevalence for each samples
# Note: for this script a complete metadata file, with all the sample IDs and IDs for paper.

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(scales)
library(ggpubr)

# Paths
setwd("/home/")
getwd()

# Create a resume table (both genus and species) 
#with row = taxa, column = samples and two column with average of relative abundance and prevalence in all samples for each taxa

# Import table (output previuos script '8a.Relative_abundance.R')
# Species table
df_sp <- read.table('Total_abundance_species.txt', sep = '\t', header = T)
head(df_sp)[1:5]
dim(df_sp)

# Genus table
df_g <- read.table('Total_abundance_genus.txt', sep = '\t', header = T)
head(df_g)[1:5]
dim(df_g)

# Calculate an average and prevalence (for all taxa)
df_sp$Average <- rowMeans(df_sp[,-1])
head(df_sp)[210:214]

col <- colnames(df_sp[-1])
threshold <- 0.0
df_sp$Prevalence <- (rowSums(df_sp[,-1] > threshold) / (ncol(df_sp) - 1)) * 100
head(df_sp)[210:215]
dim(df_sp)

df_sp$Taxonomic_level <- 'Species'
df_sp <- select(df_sp, Taxonomic_level,1:215)
head(df_sp)

df_g$Average <- rowMeans(df_g[,-1])
head(df_g)[210:214]

col <- colnames(df_g[-1])
threshold <- 0.0
df_g$Prevalence <- (rowSums(df_g[,-1] > threshold) / (ncol(df_g) - 1)) * 100
head(df_g)[210:215]

df_g$Taxonomic_level <- 'Genus'
df_g <- select(df_g, Taxonomic_level,1:215)
head(df_g)

# Create a unique table
Tab <- rbind(df_g, df_sp)
head(Tab) [1:5]
dim(Tab)

# Change samples name with  name for paper
# Import metadata file 
# example of metadata columns: Plate; ID; Code_women; Visits; ID_paper; Var1; Var2
Meta <- read.table('Metadata_complete.csv', sep = ';', header = T)
Meta <- na.omit(Meta)
head(Meta)
dim(Meta)

A <- colnames(Tab[,-c(1,2,215,216)])
length(A)
A

B <- Meta$Code_women
length(B)

identical(length(A), length(B)) # check 

colnames(Tab)[3:214] <- Meta$ID_paper
colnames(Tab)

# Save output
write.table(Tab, "Tab_summary_paper_Genus_Species.csv", sep = "\t")
write.table(Tab, "Tab_summary_paper_Genus_Species.txt", sep = "\t")

## == Prevalence of first two taxa for each weeks == ##
df <- df_sp

# divided the samples in weeks
# Follicular phase (F)
S_1 <- grep("_1$", names(df), value = TRUE)
S_1
length(S_1)
D1 <- df %>% select(all_of(S_1))
head(D1) [1:5]
dim(D1)

D1 <- D1 [1:2,]
D1$Taxa <- rownames(D1)
D1 <- select(D1, Taxa, 1:59)
dim(D1)

# How many samples have iners and crispatus abundances greater than zero? 
sum(D1[D1$Taxa == "Lactobacillus iners", -1] > 0) 
sum(D1[D1$Taxa == "Lactobacillus crispatus", -1] > 0) 

#  Calculate a prevalence and relative abundance
col <- colnames(D1[-1])
threshold <- 0.0
D1$Prevalence <- (rowSums(D1[,-1] > threshold) / (ncol(D1) - 1)) * 100
dim(D1)

D1$Average <- rowMeans(D1[,-1])
head(D1)

D1$Phases <- 'F'

D1 <- select(D1, Taxa, Prevalence, Average, Phases)
head(D1)

# Ovulatory phase (O)
S_2 <- grep("_2$", names(df), value = TRUE)
S_2
length(S_2)
D2 <- df %>% select(all_of(S_2))
head(D2) [1:5]
dim(D2)

D2 <- D2 [1:2,]
D2$Taxa <- rownames(D2)
dim(D2)
D2 <- select(D2, Taxa, 1:58)
dim(D2)

# How many samples have iners and crispatus abundances greater than zero? 
sum(D2[D2$Taxa == "Lactobacillus iners", -1] > 0) # 58/59
sum(D2[D2$Taxa == "Lactobacillus crispatus", -1] > 0) # 56/59

# Calculate a prevalence and relative abundance
col <- colnames(D2[-1])
threshold <- 0.0
D2$Prevalence <- (rowSums(D2[,-1] > threshold) / (ncol(D2) - 1)) * 100
dim(D2)

D2$Average <- rowMeans(D2[,-1])
head(D2)

D2$Phases <- 'O'

D2 <- select(D2, Taxa, Prevalence, Average, Phases)
head(D2)

# Early Luteal phase (EL)
S_3 <- grep("_3$", names(df), value = TRUE)
S_3
length(S_3)
D3 <- df %>% select(all_of(S_3))
head(D3) [1:5]
dim(D3)

D3 <- D3 [1:2,]
D3$Taxa <- rownames(D3)
dim(D3)
D3 <- select(D3, Taxa, 1:52)
dim(D3)
head(D3)
 
# How many samples have iners and crispatus abundances greater than zero? 
sum(D3[D3$Taxa == "Lactobacillus iners", -1] > 0) 
sum(D3[D3$Taxa == "Lactobacillus crispatus", -1] > 0) 
 
# Calculate a prevalence and relative abundance
col <- colnames(D3[-1])

threshold <- 0.0
D3$Prevalence <- (rowSums(D3[,-1] > threshold) / (ncol(D3) - 1)) * 100
dim(D3)

D3$Average <- rowMeans(D3[,-1])
head(D3)

D3$Phases <- 'EL'

D3 <- select(D3, Taxa, Prevalence, Average, Phases)
head(D3)

# Late Luteal phase (LL)
S_4 <- grep("_4$", names(df), value = TRUE)
S_4
length(S_4)
D4 <- df %>% select(all_of(S_4))
head(D4) [1:5]
dim(D4)

D4 <- D4 [1:2,]
D4$Taxa <- rownames(D4)
dim(D4)
D4 <- select(D4, Taxa, 1:43)
dim(D4)
head(D4)

# How many samples have iners and crispatus abundances greater than zero? 
sum(D4[D4$Taxa == "Lactobacillus iners", -1] > 0) # 43/43
sum(D4[D4$Taxa == "Lactobacillus crispatus", -1] > 0) # 43/43

# Calculate a prevalence and relative abundance
col <- colnames(D4[-1])

threshold <- 0.0
D4$Prevalence <- (rowSums(D4[,-1] > threshold) / (ncol(D4) - 1)) * 100
dim(D4)

D4$Average <- rowMeans(D4[,-1])
 
D4$Phases <- 'LL'

D4 <- select(D4, Taxa, Prevalence, Average, Phases)
head(D4)

# Create a only one table
Tab <- rbind(D1, D2)
Tab <- rbind(Tab, D3)
Tab <- rbind(Tab, D4)
head(Tab)

# Plot
Tab$Phases <- factor(Tab$Phases, levels = c('F', 'O', 'EL', 'LL'))
Tab$Taxa <- factor(Tab$Taxa, levels = c ('Lactobacillus iners', 'Lactobacillus crispatus'))

# Average
P_av <- ggplot(Tab, aes(x = factor(Phases), y = Average, fill = Taxa)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_bw() +
  geom_text(aes(label = round(Average, 2)), position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  labs(title = "", x = "", y = "Average of relative Abundance") +
  scale_fill_manual(values = c("Lactobacillus iners" = "#A6CEE3", "Lactobacillus crispatus" = "#1F78B4")) + 
  ylim(0,100)
P_av

ggsave("Plot_Average_rel_ab_weeks.png", P_av, width = 18, height = 7, dpi = 300)

# Prevalence
P_pre <- ggplot(Tab, aes(x = factor(Phases), y = Prevalence, fill = Taxa)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_bw() +
  geom_text(aes(label = round(Prevalence, 2)), position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  labs(title = "", x = "Phases", y = "Prevalence") +
  scale_fill_manual(values = c("Lactobacillus iners" = "#A6CEE3", "Lactobacillus crispatus" = "#1F78B4")) 
  ylim(0,100)
P_pre

ggsave("Plot_Prevalence_weeks.png", P_pre, width = 18, height = 7, dpi = 300)

# Panel
combined_plot <- ggarrange(P_av, P_pre, nrow = 2, labels = c("A", "B"),
                          common.legend = TRUE, legend = "right")+
   theme_bw()
combined_plot

# Save output
ggsave("Panel_average_Prevalence_weeks.png", combined_plot, width = 18, height = 7, dpi = 300)


