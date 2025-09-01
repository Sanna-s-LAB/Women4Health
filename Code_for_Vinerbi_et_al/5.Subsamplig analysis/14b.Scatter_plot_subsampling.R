# Rscript

# Script to compare a results of subsampling data with original data

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# version R: v4.4.1

 # Upload library
library(tidyverse)
library(ggplot2)
library(VennDiagram)
library(ggvenn)
library(ggpubr)
library(cowplot)

setwd('/home/')
getwd() 

# Import table (Output previous script '8a.Relative_abundance.R')
# relative abundance of original data
df1 <- read.table("Total_abundance_species.txt", sep = "\t")
head(df1) [1:5]
dim(df1)

# relative abundance of subsampling data30k
df30k <- read.table("Total_abundance_species_30k.txt", sep = "\t")
head(df30k)[1:5]
dim(df30k)

# relative abundance of subsampling data30k
df60k <- read.table("Total_abundance_species_60k.txt", sep = "\t")
head(df60k)[1:5]
dim(df60k)

# Calculate average and prevalence for each taxa
# Original data
df1$Average <- rowMeans(df1[,-1])
View(df1)

col <- colnames(df1[-1])

threshold <- 0.0
df1$Prevalence <- (rowSums(df1[,-1] > threshold) / (ncol(df1) - 1)) * 100
#View(df1)

OD <- select(df1, Taxa, Average, Prevalence)
head(OD)

# Sub-samplig data 30k
df30k$Average <- rowMeans(df30k[,-1])
dim(df30k)
head(df30k)[212:214]

col <- colnames(df30k[-1])

threshold <- 0.0
df30k$Prevalence <- (rowSums(df30k[,-1] > threshold) / (ncol(df30k) - 1)) * 100
head(df30k)[212:215]

Sub30k <- select(df30k, Taxa, Average, Prevalence)
head(Sub30k)

# Sub-samplig data 60k
df60k$Average <- rowMeans(df60k[,-1])
dim(df60k)
head(df60k)[212:214]

col <- colnames(df60k[-1])

threshold <- 0.0
df60k$Prevalence <- (rowSums(df60k[,-1] > threshold) / (ncol(df60k) - 1)) * 100
head(df60k)[212:215]

Sub60k <- select(df60k, Taxa, Average, Prevalence)
head(Sub60k)

# Take a common taxa (OD and 30k)
A <- OD$Taxa
B <- Sub30k$Taxa

Common_taxa1 <- intersect(A, B)
Common_taxa1
length(Common_taxa1)

OD_common <- subset(OD, Taxa %in% Common_taxa1)
head(OD_common)
dim(OD_common)
colnames(OD_common)[2:3] <- c('Average_OD', 'Prevalence_OD')
head(OD_common)

Sub30k_common <- subset(Sub30k, Taxa %in% Common_taxa1)
head(Sub30k_common)
dim(Sub30k_common)
colnames(Sub30k_common)[2:3] <- c('Average_Sub30k', 'Prevalence_Sub30k')
head(Sub30k_common)

df_OD_30k <- merge(OD_common, Sub30k_common, by = 'Taxa')
head(df_OD_30k)
dim(df_OD_30k)

# Take a common taxa (OD and 60k)
A <- OD$Taxa
C <- Sub60k$Taxa

Common_taxa2 <- intersect(A, C)
#Common_taxa2
length(Common_taxa2)

OD_common <- subset(OD, Taxa %in% Common_taxa2)
head(OD_common)
dim(OD_common)
colnames(OD_common)[2:3] <- c('Average_OD', 'Prevalence_OD')
head(OD_common)

Sub60k_common <- subset(Sub60k, Taxa %in% Common_taxa2)
head(Sub60k_common)
dim(Sub60k_common)
colnames(Sub60k_common)[2:3] <- c('Average_Sub60k', 'Prevalence_Sub60k')
head(Sub60k_common)

df_OD_60k <- merge(OD_common, Sub60k_common, by = 'Taxa')
head(df_OD_60k)
dim(df_OD_60k)

# Calculate an correlation test (Spearmen test) for average and prevalence (30k)
spearman_average30k <- cor.test(df_OD_30k$Average_OD, df_OD_30k$Average_Sub30k, method = "spearman", exact = F)
spearman_average30k
rho_average30k <- spearman_average30k$estimate 

spearman_prevalence30k <- cor.test(df_OD_30k$Prevalence_OD, df_OD_30k$Prevalence_Sub30k, method = "spearman", exact = F)
spearman_prevalence30k
rho_prevalence30k <- spearman_prevalence30k$estimate 

# Calculate an correlation test (Spearmen test) for average and prevalence (60k)
spearman_average60k <- cor.test(df_OD_60k$Average_OD, df_OD_60k$Average_Sub60k, method = "spearman", exact = F)
spearman_average60k
rho_average60k <- spearman_average60k$estimate 

spearman_prevalence60k <- cor.test(df_OD_60k$Prevalence_OD, df_OD_60k$Prevalence_Sub60k, method = "spearman", exact = F)
spearman_prevalence60k
rho_prevalence60k <- spearman_prevalence60k$estimate

# Plot
# Specific palette for specific taxa 
col <-c('Lactobacillus iners' = "#A6CEE3", 'Lactobacillus crispatus' = "#1F78B4",'Lactobacillus gasseri' = "#33A02C",
        'Lactobacillus jensenii' =  "#E31A1C")
col

# create a column where specifay the color for each taxa - AVERAGE 30k
df_OD_30k$Color <- ifelse(df_OD_30k$Taxa %in% names(col), col[df_OD_30k$Taxa], "black")
df_OD_30k$Legend <- ifelse(df_OD_30k$Taxa %in% names(col), df_OD_30k$Taxa, "Other species")

P_average1 <- ggplot(df_OD_30k, aes(x = Average_Sub30k, y = Average_OD, color = Legend)) +
  geom_point(size = 5) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed", size = 0.5) +
  annotate("text",  x = 50, y = 1.5, label = paste0("Spearman rho: ", round(rho_average30k, 2)), size = 10, hjust =0, color = "orange") + 
  labs(title = "", y = "Average relative abundance - main", 
       x = "Average relative abundance - subsampling 30k", color = "Taxa")+
  ylim(0, 100) +
  xlim(0, 100) +
  theme_bw() +
  theme(#axis.title.x = element_text(size = 18),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 25))+
  scale_color_manual(values = c(col, "Other species" = "black")) 
P_average1

# create a column where specifay the color for each taxa - AVERAGE 60k
df_OD_60k$Color <- ifelse(df_OD_60k$Taxa %in% names(col), col[df_OD_60k$Taxa], "black")
df_OD_60k$Legend <- ifelse(df_OD_60k$Taxa %in% names(col), df_OD_60k$Taxa, "Other species")

P_average2 <- ggplot(df_OD_60k, aes(x = Average_Sub60k, y = Average_OD, color = Legend)) +
  geom_point(size = 5) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed", size = 0.5) +
  annotate("text",  x = 50, y = 1.5, label = paste0("Spearman rho: ", round(rho_average60k, 2)), size = 10, hjust =0, color = "orange") + 
  labs(title = "", y = "", 
       x = "Average relative abundance - subsampling 60k", 
       color = "Taxa") +
  ylim(0, 100) +
  xlim(0, 100) +
  theme_bw() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20))+
  scale_color_manual(values = c(col, "Other species" = "black")) 
P_average2

# create a column where specifay the color for each taxa - PREVALENCE 30k
P_prevalence1 <- ggplot(df_OD_30k, aes(x = Prevalence_Sub30k, y = Prevalence_OD, color = Legend)) +
  geom_point(size = 5) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed", size = 0.5) +
  annotate("text",  x = 50, y = 1.5, label = paste0("Spearman rho: ", round(rho_prevalence30k, 2)), size = 10, hjust =0, color = "orange") + 
  labs(title = "", y = "Prevalence - main", 
       x = "prevalence - Subsampling 30k", 
       color = "Taxa") +  # Titolo della legenda
  ylim(0, 100) +
  xlim(0, 100) +
  theme_bw() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20))+
  scale_color_manual(values = c(col, "Other species" = "black")) 
P_prevalence1

# create a column where specifay the color for each taxa -  60k
P_prevalence2 <- ggplot(df_OD_60k, aes(x = Prevalence_Sub60k, y = Prevalence_OD, color = Legend)) +
  geom_point(size = 5) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed", size = 0.5) +
  annotate("text",  x = 50, y = 1.5, label = paste0("Spearman rho: ", round(rho_prevalence60k, 2)), size = 10, hjust =0, color = "orange") + 
  labs(title = "", y = "Prevalence - main", 
       x = "Prevalence - Subsampling 60k") +
  ylim(0, 100) +
  xlim(0, 100) +
  theme_bw() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20))+
  scale_color_manual(values = c(col, "Other species" = "black")) 
P_prevalence2

# create a pannel with all plots
# Save output 
png("10.Subsampling_30k_60K_080125/3.Plot_comparison_OD_subsampling_data_080125/Pannel_140125.png", width = 35, height = 20, units = "in", res = 300)

# Remove legend
custom_theme <- function() {
  theme_bw() +
    theme(
      axis.title = element_text(size = 25),
      axis.text = element_text(size = 20),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    )
}

P_average1 <- P_average1 + custom_theme()
P_average2 <- P_average2 + custom_theme()
P_prevalence1 <- P_prevalence1 + custom_theme()
P_prevalence2 <- P_prevalence2 + custom_theme()


P_average1 <- P_average1 + theme(legend.position = "none")
P_average2 <- P_average2 + theme(legend.position = "none")
P_prevalence1 <- P_prevalence1 + theme(legend.position = "none")
P_prevalence2 <- P_prevalence2 + theme(legend.position = "none")

legend <- get_legend(
  P_average1 + theme(legend.position = "right",
                     legend.text = element_text(face = "italic", family = 'serif', size = 25)))
o
combined_plot <- ggarrange( P_average1, P_average2, P_prevalence1, P_prevalence2,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),font.label = list(size = 25),
  label.x = 0.1, 
  common.legend = FALSE)


# Add legend
final_plot <- ggdraw() +
  draw_plot(combined_plot, 0, 0, 1,1 ) +
  draw_plot(legend, x = 0.07, y = 0.75, width = 0.2, height = 0.2)
final_plot

# Save ouput
dev.off()



