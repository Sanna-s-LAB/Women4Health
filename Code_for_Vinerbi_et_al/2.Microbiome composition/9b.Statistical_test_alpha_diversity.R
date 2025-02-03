# Rscript

# Script to calculate a statistics test for alpha diversity
# Note: for this script a complete metadata file, with all the sample IDs even when no data were present was used

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1 

# Upload library
library(rstatix)
library(ggpubr)
library(psych)
library(xtable)
library(tidyverse)
library(reshape2)
library(readxl)
library(vegan)
library(patchwork)

# path
setwd("/home/")
getwd()

# Import table (output previuos script '9.Alpha_diversity.R')
# table with results of alpha diversity
df <- read.table("Tab_alpha_Species.txt", header = TRUE, sep = "\t" )
head(df)

# Metadata file with all IDs
# example of metadata columns: Plate; ID; Code_women; Visits; ID_paper; Var1; Var2
Metadata <- read.table("Metadata_complete.txt", header = TRUE, sep = "\t" )
head(Metadata)

# Merge Metadata and results of alpha diversity
Meta_Alpha <- merge (Metadata, df, by ="ID", all= TRUE)
dim(Meta_Alpha)
head(Meta_Alpha)

# Test for Shannon Index
# Prepare atable
SH_data <- data.frame(Meta_Alpha$Visit, Meta_Alpha$Shannon_index, Meta_Alpha$Group)
colnames(SH_data) <- c("Visit", "Shannon_index","Group")
head(SH_data)
ord <- c("First", "Second", "Third", "Fourth")
SH_data$Shannon_index <- as.numeric(SH_data$Shannon_index)
head(SH_data)

# Change value of weeks with name of phases
SH_data$Weeks <- ifelse(SH_data$Visit == 'First', 'F',
                        ifelse(SH_data$Visit == 'Second', 'O',
                               ifelse(SH_data$Visit == 'Third', 'EL',
                                      ifelse(SH_data$Visit == 'Fourth', 'LL', SH_data$Visit))))
head(SH_data)




# Friendman test
FT_SHANNON <- friedman.test(SH_data$Shannon_index, groups = SH_data$Visit, blocks = SH_data$Group, data= SH_data)
FT_SHANNON
text_FT_SHANNON <- (paste ("Friedman Test","p =",round(FT_SHANNON$p.value,6)))
text_FT_SHANNON

# Paired Wilcoxon test
WT_SHANNON <- SH_data %>%
  wilcox_test(Shannon_index ~ Visit, paired=T, p.adjust.method = "bonferroni" )
head(WT_SHANNON)
WT_SHANNON$p.adj=round(WT_SHANNON$p.adj,6)

WT_sh <- as.data.frame(WT_SHANNON)
head(WT_sh)

#Save results
write.table(WT_sh, "WT_Shannon.csv", sep = "\t")
write.table(WT_sh, "WT_Shannon.txt", sep = "\t")

# Violin-Plot
# Specific palette for the four phases
col <- c( "F" = "#FFE066","O" = "#FFF5CC", "EL" = "#FFC54C", "LL" = "#FFA32B")

P_violin_Sh <- SH_data %>%
  ggplot() +
  theme_bw() +
  geom_violin(aes(x = factor(Weeks), y = Shannon_index, fill=factor(Weeks))) +
  geom_boxplot(aes(x = factor(Weeks), y = Shannon_index), width=0.1, alpha=0.6) +
  scale_x_discrete(labels = c("F", "O", "EL", "LL")) +
  labs(x = "", y = "Shannon_index (H)") +
  scale_fill_manual(values = col, labels = c("F", "O", "EL", "LL"))+
  #ylim(0,2) +
  theme(
    axis.title = element_text(size = 7, family = 'sans'),
    axis.text.x = element_text(size = 7, family = 'sans'),
    axis.text.y = element_text(size = 7, family = 'sans'),
    legend.position = 'none')
P_violin_Sh

# Save Plot
ggsave("Plot_alpha_Shannon.png", P_violin_Sh, width = 15, height = 7, dpi = 300)

# Statistics test for Simpson Index
# Prepare table

SIMPSON_data<-data.frame(Meta_Alpha$Visit, Meta_Alpha$Simpson_index, Meta_Alpha$Group)
colnames(SIMPSON_data)<-c("Visit","Simpson_index","Group")
head(SIMPSON_data)
ord <- c("First","Second","Third","Fourth")
SIMPSON_data$Visit<-factor(SIMPSON_data$Visit, levels = ord)
SIMPSON_data$Simpson_index<-as.numeric(SIMPSON_data$Simpson_index)
head(SIMPSON_data)

SIMPSON_data$Weeks <- ifelse(SIMPSON_data$Visit == 'First', 'F',
                        ifelse(SIMPSON_data$Visit == 'Second', 'O',
                               ifelse(SIMPSON_data$Visit == 'Third', 'EL',
                                      ifelse(SIMPSON_data$Visit == 'Fourth', 'LL', SIMPSON_data$Visit))))
head(SIMPSON_data)


# Friedman test
FT_SIMPSON <- friedman.test(SIMPSON_data$Simpson_index, groups = SIMPSON_data$Visit, blocks = SIMPSON_data$Group, data = SIMPSON_data)
FT_SIMPSON
text_FT_SIMPSON<- (paste("Friedman Test","p =",round(FT_SIMPSON$p.value,6)))
text_FT_SIMPSON

# paired Wilcoxon test
WT_SIMPSON <- SIMPSON_data %>%
  wilcox_test(Simpson_index ~ Visit, paired=T, p.adjust.method = "bonferroni" )
head(WT_SIMPSON)
WT_SIMPSON$p.adj=round(WT_SIMPSON$p.adj,6)

WT_sim <- as.data.frame(WT_SIMPSON)
head(WT_sim)

# Save result
write.table(WT_sim, "WT_Simpson.csv", sep = "\t")
write.table(WT_sim, "WT_Simpson.txt", sep = "\t")

# Violin Plot
# Specific palette for the four phases
col <- c( "F" = "#FFE066","O" = "#FFF5CC", "EL" = "#FFC54C", "LL" = "#FFA32B")

P_violin_simp <- SIMPSON_data %>%
  ggplot() +
  theme_bw() +
  geom_violin(aes(x = factor(Weeks), y = Simpson_index, fill=factor(Weeks))) +
  geom_boxplot(aes(x = factor(Weeks), y = Simpson_index), width=0.1, alpha=0.6) +
  scale_x_discrete(labels = c("F", "O", "EL", "LL")) +
  labs(x = "", y = "Simpson_index (D)") +
  scale_fill_manual(name = "Weeks", values = col, labels = c("F", "O", "EL", "LL"))+
  #ylim(0,1)
  theme(
    axis.title.x = element_text(size = 7, family = 'sans'),
    axis.title.y = element_text(size = 7, family = 'sans'),
    axis.text.x = element_text(size = 7, family = 'sans'),
    axis.text.y = element_text(size = 7, family = 'sans'),
    legend.position = 'none')

P_violin_simp

# Save plot
ggsave("Plot_Simpson.png", P_violin_simp, width = 15, height = 7, dpi = 300)

# Panel of alpha diversity
combined_plot <- P_violin_Sh + P_violin_simp +
  plot_annotation(tag_levels = 'A')
combined_plot

# Save output
ggsave("Combined_plot_Shannon_Simpson.png", combined_plot, width = 15, height = 7, dpi = 300)

