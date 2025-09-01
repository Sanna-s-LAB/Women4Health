# Rscript

# script to check CST statistics and plot for cst or 
# NB: for sbCSTs it is necessary to replace the CST label in the script

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 27/08/2025

# Dependence: R 4.4.1

# Upload library
library(tidyverse)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(patchwork)


# Path
setwd('/home/')
getwd()

# Import table of CST/subCST (output by '12a.Prepare_table_for_VALENCIA_CST.R')
df <- read.table('Output_Valencia.csv', sep = ',', header = T)
head(df)[1:5]

# take samples, CST, CST_subgroup and score columns
df1 <- select(df, sampleID, CST, subCST)
head(df1)

# merge subCST of IV group in only one CST (CST IV)
keep_IV <- function(x) {   
  sub("IV-[A-Z]", "IV", x)
}

df1$CST <- sapply(df1$CST, keep_IV)
head(df1)

# Create a column with visit information
df1$Visits <- sub(".*_", "", df1$sampleID)
head(df1)
dim(df1)


# Calculate counts and  percentage of CST/subCST between all visits
cst_counts <- table(df1$CST)
head(cst_counts)

cst_percent <- prop.table(cst_counts) * 100
head(cst_percent)

# Merge table
tab <- merge(cst_counts, cst_percent, by =  'Var1')
colnames(tab) <- c('CST', 'Counts', 'Percentage')
head(tab)

# trasform in data frame the table
tab <- as.data.frame(tab)
str(tab)

# create a specific palette for CST
colors <- c("#FF8C6A", "#FFDA6A", "#CAD94A", "#6AFF8A", "#6A9EFF")
colors

# Specific palette for subCST
#   colors <- c('subCST-IA' = "#FF8C6A", 'subCST-IB' = "#FFB48A", 'subCST-II' = "#FFDA6A", 'subCST-II' = "#FFE28A",  
#   'subCST-IIIA' = "#FFFF6A", 'subCST-IIIB' = "#FFFF8A",  
#   'subCST-IVB' = "#6AFF8A", 'subCST-IVC0' = "#8AFFA5", 'subCST-IVC1' = "#AAFFBF", 'subCST-IVC3' = "#C2FFDB",  
#   'subCST-V' = "#6A9EFF")

# pie chart plot
P <- ggplot(tab, aes(x = "", y = Percentage, fill = CST)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = colors) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        legend.position = "none") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.7),
            color = "black", size = 5)

P
ggsave("Pie_chart_percentage.png", P, width = 12, height = 7, dpi = 300)

#Calculate percentage of CST for each visit
head(df1)
#
# dim(subset(df1, Visits == '1'))
# dim(subset(df1, Visits == '2'))
# dim(subset(df1, Visits == '3'))
# dim(subset(df1, Visits == '4'))

# calculate counts e percentage
freq_table <- table(df1$Visit, df1$CST)
head(freq_table)

percent_table <-prop.table(freq_table, margin = 1) * 100
head(percent_table)

# trasform table in data.frame
percent_table <- as.data.frame(percent_table)
colnames(percent_table) <- c('Phases','CST', 'Percentage')
head(percent_table)

# barplot
P2 <- ggplot(percent_table, aes(x =Phases, y = Percentage, fill = CST)) +
  geom_bar(stat = "identity", position = 'dodge') +
  labs(title = "", x = "", y = "Value (%)") +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ylim(0,50) +
  scale_x_discrete(labels = c("1" = "F", "2" = "O", "3" = "EL", "4" = "LL")) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))
P2

ggsave("CST_across_phases.png", P2, width = 12, height = 7, dpi = 300)


# create a panel with two plots
combined_plot <- P + P2 +
  plot_annotation(tag_levels = 'A')
combined_plot

#Save output
ggsave("Combined_CST_Panel.png", combined_plot, width = 15, height = 7, dpi = 300)