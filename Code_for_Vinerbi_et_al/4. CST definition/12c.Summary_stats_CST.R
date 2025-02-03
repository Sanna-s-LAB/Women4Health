# Rscript
 
# Script to make a pie chart and bar plot for CST or subCST

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# Dependence: R 4.4.1

# Upload library
library(tidyverse)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(patchwork)


# Path
setwd("/home/")
getwd()

# Import table 
# Results of VALENCIA
results_valencia <- read.table("Output_Valencia_CST.csv", sep = ",", header= T)
View(results_valencia)
head(results_valencia)[1:5]

# Select specific column
df <- select(results_valencia, sampleID, CST,subCST, score)
head(df)

keep_IV <- function(x) {    # to remove subgroups of IV gropups
  sub("IV-[A-Z]", "IV", x)
}

df$CST <- sapply(df$CST, keep_IV)
head(df)

# All these steps must be repeated for the subCSt column
# Create visits column (values are 1,2 ecc)
df$Visits <- sub(".*_", "", df$sampleID)
head(df)

# Create a percentage and counts of CST or subCST for each visits
df_summary <- df %>%
  group_by(CST) %>%  # Change with column subCST
  summarise(total_score = sum(score)) %>%
  mutate(percentage = total_score / sum(total_score) * 100)
View(df_summary)

# Specific palette for CST
colors <- c('CST-I'="#FF8C6A",'CST-II'="#FFDA6A",'CST-III'="#CAD94A",'CST-IV'="#6AFF8A",'CST-V'= "#6A9EFF")
colors

# Specific palette for subCST
# subCST_colors <- c('subCST-IA' = "#FF8C6A", 'subCST-IB' = "#FFB48A", 'subCST-II' = "#FFDA6A", 'subCST-II' = "#FFE28A",  
#   'subCST-IIIA' = "#FFFF6A", 'subCST-IIIB' = "#FFFF8A",  
#   'subCST-IVB' = "#6AFF8A", 'subCST-IVC0' = "#8AFFA5", 'subCST-IVC1' = "#AAFFBF", 'subCST-IVC3' = "#C2FFDB",  
#   'subCST-V' = "#6A9EFF")

# Pie chart
P <- ggplot(df_summary, aes(x = "", y = total_score, fill = CST)) +  # Change with column subCST
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = colors) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),  
        panel.background = element_rect(fill = "white", color = NA),  
        legend.background = element_rect(fill = "white", color = NA),
        legend.position = "none") + 
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.7), 
            color = "black", size = 5) 

P

ggsave("Pie_chart_percentage.png", P, width = 12, height = 7, dpi = 300)

# Calculate for each weeks the percentage of CST or subCST
head(df)

df_percent <- df %>%
  group_by(Visits, CST) %>%  # Change with column subCST
  summarise(count = n(), .groups = 'drop') %>%
  group_by(CST) %>%
  mutate(percentage = count / sum(count) * 100)
head(df_percent)
View(df_percent)

# Plot
df_percent$CST <- ifelse(df_percent$CST == 'I', 'CST-I', df_percent$CST)  # Change with column subCST
df_percent$CST <- ifelse(df_percent$CST == 'II', 'CST-II', df_percent$CST)  # Change with column subCST
df_percent$CST <- ifelse(df_percent$CST == 'III', 'CST-III', df_percent$CST)  # Change with column subCST
df_percent$CST <- ifelse(df_percent$CST == 'IV', 'CST-IV', df_percent$CST)  # Change with column subCST
df_percent$CST <- ifelse(df_percent$CST == 'V', 'CST-V', df_percent$CST)  # Change with column subCST
head(df_percent)

# barplot
P2 <- ggplot(df_percent, aes(x = Visits, y = percentage, fill = CST)) +
      geom_bar(stat = "identity", position = 'dodge') +
      labs(title = "", x = "", y = "Value (%)") +
      theme_bw() +
      scale_fill_manual(values = colors) +
      ylim(0,40) +
  scale_x_discrete(labels = c("1" = "F", "2" = "O", "3" = "EL", "4" = "LL")) +  
  theme(axis.text = element_text(size = 7, family = 'sans'),  
        axis.title = element_text(size = 7, family = 'sans'),
        legend.text = element_text(size = 7, family = 'sans'),  
        legend.title = element_text(size = 7, family = 'sans'))
P2  

ggsave("CST_across_visits.png", P2, width = 12, height = 7, dpi = 300)

# Create a panel
combined_plot <- P + P2 +
  plot_annotation(tag_levels = 'A')
combined_plot

#Save output
ggsave("Combined_CST_distribution.png", combined_plot, width = 15, height = 7, dpi = 300)

