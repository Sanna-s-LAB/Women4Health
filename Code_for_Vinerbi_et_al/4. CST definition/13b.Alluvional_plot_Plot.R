# Rscript

# Script to make an alluvional plot to observe if women change their CST (or subCST) during the phases

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: 4.4.1

# Upload library 
library(tidyverse)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(vegan)
library(plotly)
library(ggalluvial)

# Path
setwd("/home/") 
getwd()

# Import table
# Valencia output
results_valencia <- read.table("Output_Valencia.csv", sep = ",", header= T)
dim(results_valencia)
head(results_valencia)

df <- select(results_valencia, sampleID, CST, subCST) # select specific column
head(df)

keep_IV <- function(x) {   # to remove subgroups of IV gropups
  sub("IV-[A-Z]", "IV", x)
}

df$CST <- sapply(df$CST, keep_IV)
head(df)

# Save output
# Row= samples; colum1 = Samples; colum2 = CST; column3 =subCST
# this table is used in the script '8.tsne.R'
write.table(df, "Tab_CST_SubCST.txt", sep="\t", row.names = F)
write.table(df, "Tab_CST_SubCST.csv", sep="\t", row.names = F)

# All these steps must be repeated for the subCSt column
# Add visits column (values are 1;2 ecc)
df$Visits <- sub(".*_", "", df$sampleID) # take character after "_"
head(df)

# Add women column  (specific code for women)
df$Women <- sub("_.*", "", df$sampleID) # take character before "_"
head(df)
#View(df)

# How many samples change or not during the visits
df_var <- df %>%
  group_by(Women) %>%
  mutate(CST_change = ifelse(n_distinct(CST) > 1, TRUE, FALSE)) %>% # Change with column subCST
  ungroup()
head(df_var)

# Select women that change during the visits
df_women_var <- subset (df_var, CST_change == "TRUE") # Change with column subCST
dim(df_women_var) 
head(df_women_var)

# Select women that not change during the visits
df_women_novar <- subset (df_var, CST_change == "FALSE") # Change with column subCST
dim(df_women_novar) 

# for each woman create a column with cst for each week
df_wide <- df %>%
  select(Women, Visits, CST) %>% # Change with column subCST
  pivot_wider(names_from = Visits, values_from = CST, names_prefix = "Visit_") # Change with column subCST
colnames(df_wide) <- c("Women", "CST_week_I", "CST_week_II", "CST_week_III", "CST_week_IV") # Change with column subCST
head(df_wide)

df_1 <- select(df_var, Women, CST_change) # Change with column subCST
head(df_1)

# Take unique code for women
df_unique <- df_1 %>%
  distinct(Women, .keep_all = TRUE)
head(df_unique)
dim(df_unique) 

# Merge table
Tab <- merge(df_wide, df_unique, by = "Women")
head(Tab)
#View(Tab)

# Save output
write.table(Tab, "Women_variation_Non_variation.txt", sep="\t", row.names = F)
write.table(Tab, "Women_variation_Non_variation.csv", sep="\t", row.names = F)

# How many WOMEN change during the weeks
Tab_var <- subset(Tab, CST_change == "TRUE") # Change with column subCST
dim(Tab_var) 

rows_without_na <- sum(complete.cases(Tab_var)) # information without missing value
rows_without_na 

rows_with_na <- sum(!complete.cases(Tab_var)) # Information with missing value
rows_with_na 

## How many WOMEN no change during the weeks
Tab_novar <- subset(Tab, CST_change == "FALSE") # Change with column subCST
dim(Tab_novar) 

rows_without_na <- sum(complete.cases(Tab_novar)) 
rows_without_na 

rows_with_na <- sum(!complete.cases(Tab_novar))
rows_with_na 

## Add column that sign if there is a NA value
Tab$has_NA <- apply(Tab[, c("CST_week_I", "CST_week_II", "CST_week_III", "CST_week_IV")], 1, function(x) any(is.na(x)))
head(Tab)

for (col in colnames(Tab)) {
  Tab[[col]] <- sapply(Tab[[col]], function(X) {ifelse(is.na(X),"NA", X)})
  
}
View(Tab)

#Vector to order the column of plot (for CST value)
ord <- c("NA", "I", "II", "III", "IV", "V")
ord

#Vector to order the column of plot (for subCST value)
#ord <- c("NA", "I-A", "I-B", "II", "III-A", "III-B", "IV-B", "IV-C0", "IV-C1", "IV-C3","V")
#ord

Tab$CST_week_I <- factor(Tab$CST_week_I, levels = ord) # Change with column subCST
Tab$CST_week_II <- factor(Tab$CST_week_II, levels = ord) # Change with column subCST
Tab$CST_week_III <- factor(Tab$CST_week_III, levels = ord) # Change with column subCST
Tab$CST_week_IV <- factor(Tab$CST_week_IV, levels = ord) # Change with column subCST

Tab$label_status <- with(Tab, ifelse(CST_change == "FALSE" & has_NA == "FALSE", 
                                     "NO CH (32/61)",
                                     ifelse(CST_change == "FALSE" & has_NA == "TRUE", 
                                            "NO CH + NA (18/61)",
                                            ifelse(CST_change == "TRUE" & has_NA == "FALSE", 
                                                   "CH (9/61)",
                                                   "CH + NA (2/61)"))))

head(Tab$label_status) # Change with column subCST

# Alluvial plot
A <- ggplot(Tab, aes(axis1 = CST_week_I, axis2 = CST_week_II, axis3 = CST_week_III, axis4 = CST_week_IV)) + # Change with column subCST
  geom_alluvium(aes(fill = label_status)) +  
  theme_bw() +
  geom_stratum() +
  geom_text((stat = "stratum", aes(label = after_stat(stratum),size = after_stat(n))) # set the size of the retangles proportional to the number of samples for each CST-subCST
  scale_fill_manual(values = c('NO CH (32/61)' = "#FF8C00", 'NO CH + NA (18/61)' = "#FFD700", 'CH (9/61)' = "#006400",  # Change value
                    'CH + NA (2/61)' = "#98FB98"), name = "") +
  labs(title = "", x = NULL, y = NULL) +
  annotate("text", x = 1, y = -1, label = "F", size = 3, vjust = 1, fontface = 'bold') +  
  annotate("text", x = 2, y = -1, label = "O", size = 3, vjust = 1, fontface = 'bold') +
  annotate("text", x = 3, y = -1, label = "EL", size = 3, vjust = 1, fontface = 'bold') +
  annotate("text", x = 4, y = -1, label = "LL", size = 3, vjust = 1, fontface = 'bold') +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 7, family = 'sans'),
        legend.position = c(0.95, 0.4))

A

# Save output
ggsave("Alluvial_plot.png", A , width = 10.41, height = 6.25, dpi = 300)

