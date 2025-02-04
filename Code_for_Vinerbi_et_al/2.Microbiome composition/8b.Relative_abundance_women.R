# Rscript

# Script to make a plot for each weeks and all women

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(scales)
library(patchwork)
library(ggpubr)


# files path
setwd("/home/")
getwd()

# Import table (output of previous script '8a.Relative_abundance')
# table of relative abundance of first ten taxa and others category (row = taxon; column = samples)
df <- read.table("Top10_rel_ab_species_others.txt", sep="\t", header=T)
head(df)[1:5]
dim(df)

# Calculate an average of relative abundance to order the taxa
mean_ab <- rowMeans (df[,-1])
head(df)[1:5]

Top <- names(mean_ab[order(mean_ab, decreasing=TRUE)])
Top

a <- match(Top, rownames(df))
a
Top <- df[a,]
head(Top) [1:5]
dim(Top)

colnames(Top)[1] <- "Taxa"
head(Top) [1:5]
dim(Top)
#View(Top)

# vector witha taxa name
taxa <- rownames(Top)
taxa
length(taxa)

# Palette 
palette2 <- brewer.pal(11, "Paired")  
col_taxa <- setNames(palette2, taxa)
col_taxa

# Take samples of first week
Sample_1 <- grep("_1$", names(Top), value = TRUE)
Sample_1
length(Sample_1)
df_1 <- Top %>% select(all_of(Sample_1))
head(df_1) [1:5]
dim(df_1)

First <- data.frame(rownames(df_1), df_1)
colnames(First)[1] <- "Taxa"

## Add missing women for this week
First$Code_women <- 0

dim(First)

M <- melt(First)
dim(M)

# Order with first trhee most abundant species 
A <- M %>%
  filter(Taxa == "Lactobacillus iners" & value >= 30) %>%
  arrange(desc(value))
dim(A)

B <- M %>%
  filter(Taxa == "Lactobacillus crispatus"& value >= 30) %>%
  arrange(desc(value))
dim(B)

D <- M %>%
  filter(Taxa == "Lactobacillus gasseri") %>%
  arrange(desc(value))
dim(D)

## Merge tables
Tab_def <- rbind(A, B, D)
head(Tab_def)

S<- anti_join(M, Tab_def, by = c("Taxa", "value"))
head(S)

S<- S %>%
  arrange(desc(value))
dim(M)

Tab_def <- rbind (Tab_def,S) 
dim(Tab_def)

# Create an index for order in specific mode (same samples with a same index)
vett <- NULL
i <- 1

for (code in Tab_def$variable) {
  if (!(code %in% names(vett))){
    vett[code] <- i
    i <- i+1
  }
}  
#vett

index_convention <- function(code){ 
  return(vett[[code]])
}

A <- as.vector(Tab_def$variable)
A

Tab_def$Index <- mapply(index_convention, A)
head(Tab_def)

# Plot
Tab_def <- transform(Tab_def, Taxa = reorder(Taxa, -value))
Tab_def$Species <- factor(Tab_def$Taxa, levels = taxa)

P1 <- ggplot(Tab_def, aes(x = reorder(variable, Index), y = value, fill = Taxa)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = col_taxa)+
  labs(title = "", x = 'Samples (n=59)', y = 'Relative abundance') + # change number of women
  theme(axis.text.x = element_blank(),
    axis.title = element_text(size = 7, family = 'sans'),
    legend.text = element_text(size = 7, family = 'sans'),
    legend.title = element_blank())
P1

# save output
ggsave("Plot_First_all_women.png", P1, width = 15, height = 7, dpi = 300)

# create a column with women's code
Tab_def$Women <- sub("_.*", "", Tab_def$variable)

A <- unique(Tab_def$Women)
length(A)

head(Tab_def)
dim(Tab_def)

# Save the index used to order samples
First_index <- select(Tab_def, Women, Index)
head(First_index)

# Create a table for next plot
M_first <- select(Tab_def, Women, Taxa, variable, value, Index )
head(M_first)

## Save output
write.table(Tab_def, "Bacteria_first_visit_INDEX.csv", sep = "\t")
write.table(Tab_def, "Bacteria_first_visit_INDEX.txt", sep = "\t")

# To create a plot of total women each weeks remain we use the index created for week 1
# Second week
Sample_2 <- grep("_2$", names(Top), value = TRUE)
Sample_2
length(Sample_2)
df_2 <- Top %>% select(all_of(Sample_2))
head(df_2) [1:5]
dim(df_2)

Second <- data.frame(rownames(df_2), df_2)
colnames(Second)[1] <- "Taxa"

# Add missin women
Second$code_women <- 0
dim(Second)

M <- melt(Second)
dim(M)

M$Women  <- sub("_.*", "", M$variable)
head(M)

# order this table with first index
M_second <- merge(M, First_index, by = "Women")
head(M_second)

#Plot
M_second <- transform(M_second, Taxa = reorder(Taxa, -value))
M_second$Taxa <- factor(M_second$Taxa, levels = taxa)

P2 <- ggplot(M_second, aes(x = reorder(variable, Index), y = value, fill = Taxa)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = col_taxa) +
  labs(title = "", x = 'Samples (n=58)', y = 'Relative abundance') + # change number of women
  theme(
    axis.text.x = element_blank(),
    axis.title = element_text(size = 7, family = 'sans'),
    legend.text = element_text(size = 7, family = 'sans'),
    legend.title = element_blank())
P2 

# save output
ggsave("Plot_Second_all_women.png", P2, width = 15, height = 7, dpi = 300)

# Third week
Sample_3 <- grep("_3$", names(Top), value = TRUE)
Sample_3
length(Sample_3)
df_3 <- Top %>% select(all_of(Sample_3))
head(df_3)[1:5]

# create a table for third week
Third <- data.frame(rownames(df_3), df_3)
colnames(Third)[1] <- "Taxa"
dim(Third)

Third$code_women <- 0
dim(Third)

M <- melt(Third)
dim(M)

M$Women  <- sub("_.*", "", M$variable)
head(M)

M_third <- merge(M, First_index, by = "Women", all = TRUE)
head(M_third)

# Order table
M_third <- transform(M_third, Taxa = reorder(Taxa, -value))
M_third$Taxa <- factor(M_third$Taxa, levels = taxa)

# Plot
P3 <- ggplot(M_third, aes(x = reorder(variable, Index), y = value, fill = Taxa)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = col_taxa) +
  labs(title = "", x = 'Samples (n=53)', y = 'Relative abundance') + # change number of women
  theme(axis.text.x = element_blank(),
    axis.title = element_text(size = 7, family = 'sans'),
    legend.text = element_text(size = 7, family = 'sans'),
    legend.title = element_blank())
P3 

# save output
ggsave("Plot_third_all_women.png", P3, width = 15, height = 7, dpi = 300)

# Fourth week
Sample_4 <- grep("_4$", names(Top), value = TRUE)
Sample_4
length(Sample_4)
df_4 <- Top %>% select(all_of(Sample_4))
head(df_4) [1:5]

Fourth <- data.frame(rownames(df_4), df_4)
colnames(Fourth)[1] <- "Taxa"
dim(Fourth)

# Ad Missing value
Fourth$code_women <- 0
dim(Fourth)

M <- melt(Fourth)
M$Women  <- sub("_.*", "", M$variable)
head(M)

M_fourth <- merge(M, First_index, by = "Women", all = TRUE)
head(M_fourth)

# order table
M_fourth <- transform(M_fourth, Taxa = reorder(Taxa, -value))
M_fourth$Taxa <- factor(M_fourth$Taxa, levels = taxa)

# Plot
P4 <- ggplot(M_fourth, aes(x = reorder(variable, Index), y = value, fill = Taxa)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = col_taxa) +
  labs(title = "", x = 'Samples (n=44)', y = 'Relative abundance') + # change number of women
  theme(axis.text.x = element_blank(),
    axis.title = element_text(size = 7, family = 'sans'),
    legend.text = element_text(size = 7, family = "sans"),
    legend.title = element_blank())
P4 

# Save output
ggsave("Plot_LL_all_women.png", P4, width = 15, height = 7, dpi = 300)

# Create a pannel 
# Adapt margin of legend and plots

adjust_theme <- theme(legend.position = "right",
  legend.margin = margin(0, 30, 0, 0))
  #plot.margin = unit(c(1, 2, 1, 1), "cm"))

P1 <- P1 + adjust_theme
P2 <- P2 + adjust_theme
P3 <- P3 + adjust_theme
P4 <- P4 + adjust_theme

# panel
combined_plot <- ggarrange(P1, P2, P3, P4, ncol = 1, nrow = 4, 
                           labels = c("A", "B", "C", "D"),
                           common.legend = TRUE, legend = "right")

combined_plot

# Save output
ggsave("Panel_all_weeks_all_women.png",combined_plot , width = 15, height = 7, dpi = 300)

















