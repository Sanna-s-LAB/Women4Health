# Rscript

# Script to calculate a relative abundance of taxa (genus and species taxonomic level)
# Note: for this script a complete metadata file, with all the sample IDs even when no data were present was used.
#

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(scales)

# files path
setwd("/home/")
getwd()

# Import table (output of previus script 6.Change_name_ASV_ambiguos.R)
# Counts table with samples (column) and taxon (row)
df <- read.table("ASV_&_taxa_Subgenus_correctASV.txt", sep="\t", header=T)
head(df) [7:10]
dim(df)

# Metadati files (with all ID) 
# example of metadata columns: Plate; ID; Code_women; Visits; ID_paper; Var1; Var2
Meta<- read.table ("Metadata_complete.csv", sep = ";", header= TRUE)
dim(Meta)
head(Meta)

# Repeat these steps for both the taxonomic level 'genus' and 'species' - changing names where necessary
# Aggregate taxa column 
selected_columns <- names(df)[10:ncol(df)]
Tab <- aggregate(. ~ Species, data = df[, c("Species", selected_columns)], FUN = sum) # CHANGE NAME
colnames(Tab)[1] <- "Taxa"
rownames(Tab) <- Tab$Taxa
dim(Tab) 
head(Tab) [1:5]

# Change name of Unclassified taxa
# for genus level use "Unclassified Unclassified" and "f_uncl.g_uncl",
Tab$Taxa <- ifelse(Tab$Taxa == "Unclassified Unclassified Unclassified", "f_uncl.g_uncl.sp_uncl", Tab$Taxa) 
head(Tab)
View(Tab)

# taxa with unclassified name become "g_uncl" or "sp_uncl"
Tab$Taxa <- ifelse(grepl("Unclassified", Tab$Taxa), sub("Unclassified", "sp_uncl", Tab$Taxa), Tab$Taxa)
head(Tab)

rownames(Tab) <- Tab$Taxa
head(Tab)
View(Tab)

# Calculate relative abbundance 
Ab_taxa <- colSums(Tab[, -1])  
head(Ab_taxa)

Rel_ab_taxa <- matrix (nrow=dim(Tab)[1], ncol=dim(Tab)[2])
for (i in 1:(ncol(Tab)-1) ) {  Rel_ab_taxa[, i+1] <- as.numeric(Tab[, i+1] / Ab_taxa[i])} 
head(colSums(Rel_ab_taxa[,-1]))

Rel_ab_taxa = Rel_ab_taxa[,-1]
Rel_ab_taxa = Rel_ab_taxa * 100
rownames(Rel_ab_taxa) <- Tab[,1]
colnames(Rel_ab_taxa) <- colnames(Tab)[-1]
colSums(Rel_ab_taxa)
head(Rel_ab_taxa)[,5]

# Remove unclassified taxa 
rownames(Tab)

U <- Rel_ab_taxa[rownames(Rel_ab_taxa) == 'f_uncl.g_uncl.sp_uncl', ] # Change 
head(U) [1:5]

G <- Rel_ab_taxa[rownames(Rel_ab_taxa) != 'f_uncl.g_uncl.sp_uncl', ] # Change 
head(U) [1:5]
rownames(G)
dim(G)

# Calculate mean of abundance to order the taxa
mean_rel_ab <- rowMeans(G) 
mean_rel_ab

tot_ab.rel <- names(mean_rel_ab[order(mean_rel_ab, decreasing=TRUE)])

top_ab.rel <- names(mean_rel_ab[order(mean_rel_ab, decreasing=TRUE)])[1:10] 
top_ab.rel

# Check the names of samples with specific order
a <- match(tot_ab.rel, rownames(G))
a
Tot_ab.rel<- G[a,]
head(Tot_ab.rel)
dim(Tot_ab.rel)

# Create a colum "Others" 
df <- colSums(Tot_ab.rel [11:416, ]) # CHANGE INDEX
df_t <- t(df)
rownames(df_t)<- "Others"
head(df_t)[1:5]

# Table with top taxa
b <- match(top_ab.rel , rownames(Rel_ab_taxa))
b
Top_ab.rel<- Rel_ab_taxa[b,]
head(Top_ab.rel) [1:5]

Top_tab <- rbind(Top_ab.rel, df_t)
Top_tab <- data.frame(rownames(Top_tab), Top_tab)
colnames(Top_tab)[1] <- "Taxa"
head(Top_tab)[1:5]
dim(Top_tab)

# Sum unclassified and Others
Others <- rbind (U,df_t)
Others <- colSums(Others)
Others <- data.frame(Others)
Others_t <- t(Others)
Others <- data.frame(rownames(Others_t), Others_t)
colnames(Others)[1] <- "Taxa"
head(Others)

# Remove last row (others row without unclassified value)
Top <- Top_tab [-11,] 

Top <- rbind(Top, Others)
colSums(Top[,-1])  # Check that the total sum is 100%
head(Top)[1:5]

# Create a table with total samples with all visits
Samples <- Meta$Code_complete # column with all samples
missing_cols <- setdiff(Samples, colnames(Top)) # check missing samples
length(missing_cols)

A <- as.data.frame(matrix(0, nrow = nrow(Top), ncol = length(missing_cols))) # create a dataframe
rownames(A) <- rownames(Top)
colnames(A) <- missing_cols
dim(A)

# Merge table
Top_df <- cbind(Top, A)
dim(Top_df)
View(Top_df)
dim(Top_df)

# create a table with the samples ordered
df_1 <- Top_df$Taxa
head(df_1)
dim(Top_df)

df_2 <- select(Top_df, 2:245)
head(df_2)[1:4]

col_names <- names(df_2)
split_names <- strcapture("X(\\d+)_(\\d+)", col_names, proto = list(first = numeric(), second = numeric()))
order_col <- order(split_names$first, split_names$second)
df_order <- df_2[, order_col]
colnames(df_order)

Top_df <- cbind(df_1, df_order)
colnames(Top_df)[1] <- 'Taxa'
#View(Top_df)
head(Top_df)
dim(Top_df)
length(colnames(Top_df)[2:213])

View(Top_df)

## PLOT
# Specific palette (species level)
col_taxa_top10 <-c('Lactobacillus iners' = "#A6CEE3", 'Lactobacillus crispatus' = "#1F78B4", 'Others' =  "#B2DF8A",
                   'Lactobacillus gasseri' = "#33A02C", 'Bifidobacterium sp_uncl' = "#FB9A99", 'Lactobacillus jensenii' =  "#E31A1C",
                   'Fannyhessea vaginae'  =  "#FDBF6F", 'Streptococcus sp_uncl' =  "#FF7F00", 'Streptococcus agalactiae' = "#CAB2D6",
                   'Prevotella sp_uncl' =  "#6A3D9A", 'Bifidobacterium breve' = "#FFFF99")
col_taxa_top10

# Specific palette (genus level)
# palette1 <- brewer.pal(8, "Paired")
# palette2<- brewer.pal(8, "Set2")

combined_palette <- c(col_taxa_top10)
combined_palette <- combined_palette[1:16] #CHANGE INDEX

# Melt table
Top_tab_melted = melt(Top_df , value = Tab$Taxa)

# Reorder table and legend
Top_tab_melted <- transform(Top_tab_melted, Taxa = reorder(Taxa, -value))
colnames(Top_tab_melted)[1:2] = c("Taxa", "Sample") 

# Plot 

P <- ggplot(Top_tab_melted ,aes(x =Sample, y = value, fill = Taxa)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = combined_palette) + 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) + 
  labs (title = '') + 
  ylab("Relative abundance") +
  xlab("Samples")+
  theme(axis.text.x = element_blank(),       # Remove samples name       
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 7, family = 'Sans'),  
        legend.title = element_blank())
P

# Save plot
ggsave("Top10_Species_all_samples.png", P, width = 12, height = 7, dpi = 300)

# Add Unclassified to total taxa table
Tot <- rbind (Tot_ab.rel, U)
Tot <- data.frame(rownames(Tot), Tot)
colnames(Tot)[1] <-"Taxa"
Tot$Taxa[Tot$Taxa == "U"] <- "f_uncl.g_uncl.sp_uncl"  # Change 
head(U) [1:5]
rownames(Tot) <- Tot$Taxa
head(Tot)[1:5]
colSums(Tot[,-1])
dim(Tot)

# Create a table with relative abundance for the total taxa required in the some scripts of folder '7.statistical_analysisi_pipeline'
# need: taxa in rownames (remove column 'taxa') and samples in column.

df <- Tot$Taxa <- NULL # remove column taxa
head(df)

# save output
# Change name file for genus level
write.table(df, "Relative_abundaces_species.csv", sep = "\t") # rownames = species, column = samples
write.table(df, "Relative_abundaces_species.txt", sep = "\t")

write.table(Tot, "Total_abundance_species.csv", sep = "\t") # row = species, column = samples
write.table(Tot, "Total_abundance_species.txt", sep = "\t")

write.table(Top, "Top10_rel_ab_species_others.csv", sep = "\t") # table with top 10 species (row) and column with samples
write.table(Top, "Top10_rel_ab_species_others.txt", sep = "\t")


