# Rscript

# Script to calculate a statistics test for beta diversity 
# Note: the following tests were used: 
# Friedman test  = for both untransformed and transformed data (post CLR and adjustment for the tech variable)
# Wilcoxon test unpaired = or both untransformed and transformed data
# Wilcoxon test paired = for both untrasformed and transformed data
# t-test paired  = only for transformed data
# for this script a complete metadata file, with all the sample IDs even when no data were present was used

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

# path
setwd("/home/")
getwd()

# Import table
# Metadata file
# example of metadata columns: Plate; ID; Code_women; Visits; ID_paper; Var1; Var2
Metadata <- read.table("Metadata_complete.txt", header = TRUE, sep = "\t" )
head(Metadata)

# distribution of beta diversity (output previous script '10a.beta_diversity.R')
df <- read.table("Distribution_Species_Beta.txt", header = TRUE, sep = "\t" )
head(df)

# calculate a Wilcoxon test unpaired
WT <- df %>%
  wilcox_test(beta ~ Weeks, paired = F, p.adjust.method = "bonferroni" )
WT

# Save results in dataframe
df_WT <- as.data.frame(WT)
head(df_WT)

# Save output
write.table(df_WT, "WT_unpaired_beta.csv", sep = "\t")
write.table(df_WT, "WT_unpaired_beta.txt", sep = "\t")

# Import table
# Beta diversity(Average and median) table (output previous script '10a.beta_diversity.R')
tab <- read.table("Average_median_Beta.txt", header=TRUE, sep = "\t", stringsAsFactors = TRUE)
tab <- data.frame(rownames(tab), tab)
colnames(tab)[1] <- "Code_women"
head(tab)

# Statistics test for average and median (beta diversity)
metadata_BC <- merge (Metadati, tab, by ="Code_women", all = TRUE)
head(metadata_BC)

# These analyses were first carried out for the average and then the median values 
# Create a table
Beta <- data.frame(metadata_BC$Visit, metadata_BC$Average, metadata_BC$Group) # Change median/average
colnames(Beta)<-c("Visit","Average","Group") # Change median/average
head(Beta)
ord <- c("First","Second","Third","Fourth")
Beta$Visit <- factor(Beta$Visit, levels = ord)
Beta$Average <- as.numeric(Beta$Average) # Change median/average
head(Beta)

# Friedman test (overall)
FT_test <- friedman.test(Beta$Average, groups = Beta$Visit, blocks = Beta$Group) # Change median/average
FT_test
text_FT <- paste("Friedman Test", "p =", round(FT_test$p.value, 6))

# save results in a table
df_FT <- data.frame(
  statistic = FT_test$statistic,
  p_value = FT_test$p.value,
  parameter = FT_test$parameter
)

head(df_FT)

# Save output
write.table(df_FT, "FT_beta.csv", sep = "\t")
write.table(df_FT, "FT_beta.txt", sep = "\t")

# t-test pairwaise (paired)
TT_test <- pairwise.t.test(Beta$Average, Beta$Visit, paired = TRUE, p.adjust.method = "bonferroni") # Change median/average
TT_df <- TT_test$p.value
df_TT <- as.data.frame(as.table(TT_df))
colnames(df_TT) <- c("Group1", "Group2", "p.val.adj")
df_TT$sign <- ifelse(df_TT$p.val.adj < 0.05, "*", "ns") # Ad symbol for significant value
head(df_TT)
 
# Remove Na
df_TT <- na.omit(df_TT)
head(df_TT)

# Merge table
M <- select(df_TT, p.val.adj, sign)
head(M)

Tab_TT <- cbind(df_TT, M)
colnames(Tab_TT) <- c("Group1", "Group2", "Pvalue", "Pvalue.adj", "sign")
head(Tab_TT)

# Save output
write.table(Tab_TT, "TT_beta.csv", sep = "\t")
write.table(Tab_TT, "TT_beta.txt", sep = "\t")

# Calculate Wilcoxon test (paired)
WT  <- Beta %>%
  wilcox_test(Average ~ Visit, paired=T, p.adjust.method = "bonferroni" ) # Change median/average
head(WT)
WT$p.adj=round(WT$p.adj,6)

df_WT <- as.data.frame(WT)
head(df_WT)

# Save output
write.table(df_WT, "WT_beta.csv", sep = "\t")
write.table(df_WT, "WT_beta.txt", sep = "\t")


# Violin plot 
# Specific palette for the four phases
col <- c( "First" = "#FFE066","Second" = "#FFF5CC", "Third" = "#FFC54C", "Fourth" = "#FFA32B") 
col

Plot <- Beta %>%
  ggplot() +
  theme_bw() +
  geom_violin(aes(x = factor(Visit), y = Average, fill=factor(Visit))) + # Change median/average
  geom_boxplot(aes(x = factor(Visit), y = Average), width=0.1, alpha=0.6) + # Change median/average
  scale_x_discrete(labels = c("F", "O", "EL", "LL")) +
  labs(title = "", x = "", y = "Average(Beta diversity)") + # Change median/average
  scale_fill_manual(name = "Weeks", values = col, labels = c("F", "O", "EL", "LL"))+
  theme(axis.title = element_text(size = 7, family = 'sans'),
    axis.text.x = element_text(size = 7, family = 'sans'),
    axis.text.y = element_text(size = 7, family = 'sans'),),
    legend.position = 'none') +
  geom_segment(aes(x = 1, xend = 2, y = 80, yend = 80), color = "black", size = 0.5) + # line with significant output of statistical test
  geom_text(aes(x = 1.5, y = 82, label = "*"), size = 5) +
geom_segment(aes(x = 1, xend = 3, y = 75, yend = 75), color = "black", size = 0.5) + # line with significant output of statistical test
  geom_text(aes(x = 2.5, y = 77, label = "*"), size = 5)
Plot

# Save Plot
ggsave("Plot_beta_significant.png", Plot, width = 15, height = 7, dpi = 300)

# Beta diversity between women
# Import table with same women
df_women <- read.table("Average_Beta_same_women_PostCLR.txt", header=TRUE, sep = "\t", stringsAsFactors = TRUE)
head(df_women)

# Create a loop to compare same women with different women

results_df <- data.frame(Iteration = integer(), WT_Average_p = numeric(), WT_Median_p = numeric(), stringsAsFactors = FALSE)

for (i in 1:10) {
  file_name <- paste0("Beta_diversity_random_women_rm_", i, ".txt") # import table random women
  
  # Import table with different women
  df_rm_women <- read.table(file_name, header=TRUE, sep = "", stringsAsFactors = TRUE)
  
  # Wilcoxon test dìfor average
  WT_average <- wilcox.test(df_women$Average, df_rm_women$Average, paired = FALSE, p.adjust.method = "bonferroni") # average

  # Wilcoxon test for median
  WT_Median <- wilcox.test(df_women$Median, df_rm_women$Median, paired = FALSE, p.adjust.method = "bonferroni") # women
  
  # Save results to table
  results_df <- rbind(results_df, data.frame(Iteration = i, WT_Average_p = WT_average$p.value, WT_Median_p = WT_Median$p.value))
}

head (results_df)

# Save output
write.table(results_df, "Wilcoxon_results_women.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(results_df, "Wilcoxon_results_women.csv", sep = "\t", row.names = FALSE, quote = FALSE)


