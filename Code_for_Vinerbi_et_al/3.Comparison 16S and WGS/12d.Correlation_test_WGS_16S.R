# Rscript

# Script to calculate a correlation (Spearmean test) to observe the correlation between 16S data and WGS data. We select the top 10 common taxa 
# between dataset

# Author: Elena Vinerbi (elenavinerbi@cnr.it) and Fabio Chillotti
# Last update: 03/02/2025

# R version: R 4.4.1

# upload library
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(scales)
library(stringr)
library(gridExtra)

# path
setwd("/home/")
getwd()

# Import table (output previous script '11c.Comparison_taxa_16S_WGS.R')
# Relative abundance of 16S 
df_16S <- read.table("top10_others_Species_16S.txt", sep = "\t", header = TRUE)
head(df_16S)
dim(df_16S)
#colSums(df_16S[,-1])
df_16S <- df_16S[-11,]
dim(df_16S)

# Relative abundance of WGS
df_WGS <- read.table("top10_others_Species_WGS.txt", sep = "\t", header = TRUE)
head(df_WGS)
dim(df_WGS)
#colSums(df_WGS[,-1])
df_WGS <- df_WGS[-11,]
dim(df_WGS)

# Common species between 16S - WGS
sp_16S <- df_16S$Species
sp_16S
length(sp_16S)

sp_WGS <- df_WGS$Species
sp_WGS
length(sp_WGS)

# common taxa
common_taxa <- intersect (sp_WGS, sp_16S)
common_taxa
length(common_taxa)

taxa_16S <- subset(df_16S, Species %in% common_taxa)
taxa_16S <- taxa_16S %>%
  mutate(Species = paste0(Species, "_16S"))
head(taxa_16S)

taxa_WGS <- subset(df_WGS, Species %in% common_taxa)
taxa_WGS <- taxa_WGS %>%
  mutate(Species = paste0(Species, "_WGS"))
head(taxa_WGS)

df <- rbind(taxa_WGS, taxa_16S)
View(df)

df_1 <- t(df)
colnames(df_1) <- df$Species
df_1 <- df_1[-1,]
#class(df_1)
df_1 <- data.frame(rownames(df_1), df_1)
colnames(df_1)[1] <- "Code"
View(df_1)
dim(df_1)

# remove samples without abundant more to 0 in at least 3 samples

df_1$Gemella_asaccharolytica_WGS <- NULL # change name of taxa
df_1$Gemella_asaccharolytica_16S <- NULL
dim(df_1)

# Save output
write.table(df_1, "taxa_16S_WGS.txt" ,sep="\t", row.names = T)
write.table(df_1, "taxa_16S_WGS.csv" ,sep="\t", row.names = T)

# Trasform in numeric a column
#df_1[, 2:13] <- lapply(df_1[, 2:13], as.numeric)

# Correlation test
cols_16S <- grep("_16S$", names(df_1), value = TRUE)
cols_WGS <- grep("_WGS$", names(df_1), value = TRUE)

prefixes_16S <- sub("_16S$", "", cols_16S)
prefixes_WGS <- sub("_WGS$", "", cols_WGS)

# find a common prefix (common taxa)
common_prefixes <- intersect(prefixes_16S, prefixes_WGS)
common_prefixes

# loop for cor.test
results_list <- list()

for (prefix in common_prefixes) {
  col_16S <- paste0(prefix, "_16S")
  col_WGS <- paste0(prefix, "_WGS")
  
  test_result <- cor.test(df_1[[col_16S]], df_1[[col_WGS]], method = "spearman", exact = FALSE)
  
  results_list[[prefix]] <- data.frame(
    prefix = prefix,
    rho = test_result$estimate,
    p.value = test_result$p.value
  )
}

# save a results
results_df <- bind_rows(results_list)
colnames(results_df) <- c("taxa","rho", "p.value")
View(results_df)

write.table(results_df, "Correlation_16S_WGS.csv", sep = "\t")
write.table(results_df, "Correlation_16S_WGS.txt", sep = "\t")

# Plot
# loop for plot
for (prefix in common_prefixes) {
  col_16S <- paste0(prefix, "_16S")
  col_WGS <- paste0(prefix, "_WGS")
  
  plot_data <- df_1 %>%
    select(all_of(col_16S), all_of(col_WGS)) %>%
    rename(
      `16S` = all_of(col_16S),
      `WGS` = all_of(col_WGS))
  
  plot_corr <- ggplot(plot_data, aes(x = WGS, y = `16S`)) +
    geom_point(alpha = 0.6) +  # Imposta la trasparenza dei punti per migliorare la visibilità
    geom_smooth(method = 'lm', se = FALSE, color = "blue") +  # Aggiungi linea di regressione
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Linea x = y
    labs(
      title = paste("Correlation btween", prefix, "WGS e 16S"),
      x = paste(prefix, "WGS"),
      y = paste(prefix, "16S")
    ) +
    theme(
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18)
    ) 
  
  
  # Save output
  ggsave(filename = paste0("Comparison_WGS_16S/plot_corr_", prefix, ".png"), plot = plot_corr, width = 8, height = 6)
}

# Observe correlation between Gardnerella_vaginalis(WGS) and Bifidobacterium_vaginale (16S)

bifidobacterium_16S <- subset(df_16S, Species == "Bifidobacterium_vaginale")
bifidobacterium_16S

gardnerella_WGS <- subset(df_WGS, Species == "Gardnerella_vaginalis")
gardnerella_WGS

tab <- rbind(bifidobacterium_16S, gardnerella_WGS)
#View(tab)

tab <- t(tab)
tab <- tab[-1,]
tab <- data.frame(rownames(tab), tab)
colnames(tab)[1] <- "Code"
colnames(tab)[2] <- "Bifidobacterium_vaginale_16s"
colnames(tab)[3] <- "Gardnerella_vaginalis_WGS"
#View(tab)

tab[, 2:3] <- lapply(tab[, 2:3], as.numeric)

A <- cor.test(tab$Bifidobacterium_vaginale_16s, tab$Gardnerella_vaginalis_WGS, method = "spearman", exact = FALSE)
A

class(A)
A_df <- data.frame(
  correlation = A$estimate,    # correlation coefficient
  p_value = A$p.value)
View(A_df)

A_df$Taxa <- "Bifidobacterium/Gardnerella_vaginalis"
View(A_df)
colnames(A_df) <- c("rho", "p.value", "taxa")
A_df <- select(A_df, taxa, rho, p.value)
View(A_df)

Tab <- rbind (results_df, A_df)
View(Tab)

# Save output
write.table(Tab, "Correlation_16S_WGS.csv", sep = "\t")
write.table(Tab, "Correlation_16S_WGS.txt", sep = "\t")

plot_corr <- ggplot(tab, aes(x = Gardnerella_vaginalis_WGS, y = Bifidobacterium_vaginale_16s)) +
  geom_point(alpha = 0.6) + 
  geom_smooth(method = 'lm', se = FALSE, color = "blue") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
plot_corr

ggsave(file="Gardnerella_vaginalisWGS_Bifidobacterium_vaginalis16S.png", plot_corr, width = 8, height = 6)
