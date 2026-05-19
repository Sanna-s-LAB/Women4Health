# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================


## Visualizes compositional profiles of gene families across samples and phases
## using relative abundance stacked bar plots ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(paletteer)

# ------------------------------------------------------------
### INPUT DATA
# ------------------------------------------------------------
abundanceFILE <- "genefam.MGS_FECAL_merged_abundance_table.txt"

inDF <- fread(abundanceFILE)
inDF <- data.frame(inDF, check.names = FALSE)

# ------------------------------------------------------------
### SELECT FEATURES
# ------------------------------------------------------------
colsFeatures <- colnames(inDF)[grep("genefam__", colnames(inDF))]

inDF <- inDF[, c("sample_ID", "ID", "phase", colsFeatures)]

# Convert phase to factor/character
inDF$phase <- as.character(inDF$phase)

# Clean feature names
colnames(inDF) <- gsub("genefam__", "", colnames(inDF))
colnames(inDF) <- gsub("_", " ", colnames(inDF))

# ------------------------------------------------------------
### LONG FORMAT
# ------------------------------------------------------------
long_df <- inDF %>%
  pivot_longer(
    cols = -c(sample_ID, ID, phase),
    names_to = "Feature",
    values_to = "Abundance"
  )

# ------------------------------------------------------------
### RELATIVE ABUNDANCE
# ------------------------------------------------------------
long_df <- long_df %>%
  group_by(sample_ID) %>%
  mutate(
    RelAbundance = Abundance / sum(Abundance, na.rm = TRUE)
  ) %>%
  ungroup()

# ------------------------------------------------------------
### SELECT TOP FEATURES
# ------------------------------------------------------------
top_features <- long_df %>%
  group_by(Feature) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Feature)

# Collapse others
long_df <- long_df %>%
  mutate(Feature = ifelse(Feature %in% top_features, Feature, "Other")) %>%
  group_by(sample_ID, phase, Feature) %>%
  summarise(RelAbundance = sum(RelAbundance), .groups = "drop")

# ------------------------------------------------------------
### FACTOR ORDERING
# ------------------------------------------------------------
feature_order <- c(top_features, "Other")
long_df$Feature <- factor(long_df$Feature, levels = feature_order)

# Ensure sample ordering is preserved
long_df$sample_ID <- factor(long_df$sample_ID, levels = unique(long_df$sample_ID))

# ------------------------------------------------------------
### COLORS
# ------------------------------------------------------------
paletteToPlot <- paletteer_d("ggsci::planetexpress_futurama")[-5]
colors <- setNames(paletteToPlot[1:length(top_features)], top_features)
colors <- c(colors, Other = "grey70")

# ------------------------------------------------------------
### PLOT THEME
# ------------------------------------------------------------
tm <- theme_classic() +
  theme(
    axis.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 11),

    panel.grid.major.y = element_line(color = "grey85"),

    strip.text = element_text(size = 14, face = "bold"),

    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# ------------------------------------------------------------
### GENERATE PLOTS BY PHASE
# ------------------------------------------------------------
phases <- unique(long_df$phase)

# pdf("gene_family_composition.pdf", width = 12, height = 8)

for (ph in phases) {

  df_sub <- long_df %>% filter(phase == ph)

  p <- ggplot(df_sub, aes(x = sample_ID, y = RelAbundance, fill = Feature)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(
      title = "Gene family composition",
      subtitle = paste("Phase:", ph),
      x = "Sample",
      y = "Relative abundance"
    ) +
    tm

  print(p)
}

# dev.off()



