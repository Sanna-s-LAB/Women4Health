# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================

## Calculates Shannon alpha diversity (species level) and visualizes its distribution across phases ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(data.table)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(paletteer)
library(patchwork)

# ------------------------------------------------------------
### INPUT DATA
# ------------------------------------------------------------
# Input: species-level abundance table (merged with metadata)
abundanceFILE <- "species.MGS_FECAL_merged_abundance_table.rds"

# Read data
inDF <- readRDS(abundanceFILE)

# ------------------------------------------------------------
### DEFINE FEATURE COLUMNS
# ------------------------------------------------------------
# Identify species-level features
colsFeatures <- colnames(inDF)[grep("s__", colnames(inDF))]

# ------------------------------------------------------------
### COMPUTE ALPHA DIVERSITY
# ------------------------------------------------------------
# Calculate Shannon diversity index for each sample
inDF$shannon_div <- vegan::diversity(
  inDF[, ..colsFeatures],
  index = "shannon"
)

# ------------------------------------------------------------
### DATA PREPARATION
# ------------------------------------------------------------
# Remove samples without phase annotation
inDF <- inDF[!is.na(inDF$phase), ]

# Convert phase to character (for plotting consistency)
inDF$phase <- as.character(inDF$phase)

# ------------------------------------------------------------
### PLOTTING SETTINGS
# ------------------------------------------------------------
# Color palette
paletteToPlot <- paletteer_d("ggsci::planetexpress_futurama")[-5]

# Aesthetic settings for points
size_val  <- 1
shape_val <- 21
color_val <- "gray42"

# Annotation positions
xannotPOS <- sd(inDF$shannon_div) + 0.03
yannotPOS <- mean(inDF$shannon_div, na.rm = TRUE) + 0.04

# Common theme
tm <- theme(
  axis.title   = element_text(size = 18),
  axis.text.x  = element_text(size = 13),
  axis.text.y  = element_text(size = 13),
  legend.text  = element_text(size = 13),
  legend.title = element_text(size = 14),
  panel.grid   = element_line(color = "gray86", linewidth = 0.01, linetype = "dashed"),
  panel.background = element_blank()
)

# ------------------------------------------------------------
### BOXPLOT (SHANNON BY PHASE)
# ------------------------------------------------------------
pb <- ggplot(inDF, aes(x = phase, y = shannon_div, fill = phase)) +
  geom_boxplot(
    colour = "gray66",
    width  = 0.85,
    alpha  = 0.3,
    outlier.alpha  = FALSE,
    outlier.colour = "black",
    outlier.fill   = "firebrick2",
    outlier.shape  = 23,
    outlier.size   = 2,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = paletteToPlot, name = "Phase") +
  geom_jitter(
    shape = shape_val,
    fill  = color_val,
    color = color_val,
    size  = size_val,
    height = 0,
    alpha  = 0.6
  ) +
  geom_hline(
    yintercept = mean(inDF$shannon_div, na.rm = TRUE),
    color = "red", linetype = "dashed", linewidth = 0.4
  ) +
  annotate(
    "text",
    x = -xannotPOS,
    y = yannotPOS,
    label = paste0("Mean ", round(mean(inDF$shannon_div, na.rm = TRUE), 3)),
    size = 4.5,
    color = "red"
  ) +
  labs(
    title = "Species",
    x = "Phase",
    y = "Shannon diversity"
  ) +
  coord_cartesian(clip = "off") +
  tm

# ------------------------------------------------------------
### VIOLIN + BOXPLOT
# ------------------------------------------------------------
pv <- ggplot(inDF, aes(x = phase, y = shannon_div, fill = phase)) +
  geom_violin(
    colour = "gray62",
    alpha  = 0.4
  ) +
  geom_boxplot(
    colour = "gray86",
    alpha  = 0.4,
    outlier.alpha  = FALSE,
    outlier.colour = "black",
    outlier.fill   = "firebrick2",
    outlier.shape  = 23,
    outlier.size   = 2
  ) +
  scale_fill_manual(values = paletteToPlot, name = "Phase") +
  geom_jitter(
    shape = shape_val,
    fill  = color_val,
    color = color_val,
    size  = size_val,
    height = 0,
    alpha  = 0.6
  ) +
  geom_hline(
    yintercept = mean(inDF$shannon_div, na.rm = TRUE),
    color = "red", linetype = "dashed", linewidth = 0.4
  ) +
  annotate(
    "text",
    x = -xannotPOS,
    y = yannotPOS,
    label = paste0("Mean ", round(mean(inDF$shannon_div, na.rm = TRUE), 3)),
    size = 4.5,
    color = "red"
  ) +
  labs(
    x = "Phase",
    y = "Shannon diversity"
  ) +
  coord_cartesian(clip = "off") +
  tm

# ------------------------------------------------------------
### DISTRIBUTION PLOTS (HISTOGRAM + DENSITY)
# ------------------------------------------------------------
phist <- ggplot(inDF, aes(x = shannon_div, colour = phase)) +
  geom_histogram(
    aes(y = after_stat(density)),
    fill = "white",
    alpha = 0.6
  ) +
  geom_density(
    linewidth = 1.8,
    adjust = 1.5,
    alpha = 0.3
  ) +
  scale_colour_manual(values = paletteToPlot, name = "Phase") +
  geom_vline(
    xintercept = mean(inDF$shannon_div, na.rm = TRUE),
    linetype = "dashed",
    color = "grey38",
    linewidth = 0.8
  ) +
  labs(
    x = "Alpha diversity (Shannon)",
    y = "Density"
  ) +
  coord_cartesian(clip = "off") +
  tm

# ------------------------------------------------------------
### OUTPUT PLOTS
# ------------------------------------------------------------
# Combine plots (patchwork)
print((pb | pv) / phist + plot_annotation(tag_levels = c('1','A')))

# Optional export
# pdf("Alpha_Diversity_phase.pdf", width = 14, height = 14)
# print((pb | pv) / phist)
# dev.off()
