# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================

## Visualizes CLR-transformed gene family abundances across phases
## using violin + boxplots ##

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
abundanceFILE <- "genefam.MGS_FECAL_merged_abundance_table_clr.rds"

inDF <- readRDS(abundanceFILE)
inDF <- as.data.table(inDF)

# ------------------------------------------------------------
### DEFINE FEATURE COLUMNS
# ------------------------------------------------------------
colsFeatures <- colnames(inDF)[grep("genefam__", colnames(inDF))]

# ------------------------------------------------------------
### LONG FORMAT
# ------------------------------------------------------------
df_long <- inDF %>%
  pivot_longer(
    cols = all_of(colsFeatures),
    names_to = "Feature",
    values_to = "Abundance"
  )

# Remove missing phase
df_long <- df_long[!is.na(df_long$phase), ]

# Clean feature names
df_long$Feature <- gsub("genefam__", "", df_long$Feature)
df_long$Feature <- gsub("_", " ", df_long$Feature)

# Ensure phase is factor
df_long$phase <- factor(df_long$phase)

# ------------------------------------------------------------
### PLOTTING SETTINGS
# ------------------------------------------------------------
paletteToPlot <- paletteer_d("ggsci::planetexpress_futurama")[-5]

fill_scale <- scale_fill_manual(
  values = paletteToPlot,
  name = "Phase"
)

tm <- theme(
  plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
  axis.title = element_text(size = 15),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 13),
  legend.position = "bottom",
  panel.grid = element_line(color = "gray86", linewidth = 0.01, linetype = "dashed"),
  panel.background = element_blank()
)

# ------------------------------------------------------------
### BOXPLOT + VIOLIN
# ------------------------------------------------------------
pv <- ggplot(df_long, aes(x = phase, y = Abundance, fill = phase)) +

  geom_violin(
    colour = "gray62",
    alpha  = 0.4,
    trim   = TRUE
  ) +

  geom_boxplot(
    width = 0.2,
    colour = "gray86",
    alpha = 0.4,
    outlier.alpha = FALSE
  ) +

  fill_scale +

  facet_wrap(~ Feature, scales = "free_y") +

  labs(
    title = "Gene family abundance across phases",
    x = "Phase",
    y = "CLR abundance"
  ) +

  tm

# ------------------------------------------------------------
### OUTPUT
# ------------------------------------------------------------
# pdf("genefam_boxplots.pdf", width = 12, height = 10)
print(pv)
# dev.off()
