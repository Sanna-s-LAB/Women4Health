# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================

## Generates forest plots of effect coefficients with 95% confidence intervals
## for feature–phase associations across multiple models ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(data.table)
library(ggplot2)
library(paletteer)

# ------------------------------------------------------------
### INPUT DATA
# ------------------------------------------------------------
# Select input file (association results)
inputFILE <- "phase_associations.txt"

plot_dt <- fread(inputFILE)

# ------------------------------------------------------------
### DATA PREPARATION
# ------------------------------------------------------------
# Ensure model factor ordering
plot_dt[, model := factor(
  model,
  levels = c("Base", "Nutrient three days", "Food habits")
)]

# Compute confidence intervals (95%)
plot_dt[, ci_low  := coef - 1.96 * stderr]
plot_dt[, ci_high := coef + 1.96 * stderr]

# Remove incomplete rows
plot_dt <- plot_dt[!is.na(coef) & !is.na(stderr)]

# ------------------------------------------------------------
### FEATURE LABELING AND ORDERING
# ------------------------------------------------------------
# Create combined feature label
plot_dt[, feature_label := paste(feature, metadata, value, sep = "_")]

# Order features based on Base model coefficients
base_order <- plot_dt[model == "Base"][order(coef)]$feature_label
plot_dt[, feature_label := factor(feature_label, levels = base_order)]

# Clean feature names for plotting
plot_dt[, feature_clean := gsub("_", " ", feature)]

# ------------------------------------------------------------
### FACET PREPARATION
# ------------------------------------------------------------
# Define facet variable
plot_dt[, facet_label := as.factor(metadata)]

# Order facets by frequency
facet_order <- names(sort(table(plot_dt$facet_label)))
plot_dt[, facet_label := factor(facet_label, levels = facet_order)]

# ------------------------------------------------------------
### PLOTTING SETTINGS
# ------------------------------------------------------------
palette_to_plot <- paletteer_d("ggsci::planetexpress_futurama")[-5]

# ------------------------------------------------------------
### FOREST PLOT
# ------------------------------------------------------------
p <- ggplot(plot_dt,
            aes(x = coef, y = feature_label, color = model)) +

  # Reference line at 0
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "grey50") +

  # Points
  geom_point(
    position = position_dodge(width = 0.7),
    size = 2.5
  ) +

  # Confidence intervals
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    width = 0.05,
    position = position_dodge(width = 0.7)
  ) +

  # Colors
  scale_color_manual(values = palette_to_plot) +

  # Axes
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +

  # Faceting
  facet_wrap(~ facet_label, scales = "free_y") +

  # Labels
  labs(
    x = "Effect coefficient",
    y = "Feature",
    color = "Model"
  ) +

  # Theme
  theme_bw() +
  theme(
    panel.grid.minor = element_line(color = "grey94", linetype = "dotted"),
    panel.grid.major = element_line(color = "grey92", linetype = "dotted"),
    panel.border = element_blank(),

    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.title  = element_text(size = 13),

    strip.text = element_text(size = 12, face = "bold"),

    legend.position = "top"
  )

# ------------------------------------------------------------
### OUTPUT
# ------------------------------------------------------------
# Save plot
# pdf("forest_plot.pdf", width = 12, height = 14)
print(p)
# dev.off()
