# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================


## Computes beta diversity from species-level abundances using CLR transformation,
## adjusts for technical covariates, and visualizes community differences ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(paletteer)
library(patchwork)

# ------------------------------------------------------------
### INPUT DATA
# ------------------------------------------------------------
# Input: species-level abundance table with metadata
abundanceFILE <- "species.MGS_FECAL_merged_abundance_table.rds"

# Read data
inDF <- readRDS(abundanceFILE)

# ------------------------------------------------------------
### DEFINE FEATURE COLUMNS
# ------------------------------------------------------------
colsFeatures <- colnames(inDF)[grep("s__", colnames(inDF))]

# ------------------------------------------------------------
### RELATIVE ABUNDANCE NORMALIZATION
# ------------------------------------------------------------
metabk <- inDF[, -..colsFeatures]

rel_abb <- inDF[, ..colsFeatures] /
           rowSums(inDF[, ..colsFeatures], na.rm = TRUE)

inDF <- cbind(metabk, rel_abb)

# ------------------------------------------------------------
### CLR TRANSFORMATION
# ------------------------------------------------------------
x <- as.matrix(inDF[, ..colsFeatures])

# Add pseudocount to avoid log(0)
x_pc <- x + 1e-05

# Compute CLR
log_x <- log(x_pc)
clr   <- log_x - rowMeans(log_x)

# ------------------------------------------------------------
### REMOVE TECHNICAL EFFECTS (LINEAR MODEL)
# ------------------------------------------------------------
tech_covars <- c("DNA_concentration", "readsDepth_postQC", "batch")

tech_df <- inDF[, ..tech_covars]
rownames(tech_df) <- inDF$sample_ID

# Compute residualized CLR matrix
clr_resid <- matrix(
  NA,
  nrow = nrow(clr),
  ncol = ncol(clr),
  dimnames = list(inDF$sample_ID, colnames(clr))
)

for (feat in colnames(clr)) {

  df_tmp <- data.frame(
    y = clr[, feat],
    tech_df,
    check.names = FALSE
  )

  df_tmp <- df_tmp[complete.cases(df_tmp), ]
  if (nrow(df_tmp) < 5) next

  fit <- lm(
    as.formula(paste("y ~", paste(tech_covars, collapse = " + "))),
    data = df_tmp
  )

  clr_resid[rownames(df_tmp), feat] <- resid(fit)
}

# ------------------------------------------------------------
### BETA DIVERSITY (EUCLIDEAN DISTANCE)
# ------------------------------------------------------------
btDIV <- vegan::vegdist(clr_resid, method = "euclidean")

# Convert to matrix + data.table
btDIVmat <- as.matrix(btDIV)
btDIVdf  <- data.table(btDIVmat)
colnames(btDIVdf) <- inDF$sample_ID

# ------------------------------------------------------------
### PCOA (DIMENSION REDUCTION)
# ------------------------------------------------------------
pcoaMat <- cmdscale(btDIV, eig = TRUE, k = 3)
pcoaDF  <- data.frame(pcoaMat$points)
colnames(pcoaDF) <- paste0("PCOA", 1:3)

# Variance explained
variance <- eigenvals(pcoaMat) / sum(eigenvals(pcoaMat))

# ------------------------------------------------------------
### LONG FORMAT BETA DISTANCE TABLE
# ------------------------------------------------------------
btDIVdf$sample1 <- inDF$sample_ID

bt_long <- melt(
  btDIVdf,
  id.vars = "sample1",
  variable.name = "sample2",
  value.name = "beta_distance"
)

# Add phase metadata
visit_lookup <- inDF[, c("sample_ID", "phase")]
names(visit_lookup)[2] <- "visit"

bt_long <- merge(bt_long, visit_lookup, by.x = "sample1", by.y = "sample_ID")
setnames(bt_long, "visit", "visit1")

bt_long <- merge(bt_long, visit_lookup, by.x = "sample2", by.y = "sample_ID")
setnames(bt_long, "visit", "visit2")

# Rename variables
setnames(bt_long,
  old = c("sample1","sample2","visit1","visit2"),
  new = c("query_sample","comparison_sample","query_visit","comparison_visit")
)

# ------------------------------------------------------------
### FILTER COMPARISONS
# ------------------------------------------------------------
# Extract individual IDs
bt_long$query_individual      <- sub("_[0-9]+$", "", bt_long$query_sample)
bt_long$comparison_individual <- sub("_[0-9]+$", "", bt_long$comparison_sample)

# Keep:
#   - same visit
#   - different samples
#   - different individuals
bt_filtered <- bt_long[
  query_sample != comparison_sample &
  query_visit  == comparison_visit &
  query_individual != comparison_individual
]

# ------------------------------------------------------------
### SUMMARY STATISTICS PER SAMPLE
# ------------------------------------------------------------
bt_stats <- bt_filtered[, .(
  btdiv_mean   = mean(beta_distance, na.rm = TRUE),
  btdiv_median = median(beta_distance, na.rm = TRUE)
), by = .(query_sample, query_visit)]

# ------------------------------------------------------------
### PLOTTING SETTINGS
# ------------------------------------------------------------
paletteToPlot <- paletteer_d("ggsci::planetexpress_futurama")[-5]

tm <- theme(
  plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
  axis.title = element_text(size = 15),
  axis.text  = element_text(size = 13),
  legend.position = "bottom",
  panel.grid = element_line(color = "gray86", linewidth = 0.01, linetype = "dashed"),
  panel.background = element_blank()
)

# ------------------------------------------------------------
### BETA DISTRIBUTION BY PHASE
# ------------------------------------------------------------
mean_beta <- mean(bt_filtered$beta_distance, na.rm = TRUE)

pbd <- ggplot(bt_filtered,
              aes(x = factor(query_visit), y = beta_distance, fill = query_visit)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(alpha = 0.4, outlier.alpha = FALSE) +
  scale_fill_manual(values = paletteToPlot) +
  geom_hline(yintercept = mean_beta, linetype = "dashed", color = "red") +
  labs(x = "Phase", y = "Beta distance") +
  tm

pbdhist <- ggplot(bt_filtered,
                  aes(x = beta_distance, colour = query_visit)) +
  geom_histogram(aes(y = after_stat(density)), fill = "white", alpha = 0.6) +
  geom_density(linewidth = 1.5) +
  scale_colour_manual(values = paletteToPlot) +
  labs(x = "Beta distance", y = "Density") +
  tm

# ------------------------------------------------------------
### OUTPUT BETA DISTRIBUTION
# ------------------------------------------------------------
print((pbd / pbdhist))

# ------------------------------------------------------------
### MERGE BETA STATS WITH PCOA
# ------------------------------------------------------------
pcoaDF$sample_ID <- inDF$sample_ID

inDF_plot <- merge(
  pcoaDF,
  bt_stats,
  by.x = "sample_ID",
  by.y = "query_sample"
)

colnames(inDF_plot)[5] <- "phase"

# ------------------------------------------------------------
### PCOA VISUALIZATION
# ------------------------------------------------------------
inDF_plot <- inDF_plot[!is.na(inDF_plot$phase), ]

centroids <- aggregate(
  inDF_plot[, c("PCOA1","PCOA2","PCOA3")],
  by = list(Phase = inDF_plot$phase),
  FUN = mean
)

# PCOA plots
p1 <- ggplot(inDF_plot, aes(PCOA1, PCOA2, col = phase)) +
  geom_point(alpha = 0.6) +
  geom_point(
    data = centroids,
    aes(PCOA1, PCOA2, fill = Phase),
    shape = 23,
    size = 5,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = paletteToPlot) +
  scale_fill_manual(values = paletteToPlot) +
  labs(
    x = paste0("PCOA1 (", round(variance[1]*100), "%)"),
    y = paste0("PCOA2 (", round(variance[2]*100), "%)")
  ) +
  tm

p2 <- ggplot(inDF_plot, aes(PCOA2, PCOA3, col = phase)) +
  geom_point(alpha = 0.6) +
  geom_point(
    data = centroids,
    aes(PCOA2, PCOA3, fill = Phase),
    shape = 23,
    size = 5,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = paletteToPlot) +
  scale_fill_manual(values = paletteToPlot) +
  labs(
    x = paste0("PCOA2 (", round(variance[2]*100), "%)"),
    y = paste0("PCOA3 (", round(variance[3]*100), "%)")
  ) +
  tm

# ------------------------------------------------------------
### OUTPUT PCOA
# ------------------------------------------------------------
print(p1 | p2)

