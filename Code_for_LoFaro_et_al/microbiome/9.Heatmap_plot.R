# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================

## Generates heatmaps from association results (e.g., LMM outputs),
## displaying effect sizes and highlighting significant features ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# ------------------------------------------------------------
### INPUT DATA
# ------------------------------------------------------------
# Select input file (LMM or Wilcoxon results)
inFILE <- "phase_associations.txt"

# Read table
inDF <- fread(inFILE)

# ------------------------------------------------------------
### PREPARE DATA
# ------------------------------------------------------------
# Keep only phase-related terms
inDFsub <- inDF[inDF$term %like% "^phase", ]

# Standardize columns
inDFsub$qval <- inDFsub$FDR
inDFsub$value <- inDFsub$term
inDFsub$metadata <- "phase"

# ------------------------------------------------------------
### FILTER SIGNIFICANT FEATURES
# ------------------------------------------------------------
# Select features passing threshold
feats_pass_qval <- inDFsub[qval <= 0.25]$feature

inDFsub <- inDFsub[inDFsub$feature %in% feats_pass_qval, ]

# Define significance flag
inDFsub$qval_025 <- inDFsub$qval <= 0.25

# ------------------------------------------------------------
### RESHAPE TO MATRIX FORMAT
# ------------------------------------------------------------
cols_keep <- colnames(inDFsub)[
  grep(paste(c("feature","value","metadata","coef"), collapse = "|"),
       colnames(inDFsub))
]

# Wide matrix (effect sizes)
mat_df <- pivot_wider(
  inDFsub[, ..cols_keep],
  names_from = feature,
  values_from = coef
) %>%
  data.frame(check.names = FALSE)

# Wide matrix (significance)
sig_df <- inDFsub %>%
  pivot_wider(
    id_cols = c(value, metadata),
    names_from = feature,
    values_from = qval_025
  ) %>%
  data.frame(check.names = FALSE)

# ------------------------------------------------------------
### CLEAN ROW LABELS
# ------------------------------------------------------------
if (all(mat_df$metadata == "phase")) {
  mat_df$metadata <- paste0(mat_df$metadata, "_", mat_df$value)
  sig_df$metadata <- paste0(sig_df$metadata, "_", sig_df$value)
}

# Remove helper column
mat_df$value <- NULL
sig_df$value <- NULL

# Store raw names
row_labels <- mat_df$metadata
col_labels <- setdiff(colnames(mat_df), "metadata")

# Align significance matrix
sig_df <- sig_df[match(row_labels, sig_df$metadata), ]
sig_df <- sig_df[, c("metadata", col_labels)]

# ------------------------------------------------------------
### BUILD MATRICES
# ------------------------------------------------------------
mat <- as.matrix(mat_df[, col_labels])
sig_matrix <- as.matrix(sig_df[, col_labels])

# Clean names for visualization
rownames(mat) <- gsub("_", " ", row_labels)
rownames(sig_matrix) <- gsub("_", " ", row_labels)

colnames(mat) <- gsub("_", " ", col_labels)
colnames(sig_matrix) <- gsub("_", " ", col_labels)

# Transpose (features as rows)
mat <- t(mat)
sig_matrix <- t(sig_matrix)

# ------------------------------------------------------------
### COLOR SCALE
# ------------------------------------------------------------
custom_pal <- c("#C71000FF", "white", "#008EA0FF")
col_fun <- colorRamp2(
  c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)),
  custom_pal
)

# ------------------------------------------------------------
### BUILD HEATMAP
# ------------------------------------------------------------
ht <- Heatmap(
  mat,
  name = "Effect coefficient",
  
  col = col_fun,
  
  column_title = "Species - Phase",
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    title_position = "topleft"
  ),
  
  width  = unit(ncol(mat) * 0.5, "cm"),
  height = unit(nrow(mat) * 0.5, "cm"),
  
  # Add significance markers
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(sig_matrix[i, j]) && sig_matrix[i, j]) {
      grid::grid.text("*", x, y,
        gp = gpar(fontsize = 14, col = "black")
      )
    }
  }
)

# ------------------------------------------------------------
### OUTPUT
# ------------------------------------------------------------
pdf("heatmap_phase.pdf", width = 14, height = 14)

draw(
  ht,
  heatmap_legend_side = "bottom"
)

dev.off()
