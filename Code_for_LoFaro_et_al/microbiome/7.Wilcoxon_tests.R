# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================

## Performs paired Wilcoxon signed-rank tests across phases
## for microbiome features (species or gene families) and beta diversity and alpha diversity ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(data.table)
library(dplyr)

# ------------------------------------------------------------
### INPUT DATA
# ------------------------------------------------------------
# Select input dataset:
abundanceFILE <- "species.MGS_FECAL_merged_abundance_table_clr.rds"
# abundanceFILE <- "genefam.MGS_FECAL_merged_abundance_table_clr.rds"

# Covariates for residualization
covars_for_residuals <- c("DNA_concentration","readsDepth_postQC","batch")

# Read dataset
inDF <- readRDS(abundanceFILE)

# ------------------------------------------------------------
### DEFINE FEATURE COLUMNS
# ------------------------------------------------------------
# Species OR gene families (switch pattern if needed)
colsFeatures <- colnames(inDF)[grep("s__", colnames(inDF))]
# colsFeatures <- colnames(inDF)[grep("genefam__", colnames(inDF))]

# ------------------------------------------------------------
### REMOVE TECHNICAL EFFECTS (RESIDUALIZATION)
# ------------------------------------------------------------
covariate_formula <- paste(covars_for_residuals, collapse = " + ")

for (col in colsFeatures) {

  fit <- lm(
    as.formula(paste(col, "~", covariate_formula)),
    data = inDF,
    na.action = "na.exclude"
  )

  inDF[[col]] <- resid(fit)
}

# ------------------------------------------------------------
### FUNCTION: PAIRED WILCOXON TEST
# ------------------------------------------------------------
run_wilcox_pair <- function(df, col, p1, p2) {

  tmp <- df %>%
    filter(phase %in% c(p1, p2)) %>%
    group_by(ID) %>%
    filter(n_distinct(phase) == 2 &
           all(!is.na(.data[[col]]))) %>%
    distinct(ID, phase, .keep_all = TRUE) %>%
    ungroup()

  x <- tmp[tmp$phase == p1, ][[col]]
  y <- tmp[tmp$phase == p2, ][[col]]

  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]
  y <- y[valid]

  if (length(unique(x)) > 1 && length(unique(y)) > 1) {

    test <- wilcox.test(
      x = x, y = y,
      paired = TRUE,
      exact = FALSE,
      conf.int = TRUE
    )

    return(data.table(
      comparison = paste0("phase_", p1, "_", p2),
      p.value = test$p.value,
      estimate = test$estimate,
      ci_low = test$conf.int[1],
      ci_high = test$conf.int[2],
      feature = col
    ))
  }

  return(NULL)
}

# ------------------------------------------------------------
### RUN TESTS FOR FEATURES
# ------------------------------------------------------------
phase_pairs <- list(
  c(1,2), c(1,3), c(1,4),
  c(2,3), c(2,4), c(3,4)
)

res_df <- rbindlist(
  lapply(colsFeatures, function(col) {
    rbindlist(
      lapply(phase_pairs, function(p) {
        run_wilcox_pair(inDF, col, p[1], p[2])
      }),
      fill = TRUE
    )
  }),
  fill = TRUE
)

# ------------------------------------------------------------
### MULTIPLE TESTING CORRECTION
# ------------------------------------------------------------
res_df$p.value.BH <- p.adjust(res_df$p.value, method = "BH")
res_df$p.value.Bf <- p.adjust(res_df$p.value, method = "bonferroni")

# ------------------------------------------------------------
### OUTPUT (FEATURE RESULTS)
# ------------------------------------------------------------
# write.csv(res_df, "Wilcoxon_results_features.csv", row.names = FALSE)

# ============================================================
### BETA DIVERSITY TEST (btdiv_mean)
# ============================================================
# Assumes 'btdiv_mean' is already present in inDF

res_df_beta <- rbindlist(
  lapply(phase_pairs, function(p) {
    run_wilcox_pair(inDF, "btdiv_mean", p[1], p[2])
  }),
  fill = TRUE
)

# Multiple testing correction
res_df_beta$p.value.BH <- p.adjust(res_df_beta$p.value, method = "BH")
res_df_beta$p.value.Bf <- p.adjust(res_df_beta$p.value, method = "bonferroni")

# ------------------------------------------------------------
### OUTPUT (BETA RESULTS)
# ------------------------------------------------------------
# write.csv(res_df_beta, "Wilcoxon_results_beta.csv", row.names = FALSE)


# ============================================================
### ALPHA DIVERSITY (SHANNON) – PAIRED WILCOXON TESTS
# ============================================================
# Performs paired Wilcoxon signed-rank tests across phases
# for Shannon alpha diversity

# ------------------------------------------------------------
### INITIALIZE RESULTS
# ------------------------------------------------------------
res_df <- data.table()

# ------------------------------------------------------------
### FUNCTION: PAIRED WILCOXON TEST
# ------------------------------------------------------------
run_wilcox_pair <- function(df, col, p1, p2) {

  tmp <- df %>%
    filter(phase %in% c(p1, p2)) %>%
    group_by(ID) %>%
    filter(n_distinct(phase) == 2 &
           all(!is.na(.data[[col]]))) %>%
    distinct(ID, phase, .keep_all = TRUE) %>%
    ungroup()

  x <- tmp[tmp$phase == p1, ][[col]]
  y <- tmp[tmp$phase == p2, ][[col]]

  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]
  y <- y[valid]

  # Only test if variability exists
  if (length(unique(x)) > 1 && length(unique(y)) > 1) {

    test <- wilcox.test(
      x = x, y = y,
      paired = TRUE,
      exact = FALSE,
      conf.int = TRUE
    )

    return(data.table(
      comparison = paste0("phase_", p1, "_", p2),
      p.value = test$p.value,
      estimate = test$estimate,
      ci_low = test$conf.int[1],
      ci_high = test$conf.int[2],
      feature = col
    ))
  }

  return(NULL)
}

# ------------------------------------------------------------
### RUN TESTS
# ------------------------------------------------------------
# Define all phase comparisons
phase_pairs <- list(
  c(1,2), c(1,3), c(1,4),
  c(2,3), c(2,4), c(3,4)
)

# Apply to Shannon diversity
res_df <- rbindlist(
  lapply(phase_pairs, function(p) {
    run_wilcox_pair(inDF, "shannon_div", p[1], p[2])
  }),
  fill = TRUE
)

# ------------------------------------------------------------
### MULTIPLE TESTING CORRECTION
# ------------------------------------------------------------
res_df$p.value.BH <- p.adjust(res_df$p.value, method = "BH")
res_df$p.value.Bf <- p.adjust(res_df$p.value, method = "bonferroni")

# ------------------------------------------------------------
### OUTPUT
# ------------------------------------------------------------
# write.csv(res_df, "Wilcoxon_results_alpha.csv", row.names = FALSE)
