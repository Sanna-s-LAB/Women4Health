# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================


## Performs linear mixed models to assess associations between microbiome features
## and phase, hormones, and covariates, accounting for repeated measures (ID) ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(lme4)
library(lmerTest)
library(dplyr)
library(data.table)

# ------------------------------------------------------------
### INPUT DATA
# ------------------------------------------------------------
# Select dataset
inputFILE <- "species.MGS_FECAL_merged_abundance_table_clr.rds"
# inputFILE <- "genefam.MGS_FECAL_merged_abundance_table_clr.rds"
# inputFILE <- "pathways.MGS_FECAL_merged_abundance_table_clr.rds"

inDF <- readRDS(inputFILE)

# ------------------------------------------------------------
### DEFINE FEATURE AND METADATA COLUMNS
# ------------------------------------------------------------
colsFeatures <- grep("s__", colnames(inDF), value = TRUE)

feat_df <- inDF[, colsFeatures, drop = FALSE]
meta_df <- inDF[, !colnames(inDF) %in% colsFeatures, drop = FALSE]

# ------------------------------------------------------------
### REMOVE TECHNICAL EFFECTS (RESIDUALIZATION)
# ------------------------------------------------------------
tech_covars <- c("DNA_concentration","readsDepth_postQC","batch")

resid_feat_df <- feat_df

for (feat in colsFeatures) {

  df_tmp <- data.frame(
    y = feat_df[[feat]],
    meta_df[, tech_covars, drop = FALSE],
    check.names = FALSE
  )

  df_tmp <- df_tmp[complete.cases(df_tmp), ]

  if (nrow(df_tmp) < 5) {
    resid_feat_df[[feat]] <- NA
    next
  }

  fit <- lm(
    as.formula(paste("y ~", paste(tech_covars, collapse = " + "))),
    data = df_tmp
  )

  resid_feat_df[[feat]] <- NA
  resid_feat_df[rownames(df_tmp), feat] <- resid(fit)
}

# Replace features with residualized values
feat_df <- resid_feat_df

# Align rownames
rownames(feat_df) <- inDF$sample_ID
rownames(meta_df) <- inDF$sample_ID

# ------------------------------------------------------------
### BASE MODEL: PHASE ASSOCIATIONS
# ------------------------------------------------------------
vars_base <- c(
  "phase",
  "center",
  "Bristol_stool_scale",
  "Age",
  "BMI",
  "Collection_feces_daytime"
)

results <- list()

for (feat in colnames(feat_df)) {

  df_temp <- data.frame(
    abundance = feat_df[[feat]],
    meta_df[, vars_base, drop = FALSE],
    ID = meta_df$ID
  )

  df_temp <- df_temp[complete.cases(df_temp), ]
  if (nrow(df_temp) < 5) next

  # Model formula
  fmla <- as.formula(
    paste0(
      "abundance ~ phase + ",
      paste(vars_base[vars_base != "phase"], collapse = " + "),
      " + (1 | ID)"
    )
  )

  fit <- tryCatch(
    lmer(fmla, data = df_temp, REML = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit)) next

  s <- summary(fit)

  coef_phase <- s$coefficients[
    grep("^phase", rownames(s$coefficients)),
    , drop = FALSE
  ]

  results[[feat]] <- data.frame(
    feature = feat,
    variable = rownames(coef_phase),
    estimate = coef_phase[, "Estimate"],
    std_error = coef_phase[, "Std. Error"],
    p_value = coef_phase[, "Pr(>|t|)"],
    N = nrow(df_temp)
  )
}

# Combine + adjust p-values
results_phase <- bind_rows(results)

if (nrow(results_phase) > 0) {
  results_phase$FDR <- p.adjust(results_phase$p_value, method = "BH")
}

# Output
write.table(
  results_phase,
  file = "phase_associations.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
### EXTENDED MODEL: ADD HORMONE OR DIET VARIABLE
# ============================================================
# Example: PROG (can replace with any variable)

hormone_var <- "PROG"

meta_df$phase <- as.numeric(as.character(meta_df$phase))

vars_ext <- c(
  "phase",
  hormone_var,
  "center",
  "Bristol_stool_scale",
  "Age",
  "BMI",
  "Collection_feces_daytime"
)

results <- list()

for (feat in colnames(feat_df)) {

  df_temp <- data.frame(
    abundance = feat_df[[feat]],
    meta_df[, vars_ext, drop = FALSE],
    ID = meta_df$ID
  )

  df_temp <- df_temp[complete.cases(df_temp), ]
  if (nrow(df_temp) < 5) next

  fmla <- as.formula(
    paste0(
      "abundance ~ phase + ", hormone_var, " + ",
      paste(vars_ext[!vars_ext %in% c("phase", hormone_var)], collapse = " + "),
      " + (1 | ID)"
    )
  )

  fit <- tryCatch(
    lmer(fmla, data = df_temp, REML = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit)) next

  s <- summary(fit)

  coef_interest <- s$coefficients[
    rownames(s$coefficients) %in% c(hormone_var),
    , drop = FALSE
  ]

  if (nrow(coef_interest) == 0) next

  results[[feat]] <- data.frame(
    feature = feat,
    variable = rownames(coef_interest),
    estimate = coef_interest[, "Estimate"],
    std_error = coef_interest[, "Std. Error"],
    p_value = coef_interest[, "Pr(>|t|)"],
    N = nrow(df_temp)
  )
}

# Combine results
results_ext <- bind_rows(results)

if (nrow(results_ext) > 0) {
  results_ext$FDR <- p.adjust(results_ext$p_value, method = "BH")
}

# Output
write.table(
  results_ext,
  file = paste0(hormone_var, "_associations.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
