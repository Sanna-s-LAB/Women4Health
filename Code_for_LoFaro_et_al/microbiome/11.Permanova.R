# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================

## Performs beta diversity analysis using CLR-transformed species data
## and tests associations with metadata using PERMANOVA (adonis2) ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(data.table)
library(vegan)

# ------------------------------------------------------------
### INPUT DATA
# ------------------------------------------------------------
abundanceFILE <- "species.MGS_FECAL_merged_abundance_table.rds"

# Read dataset
inDF <- readRDS(abundanceFILE)

# ------------------------------------------------------------
### FUNCTION: CLR + RESIDUALIZATION + ADONIS
# ------------------------------------------------------------
runAdonis <- function(
  data,
  adonis_vars,
  max_perm = 10000,
  nrThreads = 4,
  verbose = TRUE
) {

  data <- data.frame(data)

  # ----------------------------------------------------------
  ### DEFINE FEATURES
  # ----------------------------------------------------------
  colsFeatures <- colnames(data)[grep("s__", colnames(data))]

  # Technical covariates to remove
  tech_covars <- c("DNA_concentration", "readsDepth_postQC", "batch")

  # ----------------------------------------------------------
  ### REMOVE TECHNICAL EFFECTS
  # ----------------------------------------------------------
  if (verbose) message("Residualizing technical covariates...")

  resid_df <- data[, colsFeatures]

  for (feat in colsFeatures) {

    df_tmp <- data.frame(
      y = data[[feat]],
      data[, tech_covars, drop = FALSE],
      check.names = FALSE
    )

    df_tmp <- df_tmp[complete.cases(df_tmp), ]

    if (nrow(df_tmp) < 5) {
      resid_df[[feat]] <- NA
      next
    }

    fit <- lm(
      as.formula(paste("y ~", paste(tech_covars, collapse = " + "))),
      data = df_tmp
    )

    resid_df[[feat]] <- NA
    resid_df[rownames(df_tmp), feat] <- resid(fit)
  }

  data[, colsFeatures] <- resid_df

  # ----------------------------------------------------------
  ### FILTER COMPLETE CASES
  # ----------------------------------------------------------
  cols_to_check <- c(colsFeatures, adonis_vars, "ID")
  data <- data[complete.cases(data[, cols_to_check]), ]

  # ----------------------------------------------------------
  ### CLR TRANSFORMATION
  # ----------------------------------------------------------
  if (verbose) message("Performing CLR transformation...")

  x <- as.matrix(data[, colsFeatures])

  # Add pseudocount
  x <- x + 1e-05

  log_x <- log(x)
  clr   <- sweep(log_x, 1, rowMeans(log_x), "-")

  # Beta diversity (Euclidean on CLR)
  div <- vegan::vegdist(clr, method = "euclidean")

  # ----------------------------------------------------------
  ### RUN PERMANOVA
  # ----------------------------------------------------------
  if (verbose) message("Running PERMANOVA (adonis2)...")

  formula <- as.formula(
    paste("div ~", paste(adonis_vars, collapse = " + "))
  )

  aov_table <- adonis2(
    formula,
    data = data,
    strata = data$ID,
    permutations = max_perm,
    by = "margin",
    parallel = nrThreads
  )

  # ----------------------------------------------------------
  ### FORMAT OUTPUT
  # ----------------------------------------------------------
  results <- data.frame(
    variable   = rownames(aov_table),
    DF         = aov_table[,1],
    SumsOfSqs  = aov_table[,2],
    R2         = aov_table[,3],
    F          = aov_table[,4],
    p_value    = aov_table[,5]
  )

  if (verbose) message("PERMANOVA completed.")

  return(results)
}

# ============================================================
### UNIVARIABLE PERMANOVA
# ============================================================
vars_univ <- c(
  "phase","center","Bristol_stool_scale",
  "Age","BMI","Collection_feces_daytime"
)

out_univ <- runAdonis(
  data = inDF,
  adonis_vars = vars_univ,
  max_perm = 10000,
  nrThreads = 4
)

# write.csv(out_univ, "PERMANOVA_univariable.csv", row.names = FALSE)

# ============================================================
### MULTIVARIABLE PERMANOVA
# ============================================================
vars_multi <- c(
  "Wine","Fish","Olive_oil","Butter","Fiber",
  "PROG","LH","FSH",
  "center","Bristol_stool_scale",
  "Age","BMI","Collection_feces_daytime"
)

out_multi <- runAdonis(
  data = inDF,
  adonis_vars = vars_multi,
  max_perm = 10000,
  nrThreads = 4
)

# write.csv(out_multi, "PERMANOVA_multivariable.csv", row.names = FALSE)

