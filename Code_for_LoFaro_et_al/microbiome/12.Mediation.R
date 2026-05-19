# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================


## Performs mediation analysis to assess whether microbiome features mediate
## the relationship between exposure variables and outcomes ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(mediation)
library(lme4)
library(data.table)

# ------------------------------------------------------------
### INPUT DATA
# ------------------------------------------------------------
# Select dataset
abundanceFILE <- "species.MGS_FECAL_merged_abundance_table_clr.rds"
# abundanceFILE <- "pathways.MGS_FECAL_merged_abundance_table_clr.rds"

# Read data
inDF <- readRDS(abundanceFILE)

# ------------------------------------------------------------
### DEFINE VARIABLES
# ------------------------------------------------------------
# Exposure, mediator, outcome
exposure_var  <- "Cholesterol_mean3gg"
mediator_var  <- "PRL"
outcome_var   <- "s__Alistipes_shahii"

# Covariates
covars <- c(
  "phase",
  "readsDepth_postQC",
  "DNA_concentration",
  "batch",
  "center",
  "Bristol_stool_scale",
  "Age",
  "BMI",
  "Collection_feces_daytime"
)

covars_string <- paste(covars, collapse = " + ")

# ------------------------------------------------------------
### DATA CLEANING
# ------------------------------------------------------------
# Keep only complete cases for variables involved
vars_needed <- c(exposure_var, mediator_var, outcome_var, covars, "ID")

df <- inDF[complete.cases(inDF[, ..vars_needed]), ]

# ------------------------------------------------------------
### MEDIATOR MODEL
# ------------------------------------------------------------
# Model mediator as a function of exposure + covariates
med_formula <- as.formula(
  paste(mediator_var, "~", exposure_var, "+", covars_string, "+ (1 | ID)")
)

med.fit <- lmer(
  med_formula,
  data = df,
  REML = TRUE
)

# ------------------------------------------------------------
### OUTCOME MODEL
# ------------------------------------------------------------
# Model outcome as a function of mediator + exposure + covariates
out_formula <- as.formula(
  paste(outcome_var, "~", mediator_var, "+", exposure_var, "+", covars_string, "+ (1 | ID)")
)

out.fit <- lmer(
  out_formula,
  data = df,
  REML = TRUE
)

# ------------------------------------------------------------
### MEDIATION ANALYSIS
# ------------------------------------------------------------
med.out <- mediate(
  med.fit,
  out.fit,
  treat    = exposure_var,
  mediator = mediator_var,
  sims     = 10000,
  boot     = FALSE
)

# ------------------------------------------------------------
### FORMAT OUTPUT
# ------------------------------------------------------------
mediation_results <- data.frame(
  Outcome  = outcome_var,
  Mediator = mediator_var,
  Exposure = exposure_var,

  Effect_Type = c(
    "ACME (Indirect)",
    "ADE (Direct)",
    "Total Effect",
    "Proportion Mediated"
  ),

  Estimate = c(
    med.out$d0,
    med.out$z0,
    med.out$tau.coef,
    med.out$n0
  ),

  P_Value = c(
    med.out$d0.p,
    med.out$z0.p,
    med.out$tau.p,
    med.out$n0.p
  ),

  CI_Lower = c(
    med.out$d0.ci[1],
    med.out$z0.ci[1],
    med.out$tau.ci[1],
    med.out$n0.ci[1]
  ),

  CI_Upper = c(
    med.out$d0.ci[2],
    med.out$z0.ci[2],
    med.out$tau.ci[2],
    med.out$n0.ci[2]
  )
)

# ------------------------------------------------------------
### OUTPUT
# ------------------------------------------------------------
# write.csv(
#   mediation_results,
#   "mediation_results.csv",
#   row.names = FALSE
# )
