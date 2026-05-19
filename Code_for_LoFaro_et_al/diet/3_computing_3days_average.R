# ANDREA CARTA 15/04/2026

# =============================================================================
# Script: computing_3days_average.R
# Purpose: Compute 3-day average dietary intake per participant per week,
#          with outlier removal based on total daily caloric intake thresholds.
# Input:   A longitudinal dietary diary dataset (one row per meal per day)
# Output:  A data frame with one row per participant per week,
#          containing mean caloric and macronutrient intake over the last 3 days
# =============================================================================

library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

# Load the dietary diary dataset.
# The object must contain a data frame with at least the following columns:
#   ID, date, day, meal_type, collection_week,
#   calories, carbohydrates, fat, proteins, cholesterol,
#   sodium, sugars, fiber,
#   menstruation_day (logical), enrollment_date, BMI, height_cm, weight_kg

load("your_dataset.rda")   # <-- replace with your actual file path

# Rename to a generic working object
dietary_data <- mega_Final_TriesteCon_outlier_totali   # <-- replace with the actual object name


# -----------------------------------------------------------------------------
# 2. RENAME COLUMNS FROM SOURCE LANGUAGE TO ENGLISH
# -----------------------------------------------------------------------------
# If your dataset uses non-English column names, rename them here before
# proceeding. Adjust the left-hand side values to match your actual column names.

dietary_data <- dietary_data %>%
  rename(
    meal_type        = Pasti,
    collection_week  = weeks_raccolta,
    date             = Dates,
    enrollment_date  = data_arruolamento,
    menstruation_day = giorno_mestruazione,
    height_cm        = altezza_cm,
    weight_kg        = peso_kg,
    calories         = Calorie,
    carbohydrates    = Carboidrati,
    fat              = Grassi,
    proteins         = Proteine,
    cholesterol      = Colest,
    sodium           = Sodio,
    sugars           = Zuccheri,
    fiber            = Fibre
  ) %>%
  mutate(
    # Standardise meal type labels to English
    meal_type = recode(meal_type,
                       "Totale"          = "Total",
                       "Pranzo"          = "Lunch",
                       "Cena"            = "Dinner",
                       "Prima colazione" = "Breakfast",
                       "Colazione"       = "Breakfast",
                       "Snack"           = "Snack")
  )


# -----------------------------------------------------------------------------
# 3. SPLIT DATA BY PARTICIPANT
# -----------------------------------------------------------------------------

# Split the full dataset into a list, one element per participant (ID)
diary_list <- split(dietary_data, dietary_data$ID)


# -----------------------------------------------------------------------------
# 4. IDENTIFY OUTLIER DAYS BASED ON TOTAL CALORIC INTAKE
# -----------------------------------------------------------------------------

# Define caloric thresholds for outlier detection.
# Days with total caloric intake outside [CALORIE_MIN, CALORIE_MAX] are excluded.
CALORIE_MIN <- 500
CALORIE_MAX <- 2500

# Macronutrient variables to be set to NA when a day is flagged as an outlier
MACRO_VARS <- c("carbohydrates", "fat", "proteins", "cholesterol",
                "sodium", "sugars", "fiber")

# For each participant, flag days where the "Total" row (daily total) falls
# outside the caloric thresholds. When a day is flagged, both the "Total" row
# and all individual meal rows for that day (up to 4 meal rows before "Total")
# are set to NA for calories and all macronutrients.

diary_list_clean <- diary_list

for (i in seq_along(diary_list_clean)) {

  df_i <- diary_list_clean[[i]]

  for (j in seq_len(nrow(df_i))) {

    # Check only "Total" rows (daily totals)
    if (df_i$meal_type[j] == "Total" &&
        (is.na(df_i$calories[j]) ||
         df_i$calories[j] < CALORIE_MIN ||
         df_i$calories[j] > CALORIE_MAX)) {

      # Set to NA: the "Total" row and the preceding meal rows (up to 4)
      rows_to_na <- j - (0:4)
      rows_to_na <- rows_to_na[rows_to_na >= 1]  # avoid index out of bounds

      df_i[rows_to_na, c("calories", MACRO_VARS)] <- NA
    }
  }

  diary_list_clean[[i]] <- df_i
}

# Preserve participant names in the cleaned list
names(diary_list_clean) <- names(diary_list)


# -----------------------------------------------------------------------------
# 5. FUNCTION: compute_3day_means()
# -----------------------------------------------------------------------------
# For each participant, compute the mean dietary intake over the last 3 days
# of each collection week. This approach is used to obtain a representative
# estimate of habitual intake while limiting within-week variability.
#
# Arguments:
#   df_input  - a data frame for a single participant (one element of diary_list_clean)
#
# Returns:
#   A tibble with one row per expected week (first, second, third, fourth),
#   containing mean caloric intake by meal type and mean macronutrient intake.

compute_3day_means <- function(df_input) {

  WEEKS_EXPECTED <- c("first", "second", "third", "fourth")
  MEAL_TYPES     <- c("Total", "Lunch", "Dinner", "Breakfast", "Snack")
  MACRO_VARS     <- c("carbohydrates", "fat", "proteins", "cholesterol",
                      "sodium", "sugars", "fiber")

  # Prepare data: coerce numeric columns, parse dates, standardise factor levels
  df <- df_input %>%
    mutate(
      across(c("calories", all_of(MACRO_VARS)), ~ as.numeric(as.character(.))),
      date            = as.Date(date),
      collection_week = factor(collection_week,
                               levels  = WEEKS_EXPECTED,
                               ordered = TRUE)
    )

  # Extract participant-level demographic info (constant across rows)
  info <- df[1, c("ID", "enrollment_date", "BMI", "height_cm", "weight_kg")]

  # ------------------------------------------------------------------
  # Inner function: compute means for a single week
  # ------------------------------------------------------------------
  compute_week <- function(week_name) {

    block <- df %>% filter(collection_week == week_name)

    # Return an all-NA row if no data are available for this week
    if (nrow(block) == 0) {
      return(tibble(
        collection_week   = week_name,
        collection_date   = as.Date(NA),
        ID                = info$ID,
        enrollment_date   = info$enrollment_date,
        menstruation_date = as.Date(NA),
        BMI               = info$BMI,
        height_cm         = info$height_cm,
        weight_kg         = info$weight_kg,
        !!!setNames(rep(NA_real_, length(MEAL_TYPES)),
                    paste0("calories_", tolower(MEAL_TYPES))),
        !!!setNames(rep(NA_real_, length(MACRO_VARS)), MACRO_VARS)
      ))
    }

    # Sort by day descending and take the last 3 recorded days
    block <- block %>% arrange(desc(Day))

    # Mean caloric intake by meal type over the last 3 days
    calorie_means <- sapply(MEAL_TYPES, function(p) {
      block %>%
        filter(meal_type == p) %>%
        head(3) %>%
        pull(calories) %>%
        mean(na.rm = TRUE)
    })
    names(calorie_means) <- paste0("calories_", tolower(MEAL_TYPES))

    # Mean macronutrient intake from "Total" rows over the last 3 days
    total_rows  <- block %>% filter(meal_type == "Total") %>% head(3)
    macro_means <- sapply(MACRO_VARS, function(var) {
      mean(as.numeric(total_rows[[var]]), na.rm = TRUE)
    })

    # Collection date: day after the last recorded day of the week
    collection_date <- max(total_rows$date, na.rm = TRUE) + 1

    # Date of first menstruation day in the week (if any)
    menstruation_date <- as.Date(NA)
    if (any(block$menstruation_day == TRUE, na.rm = TRUE)) {
      menstruation_date <- block %>%
        filter(menstruation_day == TRUE) %>%
        slice(1) %>%
        pull(date)
    }

    tibble(
      collection_week   = week_name,
      collection_date   = collection_date,
      ID                = info$ID,
      enrollment_date   = info$enrollment_date,
      menstruation_date = menstruation_date,
      BMI               = info$BMI,
      height_cm         = info$height_cm,
      weight_kg         = info$weight_kg,
      !!!as.list(calorie_means),
      !!!as.list(macro_means)
    )
  }

  # Apply inner function over all expected weeks and stack results
  bind_rows(lapply(WEEKS_EXPECTED, compute_week))
}


# -----------------------------------------------------------------------------
# 6. APPLY compute_3day_means() TO ALL PARTICIPANTS
# -----------------------------------------------------------------------------

means_3days <- NULL

for (i in names(diary_list_clean)) {
  result_i    <- compute_3day_means(diary_list_clean[[i]])
  means_3days <- rbind(means_3days, result_i)
}

# Ensure week factor is ordered consistently for downstream analyses
means_3days <- means_3days %>%
  mutate(collection_week = factor(collection_week,
                                  levels = c("first", "second", "third", "fourth")))


# -----------------------------------------------------------------------------
# 7. ADD SUFFIX "_mean3d" TO DIETARY VARIABLES
# -----------------------------------------------------------------------------
# All dietary columns are suffixed with "_mean3d" to identify the aggregation
# method (3-day mean) and distinguish them from other variables in merged datasets.

diet_vars <- c("calories_total", "calories_lunch", "calories_dinner",
               "calories_breakfast", "calories_snack",
               "carbohydrates", "fat", "proteins", "cholesterol",
               "sodium", "sugars", "fiber")

positions <- match(diet_vars, names(means_3days))
names(means_3days)[positions] <- paste0(diet_vars, "_mean3d")


# -----------------------------------------------------------------------------
# 8. OPTIONAL: MANUAL EXCLUSION OF SPECIFIC PARTICIPANT-WEEK COMBINATIONS
# -----------------------------------------------------------------------------
# Use this section to flag known data quality issues.
# Set dietary variables to NA for specific participant-week pairs if needed.

cols_to_na <- paste0(diet_vars, "_mean3d")

# Example: exclude participant 117, second week
# rows_to_exclude <- which(means_3days$ID == 117 & means_3days$collection_week == "second")
# means_3days[rows_to_exclude, cols_to_na] <- NA

# Example: exclude participant 133, first week
# rows_to_exclude <- which(means_3days$ID == 133 & means_3days$collection_week == "first")
# means_3days[rows_to_exclude, cols_to_na] <- NA


# -----------------------------------------------------------------------------
# 9. STATISTICAL ANALYSIS: FRIEDMAN TEST + WILCOXON POST-HOC
# -----------------------------------------------------------------------------
# Test whether dietary intake differs significantly across the four collection
# weeks using the Friedman test (non-parametric repeated measures).
# Post-hoc pairwise comparisons use the Wilcoxon signed-rank test with
# Holm correction for multiple comparisons.

variables_to_test <- cols_to_na

plot_list <- list()

for (var in variables_to_test) {

  # Subset to relevant columns
  df_tmp <- means_3days %>%
    select(ID, collection_week, value = all_of(var))

  # Friedman test across the four collection weeks
  friedman_res <- friedman.test(
    df_tmp$value,
    groups = df_tmp$collection_week,
    blocks = df_tmp$ID
  )

  friedman_label <- if (!is.na(friedman_res$p.value)) {
    paste("Friedman test p =", round(friedman_res$p.value, 4))
  } else {
    "Friedman test: not computable"
  }

  # Pairwise Wilcoxon signed-rank tests with Holm correction
  wilcox_res <- df_tmp %>%
    wilcox_test(value ~ collection_week, paired = TRUE, p.adjust.method = "holm") %>%
    add_xy_position(x = "collection_week", step.increase = 0.5)

  wilcox_res$p.adj <- round(wilcox_res$p.adj, 6)

  # Boxplot with individual data points and significance brackets
  p <- ggboxplot(df_tmp,
                 x       = "collection_week",
                 y       = "value",
                 color   = "collection_week",
                 palette = "jco",
                 add     = "jitter") +
    labs(
      title    = paste("3-Day Mean:", var),
      subtitle = friedman_label,
      x        = "Collection week",
      y        = var
    ) +
    stat_pvalue_manual(wilcox_res, label = "p.adj", tip.length = 0.01)

  plot_list[[var]] <- p
}


# -----------------------------------------------------------------------------
# 10. SAVE OUTPUT
# -----------------------------------------------------------------------------

save(means_3days, file = "means_3days_output.rda")

# Optional: export as CSV
# write.csv(means_3days, file = "means_3days_output.csv", row.names = FALSE)
