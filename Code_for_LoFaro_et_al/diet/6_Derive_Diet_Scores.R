# PAOLA FORABOSCO 28/04/2026

# ============================================================================ #
# DESCRIPTION
# ============================================================================ #

# This R script computes four commonly used dietary indices from FFQ data:
# - Mediterranean Diet Score (aMED)
# - Healthy Eating Index, simplified version (HEI)
# - Dietary Approaches to Stop Hypertension score (DASH)
# - Alternative Healthy Eating Index (AHEI)
#
# Scores are derived using standard definitions adapted to the available FFQ
# food items. The script reads weekly FFQ data, calculates index-specific
# component scores, and outputs total dietary scores for each participant.

# INPUT FILE
# - FFQ_week.csv :  weekly intake of each food item by subject

# ============================================================================ #
# FOOD TRANSLATION (Italian / English)
# ============================================================================ #
# caffe                       - coffee
# birra                       - beer
# vino                        - wine
# superalcolici               - spirits / alcoholic drinks
# latte_intero                - whole milk
# latte_scremato              - skimmed milk
# latte_parzialmente_scremato - semi-skimmed milk
# yogurt                      - yogurt
# formaggi                    - cheese
# carne_rossa                 - red meat
# carne_bianca                - white meat
# pesce                       - fish
# insaccati                   - cured meats / sausages
# uova                        - eggs
# pasta_riso                  - pasta or rice
# pane                        - bread
# verdure                     - vegetables
# frutta                      - fruit
# olio_oliva                  - olive oil
# burro                       - butter

# ============================================================================ #
# INITIALIZATION
# ============================================================================ #

library(dplyr)
library(readr)
library(scales)

# Load FFQ data
ffq <- read_csv("FFQ_week.csv")

# Convert sample to factor
ffq <- ffq %>% mutate(sample = as.factor(sample))

# ============================================================================ #
# Dietary indices calculation: aMED, HEI, DASH, AHEI
# ============================================================================ #

#------------------------------------------#
#  aMED (Adherence to Mediterranean Diet)
#------------------------------------------#

# Calculate median intake for aMED components
medians <- ffq %>%
  summarise(
    veg = median(verdure, na.rm = TRUE),
    fruit = median(frutta, na.rm = TRUE),
    cereals = median(pasta_riso + pane, na.rm = TRUE),
    fish = median(pesce, na.rm = TRUE),
    meat = median(carne_rossa + insaccati, na.rm = TRUE),
    dairy = median(latte_intero + formaggi, na.rm = TRUE)
  )

# Assign binary scores and compute total aMED score
ffq <- ffq %>%
  mutate(
    aMED_veg = if_else(verdure >= medians$veg, 1, 0),
    aMED_fruit = if_else(frutta >= medians$fruit, 1, 0),
    aMED_cereals = if_else((pasta_riso + pane) >= medians$cereals, 1, 0),
    aMED_fish = if_else(pesce >= medians$fish, 1, 0),
    aMED_fat = if_else(olio_oliva > burro, 1, 0),
    aMED_meat = if_else((carne_rossa + insaccati) <= medians$meat, 1, 0),
    aMED_dairy = if_else((latte_intero + formaggi) <= medians$dairy, 1, 0),
    aMED_alcohol = if_else(
      (vino + birra + superalcolici) >= 1 &
        (vino + birra + superalcolici) <= 2, 1, 0
    ),
    aMED_total = rowSums(across(starts_with("aMED_")))
  )

#--------------------------------------------#
#   HEI (Healthy Eating Index simplified)
#--------------------------------------------#

ffq <- ffq %>%
  mutate(
    HEI_veg = rescale(verdure, to = c(0, 10)),
    HEI_fruit = rescale(frutta, to = c(0, 10)),
    HEI_fish = rescale(pesce, to = c(0, 10)),
    HEI_cereals = rescale(pasta_riso + pane, to = c(0, 10)),
    HEI_dairy = rescale(
      latte_scremato + latte_parzialmente_scremato + yogurt, to = c(0, 10)
    ),
    HEI_meat = 10 - rescale(carne_rossa + insaccati, to = c(0, 10)),
    HEI_saturated_fat = 10 - rescale(
      burro + formaggi + latte_intero, to = c(0, 10)
    ),
    HEI_olive_oil = rescale(olio_oliva, to = c(0, 10)),
    HEI_alcohol = case_when(
      (vino + birra + superalcolici) >= 1 &
        (vino + birra + superalcolici) <= 2 ~ 10,
      (vino + birra + superalcolici) > 2 ~ 5,
      TRUE ~ 8
    ),
    HEI_total = rowMeans(across(starts_with("HEI_")), na.rm = TRUE) * 10
  )


#----------------------------------#
#    DASH score (without sweets)
#----------------------------------#

ffq <- ffq %>%
  mutate(
    DASH_veg = ntile(verdure, 5),
    DASH_fruit = ntile(frutta, 5),
    DASH_cereals = ntile(pasta_riso + pane, 5),
    DASH_dairy = ntile(
      latte_scremato + latte_parzialmente_scremato + yogurt, 5
    ),
    DASH_meat = 6 - ntile(carne_rossa + insaccati, 5),
    DASH_total = DASH_veg + DASH_fruit + DASH_cereals +
      DASH_dairy + DASH_meat
  )

#--------------------------------------------#
#   AHEI (Alternative Healthy Eating Index)
#--------------------------------------------#

ffq <- ffq %>%
  mutate(
    AHEI_fruit = rescale(frutta, to = c(0, 10)),
    AHEI_veg = rescale(verdure, to = c(0, 10)),
    AHEI_cereals = rescale(pasta_riso + pane, to = c(0, 10)),
    AHEI_fish = rescale(pesce, to = c(0, 10)),
    AHEI_meat = 10 - rescale(carne_rossa + insaccati, to = c(0, 10)),
    AHEI_healthy_fat = rescale(olio_oliva, to = c(0, 10)),
    AHEI_alcohol = case_when(
      (vino + birra + superalcolici) >= 1 &
        (vino + birra + superalcolici) <= 2 ~ 10,
      (vino + birra + superalcolici) > 2 ~ 5,
      TRUE ~ 8
    ),
    AHEI_total = rowMeans(across(starts_with("AHEI_")), na.rm = TRUE) * 10
  )

#------------------------#
#  Final scores output
#------------------------#

scores <- ffq %>%
  select(ID, sample, aMED_total, HEI_total, DASH_total, AHEI_total)

write_csv(scores, "Diet_Scores.csv")
message("Dietary score calculation completed successfully.")

