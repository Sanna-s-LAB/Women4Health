# PAOLA FORABOSCO 28/04/2026

# ============================================================================ #
#  DESCRIPTION
# ============================================================================ #

# This script loads 3-day nutrition data from three locations (TRIESTE, CAGLIARI, 
# BOLOGNA), merges with Age and BMI data, adds visit IDs and phase information,
# calculates mean nutrient intake per participant per hormonal phase, 
# and exports the data in wide format.
#
# input files:
#   - TRIESTE_3day.rda
#   - CAGLIARI_3day.rda
#   - BOLOGNA_3day.rda
#   - PHENO.csv (Age and BMI)
#   - PHASE.csv (phase information for each ID_visit)
#
# output files:
#   - Diaries_3day_phases.csv (3day means for each hormonal phase, wide format)
#
# -----------------
# Input file format
# -----------------
# 
# PHENO.csv 
# ---------
# This CSV file contains Age and BMI information for each participant.
# Required columns:
# - ID2: participant identifier
# - Age: age of the participant
# - BMI: body mass index
# Only these three columns are used.
# ID2 is renamed to ID inside the script and must match the formatted IDs used 
# in the dietary data (e.g. X001, S023, B105).
#
# PHASE.csv (phase information)
# ---------
# This CSV file links each participant visit to a hormonal phase.
# Required columns:
# - ID_visit: unique visit identifier, created as ID_number (e.g. X001_1, S023_2)
# - phase: numeric phase label (expected values: 1, 2, 3, 4)
# This file is used to assign each visit to a hormonal phase.
#
# TRIESTE_3day.rda/CAGLIARI_3day.rda/BOLOGNA_3day.rda (3-day dietary records)
# ---------------------------------------------------
# Each .rda file must contain an object named DATA, which is a data frame with 
# 3-day averaged dietary intake for each visit 
# Required columns in DATA:
# - ID: numeric participant ID (will be reformatted with prefixes)
# - weeks_raccolta: visit label (first, second, third, fourth)
# - Nutrient variables (original Italian names):
# Calorie_medie3gg/Carboidrati_medie3gg/Grassi_medie3gg/Proteine_medie3gg
# Colest_medie3gg/Sodio_medie3gg/Zuccheri_medie3gg/Fibre_medie3gg
#
#
# ============================================================================ #
#  INITIALIZATION
# ============================================================================ #

# load libraries
library(dplyr)       # data manipulation
library(tidyr)       # data completion and reshaping

# load phenotype data
pheno <- read.csv("PHENO.csv", header=T)         # load age & bmi data
pheno <- pheno[, c("ID2", "Age","BMI")]
colnames(pheno) <- c("ID","Age", "BMI")          # standardize column names

# load phase data
phase <- read.csv("PHASE.csv", header=T) # load phase data
colnames(phase) <- c("ID_visit","phase")

# define nutrient variables
traits <- c("Calories", "Carbohydrates", "Fats", "Proteins",
            "Cholesterol", "Sodium", "Sugar", "Fiber")  

# original column names in rda files
orig.traits <- c("Calorie_medie3gg", "Carboidrati_medie3gg", "Grassi_medie3gg", 
                 "Proteine_medie3gg", "Colest_medie3gg", "Sodio_medie3gg", 
                 "Zuccheri_medie3gg", "Fibre_medie3gg")      

# ============================================================================ #
#  LOAD AND MERGE 3-day DATA
# ============================================================================ #

# TRIESTE

load("TRIESTE_3day.rda")
DATA[is.nan(as.matrix(DATA))] <- NA  # replace NaN with NA

dat.tmp <- data.frame(DATA[, c("ID", "weeks_raccolta")])
colnames(dat.tmp) <- c("ID", "visit")
dat.tmp <- dat.tmp %>%
  mutate(ID = sprintf("X%03d", ID))  # IDs for TRIESTE
dat.tmp$sample <- "TRIESTE"
dat.tmp[, traits] <- DATA[, orig.traits]
dat <- dat.tmp  # initialize main dataframe


# CAGLIARI

load("CAGLIARI_3day.rda")
DATA[is.nan(as.matrix(DATA))] <- NA
dat.tmp <- data.frame(DATA[, c("ID", "weeks_raccolta")])
colnames(dat.tmp) <- c("ID", "visit")
dat.tmp <- dat.tmp %>%
  mutate(ID = sprintf("S%03d", ID))  # IDs for CAGLIARI
dat.tmp$sample <- "CAGLIARI"
dat.tmp[, traits] <- DATA[, orig.traits]
dat <- rbind(dat, dat.tmp)  # combine with main dataframe


# BOLOGNA

load("BOLOGNA_3day.rda")
DATA[is.nan(as.matrix(DATA))] <- NA
dat.tmp <- data.frame(DATA[, c("ID", "weeks_raccolta")])
colnames(dat.tmp) <- c("ID", "visit")
dat.tmp <- dat.tmp %>%
  mutate(ID = sprintf("B%03d", ID))  # IDs for BOLOGNA
dat.tmp$sample <- "BOLOGNA"
dat.tmp[, traits] <- DATA[, orig.traits]
dat <- rbind(dat, dat.tmp)  # combine with main dataframe

# ============================================================================ #
#  Merge with phenotype data 
# ============================================================================ #

dat <- merge(dat, pheno, by = "ID")           # add Age & BMI
dat <- dat %>%
  select(ID, Age, BMI, visit, everything())   # reorder columns

# Create ID_visit
visit_order <- c("first", "second", "third", "fourth")  # order of visits
dat$ID_visit <- paste0(dat$ID, "_", match(dat$visit, visit_order))  

dat <- dat %>%
  select(ID, ID_visit, Age, BMI, everything())  # reorder columns

# ============================================================================ #
#  Merge with phase data 
# ============================================================================ #

dat <- merge(dat, phase, by = "ID_visit")  # add phase information
dat <- dat %>%
  select(ID, ID_visit, phase, Age, BMI, everything(), -visit)  # reorder columns

# Define all possible phases
all_phases <- 1:4

# Calculate mean nutrients per ID × phase and add missing phases
dat.phase <- dat %>%
  group_by(ID, phase) %>%
  summarise(
    across(all_of(traits), ~ mean(.x, na.rm = TRUE)),  # mean of nutrients
    Age = first(Age),
    BMI = first(BMI),
    sample = first(sample),
    .groups = "drop"
  ) %>%
  # Ensure each ID has all phases
  complete(ID, phase = all_phases) %>%  # don't fill Age/BMI yet
  group_by(ID) %>%
  fill(Age, BMI, sample, .direction = "downup") %>%  # propagate known values
  ungroup() %>%
  select(ID, phase, Age, BMI, sample, everything())

# ============================================================================ #
#  Transform to wide format & export file 
# ============================================================================ #

dat.phase.wide <- dat.phase %>%
  pivot_wider(
    id_cols = c(ID, Age, BMI, sample),  # identifiers
    names_from = phase,                 # each phase becomes new columns
    values_from = all_of(traits),       # nutrient values fill columns
    names_sep = "_"                     # separator: nutrient_phase
  )

# Export wide-format CSV
write.csv(dat.phase.wide, "Diaries_3day_phases.csv", row.names = FALSE) 


  