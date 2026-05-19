# PAOLA FORABOSCO 28/04/2026
#
# ============================================================================ #
#  DESCRIPTION
# ============================================================================ #

# This script processes Food Frequency Questionnaire (FFQ) data to estimate
# average daily food and nutrient intake for each participant.
#
# MAIN STEPS
# 1. Reads the original FFQ input file (cleaned_questionnaire.csv, visit=0) 
# 2. Converts FFQ consumption frequencies into daily intake frequencies
#    (servings/day) 
# 3. Transforms daily frequencies into daily amounts (grams or mL per day)
#    by multiplying by standard portion sizes
# 4. Merges the data with BDA food composition table to calculate daily 
#    nutrient intakes
# 5. Aggregates nutrient intakes by subject ID
#
# INPUT FILES
# - cleaned_questionnaire.csv (visit=0): original FFQ dataset with consumption
#                                        and frequency data
# - BDA_rid.csv (or BDA_extended.csv) : food composition table (nutrients per 100g/100mL)
#
# OUTOUT FILES
# - FFQ_die.csv :  daily intake of each food item by subject
# - FFQ.csv     :  total daily nutrient intake by subject
# - warnings.txt : log of data consistency checks and warnings

# FOOD TRASLATION
# caffe                    - coffee
# birra                    - beer
# vino                     - wine
# superalcolici            - spirits / alcoholic drinks
# latte_intero             - whole milk
# latte_scremato           - skimmed milk
# latte_parzialmente_scremato - semi-skimmed milk
# yogurt                   - yogurt
# formaggi                 - cheese
# carne_rossa              - red meat
# carne_bianca             - white meat
# pesce                    - fish
# insaccati                - cured meats / sausages
# uova                     - eggs
# pasta_riso               - pasta or rice
# pane                     - bread
# verdure                  - vegetables
# frutta                   - fruit
# olio_oliva               - olive oil
# burro                    - butter

# ============================================================================ #
# INITIALIZATION
# ============================================================================ #

library(tidyr)
library(dplyr)

FFQ.orig <- read.csv("cleaned_questionnaire.csv", header = TRUE, na.strings = "NA")


#-------------------------------------------------#
# STANDARD PORTIONS (grams/mL) for each food item
#-------------------------------------------------#

standard_portions <- c(50, 330, 125, 40, 
                       125, 125, 
                       125, 125, 
                       75, 100, 100, 
                       100, 50, 50, 100, 
                       50, 200, 150, 10, 10)
names(standard_portions) <- c("caffe", "birra", "vino", "superalcolici", 
                              "latte_intero", "latte_scremato", 
                              "latte_parzialmente_scremato", "yogurt", 
                              "formaggi", "carne_rossa", "carne_bianca", 
                              "pesce", "insaccati", "uova", "pasta_riso", 
                              "pane", "verdure", "frutta", "olio_oliva", "burro")


#----------------------------------------#
# IMPORT BDA DATA (nutrient composition)
#----------------------------------------#

# Use reduced or extended BDA table depending on need
BDA <- read.csv("BDA_rid.csv", header = TRUE)
# BDA <- read.csv("BDA_extended.csv", header = TRUE)

# Save warnings to a file
sink("warnings.txt")


# ============================================================================ #
# STEP 1 - Convert FFQ.orig into daily intake (FFQ.die)
# ============================================================================ #

rownames(FFQ.orig) <- FFQ.orig[,1]
FFQ.die <- data.frame(FFQ.orig[,1])
colnames(FFQ.die) <- "ID"
rownames(FFQ.die) <- FFQ.die$ID

#-----------------------------#
# EXTRACT COFFEE DAILY INTAKE
#-----------------------------#

alimento <- "caffe"

dat.tmp <- FFQ.orig[ , grepl(alimento, names(FFQ.orig))]

apply_logic_multiple_caffe <- function(df) {
  
    col1 <- df[["caffe"]]
    col2 <- df[["caffe_frequenza"]]
    col3 <- df[["caffe_quantita_die"]]
    
    # Check 1: caffe == 1 but caffe_frequenza not 1 or 2
    problematici1 <- which(col1 == 1 & !(col2 %in% c(1, 2)))
    # Check 1: consumption == 1 but frequency not 1 or 2
    if (length(problematici1) > 0) {
      problematici1_names <- rownames(df)[problematici1]  # rownames
      warning(sprintf("Valori non attesi in '%s': %d righe con caffe==1 ma caffe_frequenza non in (1,2) per ID: %s", 
                      alimento, length(problematici1), paste(problematici1_names, collapse = ", ")))
    }
    
    # Check 2: caffe_frequenza == 1 but caffe_quantita_die missing
    problematici2 <- which(col2 == 1 & (is.na(col3) | col3 == ""))
    if (length(problematici2) > 0) {
      problematici2_names <- rownames(df)[problematici2]  
      warning(sprintf("Valori mancanti in '%s': %d righe con frequenza==1 ma freq_sett mancante per ID: %s", 
                      alimento, length(problematici2), paste(problematici2_names, collapse = ", ")))
    }
    
    # Check 3: caffe_quantita_die > 7 (possible weekly/day confusion)
    problematici3 <- which(col3 > 7)
    if (length(problematici3) > 0) {
      problematici3_names <- rownames(df)[problematici3]  
      warning(sprintf("Valori eccessivi in '%s': %d righe con caffe_quantita_die >==7 per ID: %s", 
                      alimento, length(problematici3), paste(problematici3_names, collapse = ", ")))
    }
    
    # Main logic to convert to daily intake
    risultato <- case_when(
      col1 == 2 ~ 0,
      col1 == 1 & col2 == 2 ~ 0.5,
      col1 == 1 & col2 == 1 ~ as.numeric(col3),
      TRUE ~ NA_real_
    )
    
    return(risultato)
}

FFQ.die [,alimento] <- apply_logic_multiple_caffe(dat.tmp)

print(warnings())


#------------------------------#
# EXTRACT ALCOHOL DAILY INTAKE
#------------------------------#

alimento <- "alcol"

dat.tmp <- FFQ.orig[ , grepl(alimento, names(FFQ.orig))]

dat.tmp$birra <- FFQ.orig[,"birra_frequenza"]
dat.tmp$vino <- FFQ.orig[,"vino_frequenza"]
dat.tmp[is.na(dat.tmp)] <- 0

apply_logic_multiple_alcol <- function(df) {
  
  col1 <- df[["alcol"]]
  col2 <- df[["alcol_frequenza"]]

  # Check 1:  alcol == 1 but alcol_frequenza not 1 or 2
  problematici1 <- which(col1 == 1 & !(col2 %in% c(1, 2)))
  if (length(problematici1) > 0) {
    problematici1_names <- rownames(df)[problematici1]  
    warning(sprintf("Valori non attesi in '%s': %d righe con alcol==1 ma alcol_frequenza non in (1,2) per ID: %s", 
                    alimento, length(problematici1), paste(problematici1_names, collapse = ", ")))
  } 
  
    col3 <- df[["birra"]]
    # Main logic to convert to daily intake
    risultato1 <- case_when(
    col1 == 2 ~ 0,
    col1 == 1 & col2 == 2 ~ 0.033,  # 0.1 divided by 3 (not knowing alcohol type)
    col1 == 1 & col2 == 1 ~ as.numeric(col3)/7,
    TRUE ~ NA_real_
    )
    
    col3 <- df[["vino"]]
    # Main logic to convert to daily intake
    risultato2 <- case_when(
      col1 == 2 ~ 0,
      col1 == 1 & col2 == 2 ~ 0.033,
      col1 == 1 & col2 == 1 ~ as.numeric(col3)/7,
      TRUE ~ NA_real_
      )
    
    col3 <- df[["superalcolici_frequenza"]]
    # Main logic to convert to daily intake
    risultato3 <- case_when(
      col1 == 2 ~ 0,
      col1 == 1 & col2 == 2 ~ 0.033,
      col1 == 1 & col2 == 1 ~ as.numeric(col3)/7,
      TRUE ~ NA_real_
    )

      return(c(risultato1,risultato2,risultato3))
  
}

FFQ.die [,c("birra","vino", "superalcolici")] <- apply_logic_multiple_alcol(dat.tmp)

print(warnings())

#-----------------------------#
# EXTRACT MILK DAILY INTAKE
#-----------------------------#

alimento <- "latte"

dat.tmp <- FFQ.orig[ , grepl(alimento, names(FFQ.orig))]

# Separate milk types
tipo.latte <- data.frame(
  intero = dat.tmp[,2],
  scremato = dat.tmp[,3],
  parzialmente_scremato = dat.tmp[,4])
tipo.latte$somma <- rowSums(tipo.latte[, c("intero", "scremato", "parzialmente_scremato")])
rownames(tipo.latte) <- rownames(dat.tmp)

apply_logic_multiple_latte <- function(df) {
 
    col1 <- df[["latte_consumo"]]
    col2 <- df[["latte_frequenza"]]
    col3 <- df[["latte_freq_settimanale"]]
    
    # Check 1: consumo == 1 but frequenza not 1 or 2
    problematici1 <- which(col1 == 1 & !(col2 %in% c(1, 2)))
    if (length(problematici1) > 0) {
    problematici1_names <- rownames(df)[problematici1]  
    warning(sprintf("Valori non attesi in '%s': %d righe con consumo==1 ma frequenza non in (1,2) per ID: %s", 
                  alimento, length(problematici1), paste(problematici1_names, collapse = ", ")))
    }
   
    # Check 2: frequenza == 1 but freq_sett missing
    problematici2 <- which(col2 == 1 & (is.na(col3) | col3 == ""))
    if (length(problematici2) > 0) {
    problematici2_names <- rownames(df)[problematici2]  
     warning(sprintf("Valori mancanti in '%s': %d righe con frequenza==1 ma freq_sett mancante per ID: %s", 
                  alimento, length(problematici2), paste(problematici2_names, collapse = ", ")))
    }

    # Convert milk intake to daily frequency
    risultato <- case_when(
    col1 == 2 ~ 0,
    col1 == 1 & col2 == 2 ~ 0.1,
    col1 == 1 & col2 == 1 ~ as.numeric(col3) / 7,
    TRUE ~ NA_real_
    )

return(risultato)

}

print(warnings())

# Apply daily frequency to each milk type proportional to total

FFQ.die [,"latte_intero"] <- apply_logic_multiple_latte(dat.tmp)*tipo.latte[["intero"]]
FFQ.die [,"latte_intero"] <- ifelse(FFQ.die [,"latte_intero"]  != 0, 
                                      FFQ.die [,"latte_intero"]  / tipo.latte$somma, 0)  

FFQ.die [,"latte_scremato"] <- apply_logic_multiple_latte(dat.tmp)*tipo.latte[["scremato"]]
FFQ.die [,"latte_scremato"] <- ifelse(FFQ.die [,"latte_scremato"]  != 0, 
                                     FFQ.die [,"latte_scremato"]  / tipo.latte$somma, 0)  

FFQ.die [,"latte_parzialmente_scremato"] <- apply_logic_multiple_latte(dat.tmp)*tipo.latte[["parzialmente_scremato"]]
FFQ.die [,"latte_parzialmente_scremato"] <- ifelse(FFQ.die [,"latte_parzialmente_scremato"]  != 0, 
                                      FFQ.die [,"latte_parzialmente_scremato"]  / tipo.latte$somma, 0)  

rm(tipo.latte)

#-------------------------#
# EXTRACT ALL OTHER FOODS
#-------------------------#

apply_logic_multiple <- function(df) {
  
  if (ncol(df) %% 3 != 0) {
    stop("The number of columns must be a multiple of 3.")
  }
  col1 <- df[[paste0(alimento, "_consumo")]]
  col2 <- df[[paste0(alimento, "_freq")]]
  col3 <- df[[paste0(alimento, "_freq_sett")]]
  
  # Check 1: consumo == 1 but frequenza not 1 or 2
  problematici1 <- which(col1 == 1 & !(col2 %in% c(1, 2)))
  if (length(problematici1) > 0) {
    problematici1_names <- rownames(df)[problematici1]  # Ottieni i nomi delle righe
    warning(sprintf("Valori non attesi in '%s': %d righe con consumo==1 ma frequenza non in (1,2) per ID: %s", 
                    alimento, length(problematici1), paste(problematici1_names, collapse = ", ")))
  }
  
  # Check 2: frequenza == 1 but freq_sett missing
  problematici2 <- which(col2 == 1 & (is.na(col3) | col3 == ""))
  if (length(problematici2) > 0) {
    problematici2_names <- rownames(df)[problematici2]  # Ottieni i nomi delle righe
    warning(sprintf("Valori mancanti in '%s': %d righe con frequenza==1 ma freq_sett mancante per ID: %s", 
                    alimento, length(problematici2), paste(problematici2_names, collapse = ", ")))
  }
  
  # Main logic to convert to daily intake
  risultato <- case_when(
    col1 == 2 ~ 0,
    col1 == 1 & col2 == 2 ~ 0.1,
    col1 == 1 & col2 == 1 ~ as.numeric(col3) / 7,
    TRUE ~ NA_real_
  )
  
  return(risultato)
}

alimenti = c("yogurt","formaggi","carne_rossa","carne_bianca","pesce",
               "insaccati","uova","pasta_riso","pane","verdure","frutta",
               "olio_oliva","burro")

for (alimento in alimenti) {
 dat.tmp <- FFQ.orig[ , grepl(alimento, names(FFQ.orig))]
  colnames(dat.tmp) <- c(paste0(alimento,"_consumo"),paste0(alimento,"_freq"),
                         paste0(alimento,"_freq_sett"))

  FFQ.die [,alimento] <- apply_logic_multiple(dat.tmp)
  
}

rm(alimento)
rm(alimenti)
rm(dat.tmp)
print(warnings())
sink()

# ============================================================================ #
# STEP 2 - Compute daily g/mL by multiplying frequencies by standard portions
# ============================================================================ #

# Column check: the names of standard portions and FFQ.die must match exactly
if (!identical(colnames(FFQ.die[,-1]), names(standard_portions))) {
  message("I nomi delle colonne in FFQ e i nomi di standard portions NON sono uguali.")
} else {
  message("I nomi delle colonne in FFQ e i nomi di standard portions sono uguali. OK")
}


# Method: multiplication by column
FFQ.input<- sweep(FFQ.die[,-1], 2, standard_portions, `*`)

# ============================================================================ #
# STEP 3 - Convert FFQ.input to long format
# ============================================================================ #

FFQ.output <- FFQ.input %>%
  mutate(ID = rownames(.)) %>%
  pivot_longer(cols = -ID, names_to = "alimento", values_to = "quantita")

rm (FFQ.input)

# ============================================================================ #
# STEP 4 - Merge with BDA nutrient data and sum by ID. Save output files
# ============================================================================ #

FFQ.output.bda <- merge(FFQ.output, BDA, by = "alimento", all.x = TRUE)
kk <- length(FFQ.output.bda[1,])

rm(FFQ.output)

# Multiply nutrient content by amount and scale per 100g
FFQ.output.bda <- cbind(FFQ.output.bda[,2], 
                        (FFQ.output.bda[,c(4:kk)] * FFQ.output.bda$quantita)/100)
colnames(FFQ.output.bda) <- c("ID", "Calories", "Carbohydrates", "Fats", "Proteins",
                              "Cholesterol", "Sodium", "Sugar", "Fiber")
FFQ.sum <- aggregate(. ~ ID, data = FFQ.output.bda , FUN = sum)

rm(FFQ.output.bda)
rm(kk)

write.csv(FFQ.sum,file="FFQ.csv",row.names=F)
write.csv(FFQ.die,file="FFQ_die.csv",row.names=F)



