library(dplyr)
library(readr)
library(stringr)
library(purrr)

path <- "/mnt/sannaLAB-Temp/dasha/PRS/results_b12/Plink/"

files <- list.files(path, pattern = "5e-08\\.sscore$", full.names = TRUE)

read_sscore <- function(f) {
  
  protein <- basename(f) %>%
    str_remove("\\.sscore$") %>%
    str_remove("\\.[0-9]+e-[0-9]+$")    # remove scientific-notation suffix
  
  df <- read_tsv(f, show_col_types = FALSE)
  
  df %>%
    select(IID, PRS = SCORE1_AVG) %>%
    rename(!!protein := PRS)
}

merged <- files %>%
  map(read_sscore) %>%
  reduce(full_join, by = "IID") %>%
  arrange(IID)

write_tsv(merged, "merged_protein_PRS.tsv")
