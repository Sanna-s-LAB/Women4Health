# main directory
base_dir <- "../Final_lagged"

# list of subfolders to scan
subfolders <- c(
  "Covariates",
  "linear_on_non_linear",
  "linear_on_non_linear_with_t",
  "Nonlinear_with_time",
  "Linear_simulation_with_-0.2",
  "Linear_simulation_with_-1.2",
  "Linear_simulation_with_0",
  "Linear_simulation_with_0.2",
  "Linear_simulation_with_1.2",
  "Quadratic_simulation_with_0.2",
  "Quadratic_simulation_with_1.2",
  "Missing",
  "Quadratic_simulation_with_0",
  "Nonlinear",
  "Quadratic_simulation_with_-0.2",
  "Quadratic_simulation_with_-1.2"
)

library(dplyr)
library(readr)
library(purrr)
library(fs)
library(stringr)

read_csv_from_subfolder <- function(main_folder) {

  main_path <- file.path(base_dir, main_folder)

  # list all CSV files recursively
  csv_files <- list.files(
    main_path,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )

  map_dfr(csv_files, function(f) {

    # read the CSV file
    df <- read_csv(f, show_col_types = FALSE)

    # extract relative path and file name components
    rel_path  <- path_rel(f, main_path)
    parts     <- str_split(rel_path, .Platform$file.sep)[[1]]
    file_name <- tools::file_path_sans_ext(tail(parts, 1))

    # build a base row identifier depending on folder depth
    if (length(parts) > 1) {
      sub_folder <- parts[1]
      base_row_id <- paste(main_folder, sub_folder, file_name, sep = "_")
    } else {
      base_row_id <- paste(main_folder, file_name, sep = "_")
    }

    # if the file has more than one row, append the row index
    df %>%
      mutate(
        row_id = if (n() > 1) {
          paste0(base_row_id, "_int", row_number())
        } else {
          base_row_id
        }
      )
  })
}

# merge all subfolders
final_df <- map_dfr(subfolders, read_csv_from_subfolder)

# reorder columns and remove unnecessary one
final_df <- final_df %>% relocate(row_id)
final_df$...1 <- NULL

library(openxlsx)
write.xlsx(final_df, "../table_power.xlsx")
