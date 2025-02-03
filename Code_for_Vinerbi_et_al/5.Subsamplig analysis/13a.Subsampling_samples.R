# Rscript

# Script to make a subsamplig fastq to 60kreads and 30k reads

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# Version R: v4.4.1

# Upload library
library(ShortRead)

# paths
setwd('/home/')
getwd()

# Set paths files fastQ
input_dir <- "FastQ_file/"
output_dir <- "FastQ_subsampling/"

# List of samples R1 and R2
files_R1 <- list.files(input_dir, pattern = "_R1_paired\\.fq\\.gz$", full.names = TRUE)
files_R2 <- list.files(input_dir, pattern = "_R2_paired\\.fq\\.gz$", full.names = TRUE)

# check lenght R1 and R2
if (length(files_R1) != length(files_R2)) {
  stop("The file numbers R1 and R2 do not match")
}

# Subsampling
set.seed(23122430) # Seed 

for (i in seq_along(files_R1)) {
  reads_R1 <- readFastq(files_R1[i])
  reads_R2 <- readFastq(files_R2[i])
  
  # Ceck same lenght between R1 and R2
  if (length(reads_R1) != length(reads_R2)) {
    stop(paste("The number of readings does not match between R1 and R2 for the file:", basename(files_R1[i])))
  }
  
  FQ_sub_R1 <- sample(reads_R1, 60000, replace = FALSE) # change with 30000
  FQ_sub_R2 <- sample(reads_R2, 60000, replace = FALSE) # change with 30000
  
  # Save file subsampled
  output_R1 <- file.path(output_dir, gsub("\\.fq\\.gz$", ".fq.gz", basename(files_R1[i]))) 
  output_R2 <- file.path(output_dir, gsub("\\.fq\\.gz$", ".fq.gz", basename(files_R2[i]))) 
  
  writeFastq(FQ_sub_R1, output_R1, compress = TRUE)
  writeFastq(FQ_sub_R2, output_R2, compress = TRUE)
  
  # final message
  message("Subsampling completed for: ", basename(files_R1[i]))
}
