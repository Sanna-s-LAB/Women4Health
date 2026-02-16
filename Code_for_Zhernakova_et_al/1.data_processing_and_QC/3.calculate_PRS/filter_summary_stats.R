args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

prot = args[1]
pqtl_fname = args[2]
bim_fname <- args[3]
snp_col = args[4]
#snp_col = 'chrpos'

#pval_thres <- args[4]

setwd('/mnt/sannaLAB-Temp/dasha/PRS/data')
#prot = 'DFFA'
#pqtl_fname <- "pQTLs/invn_OID20620_UKBB_PPP_F_only_regenie.gz"
#bim_fname <- "vaginal_MGS/MGS_VAGINAL.QC.hg19.rmdup.bim"
pval_thres <- 1


pqtls <- fread(pqtl_fname, data.table = F)
bim <- fread(bim_fname, data.table = F, col.names = c("CHROM", "SNP", "v1", "POS", "A1", "A2"))
cat("\n\nFiltering base data:", pqtl_fname, "\n")
cat("Read summary statistics for", nrow(pqtls), "SNPs\n")

pqtls$MarkerName <- paste0(pqtls$chromosome, ":", pqtls$base_pair_location)
colnames(pqtls) <- gsub("MarkerName", "chrpos", colnames(pqtls))

pqtls <- pqtls[pqtls[,snp_col] %in% bim$SNP,]
cat("After removing SNPs not present in the genotype data", nrow(pqtls), "SNPs remain\n")

pqtls <- pqtls[pqtls$p_value < pval_thres,]
cat("After selecting SNPs with association p-value below", pval_thres, nrow(pqtls), "SNPs remain\n")


# Remove duplicate SNP IDs

snp_counts <- as.data.frame(table(pqtls[,snp_col]))
duplicate_ids <- snp_counts[snp_counts$Freq > 1, "Var1"]

cat("Identified", length(duplicate_ids), "duplicate SNP IDs\n")

if (length(duplicate_ids) == 0){
  cat ("No duplicates identified, writing the filtered file\n")
  write.table(pqtls, file = gsub(".gz$",".in_microarray.rmdup.txt",pqtl_fname), quote = F, sep = "\t", row.names = FALSE)
} else {
  duplicates <- pqtls[pqtls[,snp_col] %in% duplicate_ids,]
  
  # Function to check if two allele pairs match (considering strand and order)
  alleles_match <- function(a1, a2, b1, b2) {
    # Check direct match or reverse order
    direct_match <- (a1 == b1 & a2 == b2) | (a1 == b2 & a2 == b1)
    
    # Check complementary pairs (A/T, C/G)
    complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
    a1_comp <- complement[a1]
    a2_comp <- complement[a2]
    
    if (is.na(a1_comp) || is.na(a2_comp)) {
      # If alleles aren't standard nucleotides, just check direct match
      return(direct_match)
    }
    
    comp_match <- (a1_comp == b1 & a2_comp == b2) | (a1_comp == b2 & a2_comp == b1)
    
    return(direct_match | comp_match)
  }
  
  # Process duplicates
  keep_rows <- logical(nrow(pqtls))  # Initialize all FALSE
  keep_rows[!pqtls[,snp_col] %in% duplicate_ids] <- TRUE  # Keep all non-duplicates
  
  for (snp_id in duplicate_ids) {
    # Get all rows with this SNP ID in pqtls
    pqtl_rows <- which(pqtls[,snp_col] == snp_id)
    # Get the corresponding row in bim file
    bim_row <- bim[bim$SNP == snp_id, ]
    
    if (nrow(bim_row) == 0) next  # Shouldn't happen due to previous filtering
    
    # Get bim alleles (assuming columns 5 and 6 are the alleles)
    bim_allele1 <- bim_row[1, 5]
    bim_allele2 <- bim_row[1, 6]
    
    # Check which pqtl rows match the bim alleles
    matches <- sapply(pqtl_rows, function(i) {
      alleles_match(pqtls$effect_allele[i], pqtls$other_allele[i], bim_allele2, bim_allele1)
    })
    
    if (sum(matches) == 0) {
      # No matches found, just keep the first one
      warning(paste("No allele match found for SNP", snp_id, 
                    "between summary stats and bim file. Keeping first occurrence."))
      keep_rows[pqtl_rows[1]] <- TRUE
    } else {
      # Keep only the first matching row
      first_match <- pqtl_rows[which(matches)[1]]
      keep_rows[first_match] <- TRUE
      
      # If there are multiple matches, warn
      if (sum(matches) > 1) {
        warning(paste("Multiple allele matches found for SNP", snp_id, 
                      "between summary stats and bim file. Keeping first match."))
      }
    }
  }
  
  # Filter the pqtls data frame
  pqtls <- pqtls[keep_rows, ]
  cat("After removing duplicate SNPs", nrow(pqtls), "SNPs remain\n")
  
  # Verify no duplicates remain
  if (any(duplicated(pqtls$ID))) {
    warning("Some duplicate SNPs still remain after filtering")
  } else {
    cat("Successfully removed all duplicate SNPs\n\n\n")
  }
  
  write.table(pqtls, file = gsub(".gz$",".in_microarray.rmdup.txt",pqtl_fname), quote = F, sep = "\t", row.names = FALSE)
}
