# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================

## Integrates metadata, sequencing statistics, phenotypes, and dietary information ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(data.table)
library(dplyr)

# ------------------------------------------------------------
### INPUT FILES
# ------------------------------------------------------------
# Core datasets
dna_file        <- "W4H_MGS_Fecal_dna_12.txt"
h_file_phases   <- "cleaned_phenotypes_251125_uniformed_adjusted.withHOMA.log_some.phase_avg.txt"
ph_file         <- "phases_251125.csv"

# Sequencing read statistics (from preprocessing)
rs_file1 <- "MGS_FECAL_batch1_reads.summary"
rs_file2 <- "MGS_FECAL_batch2_reads.summary"

# Samples to exclude
sample_exclusions <- c("")

# ------------------------------------------------------------
### IMPORT DATA
# ------------------------------------------------------------
rs1  <- fread(rs_file1)
rs2  <- fread(rs_file2)
meta <- fread(dna_file)
horm <- fread(h_file_phases)
ph   <- fread(ph_file)

# ------------------------------------------------------------
### READ STATISTICS INTEGRATION
# ------------------------------------------------------------
# Harmonize column names
colnames(rs1)[c(1,3)] <- c("sample_ID", "clean_paired")
colnames(rs2)[c(1,3)] <- c("sample_ID", "clean_paired")

# Keep relevant columns only
rs1 <- rs1[, c("sample_ID", "clean_paired")]
rs2 <- rs2[, c("sample_ID", "clean_paired")]

# Combine batches and rename
rs <- rbind(rs1, rs2)
colnames(rs)[2] <- "total_reads_postQC"

# Merge with metadata
meta <- merge(meta, rs, by = "sample_ID", all.x = TRUE)

# ------------------------------------------------------------
### BASIC METADATA CLEANING
# ------------------------------------------------------------
# Remove samples with zero DNA concentration
meta <- meta[meta$DNA_concentration > 0,]

# Create simplified sample identifiers
meta$simplified_ID_visit <- gsub("CNRITA-001-|CNRITA-004-", "", meta$sample_ID)

# Compute sequencing depth (in millions of reads)
meta$readsDepth_postQC <- meta$total_reads_postQC / 1e6

# ------------------------------------------------------------
### PHASE ANNOTATION
# ------------------------------------------------------------
# Prepare phase table
colnames(ph)[2] <- "phase"
ph$tmp_Code <- gsub("X", "", ph$Code)

# Merge phase information
meta <- merge(meta, ph[, -"Code"],
              by.x = "simplified_ID_visit",
              by.y = "tmp_Code",
              all.x = TRUE)

# ------------------------------------------------------------
### HORMONE DATA INTEGRATION
# ------------------------------------------------------------
# Harmonize variable names
colnames(horm) <- gsub("17BES", "BES", colnames(horm))

# Create phase-compatible sample IDs
horm$tmp_Code <- gsub("X", "", horm$SampleID)
horm$simplified_ID_phase <- gsub("_F$", "_1", horm$tmp_Code)
horm$simplified_ID_phase <- gsub("_O$", "_2", horm$simplified_ID_phase)
horm$simplified_ID_phase <- gsub("_EL$", "_3", horm$simplified_ID_phase)
horm$simplified_ID_phase <- gsub("_LL$", "_4", horm$simplified_ID_phase)

# Generate matching IDs in metadata
meta$simplified_ID_phase <- paste0(
  gsub("_.*", "", meta$simplified_ID_visit), "_", meta$phase
)

# Merge hormone data
meta <- merge(meta, horm, by = "simplified_ID_phase", all.x = TRUE)

# ------------------------------------------------------------
### QUESTIONNAIRE / PHENOTYPES
# ------------------------------------------------------------
p_file <- "cleaned_questionnaire_251125.csv"

pheno <- fread(p_file)

# Keep visits only (exclude baseline if coded as 0)
pheno <- pheno[pheno$Visit_number != 0,]

# Select relevant variables
pheno <- pheno[, c("Code", "Age", "BMI",
                   "raccolta_quando", "Bristol_stool_scale")]

# Harmonize sample IDs
colnames(pheno)[1] <- "simplified_ID_visit"
pheno$simplified_ID_visit <- gsub("X", "", pheno$simplified_ID_visit)

# ------------------------------------------------------------
### SAMPLE COLLECTION TIME
# ------------------------------------------------------------
# Recode collection timing
pheno$raccolta_quando <- gsub(2, "yesterday_morning", pheno$raccolta_quando)
pheno$raccolta_quando <- gsub(1, "today_morning", pheno$raccolta_quando)
pheno$raccolta_quando <- gsub(3, "yesterday_day", pheno$raccolta_quando)

# Define grouped variable
pheno$Collection_feces_daytime <- "NA"

pheno[
  pheno$raccolta_quando == "yesterday_morning" |
  pheno$raccolta_quando == "today_morning" &
  pheno$Collection_feces_daytime == "NA"
]$Collection_feces_daytime <- "morning"

pheno[
  pheno$raccolta_quando == "yesterday_day" &
  pheno$Collection_feces_daytime != "morning"
]$Collection_feces_daytime <- "not_morning"

# Set remaining to NA
pheno[pheno$Collection_feces_daytime == "NA",]$Collection_feces_daytime <- NA

pheno$raccolta_quando <- NULL

# ------------------------------------------------------------
### BRISTOL STOOL SCALE CLEANING
# ------------------------------------------------------------
pheno$Bristol_stool_scale <- gsub(
  "Hard stools \\(constipation\\)",
  "hard_stools_constipation",
  pheno$Bristol_stool_scale
)

pheno$Bristol_stool_scale <- gsub(
  "Loose stools \\(diarrhea\\)",
  "loose_stools_diarrhea",
  pheno$Bristol_stool_scale
)

pheno$Bristol_stool_scale <- gsub(
  "Normal stools",
  "normal_stools",
  pheno$Bristol_stool_scale
)

# Merge phenotype data
meta <- merge(meta, pheno, by = "simplified_ID_visit", all.x = TRUE)

# ------------------------------------------------------------
### FOOD DATA INTEGRATION
# ------------------------------------------------------------
ffq_week_file     <- "Freeze_DataCleaned_Oct2025/FFQ_week.csv"
ffq_dscores_file  <- "Freeze_DataCleaned_Oct2025/FFQ_diet_scores.csv"
food_cluster_file <- "Freeze_DataCleaned_Oct2025/ALL_clusters.csv"
f_phases_file     <- "Freeze_DataCleaned_Oct2025/ALL_Diaries_HP_3days_corrected.csv"

ffq_week   <- fread(ffq_week_file)
ffq_dscores<- fread(ffq_dscores_file)
f_clst     <- fread(food_cluster_file)
f_means    <- fread(f_phases_file)

# ------------------------------------------------------------
### FFQ PROCESSING
# ------------------------------------------------------------
# Translate column names to English
eng_colnames <- c(
  "ID","Coffee_ffq_week","Beer_ffq_week","Wine_ffq_week",
  "Spirits_ffq_week","Whole_milk_ffq_week","Skimmed_milk_ffq_week",
  "Semi_skimmed_milk_ffq_week","Yogurt_ffq_week","Cheese_ffq_week",
  "Red_meat_ffq_week","White_meat_ffq_week","Fish_ffq_week",
  "Cold_cuts_ffq_week","Eggs_ffq_week","Pasta_or_rice_ffq_week",
  "Bread_ffq_week","Vegetables_ffq_week","Fruit_ffq_week",
  "Olive_oil_ffq_week","Butter_ffq_week"
)

colnames(ffq_week) <- eng_colnames

# Remove sparsely represented variables
ffq_week$Skimmed_milk_ffq_week <- NULL

# ------------------------------------------------------------
### DIET SCORES + CLUSTERS
# ------------------------------------------------------------
# Clean diet scores
ffq_dscores <- ffq_dscores[, -c("sample")]
colnames(ffq_dscores) <- gsub("_total", "_score", colnames(ffq_dscores))

# Clean clusters
f_clst <- f_clst[, -c("PCA1","PCA2","PCA3",
                     "cluster_extended","sample")]

colnames(f_clst) <- c("ID","cluster_FFQ_macro","cluster_diaries")

# Merge all food-related data
food_tmp <- merge(ffq_week, ffq_dscores, by = "ID", all = TRUE)
food_tmp <- merge(food_tmp, f_clst, by = "ID", all = TRUE)

# ------------------------------------------------------------
### 3-DAY DIETARY AVERAGES
# ------------------------------------------------------------
# Clean column names
colnames(f_means) <- gsub("_3days", "", colnames(f_means))

# Create matching IDs
f_means$simplified_ID_phase <- paste0(
  gsub("X", "", f_means$ID), "_", f_means$phase
)

# Remove redundant variables
f_means <- f_means[, -c("sample","Age","BMI","phase")]

# Rename columns
colnames(f_means) <- paste0(colnames(f_means), "_mean3gg")
colnames(f_means) <- gsub("ID_mean3gg", "ID", colnames(f_means))
colnames(f_means) <- gsub("simplified_ID_phase_mean3gg",
                          "simplified_ID_phase",
                          colnames(f_means))

# ------------------------------------------------------------
### MERGE FOOD WITH METADATA
# ------------------------------------------------------------
tmp <- merge(meta, f_means[, -"ID"],
             by = "simplified_ID_phase", all.x = TRUE)

tmp$tmpID <- gsub("_.*", "", tmp$simplified_ID_visit)

food_tmp$ID <- gsub("X", "", food_tmp$ID)

meta <- merge(tmp, food_tmp,
              by.x = "tmpID", by.y = "ID", all.x = TRUE)

# ------------------------------------------------------------
### FINAL METADATA CLEANING
# ------------------------------------------------------------
to_remove <- c("tmpID","simplified_ID_phase","simplified_ID_visit",
               "total_reads_postQC","SampleID","GL","AST","ALT",
               "TRI","COL","HDL","LDL","INS","TSH","FT4","TST",
               "HOMA_IR","HOMA_B","tmp_Code")

idx_remove <- match(to_remove, names(meta), nomatch = 0)
meta <- meta[, -idx_remove, with = FALSE]

# Remove excluded samples
meta_clean <- meta[!(meta$sample_ID %in% sample_exclusions),]


# ============================================================
### SPECIES ABUNDANCE PROCESSING (METAPHLAN)
# ============================================================
# This section processes MetaPhlAn species-level abundances:
#   - filtering low-prevalence/low-abundance taxa
#   - normalization to relative abundances
#   - CLR transformation
#   - averaging across visits/phases

# ------------------------------------------------------------
### FILTERING FUNCTION (Metagenomic taxa)
# ------------------------------------------------------------
# Adapted from:
# https://github.com/GRONINGEN-MICROBIOME-CENTRE/DMP

filterMetaGenomeDF <- function(
  inDF,
  presPerc = 0.1,
  minMRelAb = 0.01,
  minMedRelAb = 0.0,
  rescaleTaxa = FALSE,
  verbose = TRUE,
  keepDomains = c('Bacteria','Archaea'),
  keepLevels  = c('T','S','G','F','O','C','P')
) {

  # ----------------------------------------------------------
  # Replace NA values with 0 (microbiome columns only)
  # ----------------------------------------------------------
  tCols <- grep('k__', colnames(inDF))
  for (c in tCols) {
    inDF[,c][is.na(inDF[,c])] <- 0.0
  }

  # ----------------------------------------------------------
  # Presence filter
  # ----------------------------------------------------------
  nrRemoved <- 0
  toRemove  <- c()

  for (c in tCols) {
    nrnZ <- sum(inDF[,c] != 0.0)
    if ((nrnZ / nrow(inDF)) < presPerc) {
      nrRemoved <- nrRemoved + 1
      toRemove  <- c(toRemove, c)
    }
  }

  if (length(toRemove) > 0) {
    inDF <- inDF[, -toRemove]
  }

  if (verbose) {
    print(paste(' > presence filter: removed', nrRemoved, 'taxa'))
  }

  # ----------------------------------------------------------
  # Mean abundance filter
  # ----------------------------------------------------------
  toRemove <- c()
  for (c in grep('k__', colnames(inDF))) {
    if (mean(inDF[,c]) < minMRelAb) {
      toRemove <- c(toRemove, c)
    }
  }

  if (length(toRemove) > 0) {
    inDF <- inDF[, -toRemove]
  }

  if (verbose) {
    print(paste(' > mean abundance filter applied'))
  }

  # ----------------------------------------------------------
  # Median abundance filter
  # ----------------------------------------------------------
  toRemove <- c()
  for (c in grep('k__', colnames(inDF))) {
    if (median(inDF[,c]) < minMedRelAb) {
      toRemove <- c(toRemove, c)
    }
  }

  if (length(toRemove) > 0) {
    inDF <- inDF[, -toRemove]
  }

  if (verbose) {
    print(paste(' > median abundance filter applied'))
  }

  # ----------------------------------------------------------
  # Keep selected taxonomic levels
  # ----------------------------------------------------------
  inDFnonTaxa <- as.data.frame(inDF[, grep('k__', colnames(inDF), invert = TRUE)])

  taxaS <- inDF[, grep('s__', colnames(inDF)), drop = FALSE]

  if (rescaleTaxa) {
    taxaS <- taxaS / rowSums(taxaS)
  }

  oDF <- cbind(inDFnonTaxa, taxaS)

  if (verbose) print(' > returning filtered dataframe')

  return(oDF)
}

# ------------------------------------------------------------
### LOAD AND PREPROCESS ABUNDANCE DATA
# ------------------------------------------------------------
abundanceFILE <- "MGS_FECAL_merged_abundance_table.txt"

inDF <- read.delim(abundanceFILE, check.names = FALSE)

# Clean column names
colnames(inDF) <- gsub("_metaphlan", "", colnames(inDF))

# Keep bacterial taxa only
inDF <- inDF[(inDF$clade_name %like% "k__Bacteria"), ]

# ------------------------------------------------------------
### TRANSPOSE TO SAMPLE x FEATURES MATRIX
# ------------------------------------------------------------
tx <- as.data.table(t(inDF))
colnames(tx) <- unlist(tx[1,])
tx <- tx[-1,]

# Convert to numeric
tx <- tx %>% mutate_if(is.character, as.numeric)

# Add sample IDs
tx$sample_ID <- colnames(inDF[, -1])

# Remove excluded samples
tx <- tx[!(tx$sample_ID %in% sample_exclusions), ]

# ------------------------------------------------------------
### APPLY FILTERING THRESHOLDS
# ------------------------------------------------------------
inDF <- data.frame(tx, check.names = FALSE)

keep_feat <- colnames(
  filterMetaGenomeDF(
    inDF,
    presPerc = 0.2,
    minMRelAb = 0.1,
    keepDomains = 'Bacteria',
    keepLevels  = 'S',
    rescaleTaxa = TRUE,
    verbose     = FALSE
  )
)

# Clean feature names
keep_feat <- keep_feat[!(keep_feat %in% c("UNCLASSIFIED","sample_ID"))]
keep_feat <- sub(".*(s__[^|]+).*", "\\1", keep_feat)

# ------------------------------------------------------------
### SPECIES-LEVEL DATA CLEANING
# ------------------------------------------------------------
inDF_tmp <- tx

# Keep species only
inDF_tmp <- inDF_tmp[, grep("s__", colnames(inDF_tmp)), with = FALSE]

# Convert to numeric
inDF_tmp <- inDF_tmp %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

# Remove samples with zero total abundance
rows_keep <- rowSums(inDF_tmp, na.rm = TRUE) > 0
inDF_tmp  <- inDF_tmp[rows_keep, ]

sampleIDs <- tx[rows_keep, ]$sample_ID

# Normalize to relative abundance
inDF_tmp <- inDF_tmp / rowSums(inDF_tmp, na.rm = TRUE)

# Reassign IDs
inDF_tmp$sample_ID <- sampleIDs
inDF_tmp <- as.data.table(inDF_tmp)

# Simplify species names
colnames(inDF_tmp) <- sub(".*(s__[^|]+).*", "\\1", colnames(inDF_tmp))

# ------------------------------------------------------------
### MERGE WITH METADATA
# ------------------------------------------------------------
inDFm <- merge(inDF_tmp, meta_clean, by = "sample_ID", all.y = TRUE)

# Identify metadata vs feature columns
colsMeta     <- colnames(inDFm)[grep("s__", colnames(inDFm), invert = TRUE)]
colsFeatures <- colnames(inDFm)[grep("s__", colnames(inDFm))]

# Define factor levels
inDFm$Bristol_stool_scale <- factor(
  inDFm$Bristol_stool_scale,
  levels = c("normal_stools","loose_stools_diarrhea","hard_stools_constipation")
)

# ------------------------------------------------------------
### OUTPUT (RELATIVE ABUNDANCE)
# ------------------------------------------------------------
write.table(inDFm,
  file = "species.MGS_FECAL_merged_abundance_table.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

saveRDS(inDFm,
  file = "species.MGS_FECAL_merged_abundance_table.rds"
)

# ------------------------------------------------------------
### CLR TRANSFORMATION
# ------------------------------------------------------------
x <- data.table(inDFm[, ..colsFeatures])
x[is.na(x)] <- 0

# Add pseudocount
pseudocount <- min(as.matrix(x)[as.matrix(x) > 0]) / 2
x_pc <- x + pseudocount

# Convert to relative abundances
rel_abb <- x_pc / rowSums(x_pc)

# Compute CLR
clr_abb <- compositions::clr(rel_abb)

# Combine with metadata
metabk   <- inDFm[, ..colsMeta]
inDFm_clr <- cbind(metabk, clr_abb)

# ------------------------------------------------------------
### AVERAGE BY SUBJECT AND PHASE
# ------------------------------------------------------------
colsTechcovars <- c("readsDepth_postQC", "DNA_concentration")
cols_to_average <- c(colsTechcovars, colsFeatures)

inDFm_clr_avg <- inDFm_clr %>%
  group_by(across(all_of(c("ID","phase")))) %>%
  mutate(across(all_of(cols_to_average), ~mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  distinct(across(all_of(c("ID","phase"))), .keep_all = TRUE)

inDFm_clr_avg <- data.table(inDFm_clr_avg)

# Remove samples without phase
inDFm_clr_avg <- inDFm_clr_avg[!is.na(inDFm_clr_avg$phase), ]

# Keep selected features
inDFm_clr_avg <- inDFm_clr_avg[, c(..colsMeta, ..keep_feat)]

# Re-apply factor levels
inDFm_clr_avg$Bristol_stool_scale <- factor(
  inDFm_clr_avg$Bristol_stool_scale,
  levels = c("normal_stools","loose_stools_diarrhea","hard_stools_constipation")
)

# ------------------------------------------------------------
### OUTPUT (CLR TRANSFORMED)
# ------------------------------------------------------------
write.table(inDFm_clr_avg,
  file = "species.MGS_FECAL_merged_abundance_table_clr.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

saveRDS(inDFm_clr_avg,
  file = "species.MGS_FECAL_merged_abundance_table_clr.rds"
)



# ============================================================
### PATHWAYS ABUNDANCE PROCESSING (HUMANN / MetaCyc)
# ============================================================
# This section processes HUMAnN pathways:
#   - filtering low-prevalence / low-abundance pathways
#   - normalization to relative abundance
#   - CLR transformation
#   - aggregation at subject-phase level

# ------------------------------------------------------------
### FILTERING FUNCTION (HUMANN pathways)
# ------------------------------------------------------------
# Adapted from:
# https://github.com/GRONINGEN-MICROBIOME-CENTRE/DMP

filterHumannDF <- function(
  inDF,
  presPerc = 0.05,
  minMRelAb = 0.001,
  minMedRelAb = 0.0,
  rescale = TRUE,
  verbose = TRUE,
  type = 'MetaCyc'
) {

  # ----------------------------------------------------------
  # Separate pathways from metadata columns
  # ----------------------------------------------------------
  yesPWYdf <- inDF[, grep('PWY', colnames(inDF)), drop = FALSE]
  nonPWYdf <- inDF[, grep('PWY', colnames(inDF), invert = TRUE), drop = FALSE]

  # Replace NAs with zero
  yesPWYdf[is.na(yesPWYdf)] <- 0

  # ----------------------------------------------------------
  # Normalize to relative abundance
  # ----------------------------------------------------------
  if (rescale) {
    rsums <- rowSums(yesPWYdf)
    rsums[rsums == 0] <- 1
    yesPWYdf <- yesPWYdf / rsums
  }

  # ----------------------------------------------------------
  # Presence filter
  # ----------------------------------------------------------
  keep <- colMeans(yesPWYdf > 0) >= presPerc
  yesPWYdf <- yesPWYdf[, keep, drop = FALSE]

  if (verbose) {
    print(paste(" > presence filter:", sum(keep), "pathways retained"))
  }

  # ----------------------------------------------------------
  # Mean abundance filter
  # ----------------------------------------------------------
  keep <- colMeans(yesPWYdf) >= minMRelAb
  yesPWYdf <- yesPWYdf[, keep, drop = FALSE]

  # ----------------------------------------------------------
  # Median abundance filter
  # ----------------------------------------------------------
  keep <- apply(yesPWYdf, 2, median) >= minMedRelAb
  yesPWYdf <- yesPWYdf[, keep, drop = FALSE]

  if (verbose) {
    print(paste(" > filtering complete:", ncol(yesPWYdf), "pathways remaining"))
  }

  # Recombine metadata + filtered pathways
  return(cbind(nonPWYdf, yesPWYdf))
}

# ------------------------------------------------------------
### LOAD AND PREPROCESS HUMANN OUTPUT
# ------------------------------------------------------------
abundanceFILE <- "MGS_FECAL_humann_community.tsv.gz"

inDF <- read.delim(abundanceFILE, check.names = FALSE)

# Clean column names
colnames(inDF)[1] <- gsub("# ", "", colnames(inDF)[1])
colnames(inDF) <- gsub("_humann_Abundance.RELAB", "", colnames(inDF))

# ------------------------------------------------------------
### TRANSPOSE TO SAMPLE x PATHWAYS MATRIX
# ------------------------------------------------------------
tx <- as.data.table(t(inDF))
colnames(tx) <- unlist(tx[1,])
tx <- tx[-1,]

# Harmonize pathway names
colnames(tx) <- sub("^((?:(?!PWY)[^:|])+)(?=:|$)", "\\1-PWY", colnames(tx), perl = TRUE)
colnames(tx) <- gsub("UNMAPPED-PWY", "UNMAPPED", colnames(tx))
colnames(tx) <- gsub("UNINTEGRATED-PWY", "UNINTEGRATED", colnames(tx))
colnames(tx) <- gsub(" ", "", colnames(tx))

# Convert to numeric
tx <- tx %>% mutate(across(everything(), as.numeric))

# Add sample IDs
tx$sample_ID <- colnames(inDF[, -1])

# Remove excluded samples
tx <- tx[!(tx$sample_ID %in% sample_exclusions), ]

# ------------------------------------------------------------
### CLEAN AND FILTER PATHWAYS
# ------------------------------------------------------------
inDF <- data.frame(tx, check.names = FALSE)

# Remove non-informative features
inDF_tmp <- inDF[, !(names(inDF) %in% c("UNMAPPED", "UNINTEGRATED"))]
inDF_tmp <- inDF_tmp[, grep("\\|.*", colnames(inDF_tmp), invert = TRUE)]

# Apply filtering thresholds
keep_feat <- colnames(
  filterHumannDF(
    inDF_tmp,
    presPerc = 0.2,
    minMRelAb = 0.001,
    minMedRelAb = 0.0,
    rescale = TRUE,
    verbose = FALSE
  )
)

keep_feat <- setdiff(keep_feat, "sample_ID")

# Keep pathways only
inDF_tmp <- inDF_tmp[, grep("PWY", colnames(inDF_tmp)), drop = FALSE]
inDF_tmp <- inDF_tmp %>% mutate(across(everything(), as.numeric))

# Remove empty samples
rows_keep <- rowSums(inDF_tmp, na.rm = TRUE) > 0
inDF_tmp  <- inDF_tmp[rows_keep, ]

sampleIDs <- tx[rows_keep, ]$sample_ID

# Normalize to relative abundance
inDF_tmp <- inDF_tmp / rowSums(inDF_tmp, na.rm = TRUE)

# Assign sample IDs
inDF_tmp$sample_ID <- sampleIDs

# ------------------------------------------------------------
### MERGE WITH METADATA
# ------------------------------------------------------------
inDF_tmp <- as.data.table(inDF_tmp)
inDFm <- merge(inDF_tmp, meta_clean, by = "sample_ID")

# Identify metadata vs pathway columns
colsMeta     <- colnames(inDFm)[grep("PWY", colnames(inDFm), invert = TRUE)]
colsFeatures <- colnames(inDFm)[grep("PWY", colnames(inDFm))]

# ------------------------------------------------------------
### CLR TRANSFORMATION
# ------------------------------------------------------------
x <- data.table(inDFm[, ..colsFeatures])
x[is.na(x)] <- 0

# Add pseudocount
pseudocount <- min(as.matrix(x)[as.matrix(x) > 0]) / 2
x_pc <- x + pseudocount

# Normalize
rel_abb <- x_pc / rowSums(x_pc)

# CLR
clr_abb <- compositions::clr(rel_abb)

# Combine with metadata
metabk   <- inDFm[, ..colsMeta]
inDFm_clr <- cbind(metabk, clr_abb)

# ------------------------------------------------------------
### AVERAGE BY SUBJECT AND PHASE
# ------------------------------------------------------------
colsTechcovars <- c("readsDepth_postQC", "DNA_concentration")
cols_to_average <- c(colsTechcovars, colsFeatures)

inDFm_clr_avg <- inDFm_clr %>%
  group_by(across(all_of(c("ID", "phase")))) %>%
  mutate(across(all_of(cols_to_average), ~mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  distinct(across(all_of(c("ID", "phase"))), .keep_all = TRUE)

inDFm_clr_avg <- data.table(inDFm_clr_avg)

# Remove samples without phase
inDFm_clr_avg <- inDFm_clr_avg[!is.na(inDFm_clr_avg$phase), ]

# Keep selected pathways
inDFm_clr_avg <- inDFm_clr_avg[, c(..colsMeta, ..keep_feat)]

# ------------------------------------------------------------
### FACTOR FORMATTING
# ------------------------------------------------------------
inDFm_clr_avg$Bristol_stool_scale <- factor(
  inDFm_clr_avg$Bristol_stool_scale,
  levels = c("normal_stools","loose_stools_diarrhea","hard_stools_constipation")
)

inDFm_clr_avg$batch  <- factor(inDFm_clr_avg$batch)
inDFm_clr_avg$center <- factor(inDFm_clr_avg$center)

inDFm_clr_avg$Collection_feces_daytime <- factor(
  inDFm_clr_avg$Collection_feces_daytime,
  levels = c("morning","not_morning")
)

# ------------------------------------------------------------
### OUTPUT (CLR PATHWAYS)
# ------------------------------------------------------------
write.table(inDFm_clr_avg,
  file = "pathways.MGS_FECAL_merged_abundance_table_clr.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

saveRDS(inDFm_clr_avg,
  file = "pathways.MGS_FECAL_merged_abundance_table_clr.rds"
)


# ============================================================
### GENE FAMILIES PROCESSING (ESTRABOLOME / HUMANN)
# ============================================================
# This section processes HUMAnN gene families (UniRef90):
#   - cleaning unwanted features
#   - normalization to relative abundance
#   - prevalence filtering
#   - CLR transformation
#   - aggregation at subject-phase level

# ------------------------------------------------------------
### LOAD AND CLEAN INPUT DATA
# ------------------------------------------------------------
abundanceFILE <- "MGS_FECAL_humann_genefamilies_uniref90_level4ec.community.tsv.gz"

inDF <- read.delim(abundanceFILE, check.names = FALSE)

# Clean column names
colnames(inDF)[1] <- gsub("# ", "", colnames(inDF)[1])
colnames(inDF)   <- gsub("_humann_Abundance.RELAB", "", colnames(inDF))

# ------------------------------------------------------------
### REMOVE UNWANTED FEATURES
# ------------------------------------------------------------
# Remove taxonomy-resolved entries and unmapped signals
inDF <- inDF[-grep("\\|.*", inDF$`Gene Family`), ]
inDF <- inDF[!inDF$`Gene Family` %in% c("UNMAPPED","UNGROUPED"), ]

# Add prefix to feature names
inDF$`Gene Family` <- paste0("genefam__", inDF$`Gene Family`)

# ------------------------------------------------------------
### TRANSPOSE TO SAMPLE x FEATURES MATRIX
# ------------------------------------------------------------
tx <- as.data.table(t(inDF))
colnames(tx) <- unlist(tx[1,])
tx <- tx[-1,]

# Convert to numeric
tx <- tx %>% mutate(across(everything(), as.numeric))

# Add sample IDs
tx$sample_ID <- colnames(inDF[, -1])

# ------------------------------------------------------------
### REMOVE EMPTY SAMPLES
# ------------------------------------------------------------
inDF_tmp <- as.data.frame(tx[, -ncol(tx)])

rows_keep <- rowSums(inDF_tmp, na.rm = TRUE) > 0
inDF_tmp  <- inDF_tmp[rows_keep, ]

sampleIDs <- tx[rows_keep, ]$sample_ID

# ------------------------------------------------------------
### NORMALIZE TO RELATIVE ABUNDANCE
# ------------------------------------------------------------
inDF_tmp <- inDF_tmp / rowSums(inDF_tmp, na.rm = TRUE)
inDF_tmp$sample_ID <- sampleIDs
inDF_tmp <- as.data.table(inDF_tmp)

# ------------------------------------------------------------
### MERGE WITH METADATA
# ------------------------------------------------------------
inDFm <- merge(inDF_tmp, meta, by = "sample_ID")

# Rename selected enzymatic activities (readability)
setnames(
  inDFm,
  old = c(
    "genefam__2.8.2.22","genefam__3.1.6.1","genefam__3.2.1.21",
    "genefam__3.2.1.23","genefam__3.2.1.31","genefam__3.2.1.86"
  ),
  new = c(
    "genefam__Aryl_sulfate_sulfotransferase",
    "genefam__Arylsulfatase_type_I",
    "genefam__Beta_glucosidase",
    "genefam__Beta_galactosidase",
    "genefam__Beta_glucuronidase",
    "genefam__6_phospho_beta_glucosidase"
  )
)

# ------------------------------------------------------------
### DEFINE FEATURE AND METADATA COLUMNS
# ------------------------------------------------------------
colsFeatures <- colnames(inDFm)[grep("genefam__", colnames(inDFm))]
colsMeta     <- colnames(inDFm)[grep("genefam__", colnames(inDFm), invert = TRUE)]

# ------------------------------------------------------------
### PREVALENCE FILTERING
# ------------------------------------------------------------
feature_mat  <- as.matrix(inDFm[, ..colsFeatures])
presence_mat <- feature_mat > 0

prevalence <- colMeans(presence_mat, na.rm = TRUE)

# Keep features present in ≥20% of samples
colsFeatures_prev <- names(prevalence[prevalence >= 0.2])

# Subset dataset
inDFm_sub <- inDFm[, c(colsMeta, colsFeatures_prev), with = FALSE]

# ------------------------------------------------------------
### OUTPUT (RELATIVE ABUNDANCE)
# ------------------------------------------------------------
write.table(
  inDFm_sub,
  file = "genefam.MGS_FECAL_merged_abundance_table.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

saveRDS(
  inDFm_sub,
  file = "genefam.MGS_FECAL_merged_abundance_table.rds"
)

# ------------------------------------------------------------
### CLR TRANSFORMATION
# ------------------------------------------------------------
x <- data.table(inDFm[, ..colsFeatures])
x[is.na(x)] <- 0

# Add pseudocount
pseudocount <- min(as.matrix(x)[as.matrix(x) > 0]) / 2
x_pc <- x + pseudocount

# Normalize
rel_abb <- x_pc / rowSums(x_pc)

# Compute CLR
clr_abb <- compositions::clr(rel_abb)

# Combine metadata + CLR matrix
metabk   <- inDFm[, ..colsMeta]
inDFm_clr <- cbind(metabk, clr_abb)

# Keep only selected features
colsFeatures <- colnames(inDFm_sub)[grep("genefam__", colnames(inDFm_sub))]
colsMeta     <- colnames(inDFm_sub)[grep("genefam__", colnames(inDFm_sub), invert = TRUE)]

inDFm_clr <- inDFm_clr[, c(..colsMeta, ..colsFeatures)]

# ------------------------------------------------------------
### AVERAGE BY SUBJECT AND PHASE
# ------------------------------------------------------------
colsTechcovars  <- c("readsDepth_postQC","DNA_concentration")
cols_to_average <- c(colsTechcovars, colsFeatures)

inDFm_clr_avg <- inDFm_clr %>%
  group_by(across(all_of(c("ID","phase")))) %>%
  mutate(across(all_of(cols_to_average), ~mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  distinct(across(all_of(c("ID","phase"))), .keep_all = TRUE)

inDFm_clr_avg <- data.table(inDFm_clr_avg)

# Remove samples without phase annotation
inDFm_clr_avg <- inDFm_clr_avg[!is.na(inDFm_clr_avg$phase), ]

# ------------------------------------------------------------
### FACTOR FORMATTING
# ------------------------------------------------------------
inDFm_clr_avg$Bristol_stool_scale <- factor(
  inDFm_clr_avg$Bristol_stool_scale,
  levels = c("normal_stools","loose_stools_diarrhea","hard_stools_constipation")
)

inDFm_clr_avg$batch  <- factor(inDFm_clr_avg$batch)
inDFm_clr_avg$center <- factor(inDFm_clr_avg$center)

inDFm_clr_avg$Collection_feces_daytime <- factor(
  inDFm_clr_avg$Collection_feces_daytime,
  levels = c("morning","not_morning")
)

# ------------------------------------------------------------
### OUTPUT (CLR TRANSFORMED)
# ------------------------------------------------------------
write.table(
  inDFm_clr_avg,
  file = "genefam.MGS_FECAL_merged_abundance_table_clr.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

saveRDS(
  inDFm_clr_avg,
  file = "genefam.MGS_FECAL_merged_abundance_table_clr.rds"
)


