#!/bin/bash

# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================


## Performs functional profiling of cleaned reads using HUMAnN,
## integrating MetaPhlAn taxonomic profiles, and generating merged functional tables ##

# ------------------------------------------------------------
### INPUT PARAMETERS
# ------------------------------------------------------------
# Space-separated list of sample IDs
samplesVAR=""

# Directories
kneadmainDIR=""   # KneadData output directory (cleaned reads)
metapDIR=""       # MetaPhlAn output directory
TOOLSPATH=""      # path to external tools (e.g., usearch)
outDIR=""         # HUMAnN output directory
mergedOUTdir=""   # directory for merged tables

# ------------------------------------------------------------
### HUMANN FUNCTIONAL PROFILING
# ------------------------------------------------------------
# For each sample:
#   - perform functional profiling from cleaned reads
#   - integrate taxonomic profile from MetaPhlAn
#   - generate gene family abundance tables

for s in $samplesVAR; do

  # Create output directory
  mkdir -p "$outDIR/$s/"

  # Run HUMAnN
  humann \
    --input "$kneadmainDIR/$s/${s}_kneaddata.fastq.gz" \
    --input-format fastq.gz \
    --taxonomic-profile "$metapDIR/$s/${s}_metaphlan.txt" \
    --usearch "$TOOLSPATH/usearch6.1.544_i86linux32" \
    --output "$outDIR/$s/" \
    --output-basename "${s}_humann" \
    --o-log "$outDIR/$s/${s}_humann.log" \
    --threads 5 &

  # Small delay to avoid system overload
  sleep 10

done

# ------------------------------------------------------------
### PREPARE MERGED OUTPUT DIRECTORY
# ------------------------------------------------------------
mkdir -p "$mergedOUTdir/temp/"

# ------------------------------------------------------------
### REGROUP GENE FAMILIES (UniRef90 → EC LEVEL 4)
# ------------------------------------------------------------
# Convert gene families into functional categories (EC level 4)

for id in $samplesVAR; do

  humann_regroup_table \
    --input "$outDIR/$id/${id}_humann_genefamilies.tsv" \
    --groups uniref90_level4ec \
    --output "$outDIR/$id/${id}_humann_genefamilies.uniref90_level4ec.tsv"

done

# ------------------------------------------------------------
### RENORMALIZATION (RELATIVE ABUNDANCE)
# ------------------------------------------------------------
# Convert raw counts into relative abundances (community level)

for id in $samplesVAR; do

  humann_renorm_table \
    --input "$outDIR/$id/${id}_humann_genefamilies.uniref90_level4ec.tsv" \
    --units relab \
    --mode community \
    --update-snames \
    --output "$mergedOUTdir/temp/${id}_humann_genefamilies.uniref90_level4ec.community.tsv"

done

# ------------------------------------------------------------
### MERGE TABLES ACROSS SAMPLES
# ------------------------------------------------------------
# Combine all samples into a single abundance matrix

humann_join_tables \
  --input "$mergedOUTdir/temp/" \
  --file_name "level4ec.community" \
  --output "$mergedOUTdir/MGS_FECAL_humann_genefamilies_uniref90_level4ec.community.tsv"
