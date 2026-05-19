#!/bin/bash

# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================


## Performs taxonomic profiling of cleaned reads using MetaPhlAn ##

# ------------------------------------------------------------
### INPUT PARAMETERS
# ------------------------------------------------------------
# Space-separated list of sample IDs
samplesVAR=""

# Directories
kneadmainDIR=""   # path to KneadData output (cleaned reads)
metapDB=""        # path to MetaPhlAn database
TOOLSPATH=""      # path to additional tools (if needed)
outDIR=""         # output directory

# ------------------------------------------------------------
### METAPHLAN TAXONOMIC PROFILING
# ------------------------------------------------------------
# For each sample:
#   - use cleaned reads from KneadData
#   - perform taxonomic profiling with MetaPhlAn
#   - generate abundance tables and Bowtie2 mappings

for s in $samplesVAR; do

  # Create output and temporary directories
  mkdir -p "$outDIR/$s/temp/"

  # Run MetaPhlAn
  zcat "$kneadmainDIR/$s/${s}_kneaddata.fastq.gz" | \
  metaphlan \
    --input_type fastq \
    --sample_id "$s" \
    --bowtie2db "$metapDB" \
    --tmp_dir "$outDIR/$s/temp/" \
    --add_viruses \
    --unclassified_estimation \
    --profile_vsc \
    --vsc_out "$outDIR/$s/${s}_metaphlan.vsc.txt" \
    --bowtie2out "$outDIR/$s/${s}_metaphlan.bowtie2" \
    --nproc 5 \
    -o "$outDIR/$s/${s}_metaphlan.txt" &

  # Small delay to avoid system overload
  sleep 10

done


