#!/bin/bash

# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================


## Performs trimming, quality control, and contaminant removal
## of paired-end reads using KneadData, and summarizes read counts ##

# ------------------------------------------------------------
### INPUT PARAMETERS
# ------------------------------------------------------------
# Space-separated list of sample IDs
samplesLIST=""

# Directories
fastqDIR=""   # path to raw FASTQ files (per sample subfolders)
refPATH=""    # path to KneadData reference database
outDIR=""     # output directory

# ------------------------------------------------------------
### KNEADDATA PROCESSING
# ------------------------------------------------------------
# For each sample:
#   - perform trimming (Trimmomatic)
#   - remove contaminants (reference DB)
#   - run QC (FastQC before/after)
#   - generate cleaned paired reads

for s in $samplesLIST; do

  # Create output folder
  mkdir -p "$outDIR/$s/"

  # Run KneadData
  kneaddata \
    --input1 "$fastqDIR/$s/${s}_R1.fastq.gz" \
    --input2 "$fastqDIR/$s/${s}_R2.fastq.gz" \
    --reference-db "$refPATH" \
    --output "$outDIR/$s/" \
    --output-prefix "${s}_kneaddata" \
    --trimmomatic ./TOOLS/Trimmomatic-0.39/ \
    --run-trf \
    --run-fastqc-start \
    --run-fastqc-end \
    --cat-final-output \
    --log "$outDIR/$s/${s}_kneaddata.log" &

  # Small delay to avoid overloading system
  sleep 10

done

# ------------------------------------------------------------
### READ COUNT SUMMARY (POST-QC)
# ------------------------------------------------------------
# Counts number of reads from KneadData paired output files

kneadDIR="$outDIR"
outFILE="reads.summary"

# Initialize output file
echo -e "sample_ID\ttotal_reads_postQC" > "$outFILE"

# Loop over processed FASTQ files
for f in $kneadDIR/*/*_kneaddata_paired.fastq.gz; do

  # Extract sample name
  sample=$(basename "$f" | sed 's/_kneaddata_paired.fastq.gz//')

  # Count reads:
  # FASTQ format: 4 lines per read and divide total lines by 4
  reads=$(zcat "$f" | wc -l)
  reads=$((reads / 4))

  # Append to summary file
  echo -e "$sample\t$reads" >> "$outFILE"

done
