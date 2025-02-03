#!/bin/bash

# This script is for QC analysis

# Author: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# Tools: FastQc v 0.11.8 ; MultiQC v 1.25.2 ; Trimmomatic v0.39

# Unzip a file (if needed)
unzip file.zip -d path/file/

# Use a FastQC tool 
fastqc -t 2 -o FastQC_plate/ File_fastQ/*fastq.gz 

# Use a MultiQC tool
multiqc FastQC_plate/*_L001_R1* -o MultiQC_plate  -n multiqc_R1
multiqc FastQC_plate/*_L001_R2* -o MultiQC_plate  -n multiqc_R2

# Trimmming proccess with Trimmomatic
# V3V4 region
for I in {1..96}; do # change the range with the number of samples available
      SAMPLE_NAME=S${I}
      java -jar /home/TOOLS/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
          -phred33 \
          /home/file_fastQ/NGS*_${SAMPLE_NAME}_V3V4_L001_R1_001.fastq.gz  /home/file_fastQ/NGS*_${SAMPLE_NAME}_V3V4_L001_R2_001.fastq.gz \
          /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_V3V4_R1_paired.fq.gz  /home/file_fastQ/File_post_trimmin/${SAMPLE_NAME}_Plate_V3V4_R1_unpaired.fq.gz \
          /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_V3V4_R2_paired.fq.gz  /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_V3V4_R2_unpaired.fq.gz \
          SLIDINGWINDOW:4:25 MINLEN:50  2>/home/Quality_trimming/QT_Plate/FILE_ERR_Plate/file_${SAMPLE_NAME}_Plate_V3V4.err
done

# Merge and save ERR.file 
cat FILE_ERR/file_S*_Plate_V3V4.err > trimming.quality_V3V4_Plate.err

# V7V9 region 
for I in {1..96}; do
      SAMPLE_NAME=S${I}
      java -jar /home/TOOLS/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
          -phred33 \
          /home/file_fastQ/NGS*_${SAMPLE_NAME}_V7V9_L001_R1_001.fastq.gz   /home/file_fastQ/NGS*_${SAMPLE_NAME}_V7V9_L001_R2_001.fastq.gz \
          /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_V7V9_R1_paired.fq.gz  /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_V7V9_R1_unpaired.fq.gz \
          /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_V7V9_R2_paired.fq.gz  /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_V7V9_R2_unpaired.fq.gz \
          SLIDINGWINDOW:4:25 MINLEN:50  2>/home/Quality_trimming/QT_Plate/FILE_ERR_Plate/file_${SAMPLE_NAME}_Plate_V7V9.err
done 

# Merge and save ERR.file 
cat FILE_ERR/file_S*_Plate_V7V9.err > trimming.quality_V7V9_Plate.err

# ITS region
for I in {1..96}; do
      SAMPLE_NAME=S${I}
      java -jar /home/TOOLS/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
          -phred33 \
          /home/file_fastQ/NGS*_${SAMPLE_NAME}_ITS_L001_R1_001.fastq.gz /home/file_fastQ/NGS*_${SAMPLE_NAME}_ITS_L001_R2_001.fastq.gz \
          /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_ITS_R1_paired.fq.gz  /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_ITS_R1_unpaired.fq.gz \
          /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_ITS_R2_paired.fq.gz /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_ITS_R2_unpaired.fq.gz \
          SLIDINGWINDOW:4:25 MINLEN:50 2>/home/Quality_trimming/QT_Plate/FILE_ERR_Plate/file_${SAMPLE_NAME}_Plate_ITS.err
done 

## Merge and save ERR.file 
cat FILE_ERR_Plate/file_S*_Plate_ITS.err > trimming.quality_ITS_Plate.err

# unknown region
for I in {1..96}; do
      SAMPLE_NAME=S${I}
      java -jar /home/TOOLS/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
          -phred33 \
          /home/file_fastQ/NGS*_${SAMPLE_NAME}_unknown_L001_R1_001.fastq.gz /home/file_fastQ/NGS*_${SAMPLE_NAME}_unknown_L001_R2_001.fastq.gz \
          /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_unknown_R1_paired.fq.gz  /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_unknown_R1_unpaired.fq.gz \
          /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_unknown_R2_paired.fq.gz /home/file_fastQ/File_post_trimming/${SAMPLE_NAME}_Plate_unknown_R2_unpaired.fq.gz \
          SLIDINGWINDOW:4:25 MINLEN:50  2>/home/Quality_trimming/QT_Plate/FILE_ERR_Plate/file_${SAMPLE_NAME}_Plate_unknown.err
done 

# Merge and save ERR.file 
cat FILE_ERR_Plat/file_S*_Plate_unknown.err > trimming.quality_unknown_Plate.err


# Merge V3V4-V7V9 region after trimming ###
for I in {1..96}; do
      SAMPLE_NAME=S${I}
     cat ${SAMPLE_NAME}_Plate_V3V4_R1_paired.fq.gz ${SAMPLE_NAME}_Plate_V7V9_R1_paired.fq.gz >>  V3V4V7V9_Plate/${SAMPLE_NAME}_Plate_V3V4V7V9_R1_paired.fq.gz
     cat ${SAMPLE_NAME}_Plate_V3V4_R2_paired.fq.gz ${SAMPLE_NAME}_Plate_V7V9_R2_paired.fq.gz >>  V3V4V7V9_Plate/${SAMPLE_NAME}_Plate_V3V4V7V9_R2_paired.fq.gz
done 

# Repeate analysis QC analysis to observ if trimming proccess has been successful
# FastQC tool 
fastqc -t 2 -o FastQC_plate/ File_fastQ/*fastq.gz 

# MultiQC tool
multiqc FastQC_plate/*_R1_paired* -o MultiQC_plate  -n multiqc_R1
multiqc FastQC_plate/*_R2_paired* -o MultiQC_plate  -n multiqc_R2
