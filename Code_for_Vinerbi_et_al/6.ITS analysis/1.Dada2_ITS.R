# Rscript

# Script to use a Dada2 package for ITS reads

# Note: # we used a UNITE-DB
# The UNITE reference database from  2024-04-22 (version 10.0) was downloaded from https://dx.doi.org/10.15156/BIO/2959330/UNITE_public_21.04.2024.fasta.gz
# TO FORMAT TO DADA2 SOFTWARE USE THE FOLLOWING:
#zless UNITE_public_21.04.2024.fasta.gz| sed 's/>[^|]*|k__\(.*\);p__\(.*\);c__\(.*\);o__\(.*\);f__\(.*\);g__\(.*\);s__\(.*\)|.*/>\1;\2;\3;\4;\5;\6;\7/'  |gzip > UNITE_INS_BIO2959330_dada2format_fasta.gz


# Author: Serena Sanna
# Modified by: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: v4.4.2,
 
# update library
library(tidyverse)
library(dada2)

# path
setwd("/home/")
getwd()

# Path of files
project.folder <- "/home/"
R1.tag <- "_R1_paired.fq.gz"
R2.tag <- "_R2_paired.fq.gz"

path <- list.files(project.folder)
list.files(project.folder)

fnFs = paste(sep = "/", project.folder, grep(R1.tag, path, value = TRUE))
fnRs = paste(sep = "/", project.folder, grep(R2.tag, path, value = TRUE))
head(fnFs)
head(fnRs)

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
length(sample.names)

filtFs <- file.path(project.folder, "filtered_ITS", paste0(sample.names, "_R1_filt.fq.gz"))
filtRs <- file.path(project.folder, "filtered_ITS", paste0(sample.names, "_R2_filt.fq.gz"))
names(filtFs) = sample.names
names(filtRs) = sample.names
head(filtFs)
head(filtRs)

# Filer and trim (for ITS change parametres)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen=50, maxLen=200,  # lenght of reads
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

head(out)
dim(out)

## LEARN ERROR RATES ##
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

## Pool=TRUE ##
dadaFs <- dada(filtFs, err = errF, pool = TRUE, multithread = TRUE)
dadaFs[[1]]

dadaRs <- dada(filtRs, err = errR, pool = TRUE, multithread = TRUE)
dadaRs[[1]]

### MERGE DENOISED PAIRED READS ###
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

## TABLE OF SEQUENCE ##
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

## REMOVE CHIMERE ##
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)

## NUMBER OF LECTURE ##
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track)

## Assign taxonomy 
# This file 
taxa <- assignTaxonomy(seqtab.nochim,"Public_data/DATABASE/UniteDB/UNITE_INS_BIO2959330_dada2format_fasta.gz", multithread = TRUE, tryRC = TRUE) 
head(taxa)

## Structure of the generated dataframe 'taxa':
# - Columns = Taxonomy levels (Kingdom, Phylum, Class, Order, Family, Genus, Species)
# - Rows = Sequence variants

## Save output## 
write.table(track, "dada2_statistics_ITS.txt",sep="\t", row.names = T)
write.table(track, "dada2_statistics_ITS.csv",sep="\t", row.names = T)

write.table(seqtab.nochim, "ASV_table_ITS.txt", sep="\t", row.names = T)
write.table(seqtab.nochim, "ASV_table_ITS.csv", sep="\t", row.names = T)

write.table(taxa,  "Taxa_table_ITS.txt" ,sep="\t", row.names = T)
write.table(taxa,  "Taxa_table_ITS.csv" ,sep="\t", row.names = T)





