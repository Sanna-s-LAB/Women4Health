# Rscript

# Script to use a Dada2 package to obtain an ASVs table and table with the recognized taxa

# Author: Serena Sanna
# Modified by: Elena Vinerbi (elenavinerbi@cnr.it)
# Last update: 03/02/2025

# R version: R 4.4.1

# Upload library
library(tidyverse)
library(dada2)

setwd('/home/')

# Path of files
project.folder <- "FastQ_16S/"
R1.tag <- "R1_paired.fq.gz"
R2.tag <- "R2_paired.fq.gz"

path <- list.files(project.folder)
list.files(project.folder)

fnFs = paste(sep = "/", project.folder, grep(R1.tag, path, value = TRUE))
fnRs = paste(sep = "/", project.folder, grep(R2.tag, path, value = TRUE))
head(fnFs)
head(fnRs)

length(fnFs) == length(fnRs)

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
length(sample.names)

filtFs <- file.path(project.folder, "filtered", paste0(sample.names, "_R1_filt.fq.gz"))
filtRs <- file.path(project.folder, "filtered", paste0(sample.names, "_R2_filt.fq.gz"))
names(filtFs) = sample.names
names(filtRs) = sample.names
head(filtFs)
head(filtRs)

# Filter data 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen = 160,     # with 160 we removed a very short sequence
                     maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
head(out)

# Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Denoising (pool = TRUE)
dadaFs <- dada(filtFs, err = errF, pool = TRUE, multithread = TRUE)
dadaFs[[1]]

dadaRs <- dada(filtRs, err = errR, pool = TRUE, multithread = TRUE)
dadaRs[[1]]

# Merge denoised paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Table of sequence
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim) # row = samples, column = ASV after removal of chimeras

# Number of lecture
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track) # track: row = samples, column = counts of reads in several steps

# Assigning the taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/DATABASE/GTDB/GTDBr95-Genus.fna", multithread = TRUE, tryRC = TRUE)
taxa <- addSpecies(taxa, "/DATABASE/GTDB/GTDBr95-Species.fna", tryRC = TRUE)
head(taxa) # row = taxon, column = taxonomic levels


#Save output
write.table(track, "dada2_statistics.txt",sep="\t", row.names = T) # table with number of reads for each steps and for each samples
write.table(track, "dada2_statistics.csv",sep="\t", row.names = T)

write.table(seqtab.nochim, "ASV_table.txt", sep="\t", row.names = T) # table with counts for each ASVs (column) for each samples (row)
write.table(seqtab.nochim, "ASV_table.csv", sep="\t", row.names = T)

write.table(taxa,  "Table_taxa.txt" ,sep="\t", row.names = T) # table with taxa identify for each ASV (column = taxonomic level; row = taxon)
write.table(taxa,  "Table_taxa.csv" ,sep="\t", row.names = T)

