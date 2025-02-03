ReadME for pipeline/scripts used for vaginal microbiome data analysis through 16S method in the paper Vinerbi et al 2025 (in preparation)

***for questions please write to Elena Vinerbi (elenavinerbi@cnr.it), Fabio Chillotti (fabiochillotti@cnr.it) and Serena Sanna (serena.sanna@cnr.it)***

Date released: 03/02/2025

Running on R version 4.4.1 and Python version 3.10.12

The pipeline is divided into several folders. Please follow the numbering of the folders and internal scripts to follow the workflow.

Folder: 1.QC analysis
Contains the scripts to perform quality checks of reads, through the tools FastQC; Trimmomatic; MultiQC.

Folder: 2.Microbiome composition
Contains the scripts for obtaining taxonomic classification by means of Dada2, and the various analyses, calculation of relative abundance, alpha and beta diversity and associated statistical analysis.

Folder: 3.Comparison 16S and WGS
Contains scripts to select samples sequenced with the 16S method to repeat with the WGS method to compare results.

Folder: 4.CST definition
Contains the scripts to get the CST (and subCST) from the VALENCIA algorithm.

Folder: 5.Subsamplig analysis
Contains scripts to subsample data.

Folder: 6.ITS analysis
contains scripts for ITS analysis and classification.

Folder: 7. statistical_analyses_pipeline
contains the pipeline to run the statistical analyses. There is a specific ReadME for this part, in to the folder.

Folder: Public data
contains the additional files needed for the analysis. sub-folder VALENCIA file you find the files to run VALENCIA algorithm. sub-folder DATABSE shows the different databases used (further subfolders).
