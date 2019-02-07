#'---
#' title: "RNA seq workflow part 1 - Counts to .dds"
#' subtitle: "HIOs PMN RNAseq"
#' author: "Ryan Berger"
#' date: "2019-01-21"
#' output: 
#'   html_document:
#'      theme: flatly
#'      highlight: tango
#'      toc: true
#'      number_sections: true
#'      toc_depth: 2	
#'      toc_float:	
#'       collapsed: false	
#'       smooth_scroll: true	  
#' ---


#' # Purpose
#' This script is the first in a series of RNA seq processing scripts for converting transcript counts from a kallisto alignment into gene counts used in analyses. This is modified from a workflow by [Bioconductor](https://www.bioconductor.org/help/workflows/rnaseqGene/).

#+ section2
#' # Summary
#' **Input:**
#' 
#' * RNA seq alignment results (kallisto output)
#'     + abundance.h5 
#' * Sample key for annotating files
#'     + sample_key.csv
#' * Gene reference file (tx2gene) for converting transcript counts to genes
#'     + EnsDb.Hsapiens.v75 
#'   
#' **Output:**
#' 
#' * DESeq2 dds object

#------------------------------------------------------------
#' # Begin script

#' ## Libraries
library(DESeq2)
library(dplyr)
library(EnsDb.Hsapiens.v75)
library(here)
library(knitr)
library(rmarkdown)
library(tximport)


# Optional: install packages from Bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite("tximport")
# biocLite("rhdf5")
# biocLite("DESeq2")
# biocLite("EnsDb.Hsapiens.v75")


#+ set_cache_dir, include=F
# Set a directory for all cache files to go
dir.create(here('results/DESeq2/src_html_output/knitr_cache'), recursive = T)
opts_chunk$set(cache.path = here('results/DESeq2/src_html_output/knitr_cache/'))

#' ## Set up directories
#' We'll use the `here` package to control our working directory. This will set the root directory to `'../HIO_dualseq2/'`.

here()

# Where are the input (kallisto alignment) files located?
dir <- here('data/Run1822/kallisto-Run_1822/')


# Where should the output files go?
output_folderID <- 'DESeq2'  # Name of output folder
dir.create(here('results', output_folderID))

results.dir <- here('results/DESeq2/')


#------------------------------------------------------------
#' ## Transcript counts to gene counts

#' Need to first set up a gene reference object (tx2gene) for converting counts of ENSEMBL transcipt IDs to counts of gene IDs. We'll use the Ensembl based annotation package `EnsDb.Hsapiens.v75` to set this up.


edb <- EnsDb.Hsapiens.v75

# Create a dataframe of transcripts from ensembl human genome
Tx <- transcripts(edb,
                  columns = c(listColumns(edb , "tx"), "gene_name"),
                  return.type = "DataFrame")

# assign columns 1 (transcript ID) and 9 (gene name) to dataframe tx2gene
tx2gene <- Tx[,c('tx_id','gene_name')]



#' ## File import
#' Import the `kallisto` alignment output files (abundance.h5) and annotate them with a sample key. These files are located in the results directory.

# Load in sample key.
sample_key <- read.csv(here('data/sample_key.csv'))


#' This RNA seq run contains 12 samples that were from a different experiment (the first dual-seq experiment). The column `dualseq_expt` indicates which experiment each sample is from so we'll use this to filter these samples out of our analysis.

# Optional: filter out samples from the key that you don't want to analyze
sample_key <- dplyr::filter(sample_key, Date == '5/22/17' & hr == 8) 

# Vector of samples used in the tximport
txiSamples <- as.vector(unique(sample_key$code_name))
txiSamples

# Set up path to read files into tximport
files <- file.path(dir, sample_key$file_name, 'abundance.h5')
# Add sample IDs to files
names(files) <- sample_key$short_name


# import the abundance.h5 files using tximport

#+ tximport, cache=T
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)


# export abundance counts
write.csv(txi$abundance, file = file.path(results.dir, "complete_dataset_txi.csv"))


#' ## Create a DESeq2 object (dds)

#' Need to set the design of dds object to the experimental condition of each sample without replicate IDs. In our case it is the `code_name` from the `sample_key`. 

# Create DESeq dataset
dds <- DESeqDataSetFromTximport(txi, 
                                colData = sample_key,
                                design = ~code_name)
head(dds)


#' Now we'll eliminate all the genes that have zero counts across all conditions.

nrow(dds)
# Filter out rows with no counts
dds <- dds[rowSums(counts(dds)) > 1, ]
nrow(dds) 


# Account for transcript length
dds <- DESeq2::estimateSizeFactors(dds)
ddscounts <- DESeq2::counts(dds, normalized = TRUE)

# Save output file:
write.csv(ddscounts, file = file.path(results.dir, "complete-dataset_DESeq2-normalized-counts.csv"))

# Save dds object
saveRDS(dds, file = file.path(results.dir, 'dds_all.rds'))

#' ## Ending notes
#' RNA seq reads are now saved in a DESeq dataset (dds object). This can be loaded for downstream analysis.







#' # Remove PMN control samples
#' There are two PMN controls for the experiment: PMN RNA and a HIO:PMN 3:1 RNA mixture. These are used to compare the transcriptome similarity of PMNs vs HIOs. We'll remove these samples and save as a separate dds object since most of our analyses won't use these samples and their profile is so different it's recommended to remove them.


# Make separate sample_keys for each subset 
alt_key <- dplyr::filter(sample_key, code_name %in% c('PBS','STM','SE',
                                                   'PBS+PMNs', 'STM+PMNs', 'SE+PMNs'))


# Vector of samples used in the tximport
txiSamples <- as.vector(unique(alt_key$code_name))
txiSamples
  
# Set up path to read files into tximport
files <- file.path(dir, alt_key$file_name, 'abundance.h5')
# Add sample IDs to files
names(files) <- alt_key$short_name
  
# import the abundance.h5 files using tximport
  
#+ tximport, cache=T
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  
  
# export abundance counts
write.csv(txi$abundance, file = file.path(results.dir, "no_controls_txi.csv"))
  
  
  #' ## Create a DESeq2 object (dds)
  
  #' Need to set the design of dds object to the experimental condition of each sample without replicate IDs. In our case it is the `code_name` from the `sample_key`. 
  
# Create DESeq dataset
dds <- DESeqDataSetFromTximport(txi, 
                                  colData = alt_key,
                                  design = ~code_name)
  
  
#' Now we'll eliminate all the genes that have zero counts across all conditions.
# Filter out rows with no counts
dds <- dds[rowSums(counts(dds)) > 1, ]
  
  
# Account for transcript length
dds <- DESeq2::estimateSizeFactors(dds)
ddscounts <- DESeq2::counts(dds, normalized = TRUE)
  
# Save output file:
write.csv(ddscounts, file = file.path(results.dir, "no_controls_DESeq2-normalized-counts.csv"))
  
# Save dds object
saveRDS(dds, file = file.path(results.dir, 'dds_no_controls.rds'))



#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/01_counts_to_dds.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)