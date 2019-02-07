#'---
#' title: "RNA seq workflow part 3 - Differential expression"
#' subtitle: "HIOs PMN RNAseq"
#' author: "Ryan Berger"
#' date: "2019-01-23"
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
#' This script is designed to calculate differential gene expression between two samples. The input is a DESeq2 `dds` object and the output is a `.csv` file for each pair of samples compared.

#+ 
#' # Begin script
#' 
#' ## Libraries and directories
library(AnnotationDbi)
library(DESeq2)
library(dplyr)
library(here)
library(knitr)
library(magrittr)
library(org.Hs.eg.db)
library(rmarkdown)

## Optional: install org.Hs.eg.db
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db") 

# Create directory for diff expression .csv files
dir.create(here('results/DESeq2/diff_expression'))

#+ set_cache_dir, include=F, warning=F
# Set a directory for all cache files to go
dir.create(here('results/DESeq2/src_html_output/knitr_cache'))
opts_chunk$set(cache.path = here('results/DESeq2/src_html_output/knitr_cache/'))


#' ## Load data
#' We previously saved the `dds` object in the results folder. Let's load it now:

#+ prepdata, cache=T
dds <- readRDS(here('results/DESeq2/dds_no_controls.Rds'))

# Prep dds for differential expression
dds <- DESeq(dds)

#' ## Calculate differential expression
#' This experiment contains six injection conditions (PBS, STM, SE, PBS+PMN, STM+PMN, SE+PMN) at 8h post injection. To start, we'll calculate the differential expression of each injection over PBS at the matching time point (e.g. STM 8h over PBS 8h, ST 2h over PBS 2h). We'll save each of these as a .csv file in the `diff_expression` folder in the results directory.

# Set up the samples for diff expression
unique(colData(dds)$code_name)

# Select all but PBS
samps <- unique(colData(dds)$code_name)[2:6] %>% 
  as.character()


#' We'll set up a loop to calculate differential expression and iterate over all the samples. We'll make sure each is matched to the proper PBS control with a logical test at the start of the loop that checks the time point of each sample.

#+ loop, cache=T
# Set up loop to calculate differential expression for all samples over PBS
for(i in seq(samps)){

    # Sample for i in loop
  sample <- samps[i]
  
  
  # Calculate differential expression - sample over PBS
  res <- results(dds, 
                 contrast = c('code_name', sample, 'PBS'))
  
  
  # Add annotation - symbol and entrezID
  res$symbol <- rownames(res)
  # Add column for gene Entrez ID
  res$entrez <- mapIds(org.Hs.eg.db,
                       keys = rownames(res),
                       column = 'ENTREZID',
                       keytype = 'SYMBOL',
                       multiVals = 'first')
  # Add column for gene name
  res$name <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = 'GENENAME',
                     keytype = 'SYMBOL',
                     multiVals = 'first')
  # Make results dataframe
  
  
  res.df <- as.data.frame(res) %>% 
    arrange(padj)
  
  
  # Save output file
  file.name <- paste0(here('results/DESeq2/diff_expression//'),
                      sample,'_over_PBS_diffexpress.csv')
  write.csv(res.df, file = file.name, row.names = F)
  
}


#' ## Differential expression over PBS+PMNs control

#' Want to compare the SE+PMNs and STM+PMNs to the PBS+PMNs control and save those files in a separate folder.


# Samples to compare to PBS+PMNs
samps <- samps[-c(1:2)]
samps

# Create directory for diff expression files over PBS+PMNs
dir.create(here('results/DESeq2/diff_expression_pmn'))

#+ loop1, cache=T
# Set up loop to calculate differential expression for all samples over PBS
for(i in seq(samps)){
  
  # Sample for i in loop
  sample <- samps[i]
  sample1 <- gsub(pattern = '\\+PMNs', replacement = '', sample)
  
  
  # Calculate differential expression - sample over PBS
  res <- results(dds, 
                 contrast = c('code_name', sample, sample1))
  
  
  # Add annotation - symbol and entrezID
  res$symbol <- rownames(res)
  # Add column for gene Entrez ID
  res$entrez <- mapIds(org.Hs.eg.db,
                       keys = rownames(res),
                       column = 'ENTREZID',
                       keytype = 'SYMBOL',
                       multiVals = 'first')
  # Add column for gene name
  res$name <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = 'GENENAME',
                     keytype = 'SYMBOL',
                     multiVals = 'first')
  # Make results dataframe
  
  
  res.df <- as.data.frame(res) %>% 
    arrange(padj)
  
  
  # Save output file
  file.name <- paste0(here('results/DESeq2/diff_expression_pmn//'),
                      sample,'_over_',sample1,'_diffexpress.csv')
  write.csv(res.df, file = file.name, row.names = F)
}



#' ## Ending notes
#' The differential expression outputs are now in the results directory. These will be used to make MA and volcano plots.




#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/03_differential_expression.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)