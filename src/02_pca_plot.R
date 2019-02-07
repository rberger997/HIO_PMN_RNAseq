#'---
#' title: "RNA seq workflow part 2 - PCA plot"
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
#' This script is part of a workflow for RNA seq processing and analysis (modified from [Bioconductor](https://www.bioconductor.org/help/workflows/rnaseqGene/)). The purpose is to take in a DESeq2 object (dds) and return a PCA plot of all the samples.

#------------------------------------------------------------
#' # Begin script

#' ## Libraries and directories
# #+ load_pkgs, error=T
library(DESeq2)
library(ggplot2)
library(here)
library(knitr)
library(magrittr)
library(rmarkdown)
source(here("src/ggplot2-themes.R"))


#+ set_cache_dir, include=F
# Set a directory for all cache files to go
dir.create(here('results/DESeq2/src_html_output/knitr_cache'))
opts_chunk$set(cache.path = here('results/DESeq2/src_html_output/knitr_cache/'))


#' The output files will go to the `DESeq2` folder in the results directory.

results.dir <- here('results/DESeq2/')

# Create new directories
dir.create(here('img/'))


#' ## Load the data
#' We previously saved the `dds` object in the results folder. Let's load it now:

dds <- readRDS(here('results/DESeq2/dds_all.rds'))


#' ## Transform data using rlog
#' Before running PCA, we'll do an rlog transformation of the data. This is a variance stabilizing transformation that "transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size" ([source](https://rdrr.io/bioc/DESeq2/man/rlog.html)). 

#+ rlog, cache=T
rld <- rlog(dds, blind = FALSE)

# Compare before and after transformation
head(assay(dds)[,1:4], 3) ; head(assay(rld)[,1:4], 3)

# Save the rlog transformed data
saveRDS(rld, file = file.path(results.dir, 'rld_all.rds'))

# Load the rlog transformed data
rld <- readRDS(here('results/DESeq2/rld_all.rds'))

#' ## PCA plot
#' Now we can run PCA using the `plotPCA` function from DESeq2 and use those values to make a custom plot with `ggplot`. First we'll set the order of injection to how we want it to show up in the plot.

colData(rld)
# Set order of injection
unique(colData(rld)$code_name)

colData(rld)$Inject <- colData(rld)$code_name %>% 
  gsub(pattern = '\\+PMNs', replacement = '', x = .)
colData(rld)$Inject

  # factor(., levels(.)[c(1,6,2,5,3,4)])


#+ pcaplot
# PCA plot (use just only for PC variance estimates)
pca <- plotPCA(rld, intgroup = c('code_name','Inject','PMN'))

# Get PCA data
pca.df <- plotPCA(rld, intgroup = c('code_name', 'Inject', 'PMN'), returnData = TRUE)


# Make plot
pca_plot <- function(input){
  ggplot(data = input, aes(x = PC1, y = PC2))+
    geom_hline(yintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    geom_point(shape = 21, stroke = 1.5, 
               aes(fill = as.factor(Inject),
                   color = as.factor(PMN)), 
               size = 6) +
    theme1 + 
    #scale_y_continuous(limits = c(-25,25), breaks = seq(-25, 25, by = 10)) +
    scale_fill_brewer(palette = "Set1", name = 'Injection') +
    #scale_color_brewer(palette = "Set2", name = 'Time', direction = 1) +    
    scale_color_manual(values=c("black", 'gray'), name = 'PMNs')+
    theme(legend.position = "right") +
    
    coord_fixed(ratio = 1) +
    xlab(pca$labels$x) + #pull variance estimates from al. plotPCA call
    ylab(pca$labels$y) +
    # Move y axis
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 1, 
                                                      b = 0, l = 0))) +
    # Move x axis
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, 
                                                      b = 0, l = 0)))+
    # Shrink axis labels down
    theme(plot.caption = element_text(vjust = 1), 
          axis.title = element_text(size = 24), 
          axis.text.x = element_text(size = 18), 
          axis.text.y = element_text(size = 18), 
          plot.title = element_text(size = 30))
}

p <- pca_plot(pca.df)
p
#' ## Observations
#' 

#' ## Save plot
#' We'll save the PCA plot as a .png file in the `img` directory.

# save png of plot
#+ save, eval=F
png(filename = here("/img/pca.png"),
    width = 12, height = 5, res = 300, units = 'in')
print(p)
dev.off()



#' # PCA plot without PMN controls
#' 



#' We previously saved the `dds` object in the results folder. Let's load it now:
dds1 <- readRDS(here('results/DESeq2/dds_no_controls.rds'))

#' ## Transform data using rlog
#+ rlog, cache=T
rld1 <- rlog(dds1, blind = FALSE)

# Save the rlog transformed data
saveRDS(rld1, file = file.path(results.dir, 'rld_no_controls.rds'))


# Load the rlog transformed data
rld1 <- readRDS(here('results/DESeq2/rld_no_controls.rds'))


#' ## PCA plot

colData(rld1)$Inject <- colData(rld1)$code_name %>% 
  gsub(pattern = '\\+PMNs', replacement = '', x = .)


# PCA plot (use just only for PC variance estimates)
pca <- plotPCA(rld1, intgroup = c('code_name','Inject','PMN'))

# Get PCA data
pca.df <- plotPCA(rld1, intgroup = c('code_name', 'Inject', 'PMN'), returnData = TRUE)

plot_nocontrols <- pca_plot(pca.df)
plot_nocontrols




# save png of plot
#+ eval=F
png(filename = here("/img/no_controls_pca.png"),
    width = 8, height = 6, res = 300, units = 'in')
print(plot_nocontrols)
dev.off()


# Save ggplot objects
dir.create(here('img/ggplot_objects/'))
saveRDS(plot_nocontrols, file = here('img/ggplot_objects/gg_nocontrols_pcaplot.rds'))


#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/02_pca_plot.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)