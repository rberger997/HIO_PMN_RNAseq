#'---
#' title: "RNA seq workflow part 4 - Volcano plot"
#' subtitle: "HIOs PMNs RNAseq"
#' author: "Ryan Berger"
#' date: "2019-01-24"
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
#' This script is part of a workflow for RNA seq processing and analysis (modified from [Bioconductor](https://www.bioconductor.org/help/workflows/rnaseqGene/)). The purpose is to take in differential expression .csv files and create volcano plots showing global gene changes. For this analysis, we'll make a single figure with the volcano plots from all the differential expression files over PBS control.

#------------------------------------------------------------
#' # Begin script

#' ## Libraries and directories
#+ load_pkgs, error=T
library(dplyr)
library(ggplot2)
library(here)
library(knitr)
library(magrittr)
library(rmarkdown)
library(stringi)


# Create directory to save ggplot objects
dir.create(here('img/ggplot_objects'), recursive = T)


#+ set_cache_dir, include=F
# Set a directory for all cache files to go
opts_chunk$set(cache.path = here('results/DESeq2/src_html_output/knitr_cache/'))


#' ## Load the data
#' Start by making a single volcano plot using `ggplot`. We'll use the STM over PBS file to make this.

# Load data 
stm <- read.csv(here('results/DESeq2/diff_expression/STM_over_PBS_diffexpress.csv'), header = T)

head(stm)


#' To make things easier for coloring and faceting the plot, let's add some categories to the data. We're going to color the points based on:
#' 
#' * Non significant changes (padj > 0.05)
#' * Significant increases (padj < 0.05, log2FoldChange > 1)
#' * Significant decreases (padj < 0.05, log2FoldChange < -1)
#' 
#' We're also going to facet the plots by sample and time point. We'll make a function to define categories for each of these in our data.

# Function to prep data
prep_df <- function(df, label, pmn){
  # Set label for conditions
  df$label <- label
  df$pmn <- pmn
  
  # Set up labels for colors  
  df$colors <- 'Non significant'
  df[which(df$log2FoldChange < 0 & df$padj < 0.05), 'colors'] <- 'Decreasing'
  df[which(df$log2FoldChange > 0 & df$padj < 0.05), 'colors'] <- 'Increasing'  
  df$colors <- factor(df$colors, levels = c('Non significant', 'Decreasing', 'Increasing'))
  
  return(df)
}

stm <- prep_df(stm, label = 'STM/PBS', pmn = 'No_PMNs')

#' ## Build volcano plot


ggplot(data = stm, aes(x=log2FoldChange, y=-log10(padj), 
                        color=colors, label=symbol))+
  geom_point()+
  geom_vline(xintercept = 0, linetype = 'dotted')+
  geom_hline(yintercept = 0, linetype = 'dotted')+
  scale_color_manual(name = '',
                     values = c('Non significant' = 'gray',
                                'Decreasing' = 'blue',
                                'Increasing' = 'red'))+
  theme_bw()+ 
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), axis.text = element_text(size = 10))


#' ## Build multi-panel volcano plot
#' This individual plot will serve as the template for our multipanel facet plot with all the samples and time points. To make it, we'll load all the files and combine them into a single dataframe that we'll facet into rows and columns in a plot.

# Locate all the files
files <- list.files(here('results/DESeq2/diff_expression/'))
files

# Create an empty object to put the files into
data1 <- NULL

# Iterate over all files and combine into a single dataframe
for(i in seq(files)){
  
  file <- files[i]
  # Extract sample label and time point
  label <- paste0(gsub('_.*$',"",file),'/PBS')
  pmn <- grep(pattern = 'PMN', x = file)
  pmn <- ifelse(length(pmn) == 1, '+PMNs', '-PMNs')
  
  # Load file
  temp <- read.csv(here(paste0('results/DESeq2/diff_expression/',file)), header = T)
  
  
  temp <- prep_df(temp, label = label, pmn = pmn)
  
  data1 <- rbind(data1, temp)
  print(paste(i, 'of', length(files), 'done:', files[i]))
}



# Set order of samples for plot
data1$label <- factor(data1$label, 
                      levels = c('SE/PBS','STM/PBS','PBS+PMNs/PBS',
                                 'SE+PMNs/PBS','STM+PMNs/PBS'))



#+ fig, fig.height = 7, fig.width = 10, fig.align = 'center'
# Make plot - facet time into rows, samples into columns

p <- ggplot(data = data1, aes(x=log2FoldChange, y=-log10(padj), 
                              color=colors, label=symbol))+
  geom_point(size=.85)+
  geom_vline(xintercept = 0, linetype = 'dotted')+
  geom_hline(yintercept = 0, linetype = 'dotted')+
  scale_color_manual(name = '',
                     values = c('Non significant' = 'gray',
                                'Decreasing' = 'blue',
                                'Increasing' = 'red'),
                     labels = c('Non significant',
                                'Significant decrease',
                                'Significant increase'))+
  theme_bw()+ 
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14, face = 'bold'),
        strip.text.y = element_text(size = 14, face = 'bold'))+
  facet_grid(cols = vars(label))+
  xlim(c(-12,12))

p


#' ## Thoughts about results


#+ new
#' ## Filter out lncRNA


# Identify long non coding RNAs
lncrna <- dplyr::filter(data1, grepl("^RP([0-9]).*",symbol) & is.na(entrez) == T)
head(lncrna)

# Select all rows in data1 that are not in lncrna
data2 <- dplyr::anti_join(data1, lncrna)


#+ fig2, fig.height = 7, fig.width = 10, fig.align = 'center'
# Make plot - facet time into rows, samples into columns

volcano_plot <- function(input){
  ggplot(data = input, aes(x=log2FoldChange, y=-log10(padj), 
                           color=colors, label=symbol))+
    geom_point(size=.85)+
    geom_vline(xintercept = 0, linetype = 'dotted')+
    geom_hline(yintercept = 0, linetype = 'dotted')+
    scale_color_manual(name = '',
                       values = c('Non significant' = 'gray',
                                  'Decreasing' = 'blue',
                                  'Increasing' = 'red'),
                       labels = c('Non significant',
                                  'Significant decrease',
                                  'Significant increase'))+
    theme_bw()+ 
    theme(panel.grid.major = element_line(linetype = "blank"), 
          panel.grid.minor = element_line(linetype = "blank"), 
          axis.title = element_text(size = 14, face = "bold"), 
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          strip.text.x = element_text(size = 14, face = 'bold'),
          strip.text.y = element_text(size = 14, face = 'bold'))+
    facet_grid(cols = vars(label))+
    xlim(-12,12)
}

p <- volcano_plot(data2)
p


#' ## Save png of plot
#+ save, eval=F
png(filename = here("/img/volcano_facet.png"),
    width = 900, height = 400)
print(p)
dev.off()

# Save data
write.csv(data2, file = here('results/DESeq2/volcano_data.csv'))



#' # Split samples
#' Make a volcano plot without the PBS+PMNs sample in it.

#+
d1 <- dplyr::filter(data2, label != 'PBS+PMNs/PBS') 
d1$pmn <- factor(d1$pmn, levels = c('-PMNs', '+PMNs'))

d1$label <- gsub('\\+PMNs', '',d1$label)
unique(d1$label)



volcano_plot <- function(input, lim, plim){
  
  # Get labels for number of significant genes
  ann_decr <- input %>% 
    filter(colors == 'Decreasing') %>% 
    group_by(label, pmn) %>% 
    dplyr::count(label, colors = colors) %>% 
    mutate(padj = plim,
           log2FoldChange = -lim+1)
  
  ann_incr <- input %>% 
    filter(colors == 'Increasing') %>% 
    group_by(label, pmn) %>% 
    dplyr::count(label, colors = colors) %>% 
    mutate(padj = plim,
           log2FoldChange = lim-1)
  
  
  ggplot(data = input, aes(x=log2FoldChange, y=-log10(padj), 
                           color=colors, label=symbol))+
    geom_point(size=.85)+
    geom_vline(xintercept = 0, linetype = 'dotted')+
    geom_hline(yintercept = 0, linetype = 'dotted')+
    scale_color_manual(name = '',
                       values = c('Non significant' = 'gray',
                                  'Decreasing' = 'blue',
                                  'Increasing' = 'red'),
                       labels = c('Non significant',
                                  'Decreasing (Adj.p < 0.05)',
                                  'Increasing (Adj.p < 0.05)'))+
    theme_bw()+ 
    theme(panel.grid.major = element_line(linetype = "blank"), 
          panel.grid.minor = element_line(linetype = "blank"), 
          axis.title = element_text(size = 14, face = "bold"), 
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          strip.text.x = element_text(size = 14, face = 'bold'),
          strip.text.y = element_text(size = 14, face = 'bold'))+
    facet_grid(cols = vars(label), rows = vars(pmn))+
  geom_text(data = ann_decr,label = ann_decr$n,
            color = 'blue', fontface = 'bold', size = 5)+
    geom_text(data = ann_incr,label = ann_incr$n,
              color = 'red', fontface = 'bold', size = 5)+
    xlim(-lim,lim)+
    ylim(0, -log10(plim))
}

volcano_plot(d1, lim = 8, plim = 1e-50)


#' ## Save as png
#+ eval=F
png(filename = here("/img/volcano_facet1.png"),
    width = 8, height = 5, units = 'in', res = 300)
volcano_plot(d1, lim = 8, plim = 1e-80)
dev.off()





#-----------------------------------------------------------------

#' # Alternative plot 
#' Make a volcano plot using David's script. This version uses geom_point and jitter to visualize gene expression changes.


library(RColorBrewer)
source(here('src/ggplot2-themes.R'))


# Load data
data2 <- read.csv(here('results/DESeq2/volcano_data.csv'))



# make a new status column that will indicate statistical significance
data2$status <- 'a'
data2[which(data2$log2FoldChange > 1 & data2$padj < 0.05), 'status'] <- 'b'
data2[which(data2$log2FoldChange < -1 & data2$padj < 0.05), 'status'] <- 'c'

# set order so that blue and red are plotted on top of grey
data2 <- data2[order(data2$status),]


# Only select Serovars samples
#data2 <- dplyr::filter(data2, label %in% c('SE','STM','SE+PMNs','STM+PMNs'))



# Build plot
ggplot(data = data2,
       aes(x = log2FoldChange, y = label)) +
  geom_point(position = position_jitter(h = 0.4),
             aes(fill = status, color = status),
             shape = 21, size = 0.5)+
  scale_fill_manual(values = c("grey70", color.set[1], color.set[2])) +
  scale_color_manual(values = c("grey70", color.set[1], color.set[2])) +
  xlim(c(-6, 6)) +
  ylab("") +
  xlab(expression(paste("log"[2],"FC (HIO + bacteria / HIO + PBS)"))) +
  theme(axis.text = element_text(size = 18),
        axis.text.y = element_text(size = 24,
                                   face = 'italic'),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey30",
                                    fill = NA),
        axis.title = element_text(size = 24),
        strip.text = element_text(size = 24),
        strip.text.y = element_text(size = 24,
                                    face = "italic",
                                    angle = 0))




#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/04_volcano_plot.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)


# Copy to dropbox
# source(here('src/Hs_align_src/XX_copy_to_dropbox.R'))
