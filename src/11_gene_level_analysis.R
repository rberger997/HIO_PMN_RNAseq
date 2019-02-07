#'---
#' title: "RNA seq workflow - Gene level analysis"
#' subtitle: "HIOs PMN RNAseq"
#' author: "Ryan Berger"
#' date: "2019-02-5"
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


#+ script, include=F

# Make volcano plots of samples over STM for gene comparisons
# Use functions from volcano plot script

library(dplyr)
library(ggplot2)
library(here)
library(knitr)
library(magrittr)
library(rmarkdown)
library(stringi)

# Function to prep data
prep_df <- function(df, label){
  # Set label for conditions
  df$label <- label
  
  # Set up labels for colors  
  df$colors <- 'Non significant'
  df[which(df$log2FoldChange < -1 & df$padj < 0.05), 'colors'] <- 'Decreasing'
  df[which(df$log2FoldChange > 1 & df$padj < 0.05), 'colors'] <- 'Increasing'  
  df$colors <- factor(df$colors, levels = c('Non significant', 'Increasing', 'Decreasing'))
  
  return(df)
}


# Locate all the files
files <- list.files(here('results/DESeq2/diff_expression_pmn/'))
files

# Create an empty object to put the files into
data1 <- NULL

# Iterate over all files and combine into a single dataframe
for(i in seq(files)){
  file <- files[i]
  
  # Extract sample label and time point
  label <- gsub('_.*$',"",file)
  label1 <- gsub('\\+PMNs','',label)
  label <- paste0(label, '/', label1)
  
  
  # Load file
  temp <- read.csv(here(paste0('results/DESeq2/diff_expression_pmn/',file)), header = T)
  
  temp <- prep_df(temp, label = label)
  
  data1 <- rbind(data1, temp)
  print(paste(i, 'of', length(files), 'done:', files[i]))
}


# Set order of samples for plot
data1$label <- factor(data1$label, levels = c('PBS+PMNs/PBS','SE+PMNs/SE','STM+PMNs/STM'))



# Identify long non coding RNAs
lncrna <- dplyr::filter(data1, grepl("^RP([0-9]).*",symbol) & is.na(entrez) == T)

# Select all rows in data1 that are not in lncrna
data2 <- dplyr::anti_join(data1, lncrna)



# Make plot - facet time into rows, samples into columns

volcano_plot <- function(input, lim=5, plim = 1e-16){
  
  # Get labels for number of significant genes
  ann_decr <- input %>% 
    filter(colors == 'Decreasing') %>% 
    group_by(label) %>% 
    dplyr::count(label, colors = colors) %>% 
    mutate(padj = plim,
           log2FoldChange = -lim+.4)
  
  ann_incr <- input %>% 
    filter(colors == 'Increasing') %>% 
    group_by(label) %>% 
    dplyr::count(label, colors = colors) %>% 
    mutate(padj = plim,
           log2FoldChange = lim-.4)

  
  # Volcano plot
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
                                  'Significant increase',
                                  'Significant decrease'))+
    theme_bw()+ 
    theme(panel.grid.major = element_line(linetype = "blank"), 
          panel.grid.minor = element_line(linetype = "blank"), 
          axis.title = element_text(size = 14, face = "bold"), 
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          strip.text.x = element_text(size = 14, face = 'bold'),
          strip.text.y = element_text(size = 14, face = 'bold'))+
    #facet_grid(rows = vars(time), cols = vars(label))+
    geom_text(data = ann_decr,label = ann_decr$n,
              color = 'blue', fontface = 'bold', size = 5)+
    geom_text(data = ann_incr,label = ann_incr$n,
              color = 'red', fontface = 'bold', size = 5)+
    xlim(-lim,lim)+
    ylim(0, -log10(plim))+
    ggtitle(paste0(unique(input$label), ' 8h pi'))+
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
  
  
}

volcano_plot(data1)

#+ fig2, fig.height = 7, fig.width = 10, fig.align = 'center'

# Save volcano plot function
save_volcano <- function(gg, sample, sample1){
  # Create folder for volcano plots in img directory
  dir.create(here('img/volcano_plots'))
  
  # Save plot
  png(filename = here(paste0('img/volcano_plots/',sample,'_over_',sample1, '_volcano_plot.png')), width = 6, height = 4, units = 'in', res = 300)
  print(gg)
  dev.off()
  
  # Save ggplot object
  saveRDS(gg, 
          file = here(paste0('img/ggplot_objects/gg_',sample,'_over_',sample1,'_volcano.rds')))
}


samps <- as.character(unique(data1$label))
# Iterate through all plots
for(i in samps){

  # Set p limit for volcano plots
    #plim <- ifelse(i %in% c('SPI1', 'SPI2'), 1e-8, 1e-15)

    # Make volcano plot
    vplot <- dplyr::filter(data2, label == i) %>% 
      arrange(colors) %>% 
      volcano_plot(., lim = 5, plim = 1e-50)
  
    
    # Save volcano plot
    # Get labels for save
    sample <- gsub('\\/.*', '', i)
    sample1 <- gsub('\\+PMNs','',i) %>% 
      gsub('\\/.*', '', .)
  
    save_volcano(vplot, sample = sample, sample1 = sample1)
    
    # Update progress
    print(paste(i, 'done'))
  
}



#' # Combine plots
combine_plot <- arrange(data2, colors) %>% 
  volcano_plot(., lim = 5, plim = 1e-50)+
  facet_grid(cols = vars(label))+
  ggtitle('')
combine_plot

# Save plot
png(filename = here('img/volcano_plots/sample+PMNs_over_sample_combined_volcano.png'), width = 9, height = 4, units = 'in', res = 300)
combine_plot
dev.off()



# Barplots of top gene changes

# set up 
mut_data1 <- filter(data2, 
                      colors != 'Non significant')




# Function for making barplots of significant gene changes
gene_barplot <- function(df, samp, t, n, remove = NA, lim){
  # Change label for samples in dataframe
  df$label <- gsub('\\/.*', '', df$label)
  
  # Get label for sample
  samp1 <- gsub('\\+PMNs', '', samp)
  
  # Set colors for barplot
  labs <- c('red', 'blue')
  names(labs) <- c(paste('Increased vs.', samp1), paste('Decreased vs.', samp1))
  
  
  # Filter top genes
  df %>% 
    filter(label == samp & !symbol %in% remove) %>% 
    arrange(-abs(log2FoldChange)) %>% 
    select(symbol, log2FoldChange, padj, name, label) %>% 
    mutate(change = as.factor(ifelse(log2FoldChange < 0, 
                                     paste('Decreased vs.',samp1), 
                                     paste('Increased vs.',samp1)))) %>% 
    head(., n) %>% 
    
    # Make barplot
    ggplot(., 
           aes(x = reorder(symbol, log2FoldChange), 
               y = log2FoldChange, fill = change))+
    geom_bar(stat = 'identity', position = 'dodge')+
    coord_flip()+
    theme(panel.grid.minor = element_line(linetype = "blank"))+
    theme_bw()+
    scale_fill_manual(values = labs, name = 'Trend')+
    labs(x = '',
         title = paste0(samp, '/',samp1, ' 8pi'),
         subtitle = paste('Top', n, 'genes by fold change (padj < 0.05)'))+
    theme(plot.title = element_text(size=12, hjust=0.5, 
                                    face="bold", vjust=-2))+
    theme(plot.subtitle = element_text(size=8, hjust=0.5, vjust = -1))+
    geom_hline(yintercept = 0)+
    ylim(-lim,lim)
  
}


# Function to save barplots (png and RDS)
save_barplot <- function(p, w, h, file_label){
  # Create directory
  dir.create(here('img/gene_barplots/'), recursive = T)
 
   # Save image
  png(filename = here(paste0('img/gene_barplots/', file_label, '_gene_barplot.png')), width = w, height = h, units = 'in', res = 300)
  print(p)
  dev.off()
  
  # Save plot as ggplot object using saveRDS
  saveRDS(p, 
          file = here(paste0('img/ggplot_objects/gg_',file_label,'_gene_barplot.rds')))
  
}


# Set top n number of genes to have in barplot
n <- 15

# PBS barplot
gene_barplot(mut_data1, samp = 'PBS+PMNs', n = n, lim = 25) %>% 
  save_barplot(., w = 5, h = 4, file_label = 'PBS+PMN_over_PBS')

gene_barplot(mut_data1, samp = 'SE+PMNs', n = n, lim = 8) %>% 
  save_barplot(., w = 5, h = 4, file_label = 'SE+PMN_over_SE')

gene_barplot(mut_data1, samp = 'STM+PMNs', n = n, lim = 8) %>% 
  save_barplot(., w = 5, h = 4, file_label = 'STM+PMN_over_STM')






#+ render, include=F
# Render source file to html 
# dir.create(here('results/DESeq2_human/src_html_output'))

# render.dir <- here('results/DESeq2_human/src_html_output/')

# render(here('src/Hs_align_src/XX_volcano_plots_stm.R'), output_dir = render.dir, intermediates_dir = render.dir, clean = TRUE)

