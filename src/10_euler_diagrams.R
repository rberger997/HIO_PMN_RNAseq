#'---
#' title: "RNA seq workflow part 10 - Euler diagrams"
#' subtitle: "HIOs PMN RNAseq"
#' author: "Ryan Berger"
#' date: "2019-02-07"
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


# Making Euler diagrams with the eulerr package

# Should make area proportional venn diagrams from lists or a matrix


#' ## Libraries and directories
library(dplyr)
library(eulerr)
library(here)
library(magrittr)
library(tidyr)


# Load data from volcano plots
data2 <- read.csv(here('results/DESeq2/volcano_data1.csv'), 
                  stringsAsFactors = F)
head(data2)
unique(data2$label)


# Function to make venn diagrams
make_venn <- function(df,x1,x2,x3, padj = 1, l2fc = 0, show_labels = T){

  # Set significance filter
  df$signif <- ifelse(df$padj < padj & abs(df$log2FoldChange) > l2fc, 'Significant', 'Non significant')
  
  # want to compare all significant genes
incr <- filter(df, signif == 'Significant') %>% 
  select(c(symbol,label)) %>% 
  arrange(label)


# Split into three groups by sample
a <- filter(incr, label == x1) %>% 
  .$symbol
b <- filter(incr, label == x2) %>% 
  .$symbol
c <- filter(incr, label == x3) %>% 
  .$symbol
  
# Calculate number of shared genes in each group
ABC <- length(c[c%in%a[a%in%b]])
AB <- length(a[a%in%b]) - ABC
AC <- length(a[a%in%c]) - ABC
A <- length(a) - AB -AC - ABC
BC <- length(b[b%in%c]) - ABC
BA <- length(b[b%in%a]) - ABC
B <- length(b)-BA-BC-ABC
C <- length(c)-AC-BC-ABC


# Formatting options
eulerr_options(pointsize = 14)
options(digits = 4)


# Input in the form of a named numeric vector
fit1 <- euler(c("A" = A, "B" = B, "C" = C,
                "A&B" = AB, "A&C" = AC, "B&C" = BC,
                "A&B&C" = ABC))

diagram <- plot(fit1, 
             quantities = T,
             fill = c("lightblue", "lightcoral", "lemonchiffon"),
             lty = 1,
             labels = if(show_labels == F){
               F
             }else{
               c(x1,x2,x3)})
return(diagram)
}

head(data2)
# STM mutants venn diagrams 
# (leave labels off for better formatting in powerpoint)
m2 <- make_venn(data2, 'STM+PMNs/STM','PBS+PMNs/PBS','SE+PMNs/SE', 
                padj = 0.05, l2fc = 1, show_labels = T)
m2



# Save venn diagram objects
# saveRDS(m2, file = here(paste0('img/ggplot_objects/gg_mut_venn_2h.rds')))



# Get list of genes shared by different subsets
a <- filter(data2, label == 'PBS+PMNs/PBS') %>% 
  filter(abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  select(symbol, log2FoldChange, padj) %>% 
  set_colnames(c('symbol', 'log2_PBS','padj_PBS'))

b <- filter(data2, label == 'STM+PMNs/STM') %>% 
  filter(abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  select(symbol, log2FoldChange, padj) %>% 
  set_colnames(c('symbol', 'log2_STM','padj_STM'))

c <- filter(data2, label == 'SE+PMNs/SE') %>% 
  filter(abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  select(symbol, log2FoldChange, padj) %>% 
  set_colnames(c('symbol', 'log2_SE','padj_SE'))

shared <- inner_join(a, b, by = 'symbol') %>% 
  inner_join(c, by = 'symbol') %>% 
  select(symbol, log2_PBS, log2_STM, log2_SE) %>% 
  mutate(avg_log2 = (log2_PBS + log2_STM + log2_SE)/3) %>% 
  arrange(-avg_log2)

# Get gene names to add to dataframe
gene_names <- select(data2, symbol, name)
gene_names <- unique(gene_names)


shared <- left_join(shared, gene_names, by = 'symbol')
