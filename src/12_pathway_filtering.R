#'---
#' title: "RNA seq workflow part 12 - Filtering Reactome pathways"
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



# Find pathway similarities and differences between samples
library(here)
library(dplyr)
library(ggplot2)
library(magrittr)


# Function to make heatmap:
# 1. Filter reactome enrichment from 3 samples by pathways significant in at least 1 sample.
# 2. Make cluster heatmap of NES values from filtered pathways list

pathway_filter <- function(d1, d2, 
                            label1, label2, p = 1){
  
  # Get list of pathways - Must be significant in at least one sample
  paths <- NULL
  
  for(i in list(d1, d2)){
    a <- filter(i, pvalue <= p) %>% 
      select(Description)
    
    print(paste(nrow(a), 'significant pathways'))
    
    paths <- rbind(paths, a) %>% 
      unique()
  }
  
  print(paste(nrow(paths), 'pathways significant in at least one sample with p value cutoff at', p))
  
  
  # Get NES for each pathway from each sample
  paths1 <- paths
  
  for(i in list(d1, d2)){
    a <- filter(i, Description %in% paths$Description) %>% 
      select(c(Description, NES))
    
    paths1 <- left_join(paths1, a, by = 'Description')
    
  }
  colnames(paths1)[2:3] <- c(label1,label2)
  
  return(paths1)
}
  

# Function for making heatmaps  
path_heatmap <- function(paths1, rows = T, cutrows = 1){
  x <- as.matrix(paths1[2:3])
rownames(x) <- paths1$Description

#install.packages("matrixStats")
library(matrixStats)

x <- x[order(-rowVars(x)),]


library(pheatmap)
# Define top variance genes, extract from full data
topVarGenes <- order(rowVars(x), decreasing = TRUE)
mat <- x[topVarGenes, ]


# convert to fold over mean of all samples
#mat <- mat - rowMeans(mat)  

pheatmap(mat,
         show_rownames = rows,
         breaks = seq(-2.6, 2.7, by = .055),
         cluster_cols = F,
         cutree_rows = cutrows,
         border_color = NA,
         fontsize = 10,
         cellwidth = 25,
         cellheight = 10,
         treeheight_row = 20)

}

# function for saving heatmaps as PNG files
dir.create(here('img/reactome_heatmap_subsets/'))

save_heatmap <- function(plot1, label, subsetid, width = 12, height = 4){
png(filename = here(paste0('img/reactome_heatmap_subsets/',
                           label,
                           '_gsea_diff_heatmap_subset_',
                           subsetid,'.png')),
    width = width, height = height, units = 'in', res = 300)
print(plot1)
dev.off()
}

  
  


# Load samples
gsea.dir <- here('results/DESeq2/GSEA/reactome/')

# SE samples
se_pbs <- read.csv(file.path(gsea.dir, 'GSEA_reactome_SE_over_PBS.csv'), stringsAsFactors = F)
se_pmn <- read.csv(file.path(gsea.dir, 'GSEA_reactome_SE+PMNs_over_PBS.csv'))


# STM samples
stm_pbs <- read.csv(file.path(gsea.dir, 'GSEA_reactome_STM_over_PBS.csv'), stringsAsFactors = F)
stm_pmn <- read.csv(file.path(gsea.dir, 'GSEA_reactome_STM+PMNs_over_PBS.csv'))




# Look at specific groups from the pathway heatmaps
# SE/PBS and SE+PMN/PBS 
se1 <- pathway_filter(se_pbs, se_pmn, p = 0.05, 
                        label1 = 'SE', label2 = 'SE_PMNs')
head(se1)

# pathways up in SE but down in SE+PMNs
# Remove redundant pathways (determined manually)
se_paths <- filter(se1, SE > 0 & SE_PMNs < 0) %>% 
  filter(!Description %in% c('Influenza Life Cycle',
                             'Influenza Infection',
                             'Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.',
                             'Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal'))
p1 <- path_heatmap(se_paths)
save_heatmap(p1, 'SE', 'pmn_down')


# pathways down in STM and SE but up in ST
se_paths1 <- filter(se1, SE < 0 & SE_PMNs > 0)
p2 <- path_heatmap(se_paths1)
save_heatmap(p2,'SE','pmn_up')

# Pathways up in all three
se_up <- filter(se1, SE > 1.5 & SE_PMNs > 1.5) #%>%
p3 <- path_heatmap(se_up)
save_heatmap(p3,'SE','all_up', height = 5)



# STM heatmaps
# STM/PBS and STM+PMN/PBS 
stm1 <- pathway_filter(stm_pbs, stm_pmn, p = 0.05, 
                      label1 = 'STM', label2 = 'STM_PMNs')
head(stm1)


# Pathways up in STM and down in STM+PMNs
stm_paths <- filter(stm1, STM > 0 & STM_PMNs < 0) #%>% 
p <- path_heatmap(stm_paths)
save_heatmap(p,'STM','pmn_down')

# pathways down in STM and up in STM+PMNs
stm_paths1 <- filter(stm1, STM < 0 & STM_PMNs > 0)
p <- path_heatmap(stm_paths1)
save_heatmap(p,'STM','pmn_up')

# Pathways up in both
stm_up <- filter(stm1, STM > 1.6 & STM_PMNs > 1.6) #%>%
p <- path_heatmap(stm_up)
save_heatmap(p,'STM','all_up',height = 5)



