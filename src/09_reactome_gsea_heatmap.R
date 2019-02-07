# Find pathway similarities and differences between samples

library(here)
library(dplyr)
library(ggplot2)
library(magrittr)


# Function to make heatmap:
# 1. Filter reactome enrichment from 2 samples by pathways significant in at least 1 sample.
# 2. Make cluster heatmap of NES values from filtered pathways list

pathway_heatmap <- function(d1, d2, 
                            label1, label2, 
                            n, rows = T, p = 1, cutrows = 1){
  
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

x <- as.matrix(paths1[2:3])
rownames(x) <- paths1$Description

#install.packages("matrixStats")
library(matrixStats)

x <- x[order(-rowVars(x)),]


library(pheatmap)
# Define top variance genes, extract from full data
topVarGenes <- head(order(rowVars(x), decreasing = TRUE), n) 
mat <- x[topVarGenes, ]


# convert to fold over mean of all samples
#mat <- mat - rowMeans(mat)  

pheatmap(mat,
         show_rownames = rows,
         breaks = seq(-2.6, 2.7, by = .055),
         cluster_cols = F,
         cutree_rows = cutrows)
         #fontsize = 10,
         #cellwidth = 20,
         #cellheight = 10

}



# Load samples
gsea.dir <- here('results/DESeq2/GSEA/reactome/')


# Samples without PMNs
stm <- read.csv(file.path(gsea.dir, 'GSEA_reactome_STM_over_PBS.csv'), stringsAsFactors = F)
se <- read.csv(file.path(gsea.dir, 'GSEA_reactome_SE_over_PBS.csv'), stringsAsFactors = F)


# Samples with PMNs
pbs_pmn <- read.csv(file.path(gsea.dir, 'GSEA_reactome_PBS+PMNs_over_PBS.csv'), stringsAsFactors = F)
stm_pmn <- read.csv(file.path(gsea.dir, 'GSEA_reactome_STM+PMNs_over_PBS.csv'), stringsAsFactors = F)
se_pmn <- read.csv(file.path(gsea.dir, 'GSEA_reactome_SE+PMNs_over_PBS.csv'), stringsAsFactors = F)





# Make zoomed out heatmaps without labels
# Heatmap with SE and SE+PMN
map1 <- pathway_heatmap(se, se_pmn,
                         'SE','SE+PMN', 
                         n = 1357, rows = F, p = 1, cutrows = 1)

map2 <- pathway_heatmap(stm, stm_pmn,
                        'STM','STM+PMN', 
                        n = 1357, rows = F, p = 1, cutrows = 1)




# Save heatmaps
width = 3
height = 6
res = 500

# Save heatmap
png(filename = here('img/se_gsea_reactome_heatmap.png'),
    width = width, height = height, units = 'in', res = res)
map1
dev.off()

# Save heatmap
png(filename = here('img/stm_gsea_reactome_heatmap.png'),
    width = width, height = height, units = 'in', res = res)
map2
dev.off()



# Find pathways different in samples vs samples+PMN 
se_diff <- pathway_heatmap(se, se_pmn,
                        'SE','SE+PMN', 
                        n = 20, rows = T, p = 1, cutrows = 1)

stm_diff <- pathway_heatmap(stm, stm_pmn,
                        'STM','STM+PMN', 
                        n = 20, rows = T, p = 1, cutrows = 1)




