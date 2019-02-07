## Custom themes for Hill et al. 2017 manuscript
library(ggplot2)
## Load standard theme ----------------------------------------------------------
library(ggplot2)
library(grid)
theme1 <-  theme(axis.text.x = element_text(size = 24,
                                            angle = 0,
                                            hjust = 0.5,
                                            face = "bold"),
                 axis.text.y = element_text(size = 24,
                                            face = "bold",
                                            hjust = 1),
                 legend.position = "none",
		 legend.key = element_rect(fill = "white"),
                 panel.background = element_rect(fill = "white"),
                 ## Remove gid background
                 ## panel.grid.major=element_line(size=0.5,
                 ##                               color = "grey40",
                 ##                               linetype = "dashed"),
                 ## panel.grid.minor=element_line(size=0.5, # element_blank()
                 ##                               color = "grey70",
                 ##                               linetype = "dashed"),
                 plot.subtitle = element_text(size = 26, hjust = 0.5, face = "bold"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.title = element_text(size = 32,
                                           face = "bold"),
                 axis.title.y = element_text(vjust = 1.5),
                 axis.title.x = element_text(vjust = -0.5),
                 legend.title = element_text(size = 18, face = 'bold'),
		 
		 #legend.title = element_blank(),  ## For no legend titles
		 # panel.border = element_rect(fill = NA, color = "white"),
                 ## add black border to panel
		 panel.border = element_rect(fill = NA,
                                              color = "grey70",
                                              size = 1),
                 plot.title = element_text(size = 45,
                                           face = "bold",
                                             hjust = 0),
                 legend.text = element_text(size = 18,
                                            face = "bold"))


## Standard color palettes ------------------------------------------------------
library(RColorBrewer) 
color.set <- brewer.pal(n = 8, name = "Set1")
paired.set <- brewer.pal(n = 11, name = "Paired")
red.set <- brewer.pal(n = 8, name = "Reds")
green.set <- brewer.pal(n = 8, name = "Greens")
blue.set <- brewer.pal(n = 8, name = "Blues")
paired.set <- brewer.pal(n = 10, name = "Paired")

library(wesanderson)
la.set <- wes_palette("Zissou1", 5, type = "discrete")

## blank theme for importing images into ggplot2 -------------------------------
img.theme <- theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks  =  element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   plot.title  =  element_text(size = 45,
                                             face = "bold",
                                             hjust  =  0), 
                   legend.position = "none")
