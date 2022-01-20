#!/usr/bin/env Rscript
depPackages = c("ggplot2", "hrbrthemes", "plotly", "htmlwidgets", "tidyr", "ggdendro", "reshape2", "grid", "heatmaply")
## Check load available dependencies, else install & load
package.check <- lapply(
  depPackages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


## Read in Cell-Variant matrix
data <- read.csv("ANG148.csv")

## Flip x and y in cell data for hierarchical clustering of cells
tdata <- data.frame(t(data[5:length(data[1,])]))
data$Id <- paste(data$Chromosome, data$Position, data$Ref, data$Alt, sep=":")
colnames(tdata) <- data$Id
tdata[is.na(tdata)] <- -1

## Generate heatmap
heatmap <- heatmaply(tdata, dendrogram = "row", seriate = "mean")
heatmap

## Save as widget
saveWidget(heatmap, file="heatmap.html")






## Get cell hierarchy
#d <- dist(tdata, method = "euclidean")
#hc <- as.dendrogram(hclust(d, method = "complete"))
#dendrogram <- ggdendrogram(data = hc, rotate = TRUE)

## Pivot to long form for ggplot
#df <- data %>% 
#  pivot_longer(cols = starts_with("Cell_"),
#               names_to = "cells",
#               names_prefix = "Cell_",
#               values_to = "genotype")

## Generate heatmap
#p <- ggplot(df, aes(Id, cells, fill= genotype)) + 
#  geom_tile() +
#  theme_ipsum() + 
#  theme(axis.text.x = element_text(angle = 90))

#grid.newpage()
#print(p, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
#print(dendrogram, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))

## Generate interactive plot
#pp <- ggplotly(p, tooltip="text")
#saveWidget(pp, file="heatmap.html")


