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

data <- read.csv("roh_len_histogram.csv")
#data$Length <- factor(data$Length)
ggplot(data, aes(x=Length, y=Count_min5))+
  geom_bar(stat="identity", color="darkblue", fill="lightblue")
