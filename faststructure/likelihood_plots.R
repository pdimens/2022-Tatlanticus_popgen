library(ggplot2)
library(gridExtra)
setwd(getwd())

fs_neutral <- read.table("./bft_maf01_kinless_neutral.faststructure.fs.summary", header=TRUE, quote="\"")
fs_outlier <- read.table("./bft_maf01_kinless_outlier.faststructure.fs.summary",  header=TRUE, quote="\"")

neutral <- ggplot(fs_neutral, aes(x = k, y = likelihood)) + 
  geom_segment(aes(x= k, xend= k, y=max(likelihood)+ 0.001, yend= likelihood)) +
  geom_point(size = 7, color = "lightblue") +
  scale_x_continuous(breaks = c(1:9)) +
  ggtitle("Neutral Loci (n = 2096)")

outliers <- ggplot(fs_outlier, aes(x = k, y = likelihood)) + 
  geom_segment(aes(x= k, xend= k, y=max(likelihood)+ 0.001, yend= likelihood)) +
  geom_point(size = 7, color = "orange") +
  scale_x_continuous(breaks = c(1:9)) +
  ggtitle("Outlier loci  (n = 43)")


grid.arrange(neutral, outliers, ncol=1)
