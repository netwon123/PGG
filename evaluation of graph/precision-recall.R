library(ggplot2)
library(cowplot)
library(RColorBrewer)
setwd('G:')

stat <- read.csv('precision_indel.txt', stringsAsFactors = FALSE,sep = '\t')
head(stat)

windowsFonts(Times_New_Roman=windowsFont("Arial"))

ggplot(stat, aes(x=METRIC.Recall, y=METRIC.Precision, size=Depth,shape=Reference,color=Type)) + 
  geom_point(alpha=0.5) +theme_bw()+ scale_size(range = c(4,9))+
  theme(axis.text = element_text(family = "Arial", size = 10, color = "black"),
        axis.title = element_text(family = "Arial", size = 12, color = "black"),
        legend.title= element_text(family = "Arial", size = 10, color = "black"))
