library(ggplot2)
library(cowplot)
library(RColorBrewer)
setwd('G:')
#读取数据
stat <- read.csv('mapping.txt', stringsAsFactors = FALSE,sep = '\t')
head(stat)

windowsFonts(Times_New_Roman=windowsFont("Arial"))
pdf('mapping_ratio.pdf',width = 9,height = 3)
ggplot(data = stat,aes(x=cov,y=ratio,group = Reference,color=Reference))+
  geom_point(size=3,alpha=0.75)+
  geom_line(size=1)+
  xlab("Proportion of reads aligned perfectly")+
  ylab("Sequencing depth")+
  theme_bw()+scale_color_manual(values = c('#8ECFC9','#FA7F6F'))+
  theme(legend.position = 'bottom',axis.text = element_text(family = "Arial", size = 10, color = "black"),
        axis.title = element_text(family = "Arial", size = 12, color = "black"),
        legend.title= element_text(family = "Arial", size = 10, color = "black"),
        legend.box.background = element_rect(colour="black",size = 1))
  
dev.off()
