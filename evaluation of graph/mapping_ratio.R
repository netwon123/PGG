library(ggplot2)
library(cowplot)
library(RColorBrewer)
setwd('G:/博士/泛基因组组/PGG-project/analysis/eval_graph/final/')
library(patchwork)
stat1 <- read.csv('precision_indel.txt', stringsAsFactors = FALSE,sep = '\t')
head(stat)

stat3 <- read.csv('mapping.txt', stringsAsFactors = FALSE,sep = '\t')

stat2 <- read.csv('sv-eval.txt', stringsAsFactors = FALSE,sep = '\t')
head(stat)


p1=ggplot(stat, aes(x=METRIC.Recall, y=METRIC.Precision, size=Depth,shape=Type,color=Reference)) + 
  geom_point(alpha=0.5) +theme_bw()+ scale_size(range = c(3,8))+
  theme(axis.text = element_text( size = 10, color = "black"),
        axis.title = element_text( size = 12, color = "black"),
        legend.title= element_text( size = 10, color = "black"))+scale_shape_manual(values = c(16,15))




p2=ggplot(data = stat2,aes(x=cov,y=F1,group = Reference,color=Reference))+
  geom_point(size=3,alpha=0.75)+
  geom_line(size=1)+
  xlab("Sequencing depth")+
  ylab("F1 score")+
  theme_bw()+scale_color_manual(values = c('#8ECFC9','#FA7F6F'))+
  theme(legend.position = "none")





p3=ggplot(data = stat3,aes(x=cov,y=ratio,group = Reference,color=Reference))+
  geom_point(size=3,alpha=0.75)+
  geom_line(size=1)+
  xlab("Proportion of reads aligned perfectly")+
  ylab("Sequencing depth")+
  theme_bw()+scale_color_manual(values = c('#8ECFC9','#FA7F6F'))+
  theme(legend.position = 'none')
p2 | p3 | p1
