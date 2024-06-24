setwd('G:/博士/泛基因组组/PGG-project/analysis/population_structure/004_netview/r2_max/')
library(data.table)
library(tidyverse)
dat.ld<-fread("LD_pop_bnd.vcf.simple.vcf.gz_snp.vcf.gz.ld_onlysv.txt_3.txt.tran.txt_max.txt")
colnames(dat.ld)=c('SNP_A','SNP_B','R2')
#dat.ld%>%filter(abs(BP_A-BP_B)<=50000) -> dat.ld#论文里记得写

dat.ld %>% head() 

dat.ld.new<-dat.ld %>% 
  mutate(group=paste0(str_extract(SNP_A,pattern = "[A-Za-z]+$"),
                      "_",
                      str_extract(SNP_B,pattern = "[A-Za-z]+$")))

dat.ld.new %>% 
  mutate(new_group=case_when(
    group == "SNP_SV"|group=="SV_SNP" ~ "SNP_SV",
    group == "INDEL_SV"|group=="SV_INDEL" ~ "INDEL_SV",
    group == "SNP_SNP" ~ "SNP_SNP",
    group == "SV_SV" ~ "SV_SV",
    TRUE ~ "OTHERS"
  )) -> dat.ld.new

#write.table(dat.ld.new,file ='LD_snp_indel_chr1.ld.txt' ,quote = FALSE,row.names = FALSE)

p1=ggplot(data=dat.ld.new %>% filter(new_group=="SNP_SV"),
          aes(x=R2))+
  geom_histogram(
    bins = 150,
    color="#82b0d2",
    fill="#82b0d2")+
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "black") +  # 添加竖直虚线
  theme_classic()+
  theme(panel.grid = element_blank())+
  labs(title = "BND")

#--------------------
dat.ld2<-fread("LD_pop_del.vcf.simple.vcf.gz_snp.vcf.gz.ld_onlysv.txt_3.txt.tran.txt_max.txt")
colnames(dat.ld2)=c('SNP_A','SNP_B','R2')
dat.ld2 %>% head() 

dat.ld2.new<-dat.ld2 %>% 
  mutate(group=paste0(str_extract(SNP_A,pattern = "[A-Za-z]+$"),
                      "_",
                      str_extract(SNP_B,pattern = "[A-Za-z]+$")))

dat.ld2.new %>% 
  mutate(new_group=case_when(
    group == "SNP_SV"|group=="SV_SNP" ~ "SNP_SV",
    group == "INDEL_SV"|group=="SV_INDEL" ~ "INDEL_SV",
    group == "SNP_SNP" ~ "SNP_SNP",
    group == "SV_SV" ~ "SV_SV",
    TRUE ~ "OTHERS"
  )) -> dat.ld2.new



p2=ggplot(data=dat.ld2.new %>% filter(new_group=="SNP_SV"),
          aes(x=R2))+
  geom_histogram(
    bins = 150,
    color="#82b0d2",
    fill="#82b0d2")+
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "black") +  # 添加竖直虚线
  theme_classic()+
  theme(panel.grid = element_blank())+
  labs(title = "DEL")
#-----------------

dat.ld3<-fread("LD_pop_dup.vcf.simple.vcf.gz_snp.vcf.gz.ld_onlysv.txt_3.txt.tran.txt_max.txt")
colnames(dat.ld3)=c('SNP_A','SNP_B','R2')
dat.ld3 %>% head() 

dat.ld3.new<-dat.ld3 %>% 
  mutate(group=paste0(str_extract(SNP_A,pattern = "[A-Za-z]+$"),
                      "_",
                      str_extract(SNP_B,pattern = "[A-Za-z]+$")))

dat.ld3.new %>% 
  mutate(new_group=case_when(
    group == "SNP_SV"|group=="SV_SNP" ~ "SNP_SV",
    group == "INDEL_SV"|group=="SV_INDEL" ~ "INDEL_SV",
    group == "SNP_SNP" ~ "SNP_SNP",
    group == "SV_SV" ~ "SV_SV",
    TRUE ~ "OTHERS"
  )) -> dat.ld3.new



p3=ggplot(data=dat.ld3.new %>% filter(new_group=="SNP_SV"),
          aes(x=R2))+
  geom_histogram(
    bins = 150,
    color="#82b0d2",
    fill="#82b0d2")+
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "black") +  # 添加竖直虚线
  theme_classic()+
  theme(panel.grid = element_blank())+
  labs(title = "DUP")

#---------
dat.ld4<-fread("LD_pop_ins.vcf.simple.vcf.gz_snp.vcf.gz.ld_onlysv.tran.txt_max.txt")
colnames(dat.ld4)=c('SNP_A','SNP_B','R2')
dat.ld4 %>% head() 

dat.ld4.new<-dat.ld4 %>% 
  mutate(group=paste0(str_extract(SNP_A,pattern = "[A-Za-z]+$"),
                      "_",
                      str_extract(SNP_B,pattern = "[A-Za-z]+$")))

dat.ld4.new %>% 
  mutate(new_group=case_when(
    group == "SNP_SV"|group=="SV_SNP" ~ "SNP_SV",
    group == "INDEL_SV"|group=="SV_INDEL" ~ "INDEL_SV",
    group == "SNP_SNP" ~ "SNP_SNP",
    group == "SV_SV" ~ "SV_SV",
    TRUE ~ "OTHERS"
  )) -> dat.ld4.new



p4=ggplot(data=dat.ld4.new %>% filter(new_group=="SNP_SV"),
          aes(x=R2))+
  geom_histogram(
    bins = 150,
    color="#82b0d2",
    fill="#82b0d2")+
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "black") +  # 添加竖直虚线
  theme_classic()+
  theme(panel.grid = element_blank())+
  labs(title = "INS")

#-----------
dat.ld5<-fread("LD_pop_inv.vcf.simple.vcf.gz_snp.vcf.gz.ld_onlysv.txt_3.txt.tran.txt_max.txt")
colnames(dat.ld5)=c('SNP_A','SNP_B','R2')
dat.ld5 %>% head() 

dat.ld5.new<-dat.ld5 %>% 
  mutate(group=paste0(str_extract(SNP_A,pattern = "[A-Za-z]+$"),
                      "_",
                      str_extract(SNP_B,pattern = "[A-Za-z]+$")))

dat.ld5.new %>% 
  mutate(new_group=case_when(
    group == "SNP_SV"|group=="SV_SNP" ~ "SNP_SV",
    group == "INDEL_SV"|group=="SV_INDEL" ~ "INDEL_SV",
    group == "SNP_SNP" ~ "SNP_SNP",
    group == "SV_SV" ~ "SV_SV",
    TRUE ~ "OTHERS"
  )) -> dat.ld5.new



p5=ggplot(data=dat.ld5.new %>% filter(new_group=="SNP_SV"),
          aes(x=R2))+
  geom_histogram(
    bins = 150,
    color="#82b0d2",
    fill="#82b0d2")+
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "black") +  # 添加竖直虚线
  theme_classic()+
  theme(panel.grid = element_blank())+
  labs(title = "INV")+theme_void()
#-------------
dat.ld6<-fread("snp-snp_max.txt")
colnames(dat.ld6)=c('SNP_A','SNP_B','R2')
dat.ld6 %>% head() 

dat.ld6.new<-dat.ld6 %>% 
  mutate(group=paste0(str_extract(SNP_A,pattern = "[A-Za-z]+$"),
                      "_",
                      str_extract(SNP_B,pattern = "[A-Za-z]+$")))

dat.ld6.new %>% 
  mutate(new_group=case_when(
    group == "SNP_SV"|group=="SV_SNP" ~ "SNP_SV",
    group == "INDEL_SV"|group=="SV_INDEL" ~ "INDEL_SV",
    group == "SNP_SNP" ~ "SNP_SNP",
    group == "SV_SV" ~ "SV_SV",
    TRUE ~ "OTHERS"
  )) -> dat.ld6.new



p6=ggplot(data=dat.ld6.new %>% filter(new_group=="SNP_SNP"),
          aes(x=R2))+
  geom_histogram(
    bins = 150,
    color="#82b0d2",
    fill="#82b0d2")+
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "black") +  # 添加竖直虚
  theme_classic()+
  theme(panel.grid = element_blank())+
  labs(title = "SNP")


p6+p1+p2+p3+p4+p5
