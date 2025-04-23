setwd('G:/博士/泛基因组组/PGG-project/project_pgg/01_analysis/04_eqtl/WGS_transcritp_span/new-RNA/04_analysis_eQTL/04_af_tpm_tss_effect')
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(scales)
# 参数
gene_id <- "ENSSSCG00000000151"
lead_snp_id <- read.table("region_info.txt")[1,4]

# 读取数据
eqtl_data <- read.table("all_sig_pair_duodenum.3model_filtedbeta.txt.gz.db.gz", header = TRUE) %>% filter(pheno_id == gene_id)
ld_data <- read.table("SV-5_11747395_ld.ld.gz", header = TRUE)

# 提取 SNP 坐标
lead_pos <- as.numeric(gsub(".*:(\\d+)", "\\1", lead_snp_id))

# 合并 LD 到 eQTL 数据
eqtl_data <- eqtl_data %>%
  mutate(Pos = as.numeric(gsub(".*:(\\d+).*", "\\1", variant_id)),
         R2 = ifelse(variant_id == lead_snp_id, 1,
                     ld_data$R2[match(variant_id, ld_data$SNP_B)])) %>%
  mutate(R2.group = cut(R2, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
                        labels = c("<0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", ">0.8")),
         sigID = ifelse(variant_id == lead_snp_id, variant_id, NA))


# Manhattan plot
p1 <- ggplot(eqtl_data, aes(x = Pos/1000000, y = -log10(pval_g1), color = R2.group)) +
  geom_point(size = 2, alpha = 0.7,shape=16) +
  geom_point(data = subset(eqtl_data, !is.na(sigID)),
             aes(x = Pos/1000000, y = -log10(pval_g1)), color = "red", shape = 18, size = 4) +
  geom_text(data = subset(eqtl_data, !is.na(sigID)),
            aes(label = variant_id), vjust = -0.5, color = "red", size = 3) +
  theme_bw() +
  labs(title = paste("eQTL around", gene_id),x='Position (Mb)')+
  theme(axis.text = element_text(colour = 'black', size = 12))+
  theme(panel.grid=element_blank())+
  scale_x_continuous(labels = comma)+
  scale_color_manual( name = "LD",values = c("<0.2"="grey","0.2–0.4"="#5cb8ce","0.4–0.6"="#23a18c","0.6–0.8"="#9ea01e",">0.8"="#dc5b48")
) 
gene_data <- read.table("gene_features_plot.tsv", header = T)

# 只保留目标 gene

gene_data <- gene_data[gene_data$gene_id == gene_id, c("chrom", "feature", "start", "end")]

# Gene structure plot
gene_data <- gene_data %>%
  mutate(feature_color = case_when(
    feature == "CDS" ~ "#4d628e",
    feature == "five_prime_utr" ~ "#4d628e",
    feature == "three_prime_utr" ~ "#4d628e",
    TRUE ~ "#CCCCCC"
  ))

gene_data=unique(gene_data)
gene_data <- gene_data %>%
  arrange(start, end)

 

p2=ggplot(gene_data, aes(xmin = start/1000000, xmax = end/1000000, y = chrom, fill = feature_color)) +
  geom_rect(aes(xmin = start/1000000, xmax = end/1000000 ,
                ymin = as.numeric(1) - 0.5, ymax = as.numeric(1) + 0.5)) +
  geom_hline(yintercept = 1, linetype = "solid", color = "#1c4593", linewidth = 0.5) +  # 添加横线
  scale_fill_identity() +  # 使用 feature_color 列的颜色
  labs(
    title = gene_id,
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  # 隐藏 y 轴标签
    axis.ticks.y = element_blank(),  # 隐藏 y 轴刻度
    panel.grid.major.y = element_blank(),  # 隐藏 y 轴网格线
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank()
  )

 
 
pdf(paste( gene_id,'locuzoom.pdf'),height = 4,width = 6)
p1 / p2 + plot_layout(heights = c(31, 1))
dev.off()

