setwd('G:/博士/泛基因组组/PGG-project/project_pgg/01_analysis/04_eqtl/WGS_transcritp_span/new-RNA/04_analysis_eQTL/07_locuzoom_plot/')
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(scales)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

genotype_file <- 'region_SV-1_181662559.vcf'
expression_file <- '../../01_raw_expression/TPM_Liver.csv'
gene_id <- 'ENSSSCG00000005029'
snp_id <- 'SV-1:181662559'

lead_snp_id <- read.table("ENSSSCG00000005029.region_info.txt")[1,4]

# 读取数据
eqtl_data <- read.table("all_sig_pair_liver.3model_filtedbeta.txt.gz.db.gz", header = TRUE) %>% filter(pheno_id == gene_id)
ld_data <- read.table("SV-1_181662559_ld.ld.gz", header = TRUE)

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
eqtl_data <- eqtl_data[!is.na(eqtl_data$R2.group), ]


# Manhattan plot
p1 <- ggplot(eqtl_data, aes(x = Pos/1000000, y = -log10(pval_g1), color = R2.group)) +
  geom_point(size = 1, alpha = 0.7,shape=16) +
  geom_point(data = subset(eqtl_data, !is.na(sigID)),
             aes(x = Pos/1000000, y = -log10(pval_g1)), color = "red", shape = 18, size = 2) +
  geom_text(data = subset(eqtl_data, !is.na(sigID)),
            aes(label = variant_id), vjust = -0.5, color = "red", size = 3) +
  theme_bw() +
  labs(title = paste("eQTL around", gene_id),x='Position (Mb)')+
  theme(axis.text = element_text(colour = 'black', size = 12))+
  theme(panel.grid=element_blank())+
  scale_x_continuous(labels = comma)+
  scale_color_manual( name = "LD",values = c("<0.2"="#4d628e","0.2–0.4"="#5688b3","0.4–0.6"="#5cb8ce","0.6–0.8"="#f5d151",">0.8"="#dc5b48")
)+theme(legend.position = 'none') 


gene_data <- read.table("gene_features_plot.tsv.gz", header = T)

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

 #===========================================================================

# 原始读取
vcf <- fread(genotype_file, header = TRUE)

# 标准 VCF 的前几列是固定的
vcf_header <- readLines( genotype_file, n = 100)
vcf_header <- vcf_header[grep("^#CHROM", vcf_header)]
sample_names <- unlist(strsplit(vcf_header, "\t"))[-(1:9)]  # 第10列之后是样本名

# 找到目标 SNP
setDF(vcf) 
snp_row <- vcf[vcf[,3] == snp_id, ]

if (nrow(snp_row) == 0) stop(paste("SNP", snp_id, "not found."))

# 提取基因型信息（第10列开始）
genotypes <- snp_row[ , 10:ncol(snp_row)]
genotype_df <- data.frame(Sample = sample_names,
                          GT = as.character(genotypes[1, ]),
                          stringsAsFactors = FALSE)

# 提取真正的 GT 字段，如 "0/1", "1|0", etc.
genotype_df$GT <- sub(":.*", "", genotype_df$GT)
genotype_df$GT <- gsub("[|]", "/", genotype_df$GT)
genotype_df$Genotype <- recode(genotype_df$GT,
                               "0/0" = "AA",
                               "0/1" = "AB", "1/0" = "AB",
                               "1/1" = "BB")


# 读取表达矩阵（支持 .gz）
expr <- fread(expression_file)
expr <- as.data.frame(expr)
rownames(expr) <- expr[[1]]
expr <- expr[,-1]

# 提取指定基因表达
if (!(gene_id %in% rownames(expr))) stop(paste("Gene", gene_id, "not found in expression file"))
expr_df <- data.frame(Sample = colnames(expr), Expression = as.numeric(expr[gene_id, ]))

# 合并
plot_df <- inner_join(expr_df, genotype_df, by = "Sample")

# 绘图
comparisons = list(c("AA", "AB"), c("AA", "BB"),c("AB", "BB"))

counts <- plot_df %>%
  group_by(Genotype) %>%
  summarise(n = n())

# 箱线图
p3 <- ggplot(plot_df, aes(x = Genotype, y = Expression, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  theme_bw() +
  labs(x = "Genotype", y = "TPM") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = comparisons,
    label = "p.signif",
    method.args = list(p.adjust.method = "fdr"),
    hide.ns = FALSE
  ) +
  theme(axis.text = element_text(size = 12, color = 'black')) +
  scale_fill_manual(values = c('AA' = '#f5d151', 'AB' = '#4d628e', 'BB' = '#dc5b48')) +
  theme(legend.position = 'none') +
  geom_text(
    data = counts,
    aes(x = Genotype, y = min(plot_df$Expression) - 0.05 * diff(range(plot_df$Expression)), label = paste0("n=", n)),
    inherit.aes = FALSE,
    size = 4
  )
 
pdf(paste( gene_id,'locuzoom.pdf'),height = 3,width = 8)
(p1 / p2 + plot_layout(heights = c(31, 1))) |p3+ plot_layout(widths = c(0.6, 1),heights = c(222, 1)) 

dev.off()

