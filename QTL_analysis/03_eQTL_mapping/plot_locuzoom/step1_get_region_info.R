#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(readr))

# 命令行参数：gene_id eqtl_file gene_anno_file [window] [margin]
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript get_region_info.R <gene_id> <eqtl_file> <gene_anno_file> [window=500000] [margin=10000]")
}

gene_id <- args[1]
eqtl_file <- args[2]
gene_anno_file <- args[3]
window <- ifelse(length(args) >= 4, as.numeric(args[4]), 1000000)
margin <- ifelse(length(args) >= 5, as.numeric(args[5]), 10000)

# 读取数据
eqtl_data <- read.table(eqtl_file, header = TRUE)
gene_anno <- read.table(gene_anno_file, header = TRUE)

# 找出 lead SNP
lead_snp <- eqtl_data %>%
  filter(pheno_id == gene_id) %>%
  slice_min(pval_g1, with_ties = FALSE)

if (nrow(lead_snp) == 0) {
  stop(paste("No SNP found for gene_id:", gene_id))
}

lead_snp_id <- lead_snp$variant_id
lead_chr <- gsub(".*-(.*):.*", "\\1", lead_snp$variant_id)
lead_pos <- as.numeric(gsub(".*:(\\d+)", "\\1", lead_snp$variant_id))

# 获取基因坐标
gene_info <- gene_anno %>% filter(id == gene_id)

if (nrow(gene_info) == 0) {
  stop(paste("No annotation found for gene_id:", gene_id))
}

gene_start <- gene_info$gene_start
gene_end <- gene_info$gene_end
chr <- gene_info$chrom

# 定义区域
region_start <- max(1, min(lead_pos - window, gene_start - margin))
region_end <- max(lead_pos + window, gene_end + margin)

# 输出区域信息
write.table(data.frame(chr, region_start, region_end, lead_snp_id),
            file=paste0(gene_id,".region_info.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

