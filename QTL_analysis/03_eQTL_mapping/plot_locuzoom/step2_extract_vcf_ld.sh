#!/bin/bash

# 用法: ./run_ld_analysis.sh region_info.txt

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 region_info.txt"
  exit 1
fi

input_file=$1

# 循环读取每一行
while read -r chr start end lead_id; do
  echo "Processing: chr=$chr, start=$start, end=$end, lead_id=$lead_id"

  region_vcf="region_${lead_id}.vcf"

  # 提取区域 VCF
  bcftools view -r ${chr}:${start}-${end} ../../../../../02_genotype/joint.vcf.gz \
    -o ${region_vcf} -Ov --threads 8

  # 计算 LD
  plink --vcf ${region_vcf} \
        --ld-snp ${lead_id} \
        --ld-window-kb 500 \
        --ld-window 99999 \
        --ld-window-r2 0 \
        --out ${lead_id}_ld \
        --r2 gz \
        --threads 16

  echo "Finished processing ${lead_id}"
done < "$input_file"
