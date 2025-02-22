# 设置最大并行任务数（根据系统资源调整）
MAX_JOBS=16
BACKGROUND_JOBS=0

for CHR in {1..22}; do
    # 检查输入文件是否存在
    INPUT_FILE="chr${CHR}.phastCons100way.wigFix.bed"
    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "Error: Input file ${INPUT_FILE} not found!"
        continue
    fi

    # 解压并分割文件
    {
        cat "${INPUT_FILE}" | split -l 1000000 -d -a 3 - "chr${CHR}.phastCons100way.hg38.chunk_"

        # 重命名文件为 .bed 后缀
        for FILE in ./chr${CHR}.phastCons100way.hg38.chunk_*; do
            if [[ ! "${FILE}" =~ \.bed$ ]]; then
                mv "${FILE}" "${FILE}.bed"
            fi
        done

        # 压缩文件为 .bed.gz
        bgzip ./chr${CHR}.phastCons100way.hg38.chunk_*.bed
    } &

    # 控制并行任务数
    BACKGROUND_JOBS=$((BACKGROUND_JOBS + 1))
    if [[ "${BACKGROUND_JOBS}" -ge "${MAX_JOBS}" ]]; then
        wait -n  # 等待任意一个后台任务完成
        BACKGROUND_JOBS=$((BACKGROUND_JOBS - 1))
    fi
done

# 等待所有后台任务完成
wait
echo "All jobs completed!"

