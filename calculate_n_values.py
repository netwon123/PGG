import argparse

def calculate_n_values(contig_lengths):
    # 计算总长度
    total_length = sum(contig_lengths)

    # 计算累计长度
    cumulative_length = 0
    n_values = {}

    # 计算 N0 到 N100 的值
    thresholds = [0.01 * i for i in range(101)]  # 从 0 到 1 的百分比
    threshold_index = 0

    for length in contig_lengths:
        cumulative_length += length

        while threshold_index < len(thresholds) and cumulative_length >= thresholds[threshold_index] * total_length:
            n_values[int(thresholds[threshold_index] * 100)] = length
            threshold_index += 1

        if threshold_index >= len(thresholds):
            break

    # 确保所有百分比都有值
    for i in range(101):
        if i not in n_values:
            n_values[i] = contig_lengths[-1]

    return n_values

def main(input_file, output_file):
    with open(input_file, 'r') as f:
        contig_lengths = [int(line.strip()) for line in f]

    # 从大到小排序
    contig_lengths.sort(reverse=True)

    # 计算 N0 到 N100 的值
    n_values = calculate_n_values(contig_lengths)

    # 写入结果到新文件
    with open(output_file, 'w') as f:
        for i in range(101):
            f.write(f"N{i}: {n_values[i]}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate N0 to N100 values for contig lengths.")
    parser.add_argument('input_file', type=str, help='Path to the input file containing contig lengths.')
    parser.add_argument('output_file', type=str, help='Path to the output file to save the results.')

    args = parser.parse_args()

    main(args.input_file, args.output_file)
