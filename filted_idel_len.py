

# 获取命令行参数
import sys
if len(sys.argv) != 3:
    sys.exit("Please provide input and output file paths as command line arguments.")

input_file = sys.argv[1]
output_file = sys.argv[2]

# 读取VCF文件并过滤结构变异
with open(input_file, 'r') as vcf_file:
    with open(output_file, 'w') as filtered_file:
        for line in vcf_file:
            if line.startswith('#'):
                # 写入注释行
                filtered_file.write(line)
            else:
                fields = line.strip().split('\t')
                info_field = fields[7].split(';')
                
                # 解析INFO字段，获取SV类型和SV长度
                svtype = None
                svlen = None
                for field in info_field:
                    if field.startswith('SVTYPE='):
                        svtype = field.split('=')[1]
                    elif field.startswith('SVLEN='):
                        svlen = int(field.split('=')[1])
                
                if svtype is not None and svlen is not None and fields[6] == 'PASS' and abs(svlen) <= 100000:
                    # 写入满足条件的变异行
                    filtered_file.write(line)
