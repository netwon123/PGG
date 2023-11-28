import argparse

def count_sv_types(vcf_file):
    categories = ['INS', 'DEL', 'DUP', 'INV', 'BND']
    sv_types = {'unknown': 0}

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line_list = line.strip().split('\t')
            sv_type = line_list[7].split(';')[1].split('=')[1]  # 根据vcf文件的具体格式进行修改
            if sv_type not in categories:
                sv_type = 'unknown'
            if sv_type in sv_types:
                sv_types[sv_type] += 1
            else:
                sv_types[sv_type] = 1

    return sv_types

def main():
    parser = argparse.ArgumentParser(description='Count structural variant types in a VCF file')
    parser.add_argument('vcf_file', help='Input VCF file')
    args = parser.parse_args()
    
    result = count_sv_types(args.vcf_file)
    
    # 输出结果
    output = args.vcf_file + "_output.txt"
    with open(output, 'w') as f:
        f.write("File Name\tunknown\tINV\tDEL\tINS\tDUP\tBND\n")
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(args.vcf_file, result['unknown'], result['INV'], result['DEL'], result['INS'], result['DUP'], result['BND']))

if __name__ == "__main__":
    main()
