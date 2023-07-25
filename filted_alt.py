#filted INS without sequences in sniffles2
import argparse

def filter_ins_vcf(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
                continue

            fields = line.strip().split('\t')
            if fields[4] == '<INS>':
                continue

            f_out.write(line)

def main():
    parser = argparse.ArgumentParser(description='Filter VCF file by removing rows with <INS> in the fifth column')
    parser.add_argument('input', help='input VCF file path')
    parser.add_argument('output', help='output filtered VCF file path')

    args = parser.parse_args()
    filter_ins_vcf(args.input, args.output)

if __name__ == '__main__':
    main()
