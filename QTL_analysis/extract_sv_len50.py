import sys

def main(vcf_file):
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                # Print header lines
                print(line.strip())
            elif line.startswith('#'):
                # Print column headers
                print(line.strip())
            else:
                fields = line.strip().split('\t')
                ref = fields[3]
                alt = fields[4]
                indel_length = abs(len(ref) - len(alt))
                # Check if it's a sv
                if len(ref) != len(alt) and indel_length > 49:
                    print(line.strip())


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <vcf_file>")
        sys.exit(1)
    vcf_file = sys.argv[1]
    main(vcf_file)
