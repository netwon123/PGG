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

                # Skip SNP (ref and alt of length 1)
                if len(ref) == 1 and len(alt) == 1:
                    continue
                else:
                    # Calculate the length of the indel
                    indel_length = abs(len(ref) - len(alt))

                    # Skip indels that are either less than 50 bp or greater than 100,000 bp
                    if indel_length < 50 or indel_length > 100000:
                        continue

                    # Print the indels with length between 50 and 100,000 bp
                    print(line.strip())

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <vcf_file>")
        sys.exit(1)
    vcf_file = sys.argv[1]
    main(vcf_file)
