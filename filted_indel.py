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
                
                # Check if it's a SNP
                if len(ref) == 1 and len(alt) == 1:
                    print(line.strip())
                else:
                    # Calculate the length of the indel
                    indel_length = abs(len(ref) - len(alt))

                    # Check if it's an indel and its length is less than 50 bp
                    if len(ref) != len(alt) and indel_length < 50:
                        print(line.strip())

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <vcf_file>")
        sys.exit(1)
    vcf_file = sys.argv[1]
    main(vcf_file)
