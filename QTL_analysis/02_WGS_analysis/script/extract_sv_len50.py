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

                if len(ref) == 1 and len(alt) > 49:
                    # Insertion (ref is length 1, alt length > 49)
                    print(line.strip())
                elif len(alt) == 1 and len(ref) > 49:
                    # Deletion (alt is length 1, ref length > 49)
                    print(line.strip())

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <vcf_file>")
        sys.exit(1)
    vcf_file = sys.argv[1]
    main(vcf_file)
