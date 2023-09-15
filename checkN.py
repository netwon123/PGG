def extract_rows_with_single_n(vcf_file):
    rows_with_single_n = []
    with open(vcf_file, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                col3 = fields[2]
                col4 = fields[3]
                if col3.count("N") == 1 or col4.count("N") == 1:
                    rows_with_single_n.append(line)
    return rows_with_single_n

vcf_file = "MZ_cas.assembly.vcf.filted.vcf.filtedins.vcf_pass.vcf.out.vcf"
rows_with_single_n = extract_rows_with_single_n(vcf_file)
for row in rows_with_single_n:
    print(row)
