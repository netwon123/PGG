with open('filtered_pop.vcf', 'r') as f, open('output.txt', 'w') as out_file:
    for line in f:
        if line.startswith('#'):
            continue # skip header lines
        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        info_items = fields[7].split(';')
        sv_type = info_items[1].split('=')
        sv_len = info_items[2].split('=')
        sv_type_value=sv_type[1]
        sv_len_value=sv_len[1]

        if sv_type_value == 'BND' :
            sv_len_value = 0

        out_file.write(f"{chrom}\t{pos}\t{sv_len_value}\n")

print("Output written to output.txt")
