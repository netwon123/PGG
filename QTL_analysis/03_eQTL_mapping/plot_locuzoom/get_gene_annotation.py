import pandas as pd
import numpy as np
from collections import defaultdict
import gzip

def gtf_to_filtered_features(annotation_gtf, features=('CDS', 'five_prime_utr', 'three_prime_utr'), exclude_chrs=[]):
    """Extract specified GTF features and output gene_id, chrom, feature, start, end"""
    gene_id_list = []
    chrom = []
    feature_list = []
    start = []
    end = []

    opener = gzip.open(annotation_gtf, 'rt') if annotation_gtf.endswith('.gz') else open(annotation_gtf, 'r')

    with opener as gtf:
        for row in gtf:
            row = row.strip().split('\t')
            if row[0].startswith('#') or row[2] not in features:
                continue  # skip headers and irrelevant features

            chr_name = row[0]
            if chr_name in exclude_chrs:
                continue

            feature_type = row[2]
            start_pos = int(row[3]) - 1  # BED-style start (0-based)
            end_pos = int(row[4])

            # parse attributes
            attributes = defaultdict(str)
            for attr in row[8].strip().replace('"', '').split(';'):
                if attr.strip():
                    kv = attr.strip().split(' ', 1)
                    if len(kv) == 2:
                        attributes[kv[0]] = kv[1]

            gene_id = attributes.get('gene_id', 'NA')

            # append to lists
            gene_id_list.append(gene_id)
            chrom.append(chr_name)
            feature_list.append(feature_type)
            start.append(start_pos)
            end.append(end_pos)

    # create DataFrame
    bed_df = pd.DataFrame({
        'gene_id': gene_id_list,
        'chrom': chrom,
        'feature': feature_list,
        'start': start,
        'end': end
    })

    return bed_df
df = gtf_to_filtered_features('Sus_scrofa.Sscrofa11.1.113.chr.gtf')
df.to_csv('gene_features_plot.tsv', sep='\t', index=False)
