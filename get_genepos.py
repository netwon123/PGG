import pandas as pd
import numpy as np
from collections import defaultdict
import subprocess
import gzip

def gtf_to_tss_bed(annotation_gtf, feature='gene', exclude_chrs=[], phenotype_id='gene_id'):
    """Parse genes and TSSs from GTF and return DataFrame for BED output"""
    chrom = []
    start = []
    end = []
    gene_id = []
    gene_name = []

    if annotation_gtf.endswith('.gz'):
        opener = gzip.open(annotation_gtf, 'rt')
    else:
        opener = open(annotation_gtf, 'r')

    with opener as gtf:
        for row in gtf:
            row = row.strip().split('\t')
            if row[0][0] == '#' or row[2] != feature: continue # skip header
            chrom.append(row[0])

            # TSS: gene start (0-based coordinates for BED)
            if row[6] == '+':
                start.append(np.int64(row[3]))
                end.append(np.int64(row[4]))
            elif row[6] == '-':
                start.append(np.int64(row[3]))  # last base of gene
                end.append(np.int64(row[4]))
            else:
                raise ValueError('Strand not specified.')

            attributes = defaultdict()
            for a in row[8].replace('"', '').split(';')[:-1]:
                kv = a.strip().split(' ')
                if kv[0]!='tag':
                    attributes[kv[0]] = kv[1]
                else:
                    attributes.setdefault('tags', []).append(kv[1])

            gene_id.append(attributes['gene_id'])
            

    if phenotype_id == 'gene_id':
        bed_df = pd.DataFrame(data={'chr':chrom, 'start':start, 'end':end, 'gene_id':gene_id}, columns=['chr', 'start', 'end', 'gene_id'], index=gene_id)
    elif phenotype_id == 'gene_name':
        bed_df = pd.DataFrame(data={'chr':chrom, 'start':start, 'end':end, 'gene_id':gene_name}, columns=['chr', 'start', 'end', 'gene_id'], index=gene_name)
    # drop rows corresponding to excluded chromosomes
    mask = np.ones(len(chrom), dtype=bool)
    for k in exclude_chrs:
        mask = mask & (bed_df['chr']!=k)
    bed_df = bed_df[mask]

    # sort by start position
    bed_df = bed_df.groupby('chr', sort=False, group_keys=False).apply(lambda x: x.sort_values('start'))

    return bed_df

tss=gtf_to_tss_bed('Sus_scrofa.Sscrofa11.1.100.gtf', feature='gene', exclude_chrs=[], phenotype_id='gene_id')
tss.to_csv('gene.bed')
