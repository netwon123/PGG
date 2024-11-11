import pysam
import numpy as np

def calculate_average_depth(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    depths = []

    for record in vcf.fetch():
        dp = record.info.get('DP')
        if dp is not None:
            depths.append(dp)


    if depths:
        average_depth = np.mean(depths)
        print(" (INFO/DP):", average_depth)
    else:
        print("no (INFO/DP)")


vcf_file = "panel_2991.vcf"
calculate_average_depth(vcf_file)
