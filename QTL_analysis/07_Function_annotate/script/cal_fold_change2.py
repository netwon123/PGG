import pandas as pd
import pybedtools
import os
import argparse

def calculate_enrichment(sv_file, repeat_folder, genome_size=2391714235):
    sv_bed = pybedtools.BedTool(sv_file)

    sv_total_length = sum([interval.length for interval in sv_bed])

    enrichment_results = pd.DataFrame(columns=["Repeat Type", "Enrichment Factor"])

    for repeat_file in os.listdir(repeat_folder):
        repeat_path = os.path.join(repeat_folder, repeat_file)

        repeat_bed = pybedtools.BedTool(repeat_path)

        overlap = sv_bed.intersect(repeat_bed, wa=True,f=0.1)
        overlap_length = sum([interval.length for interval in overlap])

        repeat_total_length = sum([interval.length for interval in repeat_bed])

        genome_fraction = repeat_total_length / genome_size
        enrichment_factor = (overlap_length / sv_total_length) / genome_fraction if genome_fraction > 0 else 0

        enrichment_results = enrichment_results.append({
            "Repeat Type": repeat_file,
            "Enrichment Factor": enrichment_factor
        }, ignore_index=True)

    return enrichment_results

def main():
    parser = argparse.ArgumentParser(description="Calculate enrichment factors for repeat regions in SVs.")
    parser.add_argument("sv_file", type=str, help="Path to the SV BED file")
    parser.add_argument("repeat_folder", type=str, help="Path to the folder containing repeat BED files")
    parser.add_argument("--genome_size", type=int, default=2391714235, help="Genome size (default: 2391714235)")

    args = parser.parse_args()

    results = calculate_enrichment(args.sv_file, args.repeat_folder, args.genome_size)
    print(results)
    results.to_csv('output_enrichment.csv',index=False)

if __name__ == "__main__":
    main()
