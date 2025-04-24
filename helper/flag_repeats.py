from collections import defaultdict

import argparse
import pandas as pd

# Load BED file
def load_bed(filepath):
    bed_data = defaultdict(list)
    with open(filepath) as f:
        for line in f:
            fields = line.strip().split('\t')
            read_id, start, end, name, strand = fields
            start, end = int(start), int(end)
            d4z4_label = name.split(":")[1]
            bed_data[(read_id, strand)].append((start, end, d4z4_label))
    return bed_data

# Analyze for gaps and overlaps
def analyse_bed(bed_data):
    results = {}
    for (read_id, strand), entries in bed_data.items():
        entries.sort()  # sort by start
        gap_found = False
        gap_distance = "NA"

        # Gaps between consecutive
        for i in range(1, len(entries)):
            prev_end = entries[i-1][1]
            curr_start = entries[i][0]
            gap = curr_start - prev_end
            if gap > 400:
                gap_found = True
                gap_distance = gap
                break

        # Overlaps between any pair (non-consecutive)
        overlap_found = False
        overlapping_labels = []
        overlapping_coords = []

        for i in range(len(entries)):
            for j in range(i+1, len(entries)):
                s1, e1, label1 = entries[i]
                s2, e2, label2 = entries[j]
                if s1 < e2 and s2 < e1:
                    overlap_found = True
                    overlapping_labels.extend([label1, label2])
                    overlapping_coords.extend([f"{s1}-{e1}", f"{s2}-{e2}"])

        # remove duplicates in overlaps
        overlapping_labels = list(set(overlapping_labels))
        overlapping_coords = list(set(overlapping_coords))

        results[(read_id, strand)] = {
            'gaps': gap_found,
            'gap_distance': gap_distance,
            'overlaps': overlap_found,
            'overlapping_repeats': overlapping_labels if overlap_found else "NA",
            'overlapping_repeats_coords': overlapping_coords if overlap_found else "NA"
        }
    return results

# Merge with main TSV
def merge_with_main(main_tsv_path, bed_results):
    main_df = pd.read_csv(main_tsv_path, sep='\t')

    def get_val(read_id, strand, key):
        return bed_results.get((read_id, strand), {}).get(key, "NA")

    main_df['gaps'] = main_df.apply(lambda row: get_val(row['ReadID'], row['strand'], 'gaps'), axis=1)
    main_df['gap_distance'] = main_df.apply(lambda row: get_val(row['ReadID'], row['strand'], 'gap_distance'), axis=1)
    main_df['overlaps'] = main_df.apply(lambda row: get_val(row['ReadID'], row['strand'], 'overlaps'), axis=1)
    main_df['overlapping_repeats'] = main_df.apply(lambda row: get_val(row['ReadID'], row['strand'], 'overlapping_repeats'), axis=1)
    main_df['overlapping_repeats_coords'] = main_df.apply(lambda row: get_val(row['ReadID'], row['strand'], 'overlapping_repeats_coords'), axis=1)

    main_df = main_df.astype(str).replace("nan", "NA")

    main_df.to_csv(main_tsv_path, sep='\t', index=False)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Flag reads with overlapping reads or big gaps.")

    parser.add_argument("--main_tsv", required=True, help="Main output TSV file to update flags.")
    parser.add_argument("--repeats_bed", required=True, help="Output BED file for d4z4 repeats.")

    args = parser.parse_args()
    
    bed_data = load_bed(args.repeats_bed)
    bed_results = analyse_bed(bed_data)
    merge_with_main(args.main_tsv, bed_results)
