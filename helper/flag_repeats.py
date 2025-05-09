from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

import argparse
import pandas as pd

from alignment_order_clipping import get_original_sequence

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



def check_re_vs_haplotype(row,  re_threshold=0.7, ratio_threshold=50):
    haplotype = row['Haplotype']
    num_xapi = row['XapI_Sensitive_Repeats']
    num_blni = row['BlnI_Sensitive_Repeats']
    xapi_ratio = row['XapI_RE_Ratio_(%)']
    blni_ratio = row['BlnI_RE_Ratio_(%)']
    copy_number = row['MappedEstimatedCopies']

    # print(type(chr), type(haplotype), type(xapi_ratio), type(blni_ratio))

    if type(haplotype) == float and pd.isna(haplotype):
        return "NA"  # Skip if haplotype is NaN
    
    flag = []

    if "4" in haplotype and blni_ratio > ratio_threshold:
        print(f"{row['ReadID']} classified as haplotype 4q but blni is higher than xapi, have {copy_number} repeats")
        flag.append("wrong_re")
    elif "10" in haplotype and xapi_ratio > ratio_threshold:
        print(f"{row['ReadID']} classified as haplotype 10q but xapi is higher than blni, have {copy_number} repeats")
        flag.append("wrong_re")
    
    if copy_number > 4:
        dynamic_threshold = int(re_threshold * copy_number)
        if ("4" in haplotype) and (num_xapi < dynamic_threshold):
            flag.append("low_xapi")
            print(f"{row['ReadID']} is {haplotype} but only {num_xapi} xapi for {copy_number} repeats")
        elif ("10" in haplotype) and (num_blni < dynamic_threshold):
            flag.append("low_blni")
            print(f"{row['ReadID']} is {haplotype} but only {num_blni} blni for {copy_number} repeats")

    return ", ".join(flag) if flag else "NA"

# def get_original_sequence(fasta_file, read_id):
    # """
    # Get the original sequence of a read from a FASTA file.
    # """
    # # Parse the input FASTA file to get the original read sequences
    # original_reads = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}

    # if read_id not in original_reads:
    #     print(f"Warning: Read {read_id} not found in input FASTA.")
    #     return

    # return original_reads[read_id]

# Function to get the start and end based on strand
def get_distal_d4z4_coords(df):
    if df.empty:
        return None, None, None

    repeat_threshold = 0.8 * 3300  # 80% of 3300 bp

    # Extract numeric values safely as a Series
    df = df.copy()  # Avoid modifying the original dataframe
    df["d4z4_num"] = df["name"].str.extract(r"d4z4_(\d+)$")[0].astype(float)

    strand = df["strand"].iloc[0]

    if strand == "+":  
        # Get the row with max d4z4_num
        distal_d4z4 = df.loc[df["d4z4_num"].idxmax()]

        # skip small d4z4 copies
        # size = distal_d4z4["BlockSizes"].split(",")[0]
        size = distal_d4z4["end"] - distal_d4z4["start"]
        if int(size) < repeat_threshold:
            # If the first block size is 0, use the second largest
            distal_copy = distal_d4z4["d4z4_num"] - 1
            distal_d4z4 = df[df["d4z4_num"] == distal_copy].iloc[0] if (df["d4z4_num"] == distal_copy).any() else df.loc[df["d4z4_num"].idxmax()]
    else:  
        # Get d4z4_1, or fallback to min if missing
        distal_d4z4 = df[df["d4z4_num"] == 1].iloc[0] if (df["d4z4_num"] == 1).any() else df.loc[df["d4z4_num"].idxmin()]

        # size = distal_d4z4["BlockSizes"].split(",")[0]
        size = distal_d4z4["end"] - distal_d4z4["start"]
        if int(size) < repeat_threshold:
            # If the first block size is 0, use the second largest
            distal_copy = distal_d4z4["d4z4_num"] + 1
            distal_d4z4 = df[df["d4z4_num"] == distal_copy].iloc[0] if (df["d4z4_num"] == distal_copy).any() else df.loc[df["d4z4_num"].idxmin()]

    # print(df["chrom"].iloc[0], strand, distal_d4z4["start"], distal_d4z4["end"])
    return strand, distal_d4z4["start"], distal_d4z4["end"]

def flag_read_label(main_tsv, repeats_bed, fasta_file):
    """
    1. get unclassified read id
    2. get original seq
    3. for each read, get the repeats from bed file
    4. Get the distal coords of read by taking the copy with the highest index for + strand and lowest index for - strand
    5. Check if coordinates in read from bed file to the end of the read is more than 3.3kb
    6. If yes, change the read label to "Partial Distal 4qB or 10qB", else keep it unclassified
    """

    inconsistent_flag_dict = {}
    
    # put main tsv to df and filter only lines with "Unclassified" in read label
    main_tsv_df = pd.read_csv(main_tsv, sep="\t")
    main_tsv_df_unclassified = main_tsv_df[main_tsv_df["ReadLabel"].str.contains("Unclassified")]
    # put read ids of unclassified reads to a list
    unclassified_read_ids = set(main_tsv_df_unclassified["ReadID"])

    for read_id in unclassified_read_ids:
        # get the original sequence of the read
        original_sequence = get_original_sequence(fasta_file, read_id)

        # main tsv line of read id with the same strand as strand from d4z4_bed_df
        # curr_read_df = main_tsv_df[(main_tsv_df["ReadID"] == read_id) & (main_tsv_df["strand"] == strand)]
        curr_read_df = main_tsv_df[(main_tsv_df["ReadID"] == read_id)]
        if curr_read_df["duplex"].values[0] == True:
            continue

        # get the d4z4 repeat coordinates from the bed file
        d4z4_bed_df = pd.read_csv(repeats_bed, sep="\t", header=None)
        d4z4_bed_df.columns = ["chrom", "start", "end", "name", "strand"]
        d4z4_bed_df = d4z4_bed_df[d4z4_bed_df["name"].str.contains(read_id)]

        if d4z4_bed_df.empty:
            continue

        strand, distal_start, distal_end = get_distal_d4z4_coords(d4z4_bed_df)

        if strand == "+":
            # check if the distance from distal start to the end of the read is more than 3.3kb
            # if chr 4 then 4qB, if chr 10 then 10qB
            if (len(original_sequence) - distal_end) > 3300:
                inconsistent_flag_dict.setdefault((read_id, strand), []).append("inconsistent_reads")

        else:
            # check if the distance from distal start to the start of the read is more than 3.3kb
            if distal_start > 3300:
                inconsistent_flag_dict.setdefault((read_id, strand), []).append("inconsistent_reads")

    return inconsistent_flag_dict

# Merge with main TSV
def merge_with_main(main_tsv_path, bed_results, hybrid_tsv_path, inconsistent_flag_dict):
    main_df = pd.read_csv(main_tsv_path, sep='\t')

    def get_val(read_id, strand, key):
        return bed_results.get((read_id, strand), {}).get(key, "NA")

    main_df['gaps'] = main_df.apply(lambda row: get_val(row['ReadID'], row['strand'], 'gaps'), axis=1)
    main_df['gap_distance'] = main_df.apply(lambda row: get_val(row['ReadID'], row['strand'], 'gap_distance'), axis=1)
    main_df['overlaps'] = main_df.apply(lambda row: get_val(row['ReadID'], row['strand'], 'overlaps'), axis=1)
    main_df['overlapping_repeats'] = main_df.apply(lambda row: get_val(row['ReadID'], row['strand'], 'overlapping_repeats'), axis=1)
    main_df['overlapping_repeats_coords'] = main_df.apply(lambda row: get_val(row['ReadID'], row['strand'], 'overlapping_repeats_coords'), axis=1)

    # Apply RE vs haplotype check
    main_df['misclassification_flags'] = main_df.apply(check_re_vs_haplotype, axis=1)

    # Load hybrid tsv misclassification file
    hybrid_df = pd.read_csv(hybrid_tsv_path, sep=' ', header=None, names=["ReadID", "PrimaryChr", "Supplementary"])

    # Build a mapping of ReadID to sa_flag
    sa_flag_dict = {}
    for _, row in hybrid_df.iterrows():
        read_id = row["ReadID"]
        primary_chr = row["PrimaryChr"]
        # print(row)
        # print(row["Supplementary"], type(row["Supplementary"]))
        supp_chr = row["Supplementary"].split(",")[0]  # grab first part of supplementary field
        sa_flag = f"sa_{primary_chr}_{supp_chr}"
        sa_flag_dict.setdefault(read_id, []).append(sa_flag)  # in case multiple per read

    # Combine misclassification flags
    def combine_flags(existing_flag, new_flags):
        if not new_flags:
            return existing_flag
        if existing_flag == "NA":
            return ", ".join(new_flags)
        return existing_flag + ", " + ", ".join(new_flags)

    # Apply SA misclassification flags
    main_df["misclassification_flags"] = main_df.apply(
        lambda row: combine_flags(row["misclassification_flags"], sa_flag_dict.get(row["ReadID"], [])),
        axis=1
    )

    main_df["misclassification_flags"] = main_df.apply(
        lambda row: combine_flags(row["misclassification_flags"], inconsistent_flag_dict.get((row["ReadID"], row["strand"]), [])),
        axis=1
    )

    main_df = main_df.astype(str).replace("nan", "NA")

    main_df.to_csv(main_tsv_path, sep='\t', index=False)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Flag reads with overlapping reads or big gaps.")

    parser.add_argument("--main_tsv", required=True, help="Main output TSV file to update flags.")
    parser.add_argument("--repeats_bed", required=True, help="Output BED file for d4z4 repeats.")
    parser.add_argument("--hybrid", required=True, help="Potential reads with hybrid file.")
    parser.add_argument("--fasta", required=True, help="Input FASTA file with full reads.")

    args = parser.parse_args()
    
    bed_data = load_bed(args.repeats_bed)
    bed_results = analyse_bed(bed_data)

    inconsistent_flag_dict = flag_read_label(args.main_tsv, args.repeats_bed, args.fasta)

    merge_with_main(args.main_tsv, bed_results, args.hybrid, inconsistent_flag_dict)
